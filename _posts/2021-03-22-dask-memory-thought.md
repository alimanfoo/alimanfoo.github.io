---
layout: post
title:  "Dask thought: distributed clusters and worker memory/storage"
description: "Some thoughts about Dask clusters and memory management."
image: /assets/2021-03-22-dask-memory-thought.png
date:   2021-03-22
categories: dask
---

*When a Dask distributed cluster runs out of memory, it can be
tricky to figure out why, and even harder to resolve. What if the
system that cluster workers use to store intermediate results was
pluggable and allowed different implementations? And what if you could
choose to plug in a storage implementation that allowed persistence of
results between cluster restarts?*

<hr/>

One of the current challenges working with Dask and using the
distributed scheduler on a cluster is that the intermediate results in
a task graph are held in memory on the worker nodes in the
cluster. This can present a challenge to the user because if workers
start running out of memory during a computation, it can be for any of
multiple reasons. E.g., it could be because some individual tasks are
consuming too much memory during execution, perhaps because chunks are
too large or the wrong shape or the task is allocating memory
unnecessarily. It could also be because the ordering of task execution
is sub-optimal, and workers are being asked to hold onto lots of
intermediate results. Or some combination of these and/or other
factors.

If workers start running out of memory, then they can start "spilling
to disk", which I understand to mean that they write out intermediate
results to their local file system. You can see this in the status
dashboard, the blue bars in the memory chart start turning
orange. This sounds like a good idea in principle, and should allow
computations to proceed, but in my experience when you start seeing
these bars turning orange, it is a precursor to everything grinding to
a halt. I.e., distributed clusters start behaving in a very
pathological way when workers start spilling to disk.

Ordering task execution in general is a difficult problem, and Dask
generally does an admirable job of coming up with an order that
manages memory and balances the other considerations. However,
sometimes it does get it wrong, and in some surprisingly simple
cases. One case that has come up is rechunking, i.e., when you want to
transform the chunks of some array to a different shape. This can come
up when you may have stored data natively in one chunk shape (e.g.,
tall and skinny) but a computation requires them to be in a different
shape (e.g., short and fat). Several people have hit situations where
rechunking causes a distributed cluster to run out of memory (e.g.,
[here](https://discourse.pangeo.io/t/best-practices-to-go-from-1000s-of-netcdf-files-to-analyses-on-a-hpc-cluster/588/24),
[here](https://github.com/dask/dask/issues/5105),
[here](https://github.com/dask/dask/issues/6745)), and the underlying
cause seems to be sub-optimal task ordering. Some of this has been
resolved with [improved heuristics around task
ordering](https://github.com/dask/dask/pull/6779), but the solution is
hard to generalise and I suspect this kind of issue will come up
again.

Also, cases like this reveal the fact that, as I understand it, Dask
does not schedule work in a way that is sensitive to memory
requirements or availability on workers. This is quite different from
the situation that those of us who previously worked with HPC
scheduling systems such as SGE or SLURM may be used to, where you have
to give each task a memory limit, tasks are scheduled to workers only
if sufficient memory is available, and killed if they exceed their
limit. This can be a pain when you get the memory limit wrong, but it
can help to avoid clusters grinding to a halt because memory is
overcommitted.

Some of these issues have led the Pangeo folks to [develop a rechunker
package](https://medium.com/pangeo/rechunker-the-missing-link-for-chunked-array-analytics-5b2359e9dc11),
which is specifically designed to deal with the rechunking use case
but in a way that strictly deals with memory availability. I.e.,
memory is carefully managed and intermediate results are written to
some kind of persistent storage when they cannot be held in RAM. This
resolves the use case where you want to perform a one-off rechunking,
but doesn't deal with other scenarios where you need to rechunk in the
middle of a workflow, or other cases where task ordering can lead to
excessive memory use.

This got me thinking that perhaps a useful modification to the dask
distributed architecture would be to allow for storage of intermediate
results in systems other than worker memory.

For example, if your cluster is running in a public cloud, what if
workers could write intermediate results to object storage. This would
of course potentially slow down a computation a lot, because writing
to and reading from object storage would add overhead. This waste
would be particularly acute where dependent tasks can be run on the
same worker, i.e., where no data transfer between workers was
required. However, in cases where otherwise the cluster would start
moving intermediate results between workers, this might not be too
bad, because the overhead of going to/from object storage might be of
the same order as the overhead of network communication between
workers.

If this cluster storage layer had an abstract interface that allowed
different implementations, then one could imagine other
approaches. E.g., a cluster could use memcached, or a distributed file
system, or any other system that could provide shared storage to a
cluster. The current approach of using worker RAM for storage of
intermediate results could then become an instance of this interface,
that could be changed for a different implementation if desired.

One of the original reasons as I understand it for using worker memory
like this is it allows peer-to-peer communication between
workers. I.e., when data needs to be shared between workers, there is
no single bottleneck through which data must flow, rather workers can
request data directly from each other. Avoiding a bottleneck is
obviously desirable, but there are plenty of approaches to providing
shared storage to a cluster that would also achieve this. E.g., using
cloud object storage, or memcached, or a distributed file system, all
would avoid a single bottleneck, and could possibly take advantage of
the network efficiently.

If different storage implementations could be used, then another
benefit might be that a call to `persist()` could allow convenient
persistence of key results that lived between cluster
incarnations. I.e., calling `persist()` currently causes a result to
be materialised in cluster workers' memory, but if the cluster is
shutdown then that disappears. If instead one separated the cluster
workers from the shared cluster storage, and allowed for shared
storage that was backed by other systems, then the shared storage
could persist after cluster shutdown, and be reused the next time a
cluster was started.

I.e., one might imagine an API something like this:

```python
# set up some kind of shared storage for the cluster, e.g.:
from dask.distributed.memory import FSMemory
memory = FSMemory("gs://my-cluster-bucket/")
 
# launch some kind of cluster, passing in reference to memory, e.g.:
from dask_kubernetes import KubeCluster
cluster = KubeCluster(..., memory=memory)
cluster.scale(10)

# connect a client to the cluster
from dask.distributed import Client
client = Client(cluster)

# set up some big computation
x = ...

# persist to the cluster's shared storage
x.persist()

# do some more work...
```

The crucial point here is that, first time `x.persist()` is called it
runs the full computation and persists the result in the cluster's
shared storage. But if cluster is restarted and you come back later,
next time it's run it realises that the result already exists and
won't recompute it. Of course you could call `x.store()` to manually
store the result somewhere, which is what we currently do, but that
can be inconvenient because your code ends up scattered with lots of
if/else statements of the form "if the result has been stored
previously, load it from there, otherwise run a computation".

Another possible advantage to abstracting out the cluster's shared
storage system and allowing for more persistent storage, is that it
could open the door to supporting resumable computations. I.e.,
currently if a computation crashes or is interrupted part-way through,
you have to restart it from the beginning. Mostly dask is targeting
interactive workloads where each computation is relatively short, and
the distributed scheduler is resilient to things like individual
workers failing, so supporting this kind of thing has not been a
priority. But I can imagine that some use cases involving
longer-running computation might benefit from the ability to resume a
computation from the point where it was interrupted, even if the
entire cluster is shut down and restarted.

Coming back to the original discussion of debugging situations where a
cluster runs out of memory, separating the cluster workers from the
cluster storage, at least architecturally if not physically, would
seem to be useful, because it would allow separate reporting and
diagnostics of how much memory (RAM) was being used by tasks during
task execution, versus how much memory (storage) was being used to
hold onto intermediate results during computation of a task graph.

Of course I have no idea how easy this would be to implement, but
thought I'd float the idea in case it was useful for discussion.
