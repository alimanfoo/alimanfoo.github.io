---
layout: post
title: Zarr 2 - groups and filters
---


Recently I've been working on [Zarr](http://zarr.readthedocs.io/en/latest/), a Python package providing chunked, compressed storage for numerical arrays. I've just released [Zarr version 2](http://zarr.readthedocs.io/en/latest/release.html) which adds two new major features: [groups](http://zarr.readthedocs.io/en/latest/tutorial.html#groups) and [filters](http://zarr.readthedocs.io/en/latest/tutorial.html#filters). This post provides a brief tour of what's new.


{% highlight python %}
import numpy as np
import zarr
zarr.__version__
{% endhighlight %}




    '2.1.0'



## Groups

Zarr is heavily inspired by [HDF5](https://www.hdfgroup.org/HDF5/), which provides the capability to organize arrays into groups. Groups can also contain other groups, forming a hierarchy.

Here's how to create a group with Zarr:


{% highlight python %}
root_group = zarr.group()
root_group
{% endhighlight %}




    Group(/, 0)
      store: DictStore



If you are familiar with [h5py](http://www.h5py.org/), the API for the Zarr [Group](http://zarr.readthedocs.io/en/latest/api/hierarchy.html#zarr.hierarchy.Group) class is very similar. For example, create a sub-group:


{% highlight python %}
foo_group = root_group.create_group('foo')
foo_group
{% endhighlight %}




    Group(/foo, 0)
      store: DictStore



The [Group](http://zarr.readthedocs.io/en/latest/api/hierarchy.html#zarr.hierarchy.Group) class has various methods for creating arrays within a group, e.g.:


{% highlight python %}
a1 = foo_group.zeros('bar', shape=(10000, 10000), chunks=(1000, 1000))
a1
{% endhighlight %}




    Array(/foo/bar, (10000, 10000), int32, chunks=(1000, 1000), order=C)
      nbytes: 381.5M; nbytes_stored: 323; ratio: 1238390.1; initialized: 0/100
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: DictStore



Group members can be accessed via the square bracket notation, e.g:


{% highlight python %}
root_group['foo']
{% endhighlight %}




    Group(/foo, 1)
      arrays: 1; bar
      store: DictStore




{% highlight python %}
foo_group['bar']
{% endhighlight %}




    Array(/foo/bar, (10000, 10000), int32, chunks=(1000, 1000), order=C)
      nbytes: 381.5M; nbytes_stored: 323; ratio: 1238390.1; initialized: 0/100
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: DictStore



Multiple hierarchy levels can also be traversed, e.g.:


{% highlight python %}
root_group['foo/bar']
{% endhighlight %}




    Array(/foo/bar, (10000, 10000), int32, chunks=(1000, 1000), order=C)
      nbytes: 381.5M; nbytes_stored: 323; ratio: 1238390.1; initialized: 0/100
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: DictStore



In the examples above, all data will be stored in memory. However, Zarr can use a variety of other storage layers. For example, data can also be stored on the file system, e.g.:


{% highlight python %}
store = zarr.DirectoryStore('example')
persistent_group = zarr.group(store)
a2 = persistent_group.require_dataset('foo/bar', shape=(10000, 10000), chunks=(1000, 1000))
a2
{% endhighlight %}




    Array(/foo/bar, (10000, 10000), float64, chunks=(1000, 1000), order=C)
      nbytes: 762.9M; nbytes_stored: 323; ratio: 2476780.2; initialized: 0/100
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: DirectoryStore



Data can also be stored in a [Zip file](http://zarr.readthedocs.io/en/latest/api/storage.html#zarr.storage.ZipStore) (with some limitations), on [S3](http://s3fs.readthedocs.io/en/latest/api.html#s3fs.mapping.S3Map), on [HDFS](http://hdfs3.readthedocs.io/en/latest/api.html#hdfs3.mapping.HDFSMap), or via any storage system that can expose a [MutableMapping](https://docs.python.org/3/library/collections.abc.html) interface.

For more information about groups, see the [groups section of the Zarr tutorial](http://zarr.readthedocs.io/en/latest/tutorial.html#groups) and the [zarr.hierarchy API docs](http://zarr.readthedocs.io/en/latest/api/hierarchy.html).

## Filters


{% highlight python %}

{% endhighlight %}
