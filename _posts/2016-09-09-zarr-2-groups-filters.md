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




    Array(/foo/bar, (10000, 10000), float64, chunks=(1000, 1000), order=C)
      nbytes: 762.9M; nbytes_stored: 323; ratio: 2476780.2; initialized: 0/100
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




    Array(/foo/bar, (10000, 10000), float64, chunks=(1000, 1000), order=C)
      nbytes: 762.9M; nbytes_stored: 323; ratio: 2476780.2; initialized: 0/100
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: DictStore



Multiple hierarchy levels can also be traversed, e.g.:


{% highlight python %}
root_group['foo/bar'] == foo_group['bar']
{% endhighlight %}




    True



In the examples above, all data will be stored in memory. However, Zarr can use a variety of other storage layers. For example, data can be stored on the file system, e.g.:


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



Data can also be stored in a [Zip file](http://zarr.readthedocs.io/en/latest/api/storage.html#zarr.storage.ZipStore) (with some limitations), on [S3](http://s3fs.readthedocs.io/en/latest/api.html#s3fs.mapping.S3Map), [HDFS](http://hdfs3.readthedocs.io/en/latest/api.html#hdfs3.mapping.HDFSMap), or any storage system that can by accessed via the [MutableMapping](https://docs.python.org/3/library/collections.abc.html) interface.

For more information about groups, see the [groups section of the Zarr tutorial](http://zarr.readthedocs.io/en/latest/tutorial.html#groups) and the [`zarr.hierarchy` API docs](http://zarr.readthedocs.io/en/latest/api/hierarchy.html).

## Filters

Zarr 2 also adds support for filters. Filters are arbitrary data transformations that can be applied to data chunks prior to compression and storage. The idea is that, for some kinds of data, certain transformations can help to improve the compression ratio, or provide other useful features such as error checking.

There are a few [built-in filter classes](http://zarr.readthedocs.io/en/latest/api/codecs.html) available in Zarr, including delta, scale-offset, quantize, packbits and categorize transformations. Here's a trivial example using the delta filter:


{% highlight python %}
data = np.arange(100000000)
z1a = zarr.array(data)  # no filters
z1a
{% endhighlight %}




    Array((100000000,), int64, chunks=(48829,), order=C)
      nbytes: 762.9M; nbytes_stored: 11.3M; ratio: 67.4; initialized: 2048/2048
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict




{% highlight python %}
z1b = zarr.array(data, filters=[zarr.Delta(dtype=data.dtype)])
z1b
{% endhighlight %}




    Array((100000000,), int64, chunks=(48829,), order=C)
      nbytes: 762.9M; nbytes_stored: 3.6M; ratio: 213.3; initialized: 2048/2048
      filters: Delta(dtype=int64)
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict



Note that the delta filter improves the compression ratio in this case. 

For floating point data you can try the quantize or scale-offset filters. Both are lossy for floating point data and allow you to store data with a given precision. E.g.:


{% highlight python %}
data = np.random.normal(loc=10, scale=2, size=10000000)
data
{% endhighlight %}




    array([  8.74674205,  12.21517843,   5.64209159, ...,   9.15225917,
            10.78280929,  10.8855431 ])




{% highlight python %}
z2a = zarr.array(data)  # no filters
z2a
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 67.2M; ratio: 1.1; initialized: 256/256
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict




{% highlight python %}
z2b = zarr.array(data, filters=[zarr.Quantize(digits=1, dtype=data.dtype)])
z2b
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 19.4M; ratio: 3.9; initialized: 256/256
      filters: Quantize(digits=1, dtype=float64)
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict




{% highlight python %}
z2b[:]
{% endhighlight %}




    array([  8.75  ,  12.1875,   5.625 , ...,   9.125 ,  10.8125,  10.875 ])




{% highlight python %}
z2c = zarr.array(data, filters=[zarr.FixedScaleOffset(offset=10, scale=10, dtype=data.dtype)])
z2c
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 17.3M; ratio: 4.4; initialized: 256/256
      filters: FixedScaleOffset(scale=10, offset=10, dtype=float64)
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict




{% highlight python %}
z2c[:]
{% endhighlight %}




    array([  8.7,  12.2,   5.6, ...,   9.2,  10.8,  10.9])



Here's an example using the packbits filter with a Boolean array:


{% highlight python %}
data = np.random.randint(0, 2, size=10000000, dtype=bool)
data
{% endhighlight %}




    array([False,  True,  True, ...,  True,  True, False], dtype=bool)




{% highlight python %}
z3a = zarr.array(data)  # no filters
z3a
{% endhighlight %}




    Array((10000000,), bool, chunks=(156250,), order=C)
      nbytes: 9.5M; nbytes_stored: 4.8M; ratio: 2.0; initialized: 64/64
      compressor: Blosc(cname='lz4', clevel=5, shuffle=1)
      store: dict




{% highlight python %}
z3b = zarr.array(data, filters=[zarr.PackBits()], compressor=None)
z3b
{% endhighlight %}




    Array((10000000,), bool, chunks=(156250,), order=C)
      nbytes: 9.5M; nbytes_stored: 1.2M; ratio: 8.0; initialized: 64/64
      filters: PackBits()
      store: dict



The Zarr packbits filter packs boolean values into single bits, hence the compression ratio of 8.0 on some Random boolean data.

More than one filter can be provided, and compressors are filters too, so you can do some fairly zany things if you want to, e.g.:


{% highlight python %}
data = np.random.normal(loc=10, scale=2, size=10000000)
zany = zarr.array(data, 
                  filters=[zarr.Quantize(digits=1, dtype=data.dtype),
                           zarr.Blosc(clevel=0, shuffle=1),
                           zarr.BZ2(level=9)],
                  compression=None)
zany
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 9.1M; ratio: 8.4; initialized: 256/256
      filters: Quantize(digits=1, dtype=float64)
               Blosc(cname='lz4', clevel=0, shuffle=1)
               BZ2(level=9)
      store: dict



Please note that the built-in filters in Zarr have not been optimized at all yet, and I am sure there is much room for performance improvement. The main idea in this release is to establish a simple API for developing and integrating new filters, so it is easier to explore different options for new data.

For more information about filters, see the [filters section of the Zarr tutorial](TODO) and the [zarr.codecs API docs](TODO).

## Acknowledgments and further reading

I hope this new release of Zarr is useful, any [feedback or suggestions]() are very welcome as always. The latest version of Zarr is available from [PyPI]() and [conda-forge](), see the [installation instructions]() for more information. Here are some other resources you might find interesting:

* Zarr documentation
* To HDF5 and beyond
* CPU blues
* ...

Development of Zarr is motivated by our work on the [genomic epidemiology of malaria](), and supported by the [MRC Centre for Genomics and Global Health]().

Thanks to Matthew Rocklin, Stephan Hoyer and Francesc Alted for much advice and inspiration. As Francesc would say, enjoy data!


{% highlight python %}

{% endhighlight %}
