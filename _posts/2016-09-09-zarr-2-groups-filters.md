---
layout: post
title: Zarr 2 - groups, filters and Zstandard
---


Recently I've been working on [Zarr](http://zarr.readthedocs.io/en/latest/), a Python package providing chunked, compressed storage for numerical arrays. I've just released [Zarr version 2](http://zarr.readthedocs.io/en/latest/release.html) which adds two new major features: [groups](http://zarr.readthedocs.io/en/latest/tutorial.html#groups) and [filters](http://zarr.readthedocs.io/en/latest/tutorial.html#filters). It also brings support for [Zstandard](http://facebook.github.io/zstd/) compression via [Blosc](http://www.blosc.org/). This post provides a brief tour of what's new.


{% highlight python %}
import numpy as np
import zarr; print('zarr', zarr.__version__)
from zarr import blosc; print('blosc', blosc.__version__)
{% endhighlight %}

    zarr 2.1.1
    blosc 1.11.1


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
a1 = foo_group.zeros('bar', shape=(10000, 10000))
a1
{% endhighlight %}




    Array(/foo/bar, (10000, 10000), float64, chunks=(157, 313), order=C)
      nbytes: 762.9M; nbytes_stored: 321; ratio: 2492211.8; initialized: 0/2048
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




    Array(/foo/bar, (10000, 10000), float64, chunks=(157, 313), order=C)
      nbytes: 762.9M; nbytes_stored: 321; ratio: 2492211.8; initialized: 0/2048
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

Zarr 2 also adds support for filters. Filters are arbitrary data transformations that can be applied to encode data chunks prior to storage, and to decode data on retrieval. The idea is that, for some kinds of data, certain transformations can help to improve the compression ratio, or provide other useful features such as error checking.

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

For floating point data you can try the quantize or scale-offset filters, which allow you to store data with a given precision. E.g.:


{% highlight python %}
data = np.random.normal(loc=10, scale=2, size=10000000)
data
{% endhighlight %}




    array([  6.72635296,   9.91028746,   7.83297722, ...,  10.12890375,
            10.54766939,   9.78125015])




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




    array([  6.75  ,   9.9375,   7.8125, ...,  10.125 ,  10.5625,   9.8125])




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




    array([  6.7,   9.9,   7.8, ...,  10.1,  10.5,   9.8])



More than one filter can be provided, and compressors are filters too, so you can do some fairly zany things if you want to, e.g.:


{% highlight python %}
zany = zarr.array(data, 
                  filters=[zarr.Quantize(digits=1, dtype=data.dtype),
                           zarr.Blosc(clevel=0, shuffle=blosc.SHUFFLE)],
                  compressor=zarr.BZ2(level=9))
zany
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 9.1M; ratio: 8.4; initialized: 256/256
      filters: Quantize(digits=1, dtype=float64)
               Blosc(cname='lz4', clevel=0, shuffle=1)
      compressor: BZ2(level=9)
      store: dict



Please note that the built-in filters in Zarr have not been optimized at all yet, I am sure there is much room for performance improvement. The main idea in the Zarr version 2 release is to establish a simple API for developing and integrating new filters, so it is easier to explore different options for new data.

For more information about filters, see the [filters section of the Zarr tutorial](http://zarr.readthedocs.io/en/latest/tutorial.html#filters) and the [`zarr.codecs` API docs](http://zarr.readthedocs.io/en/latest/api/codecs.html).

## Zstandard

One other thing I wanted to mention is that Zarr now supports compression with [Zstandard](http://facebook.github.io/zstd/). I'm not going to say much here because I'm hoping to write up some detailed benchmark data in a separate blog post soon. But from what I have seen so far, Zstandard is a superb codec, providing high compression ratios while maintaining excellent speed for both compression and decompression, although as always mileage will vary depending on the nature of the data.

Here's how to create a Zarr array using Zstandard compression via Blosc:


{% highlight python %}
z4 = zarr.array(data, compressor=zarr.Blosc(cname='zstd', clevel=5, shuffle=blosc.SHUFFLE))
z4
{% endhighlight %}




    Array((10000000,), float64, chunks=(39063,), order=C)
      nbytes: 76.3M; nbytes_stored: 66.2M; ratio: 1.2; initialized: 256/256
      compressor: Blosc(cname='zstd', clevel=5, shuffle=1)
      store: dict



## Acknowledgments and further reading

I hope this new release of Zarr is useful, any [feedback or suggestions](https://github.com/alimanfoo/zarr/issues) are very welcome as always.

The latest version of Zarr is available from [PyPI](https://pypi.python.org/pypi/zarr) and [conda-forge](https://github.com/conda-forge/zarr-feedstock), see the [installation instructions](http://zarr.readthedocs.io/en/latest/index.html#installation) for more information.  

If this is the first time you are reading about Zarr, you might like to take a look at these previous posts: [To HDF5 and beyond](http://alimanfoo.github.io/2016/04/14/to-hdf5-and-beyond.html), and [CPU blues](http://alimanfoo.github.io/2016/05/16/cpu-blues.html).

Development of Zarr is motivated by our work on the [genomic epidemiology of malaria](https://www.malariagen.net/) and supported by the [MRC Centre for Genomics and Global Health](http://www.cggh.org/).

Thanks to [Matthew Rocklin](https://github.com/mrocklin), [Stephan Hoyer](https://github.com/shoyer) and [Francesc Alted](https://github.com/FrancescAlted) for much good advice and inspiration. As Francesc would say, enjoy data!
