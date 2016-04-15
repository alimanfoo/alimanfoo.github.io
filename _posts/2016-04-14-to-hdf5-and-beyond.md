---
layout: post
title: To HDF5 and beyond
---


This post contains some notes about 3 Python libraries for working with numerical data too large to fit into main memory: [``h5py``](http://www.h5py.org/), [``bcolz``](http://bcolz.blosc.org/) and [``zarr``](https://github.com/alimanfoo/zarr).

## HDF5 (``h5py``)

When I first discovered the [HDF5 file format](https://www.hdfgroup.org/HDF5/) a few years ago it was pretty transformative. I'd been struggling to find efficient ways of exploring and analysing data coming from [our research](http://www.malariagen.net/ag1000g) on genetic variation in the mosquitoes that carry malaria. The data are not enormous - typically integer arrays with around 20 billion elements - but they are too large to work with in memory on a typical laptop or desktop machine. The file formats traditionally used for these data are very slow to parse, so I went looking for an alternative.

HDF5 files provide a great solution for storing multi-dimensional arrays of numerical data. The arrays are divided up into chunks and each chunk is compressed, enabling data to be stored efficiently in memory or on disk. Depending on how the chunks are configured, usually a good compromise can be achieved that means data can be read very quickly even when using different access patterns, e.g., taking horizontal or vertical slices of a matrix. Also, the [``h5py``](http://www.h5py.org/) Python library provides a very convenient API for working with HDF5 files. I found there was a whole range of analyses I could happily get done on my laptop on the train home from work.

## ``bcolz.carray``

A bit later on I discovered the [``bcolz``](http://bcolz.blosc.org/) library. ``bcolz`` is primarily intended for storing and querying large tables of data, but it does provide a [``carray``](http://bcolz.blosc.org/reference.html#the-carray-class) class that is roughly analogous to an HDF5 dataset in that it can store numerical arrays in a chunked, compressed form, either in memory or on disk. 

Reading and writing data to a ``carray`` is typically a lot faster than HDF5. For example:


{% highlight python %}
import numpy as np
import h5py
import bcolz
import tempfile


def h5fmem(**kwargs):
    """Convenience function to create an in-memory HDF5 file."""

    # need a file name even tho nothing is ever written
    fn = tempfile.mktemp()

    # file creation args
    kwargs['mode'] = 'w'
    kwargs['driver'] = 'core'
    kwargs['backing_store'] = False

    # open HDF5 file
    h5f = h5py.File(fn, **kwargs)

    return h5f
{% endhighlight %}

Setup a simple array of integer data to store.


{% highlight python %}
a1 = np.arange(1e8, dtype='i4')
a1
{% endhighlight %}




    array([       0,        1,        2, ..., 99999997, 99999998, 99999999], dtype=int32)



Time how long it takes to store in an HDF5 dataset.


{% highlight python %}
%timeit h5fmem().create_dataset('arange', data=a1, chunks=(2**18,), compression='gzip', compression_opts=1, shuffle=True)
{% endhighlight %}

    1 loop, best of 3: 1.41 s per loop


Time how long it takes to store in a ``carray``.


{% highlight python %}
%timeit bcolz.carray(a1, chunklen=2**18, cparams=bcolz.cparams(cname='lz4', clevel=5, shuffle=1))
{% endhighlight %}

    10 loops, best of 3: 95.6 ms per loop


In the example above, ``bcolz`` is more than 10 times faster at storing (compressing) the data than HDF5. As I understand it, this performance gain comes from several factors. ``bcolz`` uses a C library called [``blosc``](https://github.com/blosc/c-blosc) internally to perform compression and decompression operations. ``blosc`` can use multiple threads, so some of the work is done in parallel. ``blosc`` also splits data up in a way that is designed to work well with the CPU cache architecture. Finally, ``blosc`` is a meta-compressor and several different compression libraries can be used - above I used the ``lz4`` compressor, which does not achieve quite the same compression ratios as ``gzip`` (``zlib``) but is much faster with numerical data.

## ``zarr``

Speed really makes a difference when working interactively with data, so I started using the ``carray`` class where possible in my analyses, especially for storing intermediate data. However, it does have some limitations. A ``carray`` can be multidimensional, but because ``bcolz`` is not really designed for multi-dimensional data, a ``carray`` can only be chunked along the first dimension. This means taking slices of the first dimension is efficient, but slicing any other dimension will be very inefficient, because the entire array will need to be read and decompressed to access even a single column of a matrix.

To explore better ways of working with large multi-dimensional data, I recently created a new library called [``zarr``](https://github.com/alimanfoo/zarr). ``zarr`` borrows code heavily from ``bcolz``, and in particular it also uses ``blosc`` internally to handle all compression and decompression operations. However, ``zarr`` supports chunking of arrays along multiple dimensions, enabling good performance for multiple data access patterns. For example:


{% highlight python %}
import zarr
{% endhighlight %}

Setup a 2-dimensional array of integer data.


{% highlight python %}
a2 = np.arange(1e8, dtype='i4').reshape(10000, 10000)
{% endhighlight %}

Store the data in a ``carray``.


{% highlight python %}
c2 = bcolz.carray(a2, chunklen=100)
c2
{% endhighlight %}




    carray((10000, 10000), int32)
      nbytes: 381.47 MB; cbytes: 10.63 MB; ratio: 35.87
      cparams := cparams(clevel=5, shuffle=1, cname='blosclz')
    [[       0        1        2 ...,     9997     9998     9999]
     [   10000    10001    10002 ...,    19997    19998    19999]
     [   20000    20001    20002 ...,    29997    29998    29999]
     ..., 
     [99970000 99970001 99970002 ..., 99979997 99979998 99979999]
     [99980000 99980001 99980002 ..., 99989997 99989998 99989999]
     [99990000 99990001 99990002 ..., 99999997 99999998 99999999]]



Store the data in a ``zarr`` array.


{% highlight python %}
z = zarr.array(a2, chunks=(1000, 1000))
z
{% endhighlight %}




    zarr.ext.SynchronizedArray((10000, 10000), int32, chunks=(1000, 1000))
      cname: blosclz; clevel: 5; shuffle: 1 (BYTESHUFFLE)
      nbytes: 381.5M; cbytes: 10.0M; ratio: 38.0; initialized: 100/100



Time how long it takes to access a slice along the first dimension.


{% highlight python %}
%timeit c2[:1000]
{% endhighlight %}

    10 loops, best of 3: 22.1 ms per loop



{% highlight python %}
%timeit z[:1000]
{% endhighlight %}

    10 loops, best of 3: 39.3 ms per loop


Time how long it takes to access a slice along the second dimension.


{% highlight python %}
%timeit c2[:, :1000]
{% endhighlight %}

    1 loop, best of 3: 250 ms per loop



{% highlight python %}
%timeit z[:, :1000]
{% endhighlight %}

    10 loops, best of 3: 22.5 ms per loop


By using ``zarr`` and chunking along both dimensions of the array, we have forfeited a small amount of speed when slicing the first dimension to gain a lot of speed when accessing the second dimension.

Like ``h5py`` and ``bcolz``, ``zarr`` can store data either in memory or on disk. However, ``zarr`` has some unique features too. For example, multi-dimensional arrays can be resized along any dimension, allowing an array to be grown by appending new data in a flexible way. Also, ``zarr`` arrays can be used in parallel computations, supporting concurrent reads and writes in either a multi-threaded or multi-process context. That is something I am just beginning to explore, and hope to follow up in a separate post.

``zarr`` is still in an experimental phase, but if you do try it out, any feedback is very welcome.

## Further reading

* [HDF5](https://www.hdfgroup.org/HDF5/)
* [h5py](http://www.h5py.org/)
* [bcolz](http://bcolz.blosc.org/)
* [blosc](http://blosc.org/)
* [zarr](https://github.com/alimanfoo/zarr)
