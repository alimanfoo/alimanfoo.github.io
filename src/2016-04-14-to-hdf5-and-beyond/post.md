---
layout: post
title: To HDF5 and beyond
---


*This post contains some notes about three Python libraries for working with numerical data too large to fit into main memory: [h5py](http://www.h5py.org/), [Bcolz](http://bcolz.blosc.org/) and [Zarr](https://github.com/alimanfoo/zarr).*

*2016-05-18: Updated to use the new [1.0.0 release of Zarr](http://zarr.readthedocs.io/en/latest/index.html).*

## HDF5 (``h5py``)

When I first discovered the [HDF5 file format](https://www.hdfgroup.org/HDF5/) a few years ago it was pretty transformative. I'd been struggling to find efficient ways of exploring and analysing data coming from [our research](http://www.malariagen.net/ag1000g) on genetic variation in the mosquitoes that carry malaria. The data are not enormous - typically integer arrays with around 20 billion elements - but they are too large to work with in memory on a typical laptop or desktop machine. The file formats traditionally used for these data are very slow to parse, so I went looking for an alternative.

HDF5 files provide a great solution for storing multi-dimensional arrays of numerical data. The arrays are divided up into chunks and each chunk is compressed, enabling data to be stored efficiently in memory or on disk. Depending on how the chunks are configured, usually a good compromise can be achieved that means data can be read very quickly even when using different access patterns, e.g., taking horizontal or vertical slices of a matrix. Also, the [``h5py``](http://www.h5py.org/) Python library provides a very convenient API for working with HDF5 files. I found there was a whole range of analyses I could happily get done on my laptop on the train home from work.

## Bcolz

A bit later on I discovered the [Bcolz](http://bcolz.blosc.org/) library. Bcolz is primarily intended for storing and querying large tables of data, but it does provide a [``carray``](http://bcolz.blosc.org/reference.html#the-carray-class) class that is roughly analogous to an HDF5 dataset in that it can store numerical arrays in a chunked, compressed form, either in memory or on disk. 

Reading and writing data to a ``bcolz.carray`` is typically a lot faster than HDF5. For example:


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

    1 loop, best of 3: 1.39 s per loop


Time how long it takes to store in a ``bcolz.carray``.


{% highlight python %}
%timeit bcolz.carray(a1, chunklen=2**18, cparams=bcolz.cparams(cname='lz4', clevel=5, shuffle=1))
{% endhighlight %}

    10 loops, best of 3: 81.5 ms per loop


In the example above, Bcolz is more than 10 times faster at storing (compressing) the data than HDF5. As I understand it, this performance gain comes from several factors. Bcolz uses a C library called [Blosc](https://github.com/blosc/c-blosc) to perform compression and decompression operations. Blosc can use multiple threads internally, so some of the work is done in parallel. Blosc also splits data up in a way that is designed to work well with the CPU cache architecture. Finally, Blosc is a meta-compressor and several different compression libraries can be used internally - above I used the LZ4 compressor, which does not achieve quite the same compression ratios as gzip (zlib) but is much faster with numerical data.

## Zarr

Speed really makes a difference when working interactively with data, so I started using the ``bcolz.carray`` class where possible in my analyses, especially for storing intermediate data. However, it does have some limitations. A ``bcolz.carray`` can be multidimensional, but because Bcolz is not really designed for multi-dimensional data, a ``bcolz.carray`` can only be chunked along the first dimension. This means taking slices of the first dimension is efficient, but slicing any other dimension will be very inefficient, because the entire array will need to be read and decompressed to access even a single column of a matrix.

To explore better ways of working with large multi-dimensional data, I recently created a new library called [Zarr](https://github.com/alimanfoo/zarr). Zarr like Bcolz uses Blosc internally to handle all compression and decompression operations. However, Zarr supports chunking of arrays along multiple dimensions, enabling good performance for multiple data access patterns. For example:


{% highlight python %}
import zarr
zarr.__version__
{% endhighlight %}




    '1.0.0'



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




    zarr.core.Array((10000, 10000), int32, chunks=(1000, 1000), order=C)
      compression: blosc; compression_opts: {'cname': 'blosclz', 'clevel': 5, 'shuffle': 1}
      nbytes: 381.5M; nbytes_stored: 10.0M; ratio: 38.0; initialized: 100/100
      store: builtins.dict



Time how long it takes to access a slice along the first dimension.


{% highlight python %}
%timeit c2[:1000]
{% endhighlight %}

    100 loops, best of 3: 11.5 ms per loop



{% highlight python %}
%timeit z[:1000]
{% endhighlight %}

    10 loops, best of 3: 19 ms per loop


Time how long it takes to access a slice along the second dimension.


{% highlight python %}
%timeit c2[:, :1000]
{% endhighlight %}

    1 loop, best of 3: 148 ms per loop



{% highlight python %}
%timeit z[:, :1000]
{% endhighlight %}

    100 loops, best of 3: 12.9 ms per loop


By using Zarr and chunking along both dimensions of the array, we have forfeited a small amount of speed when slicing the first dimension to gain a lot of speed when accessing the second dimension.

Like h5py and Bcolz, Zarr can store data either in memory or on disk. Zarr has some other notable features too. For example, multi-dimensional arrays can be resized along any dimension, allowing an array to be grown by appending new data in a flexible way. Also, [Zarr arrays can be used in parallel computations](http://alimanfoo.github.io/2016/05/16/cpu-blues.html), supporting concurrent reads and writes in either a multi-threaded or multi-process context. 

Zarr is still in an experimental phase, but if you do try it out, any feedback is very welcome.

## Further reading

* [HDF5](https://www.hdfgroup.org/HDF5/)
* [h5py](http://www.h5py.org/)
* [Bcolz](http://bcolz.blosc.org/)
* [Blosc](http://blosc.org/)
* [Zarr](https://github.com/alimanfoo/zarr)

## Post-script: Performance with real genotype data

Here are benchmarks with some real data.


{% highlight python %}
import operator
from functools import reduce


def human_readable_size(size):
    if size < 2**10:
        return "%s" % size
    elif size < 2**20:
        return "%.1fK" % (size / float(2**10))
    elif size < 2**30:
        return "%.1fM" % (size / float(2**20))
    elif size < 2**40:
        return "%.1fG" % (size / float(2**30))
    else:
        return "%.1fT" % (size / float(2**40))

    
def h5d_diagnostics(d):
    """Print some diagnostics on an HDF5 dataset."""
    
    print(d)
    nbytes = reduce(operator.mul, d.shape) * d.dtype.itemsize
    cbytes = d._id.get_storage_size()
    if cbytes > 0:
        ratio = nbytes / cbytes
    else:
        ratio = np.inf
    r = '  cname=%s' % d.compression
    r += ', clevel=%s' % d.compression_opts
    r += ', shuffle=%s' % d.shuffle
    r += '\n  nbytes=%s' % human_readable_size(nbytes)
    r += ', cbytes=%s' % human_readable_size(cbytes)
    r += ', ratio=%.1f' % ratio
    r += ', chunks=%s' % str(d.chunks)
    print(r)
    
{% endhighlight %}

Locate a genotype dataset within an HDF5 file from the [Ag1000G project](http://www.malariagen.net/data/ag1000g-phase1-AR3).


{% highlight python %}
callset = h5py.File('/data/coluzzi/ag1000g/data/phase1/release/AR3/variation/main/hdf5/ag1000g.phase1.ar3.pass.h5',
                    mode='r')
genotype = callset['3R/calldata/genotype']
genotype
{% endhighlight %}




    <HDF5 dataset "genotype": shape (13167162, 765, 2), type "|i1">



Extract the first million rows of the dataset to use for benchmarking.


{% highlight python %}
a3 = genotype[:1000000]
{% endhighlight %}

Benchmark compression performance.


{% highlight python %}
%%time
genotype_hdf5_gzip = h5fmem().create_dataset('genotype', data=a3, 
                                             compression='gzip', compression_opts=1, shuffle=False,
                                             chunks=(10000, 100, 2))
{% endhighlight %}

    CPU times: user 6.58 s, sys: 72 ms, total: 6.66 s
    Wall time: 6.62 s



{% highlight python %}
h5d_diagnostics(genotype_hdf5_gzip)
{% endhighlight %}

    <HDF5 dataset "genotype": shape (1000000, 765, 2), type "|i1">
      cname=gzip, clevel=1, shuffle=False
      nbytes=1.4G, cbytes=51.1M, ratio=28.5, chunks=(10000, 100, 2)



{% highlight python %}
%%time
genotype_hdf5_lzf = h5fmem().create_dataset('genotype', data=a3, 
                                             compression='lzf', shuffle=False,
                                             chunks=(10000, 100, 2))
{% endhighlight %}

    CPU times: user 2.01 s, sys: 100 ms, total: 2.11 s
    Wall time: 2.1 s



{% highlight python %}
h5d_diagnostics(genotype_hdf5_lzf)
{% endhighlight %}

    <HDF5 dataset "genotype": shape (1000000, 765, 2), type "|i1">
      cname=lzf, clevel=None, shuffle=False
      nbytes=1.4G, cbytes=71.7M, ratio=20.3, chunks=(10000, 100, 2)



{% highlight python %}
%%time
genotype_carray = bcolz.carray(a3, cparams=bcolz.cparams(cname='lz4', clevel=1, shuffle=2))
{% endhighlight %}

    CPU times: user 2.29 s, sys: 0 ns, total: 2.29 s
    Wall time: 669 ms



{% highlight python %}
genotype_carray
{% endhighlight %}




    carray((1000000, 765, 2), int8)
      nbytes: 1.42 GB; cbytes: 48.70 MB; ratio: 29.96
      cparams := cparams(clevel=1, shuffle=2, cname='lz4')
    [[[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]
    
     [[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]
    
     [[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]
    
     ..., 
     [[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]
    
     [[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]
    
     [[0 0]
      [0 0]
      [0 0]
      ..., 
      [0 0]
      [0 0]
      [0 0]]]




{% highlight python %}
%%time
genotype_zarr = zarr.array(a3, chunks=(10000, 100, 2), compression='blosc', 
                           compression_opts=dict(cname='lz4', clevel=1, shuffle=2))
{% endhighlight %}

    CPU times: user 2.85 s, sys: 56 ms, total: 2.9 s
    Wall time: 1.26 s



{% highlight python %}
genotype_zarr
{% endhighlight %}




    zarr.core.Array((1000000, 765, 2), int8, chunks=(10000, 100, 2), order=C)
      compression: blosc; compression_opts: {'cname': 'lz4', 'clevel': 1, 'shuffle': 2}
      nbytes: 1.4G; nbytes_stored: 50.2M; ratio: 29.1; initialized: 800/800
      store: builtins.dict



Note that although I've used the LZ4 compression library with Bcolz and Zarr, the compression ratio is actually better than when using gzip (zlib) with HDF5. This is due to the bitshuffle filter, which comes bundled with Blosc. The bitshuffle filter can also be used with HDF5 with some configuration I believe.

Compression with Zarr is slightly slower than Bcolz above, but this is entirely due to the choice of chunk shape and the correlation structure in the data. If we use the same chunking for both, compression speed is similar...


{% highlight python %}
%%time 
_ = zarr.array(a3, chunks=(genotype_carray.chunklen, 765, 2), compression='blosc',
               compression_opts=dict(cname='lz4', clevel=1, shuffle=2))
{% endhighlight %}

    CPU times: user 2.21 s, sys: 56 ms, total: 2.26 s
    Wall time: 675 ms


Benchmark data access via slices along first and second dimensions. 


{% highlight python %}
%timeit genotype_hdf5_gzip[:10000]
{% endhighlight %}

    10 loops, best of 3: 24.4 ms per loop



{% highlight python %}
%timeit genotype_hdf5_lzf[:10000]
{% endhighlight %}

    10 loops, best of 3: 27.8 ms per loop



{% highlight python %}
%timeit genotype_carray[:10000]
{% endhighlight %}

    100 loops, best of 3: 9.32 ms per loop



{% highlight python %}
%timeit genotype_zarr[:10000]
{% endhighlight %}

    100 loops, best of 3: 13 ms per loop



{% highlight python %}
%timeit genotype_hdf5_gzip[:, :10]
{% endhighlight %}

    1 loop, best of 3: 279 ms per loop



{% highlight python %}
%timeit genotype_hdf5_lzf[:, :10]
{% endhighlight %}

    1 loop, best of 3: 290 ms per loop



{% highlight python %}
%timeit genotype_carray[:, :10]
{% endhighlight %}

    1 loop, best of 3: 1.09 s per loop



{% highlight python %}
%timeit genotype_zarr[:, :10]
{% endhighlight %}

    10 loops, best of 3: 123 ms per loop


## Setup


{% highlight python %}
# run on my laptop, which doesn't have AVX2
import cpuinfo
cpuinfo.main()
{% endhighlight %}

    Vendor ID: GenuineIntel
    Hardware Raw: 
    Brand: Intel(R) Core(TM) i7-3667U CPU @ 2.00GHz
    Hz Advertised: 2.0000 GHz
    Hz Actual: 3.1355 GHz
    Hz Advertised Raw: (2000000000, 0)
    Hz Actual Raw: (3135546000, 0)
    Arch: X86_64
    Bits: 64
    Count: 4
    Raw Arch String: x86_64
    L2 Cache Size: 4096 KB
    L2 Cache Line Size: 0
    L2 Cache Associativity: 0
    Stepping: 9
    Model: 58
    Family: 6
    Processor Type: 0
    Extended Model: 0
    Extended Family: 0
    Flags: acpi, aes, aperfmperf, apic, arat, arch_perfmon, avx, bts, clflush, cmov, constant_tsc, cx16, cx8, de, ds_cpl, dtes64, dtherm, dts, eagerfpu, epb, ept, erms, est, f16c, flexpriority, fpu, fsgsbase, fxsr, ht, ida, lahf_lm, lm, mca, mce, mmx, monitor, msr, mtrr, nonstop_tsc, nopl, nx, pae, pat, pbe, pcid, pclmulqdq, pdcm, pebs, pge, pln, pni, popcnt, pse, pse36, pts, rdrand, rdtscp, rep_good, sep, smep, smx, ss, sse, sse2, sse4_1, sse4_2, ssse3, syscall, tm, tm2, tpr_shadow, tsc, tsc_deadline_timer, vme, vmx, vnmi, vpid, x2apic, xsave, xsaveopt, xtopology, xtpr

