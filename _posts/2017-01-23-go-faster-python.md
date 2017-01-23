---
layout: post
title: Go faster Python
---


## Preamble

This blog post gives an introduction to some techniques for benchmarking, profiling and optimising Python code. If you would like to try the code examples for yourself, you can [download the Jupyter notebook](@@TODO) that this blog post was generated from. To run the notebook, you will need a working Python 3 installation, and will also need to install a couple of Python packages. The way I did that was to first [install miniconda]() into my home directory. I then ran the following commands from the command line:

<pre>
user@host:~$ export PATH=~/miniconda3/bin/:$PATH
user@host:~$ conda create -n go_faster_python python=3.5
user@host:~$ source activate go_faster_python
(go_faster_python) user@host:~$ conda config --add channels conda-forge
(go_faster_python) user@host:~$ conda install cython numpy jupyter line_profiler
(go_faster_python) user@host:~$ jupyter notebook &
</pre>

## Introduction

I use Python both for writing software libraries and for interactive data analysis. Python is an **interpreted**, **dynamically-typed** language, with lots of convenient **high-level data structures** like lists, sets, dicts, etc. It's great for getting things done, because no compile step means no delays when developing and testing code or when exploring data. No type declarations means lots of flexibility and less typing (on the keyboard I mean - sorry, bad joke). The high-level data structures mean you can focus on solving the problem rather than low-level nuts and bolts.

But the down-side is that Python can be slow. If you have a Python program that's running slowly, what are your options?

## Benchmarking and profiling

You've probably heard someone say that **premature optimisation is the root of all evil**. That's a pretty extreme statement, but I think it doesn't hurt to emphasise the experience that many people have had, which is that my intuitions for why a piece of code is running slowly are nearly always wrong. If you're going to try and optimise something, you need to do some benchmarking and profiling first, to find out (1) exactly how slow it goes, and (2) where the bottleneck is.

To introduce some basic Python benchmarking and profiling tools, let's look at a toy example: computing the sum of 2-dimensional array of integers. Here's some example data:


{% highlight python %}
data = [[90, 62, 33, 78, 82],
        [37, 31, 0, 72, 32],
        [7, 71, 79, 81, 100],
        [33, 50, 66, 81, 71],
        [87, 26, 54, 78, 81],
        [37, 22, 96, 79, 41],
        [88, 75, 100, 19, 88],
        [24, 72, 59, 33, 92],
        [71, 6, 59, 8, 11],
        [89, 76, 65, 12, 13]]
{% endhighlight %}

Strictly speaking this isn't an array, it's a list of lists. But using lists is a common and natural way to store data in Python. 

Here is a naive implementation of a function called `sum2d` to compute the overall sum of a 2-dimensional data structure:


{% highlight python %}
def sum1d(l):
    """Compute the sum of a list of numbers."""
    s = 0
    for x in l:
        s += x
    return s


def sum2d(ll):
    """Compute the sum of a list of lists of numbers."""
    s = 0
    for l in ll:
        s += sum1d(l)
    return s

{% endhighlight %}

Run the implementation to check it works:


{% highlight python %}
sum2d(data)
{% endhighlight %}




    2817



We need a bigger dataset to illustrate slow performance. To make a bigger test dataset I'm going to make use of the multiplication operator ('\*') which when applied to a Python list will create a new list by repeating the elements of the original list. E.g., here's the original list repeated twice: 


{% highlight python %}
data * 2
{% endhighlight %}




    [[90, 62, 33, 78, 82],
     [37, 31, 0, 72, 32],
     [7, 71, 79, 81, 100],
     [33, 50, 66, 81, 71],
     [87, 26, 54, 78, 81],
     [37, 22, 96, 79, 41],
     [88, 75, 100, 19, 88],
     [24, 72, 59, 33, 92],
     [71, 6, 59, 8, 11],
     [89, 76, 65, 12, 13],
     [90, 62, 33, 78, 82],
     [37, 31, 0, 72, 32],
     [7, 71, 79, 81, 100],
     [33, 50, 66, 81, 71],
     [87, 26, 54, 78, 81],
     [37, 22, 96, 79, 41],
     [88, 75, 100, 19, 88],
     [24, 72, 59, 33, 92],
     [71, 6, 59, 8, 11],
     [89, 76, 65, 12, 13]]



Make a bigger dataset by repeating the original data a million times:


{% highlight python %}
big_data = data * 1000000
{% endhighlight %}

Now we have a dataset that is 10,000,000 rows by 5 columns:


{% highlight python %}
len(big_data)
{% endhighlight %}




    10000000




{% highlight python %}
len(big_data[0])
{% endhighlight %}




    5



Try running the function on these data:


{% highlight python %}
sum2d(big_data)
{% endhighlight %}




    2817000000



On my laptop this takes a few seconds to run.

### Benchmarking:  `%time`, `%timeit`, `timeit`

Before you start optimising, you need a good estimate of performance as a place to start from, so you know when you've improved something. If you're working in a Jupyter notebook there are a couple of magic commands available which are very helpful for benchmarking: [`%time`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25time#magic-time) and [`%timeit`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25timeit#magic-timeit). If you're writing a Python script to do the benchmarking, you can use the [`timeit`](https://docs.python.org/3/library/timeit.html) module from the Python standard library.

Let's look at the output from [`%time`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25time#magic-time):


{% highlight python %}
%time sum2d(big_data)
{% endhighlight %}

    CPU times: user 2.85 s, sys: 0 ns, total: 2.85 s
    Wall time: 2.85 s





    2817000000



The first line of the output gives the amount of CPU time, broken down into 'user' (your code) and 'sys' (operating system code). The second line gives the wall time, which is the actual amount of time elapsed. Generally the total CPU time and the wall time will be the same, but sometimes not. E.g., if you are benchmarking a multi-threaded program, then wall time may be less than CPU time, because CPU time counts time spent by each CPU separately and adds them together, but the CPUs may actually be working in parallel.

One thing to watch out for when benchmarking is that performance can be variable, and may be affected by other processes running on your computer. To see this happen, try running the cell above again, but while it's running, give your computer something else to do at the same time, e.g., play some music, or a video, or just scroll the page up and down a bit.

To control for this variation, it's a good idea to benchmark several runs (and avoid the temptation to check your email while it's running). The [`%timeit`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25timeit#magic-timeit) magic will automatically benchmark a piece of code several times: 


{% highlight python %}
%timeit sum2d(big_data)
{% endhighlight %}

    1 loop, best of 3: 2.82 s per loop


Alternatively, using the [`timeit`](https://docs.python.org/3/library/timeit.html) module:


{% highlight python %}
import timeit
timeit.repeat(stmt='sum2d(big_data)', repeat=3, number=1, globals=globals())
{% endhighlight %}




    [2.8302519290009513, 2.9616496020007617, 2.8006342520002363]



### Function profiling: `%prun`, `cProfile`

The next thing to do is investigate which part of the code is taking the most time. If you're working in a Jupyter notebook, the [`%prun`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25prun#magic-prun) command is a very convenient way to profile some code. Use it like this:


{% highlight python %}
%prun sum2d(big_data)
{% endhighlight %}

     

The output from [`%prun`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25prun#magic-prun) pops up in a separate panel, but for this blog post I need to get the output inline, so I'm going to use the [`cProfile`](https://docs.python.org/3/library/profile.html?highlight=cprofile) module from the Python standard library directly, which does the same thing:


{% highlight python %}
import cProfile
cProfile.run('sum2d(big_data)', sort='time')
{% endhighlight %}

             10000004 function calls in 3.766 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     10000000    2.321    0.000    2.321    0.000 <ipython-input-2-12b138a62a96>:1(sum1d)
            1    1.445    1.445    3.766    3.766 <ipython-input-2-12b138a62a96>:9(sum2d)
            1    0.000    0.000    3.766    3.766 {built-in method builtins.exec}
            1    0.000    0.000    3.766    3.766 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    


There are a couple of things to notice here. 

First, the time taken to execute the profiling run is quite a bit longer than the time we got when benchmarking earlier. This is because profiling adds some overhead. Generally this doesn't affect the conclusions you would draw about which functions take the most time, but it's something to be aware of.

Second, most of the time is being taken up in the `sum1d` function, although a decent amount of time is also being spent within the `sum2d` function. You can see this from the 'tottime' column, which shows the total amount of time spent within a function, **not** including any calls made to other functions. The 'cumtime' column shows the total amount of time spent in a function, including any function calls.

Also, you'll notice that there were 10,000,004 function calls. Calling a Python function has some overhead. Maybe the code would go faster if we reduced the number of function calls? Here's a new implementation, combining everything into a single function:


{% highlight python %}
def sum2d_v2(ll):
    """Compute the sum of a list of lists of numbers."""
    s = 0
    for l in ll:
        for x in l:
            s += x
    return s

{% endhighlight %}


{% highlight python %}
%timeit sum2d_v2(big_data)
{% endhighlight %}

    1 loop, best of 3: 2.17 s per loop


This is a bit faster. What does the profiler tell us?


{% highlight python %}
cProfile.run('sum2d_v2(big_data)', sort='time')
{% endhighlight %}

             4 function calls in 2.212 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    2.212    2.212    2.212    2.212 <ipython-input-14-81d66843d00e>:1(sum2d_v2)
            1    0.000    0.000    2.212    2.212 {built-in method builtins.exec}
            1    0.000    0.000    2.212    2.212 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    


In fact we've hit a dead end here, because there is only a single function to profile, and function profiling cannot tell us which lines of code within a function are taking up most time. To get further we need to do some...

### Line profiling: `%lprun`, `line_profiler`

You can do line profiling with a Python module called [`line_profiler`](https://github.com/rkern/line_profiler). This is not part of the Python standard library, so you have to install it separately, e.g., via pip or conda.

For convenience, the `line_profiler` module provides a `%lprun` magic command for use in a Jupyter notebook, which can be used as follows: 


{% highlight python %}
%load_ext line_profiler
%lprun -f sum2d_v2 sum2d_v2(big_data)
{% endhighlight %}

You can also do the same thing via regular Python code:


{% highlight python %}
import line_profiler
l = line_profiler.LineProfiler()
l.add_function(sum2d_v2)
l.run('sum2d_v2(big_data)')
l.print_stats()
{% endhighlight %}

    Timer unit: 1e-06 s
    
    Total time: 30.7352 s
    File: <ipython-input-14-81d66843d00e>
    Function: sum2d_v2 at line 1
    
    Line #      Hits         Time  Per Hit   % Time  Line Contents
    ==============================================================
         1                                           def sum2d_v2(ll):
         2                                               """Compute the sum of a list of lists of numbers."""
         3         1            2      2.0      0.0      s = 0
         4  10000001      2352611      0.2      7.7      for l in ll:
         5  60000000     14782948      0.2     48.1          for x in l:
         6  50000000     13599659      0.3     44.2              s += x
         7         1            0      0.0      0.0      return s
    


Notice that this takes *a lot* longer than with just function profiling or without any profiling. Line profiling adds *a lot* more overhead, and this really can skew benchmarking results sometimes, so it's a good idea after each optimisation you make to benchmark without any profiling at all, as well as running function and line profiling.

If you are getting bored waiting for the line profiler to finish, you can interrupt it and it will still output some useful statistics.

Note that you have to explicitly tell `line_profiler` which functions to do line profiling within. When using the `%lprun` magic this is done via the `-f` option. 

Most of the time is spent in the inner for loop, iterating over the inner lists, and performing the addition. Python has a built-in `sum` function, maybe we could try that? ...


{% highlight python %}
def sum2d_v3(ll):
    """Compute the sum of a list of lists of numbers."""
    s = 0
    for l in ll:
        x = sum(l)
        s += x
    return s

{% endhighlight %}


{% highlight python %}
%timeit sum2d_v3(big_data)
{% endhighlight %}

    1 loop, best of 3: 1.86 s per loop


We've shaved off a bit more time. What do the profiling results tell us?


{% highlight python %}
cProfile.run('sum2d_v3(big_data)', sort='time')
{% endhighlight %}

             10000004 function calls in 2.657 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     10000000    1.338    0.000    1.338    0.000 {built-in method builtins.sum}
            1    1.320    1.320    2.657    2.657 <ipython-input-18-baa7cce51590>:1(sum2d_v3)
            1    0.000    0.000    2.657    2.657 {built-in method builtins.exec}
            1    0.000    0.000    2.657    2.657 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    



{% highlight python %}
import line_profiler
l = line_profiler.LineProfiler()
l.add_function(sum2d_v3)
l.run('sum2d_v3(big_data)')
l.print_stats()
{% endhighlight %}

    Timer unit: 1e-06 s
    
    Total time: 9.35954 s
    File: <ipython-input-18-baa7cce51590>
    Function: sum2d_v3 at line 1
    
    Line #      Hits         Time  Per Hit   % Time  Line Contents
    ==============================================================
         1                                           def sum2d_v3(ll):
         2                                               """Compute the sum of a list of lists of numbers."""
         3         1            3      3.0      0.0      s = 0
         4  10000001      2429112      0.2     26.0      for l in ll:
         5  10000000      4091269      0.4     43.7          x = sum(l)
         6  10000000      2839153      0.3     30.3          s += x
         7         1            1      1.0      0.0      return s
    


Now a decent amount of time is being spent inside the built-in `sum` function, and there's not much we can do about that. But there's also time being spent in the for loop, and in arithmetic. To go further, we need to try...

## NumPy

For numerical problems, the first port of call is [NumPy](http://www.numpy.org/). Let's use it to solve the sum2d problem. First, let's create a new test dataset, of the same size (10,000,000 rows, 5 columns), but using the `np.random.randint` function:


{% highlight python %}
import numpy as np
big_array = np.random.randint(0, 100, size=(10000000, 5))
big_array
{% endhighlight %}




    array([[56, 28, 51, 42, 25],
           [24, 71, 30, 56,  4],
           [35, 48, 50, 91, 17],
           ..., 
           [30, 78, 50, 97, 55],
           [71, 42, 19, 38, 89],
           [71, 40, 45, 92, 55]])



The `big_array` variable is a NumPy array. Here's a few useful properties:


{% highlight python %}
# number of dimensions
big_array.ndim
{% endhighlight %}




    2




{% highlight python %}
# size of each dimension
big_array.shape
{% endhighlight %}




    (10000000, 5)




{% highlight python %}
# data type of each array element
big_array.dtype
{% endhighlight %}




    dtype('int64')




{% highlight python %}
# number of bytes of memory used to store the data
big_array.nbytes
{% endhighlight %}




    400000000




{% highlight python %}
# some other features of the array
big_array.flags
{% endhighlight %}




      C_CONTIGUOUS : True
      F_CONTIGUOUS : False
      OWNDATA : True
      WRITEABLE : True
      ALIGNED : True
      UPDATEIFCOPY : False



NumPy also has its own `np.sum()` function which can operate on N-dimensional arrays. Let's try it:


{% highlight python %}
%timeit np.sum(big_array)
{% endhighlight %}

    10 loops, best of 3: 30.6 ms per loop


So using NumPy is almost two orders of magnitude faster than our own Python implementation. Where does the speed come from? NumPy's functions are implemented in C, so all of the looping and arithmetic is done in native C code. Also, a NumPy array stores it's data in a single, contiguous block of memory, which can be accessed very quickly and efficiently.

### Aside: array arithmetic

There are lots of things you can do with NumPy arrays, without ever having to write a for loop. E.g.:


{% highlight python %}
# column sum
np.sum(big_array, axis=0)
{% endhighlight %}




    array([494888839, 494827505, 495034687, 495112245, 494940377])




{% highlight python %}
# row sum
np.sum(big_array, axis=1)
{% endhighlight %}




    array([202, 185, 241, ..., 310, 259, 303])




{% highlight python %}
# add 2 to every array element
big_array + 2
{% endhighlight %}




    array([[58, 30, 53, 44, 27],
           [26, 73, 32, 58,  6],
           [37, 50, 52, 93, 19],
           ..., 
           [32, 80, 52, 99, 57],
           [73, 44, 21, 40, 91],
           [73, 42, 47, 94, 57]])




{% highlight python %}
# multiply every array element by 2
big_array * 2
{% endhighlight %}




    array([[112,  56, 102,  84,  50],
           [ 48, 142,  60, 112,   8],
           [ 70,  96, 100, 182,  34],
           ..., 
           [ 60, 156, 100, 194, 110],
           [142,  84,  38,  76, 178],
           [142,  80,  90, 184, 110]])




{% highlight python %}
# add two arrays element-by-element
big_array + big_array
{% endhighlight %}




    array([[112,  56, 102,  84,  50],
           [ 48, 142,  60, 112,   8],
           [ 70,  96, 100, 182,  34],
           ..., 
           [ 60, 156, 100, 194, 110],
           [142,  84,  38,  76, 178],
           [142,  80,  90, 184, 110]])




{% highlight python %}
# more complicated expressions
t = (big_array * 2) == (big_array + big_array)
t
{% endhighlight %}




    array([[ True,  True,  True,  True,  True],
           [ True,  True,  True,  True,  True],
           [ True,  True,  True,  True,  True],
           ..., 
           [ True,  True,  True,  True,  True],
           [ True,  True,  True,  True,  True],
           [ True,  True,  True,  True,  True]], dtype=bool)




{% highlight python %}
np.all(t)
{% endhighlight %}




    True



## Cython

For when you can't solve a problem with NumPy...

<img src='http://docs.cython.org/en/latest/_images/math/85505aa54782f6e6e6c113b9c562478082c1bbac.png'>


{% highlight python %}
def numpy_approx_pi(n):
    pi = np.sqrt(6 * np.sum(1/(np.arange(1, n+1)**2)))
    return pi

{% endhighlight %}


{% highlight python %}
%timeit numpy_approx_pi(10000000)
{% endhighlight %}


{% highlight python %}
def recip_square(i):
    s = 1. / i**2
    return s


def approx_pi(n):
    """Compute an approximate value of pi."""
    val = 0
    for k in range(1, n+1):
        x = recip_square(k)
        val += x
    pi = (6 * val)**.5
    return pi

{% endhighlight %}

From the Cython tutorial: **... remember the golden rule of optimization: Never optimize without having profiled. Let me repeat this: Never optimize without having profiled your code. Your thoughts about which part of your code takes too much time are wrong. At least, mine are always wrong.**


{% highlight python %}
%timeit approx_pi(10000000)
{% endhighlight %}


{% highlight python %}
%prun approx_pi(10000000)
{% endhighlight %}


{% highlight python %}
%load_ext cython
{% endhighlight %}


{% highlight python %}
%%cython
# cython: profile=True


def recip_square(i):
    s = 1. / i**2
    return s


def approx_pi_cy1(n):
    """Compute an approximate value of pi."""
    val = 0
    for k in range(1, n+1):
        x = recip_square(k)
        val += x
    pi = (6 * val)**.5
    return pi

{% endhighlight %}


{% highlight python %}
%timeit approx_pi_cy1(10000000)
{% endhighlight %}


{% highlight python %}
%prun approx_pi_cy1(10000000)
{% endhighlight %}

* Copy-paste code from above, rename `approx_pi_cy2`
* Add -a to inspect Cython's generated code.
* Explain yellow.
* Break up `recip_square`
* Add type to argument `i` in `recip_square`
* Add type arguments to other variables in `recip_square` and add `@cython.cdivision(True)` and `cimport cython`
* Add line profiling support:

<pre>
# cython: linetrace=True
# cython: binding=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
</pre>

* Add `int` type to k and n in `approx_pi_cy2` to optimise for loop
* Make `recip_square` a `cdef` function
* Add `double` declaration to `val` and return type from `recip_square`
* Remove line profiling support
* Add `@cython.profile(False)` to `recip_square`
* Make `recip_square` an `inline` function


{% highlight python %}
%%cython -a


cimport cython


@cython.cdivision(True)
cdef inline double recip_square(int i):
    cdef:
        double x, s
    x = i**2
    s = 1./x
    return s


def approx_pi_cy2(int n):
    """Compute an approximate value of pi."""
    cdef:
        long k
        double val
    val = 0
    for k in range(1, n+1):
        x = recip_square(k)
        val += x
    pi = (6 * val)**.5
    return pi

{% endhighlight %}


{% highlight python %}
%timeit approx_pi_cy2(10000000)
{% endhighlight %}


{% highlight python %}
%prun approx_pi_cy2(10000000)
{% endhighlight %}


{% highlight python %}
%lprun -f approx_pi_cy2 approx_pi_cy2(10000000)
{% endhighlight %}

## Cython and NumPy

Gain efficient, low-level access to NumPy arrays...

* Copy-paste sum2d_v2
* Change name, add -a
* Add numpy import
* Add type declaration to `ll`
* Change for loops, introduce typed variables `i` and `j`
* Type `s`
* Add cython import
* Add `boundscheck` and `wraparound` annotations


{% highlight python %}
%%cython -a


cimport numpy as np
cimport cython


@cython.wraparound(False)
@cython.boundscheck(False)
def sum2d_cy(np.int64_t[:, :] ll):
    """Compute the sum of a list of lists of numbers."""
    cdef:
        int i, j
        np.int64_t s
    s = 0
    for i in range(ll.shape[0]):
        for j in range(ll.shape[1]):
            s += ll[i, j]
    return s

{% endhighlight %}


{% highlight python %}
%time sum2d_cy(big_array)
{% endhighlight %}


{% highlight python %}
%time np.sum(big_array)
{% endhighlight %}

## Further reading...


{% highlight python %}

{% endhighlight %}
