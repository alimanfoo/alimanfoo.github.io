---
layout: post
title: Go faster Python
---


## Preamble

This blog post gives an introduction to some techniques for benchmarking, profiling and optimising Python code. If you would like to try the code examples for yourself, you can [download the Jupyter notebook](https://github.com/alimanfoo/alimanfoo.github.io/blob/master/_posts/2017-01-23-go-faster-python.ipynb) (right click the "Raw" button, save link as...) that this blog post was generated from. To run the notebook, you will need a working Python 3 installation, and will also need to install a couple of Python packages. The way I did that was to first [download](http://conda.pydata.org/miniconda.html) and [install](http://conda.pydata.org/docs/install/quick.html) miniconda into my home directory. I then ran the following from the command line:

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

My intuitions for why a piece of code is running slowly are very often wrong, and many other programmers more experienced than me [say the same thing](http://cython.readthedocs.io/en/latest/src/tutorial/profiling_tutorial.html#profiling-tutorial). If you're going to try and optimise something, you need to do some benchmarking and profiling first, to find out (1) exactly how slow it goes, and (2) where the bottleneck is.

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

    CPU times: user 2.89 s, sys: 0 ns, total: 2.89 s
    Wall time: 2.89 s
    
    2817000000



The first line of the output gives the amount of CPU time, broken down into 'user' (your code) and 'sys' (operating system code). The second line gives the wall time, which is the actual amount of time elapsed. Generally the total CPU time and the wall time will be the same, but sometimes not. E.g., if you are benchmarking a multi-threaded program, then wall time may be less than CPU time, because CPU time counts time spent by each CPU separately and adds them together, but the CPUs may actually be working in parallel.

One thing to watch out for when benchmarking is that performance can be variable, and may be affected by other processes running on your computer. To see this happen, try running the cell above again, but while it's running, give your computer something else to do at the same time, e.g., play some music, or a video, or just scroll the page up and down a bit.

To control for this variation, it's a good idea to benchmark several runs (and avoid the temptation to check your email while it's running). The [`%timeit`](http://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=%25timeit#magic-timeit) magic will automatically benchmark a piece of code several times: 


{% highlight python %}
%timeit sum2d(big_data)
{% endhighlight %}

    1 loop, best of 3: 2.84 s per loop


Alternatively, using the [`timeit`](https://docs.python.org/3/library/timeit.html) module:


{% highlight python %}
import timeit
timeit.repeat(stmt='sum2d(big_data)', repeat=3, number=1, globals=globals())
{% endhighlight %}




    [2.892044251999323, 2.8372464299973217, 2.8505056829999376]



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

             10000004 function calls in 3.883 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     10000000    2.378    0.000    2.378    0.000 <ipython-input-2-12b138a62a96>:1(sum1d)
            1    1.505    1.505    3.883    3.883 <ipython-input-2-12b138a62a96>:9(sum2d)
            1    0.000    0.000    3.883    3.883 {built-in method builtins.exec}
            1    0.000    0.000    3.883    3.883 <string>:1(<module>)
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

    1 loop, best of 3: 2.22 s per loop


This is a bit faster. What does the profiler tell us?


{% highlight python %}
cProfile.run('sum2d_v2(big_data)', sort='time')
{% endhighlight %}

             4 function calls in 2.223 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    2.223    2.223    2.223    2.223 <ipython-input-14-81d66843d00e>:1(sum2d_v2)
            1    0.000    0.000    2.223    2.223 {built-in method builtins.exec}
            1    0.000    0.000    2.223    2.223 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    


In fact we've hit a dead end here, because there is only a single function to profile, and function profiling cannot tell us which lines of code within a function are taking up most time. To get further we need to do some...

### Line profiling: `%lprun`, `line_profiler`

You can do line profiling with a Python module called [`line_profiler`](https://github.com/rkern/line_profiler). This is not part of the Python standard library, so you have to install it separately, e.g., via pip or conda.

For convenience, the `line_profiler` module provides a `%lprun` magic command for use in a Jupyter notebook, which can be used as follows: 


{% highlight python %}
%load_ext line_profiler
{% endhighlight %}


{% highlight python %}
%lprun -f sum2d_v2 sum2d_v2(big_data)
{% endhighlight %}

You can also do the same thing via regular Python code:


{% highlight python %}
import line_profiler
l = line_profiler.LineProfiler()
l.add_function(sum2d_v2)
l.run('sum2d_v2(big_data)')
{% endhighlight %}




    <line_profiler.LineProfiler at 0x7fcf947c4d68>




{% highlight python %}
l.print_stats()
{% endhighlight %}

    Timer unit: 1e-06 s
    
    Total time: 31.6946 s
    File: <ipython-input-14-81d66843d00e>
    Function: sum2d_v2 at line 1
    
    Line #      Hits         Time  Per Hit   % Time  Line Contents
    ==============================================================
         1                                           def sum2d_v2(ll):
         2                                               """Compute the sum of a list of lists of numbers."""
         3         1            2      2.0      0.0      s = 0
         4  10000001      2338393      0.2      7.4      for l in ll:
         5  60000000     15348131      0.3     48.4          for x in l:
         6  50000000     14008119      0.3     44.2              s += x
         7         1            0      0.0      0.0      return s
    


Notice that this takes *a lot* longer than with just function profiling or without any profiling. Line profiling adds *a lot* more overhead, and this can skew benchmarking results sometimes, so it's a good idea after each optimisation you make to benchmark without any profiling at all, as well as running function and line profiling.

If you are getting bored waiting for the line profiler to finish, you can interrupt it and it will still output some useful statistics.

Note that you have to explicitly tell `line_profiler` which functions to do line profiling within. When using the `%lprun` magic this is done via the `-f` option. 

Most of the time is spent in the inner for loop, iterating over the inner lists, and performing the addition. Python has a built-in `sum()` function, maybe we could try that? ...


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

    1 loop, best of 3: 1.91 s per loop


We've shaved off a bit more time. What do the profiling results tell us?


{% highlight python %}
cProfile.run('sum2d_v3(big_data)', sort='time')
{% endhighlight %}

             10000004 function calls in 2.714 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     10000000    1.358    0.000    1.358    0.000 {built-in method builtins.sum}
            1    1.356    1.356    2.714    2.714 <ipython-input-21-baa7cce51590>:1(sum2d_v3)
            1    0.000    0.000    2.714    2.714 {built-in method builtins.exec}
            1    0.000    0.000    2.714    2.714 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    



{% highlight python %}
import line_profiler
l = line_profiler.LineProfiler()
l.add_function(sum2d_v3)
l.run('sum2d_v3(big_data)')
l.print_stats()
{% endhighlight %}

    Timer unit: 1e-06 s
    
    Total time: 9.65201 s
    File: <ipython-input-21-baa7cce51590>
    Function: sum2d_v3 at line 1
    
    Line #      Hits         Time  Per Hit   % Time  Line Contents
    ==============================================================
         1                                           def sum2d_v3(ll):
         2                                               """Compute the sum of a list of lists of numbers."""
         3         1            1      1.0      0.0      s = 0
         4  10000001      2536840      0.3     26.3      for l in ll:
         5  10000000      4125789      0.4     42.7          x = sum(l)
         6  10000000      2989380      0.3     31.0          s += x
         7         1            1      1.0      0.0      return s
    


Now a decent amount of time is being spent inside the built-in `sum()` function, and there's not much we can do about that. But there's also still time being spent in the for loop, and in arithmetic. 

You've probably realised by now that this is a contrived example. We've actually been going around in circles a bit, and our attempts at optimisation haven't got us very far. To optimise further, we need to try something different...

## NumPy

For numerical problems, the first port of call is [NumPy](http://www.numpy.org/). Let's use it to solve the sum2d problem. First, let's create a new test dataset, of the same size (10,000,000 rows, 5 columns), but this time using the `np.random.randint()` function:


{% highlight python %}
import numpy as np
big_array = np.random.randint(0, 100, size=(10000000, 5))
big_array
{% endhighlight %}




    array([[34, 55, 45, 16, 83],
           [20, 34, 77, 55, 64],
           [45, 61, 65,  8, 37],
           ..., 
           [73, 86, 73,  5, 88],
           [70, 57, 92, 40, 21],
           [50, 46, 93, 38, 43]])



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

    10 loops, best of 3: 30.5 ms per loop


So using NumPy is almost two orders of magnitude faster than our own Python sum2d implementation. Where does the speed come from? NumPy's functions are implemented in C, so all of the looping and arithmetic is done in native C code. Also, a NumPy array stores it's data in a single, contiguous block of memory, which can be accessed quickly and efficiently.

### Aside: array arithmetic

There are lots of things you can do with NumPy arrays, without ever having to write a for loop. E.g.:


{% highlight python %}
# column sum
np.sum(big_array, axis=0)
{% endhighlight %}




    array([494856187, 495018159, 495000327, 494824068, 494760222])




{% highlight python %}
# row sum
np.sum(big_array, axis=1)
{% endhighlight %}




    array([233, 250, 216, ..., 325, 280, 270])




{% highlight python %}
# add 2 to every array element
big_array + 2
{% endhighlight %}




    array([[36, 57, 47, 18, 85],
           [22, 36, 79, 57, 66],
           [47, 63, 67, 10, 39],
           ..., 
           [75, 88, 75,  7, 90],
           [72, 59, 94, 42, 23],
           [52, 48, 95, 40, 45]])




{% highlight python %}
# multiply every array element by 2
big_array * 2
{% endhighlight %}




    array([[ 68, 110,  90,  32, 166],
           [ 40,  68, 154, 110, 128],
           [ 90, 122, 130,  16,  74],
           ..., 
           [146, 172, 146,  10, 176],
           [140, 114, 184,  80,  42],
           [100,  92, 186,  76,  86]])




{% highlight python %}
# add two arrays element-by-element
big_array + big_array
{% endhighlight %}




    array([[ 68, 110,  90,  32, 166],
           [ 40,  68, 154, 110, 128],
           [ 90, 122, 130,  16,  74],
           ..., 
           [146, 172, 146,  10, 176],
           [140, 114, 184,  80,  42],
           [100,  92, 186,  76,  86]])




{% highlight python %}
# more complicated operations
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

For some problems, there is no convenient way to solve them with NumPy alone. Or in some cases, the problem can be solved with NumPy, but the solution requires too much memory. If this is the case, another option is to try optimising with [Cython](http://docs.cython.org/en/latest/). Cython takes your Python code and transforms it into C code, which can then be compiled just like any other C code. A key feature is that Cython allows you to add static type declarations to certain variables, which then enables it to generate highly efficient native C code for certain critical sections of your code.

To illustrate the use of Cython I'm going to walk through an example from the Cython profiling tutorial. The task is to compute an approximate value for pi, using the formula below:

<img src='/assets/approx_pi.png'>

In fact there is a way to solve this problem using NumPy:


{% highlight python %}
def approx_pi_numpy(n):
    pi = (6 * np.sum(1 / (np.arange(1, n+1)**2)))**.5
    return pi

{% endhighlight %}


{% highlight python %}
approx_pi_numpy(10000000)
{% endhighlight %}




    3.1415925580968325




{% highlight python %}
%timeit approx_pi_numpy(10000000)
{% endhighlight %}

    10 loops, best of 3: 84.6 ms per loop


The NumPy solution is pretty quick, but it does require creating several reasonably large arrays in memory. If we wanted a higher precision estimate, we might run out of memory. Also, allocating memory has some overhead, and so we should be able to find a Cython solution that is even faster.

Let's start from a pure Python solution:


{% highlight python %}
def recip_square(i):
    x = i**2
    s = 1 / x
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


{% highlight python %}
approx_pi(10000000)
{% endhighlight %}




    3.1415925580959025



Benchmark and profile:


{% highlight python %}
%timeit approx_pi(10000000)
{% endhighlight %}

    1 loop, best of 3: 3.79 s per loop



{% highlight python %}
cProfile.run('approx_pi(10000000)', sort='time')
{% endhighlight %}

             10000004 function calls in 4.660 seconds
    
       Ordered by: internal time
    
       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     10000000    3.092    0.000    3.092    0.000 <ipython-input-42-2f54265e4169>:1(recip_square)
            1    1.568    1.568    4.660    4.660 <ipython-input-42-2f54265e4169>:7(approx_pi)
            1    0.000    0.000    4.660    4.660 {built-in method builtins.exec}
            1    0.000    0.000    4.660    4.660 <string>:1(<module>)
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    
    



{% highlight python %}
l = line_profiler.LineProfiler()
# l.add_function(recip_square)
l.add_function(approx_pi)
# use a smaller value of n, otherwise line profiling takes ages
l.run('approx_pi(2000000)')
l.print_stats()
{% endhighlight %}

    Timer unit: 1e-06 s
    
    Total time: 3.36972 s
    File: <ipython-input-42-2f54265e4169>
    Function: approx_pi at line 7
    
    Line #      Hits         Time  Per Hit   % Time  Line Contents
    ==============================================================
         7                                           def approx_pi(n):
         8                                               """Compute an approximate value of pi."""
         9         1            2      2.0      0.0      val = 0
        10   2000001       676782      0.3     20.1      for k in range(1, n+1):
        11   2000000      2007925      1.0     59.6          x = recip_square(k)
        12   2000000       685011      0.3     20.3          val += x
        13         1            5      5.0      0.0      pi = (6 * val)**.5
        14         1            0      0.0      0.0      return pi
    


Notice that I commented out line profiling for the `recip_square()` function. If you're running this yourself, I recommend running line profiling with and without including `recip_square()`, and each time take a look at the results for the outer `approx_pi()` function. This is a good example of how the overhead of line profiling can skew timings, so again, something to beware of. 

Now let's start constructing a Cython implementation. If you're working in a Jupyter notebook, there is a very convenient `%%cython` magic which enables you to write a Cython module within the notebook. When you execute a `%%cython` code cell, behind the scenes Cython will generate and compile your module, and import the functions into the current session so they can be called from other code cells. To make the `%%cython` magic available, we need to load the Cython notebook extension:


{% highlight python %}
%load_ext cython
{% endhighlight %}

To begin with, all I'll do is copy-paste the pure Python implementation into the Cython module, then change the function names so we can benchmark and profile separately:


{% highlight python %}
%%cython


def recip_square_cy1(i):
    x = i**2
    s = 1 / x
    return s


def approx_pi_cy1(n):
    """Compute an approximate value of pi."""
    val = 0
    for k in range(1, n+1):
        x = recip_square_cy1(k)
        val += x
    pi = (6 * val)**.5
    return pi

{% endhighlight %}


{% highlight python %}
approx_pi_cy1(10000000)
{% endhighlight %}




    3.1415925580959025




{% highlight python %}
%timeit approx_pi_cy1(10000000)
{% endhighlight %}

    1 loop, best of 3: 2.74 s per loop


Notice that the Cython function is a bit faster than the pure-Python function, even though the code is identical. However, we're still a long way short of the NumPy implementation.

What we're going to do next is make a number of modifications to the Cython implementation. If you're working through this in a notebook, I recommend copy-pasting the Cython code cell from above into a new empty code cell below, and changing the function names to `approx_py_cy2()` and `recip_square_cy2()`. Then, work through the optimisation steps below, one-by-one. After each step, run benchmarking and profiling to see where and how much speed has been gained, and examine the HTML diagnostics generated by Cython to see how much yellow (Python interaction) there is left. The goal is to remove all yellow from critical sections of the code. 

* Copy-paste Cython code cell from above into a new cell below; rename `approx_pi_cy1` to `approx_pi_cy2` and `recip_square_cy1` to `recip_square_cy2`. Be careful to rename all occurrences of the function names.
* Change `%%cython` to `%%cython -a` at the top of the cell. Adding the `-a` flag causes Cython to generate an HTML document with some diagnostics. This includes colouring of every line of code according to how much interaction is happening with Python. The more you can remove interaction with Python, the more Cython can optimise and generate efficient native C code.
* Add support for function profiling by adding the following comment (special comments like this are interpreted by Cython as compiler directives): 

<pre># cython: profile=True
</pre>

* Add support for line profiling by adding the following comments:

<pre>
# cython: linetrace=True
# cython: binding=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
</pre>

The following steps optimise the `recip_square` function:

* Add a static type declaration for the `i` argument to the `recip_square` function, by changing the function signature from `recip_square(i)` to `recip_square(int i)`.
* Add static type declarations for the `x` and `s` variables within the `recip_square` function, via a `cdef` section at the start of the function.
* Tell Cython to use C division instead of Python division within the `recip_square` function by adding the `@cython.cdivision(True)` annotation immediately above the function definition. This also requires adding the import statement `cimport cython`.

Now the `recip_square` function should be fully optimised. The following steps optimise the outer `approx_pi_cy2` function:

* Optimise the `for` loop by adding static type declarations for the `n` argument and the `k` variable.
* Make `recip_square` a `cdef` function and add a return type.
* Add a static type declaration for the `val` variable.
* Remove line profiling support.
* Remove function profiling support.
* Make `recip_square` an `inline` function.

Here is the end result...


{% highlight python %}
%%cython -a


cimport cython


@cython.cdivision(True)
cdef inline double recip_square_cy2(int i):
    cdef:
        double x, s
    x = i**2
    s = 1 / x
    return s


def approx_pi_cy2(int n):
    """Compute an approximate value of pi."""
    cdef:
        int k
        double val
    val = 0
    for k in range(1, n+1):
        x = recip_square_cy2(k)
        val += x
    pi = (6 * val)**.5
    return pi

{% endhighlight %}




<style type="text/css">
    
.cython { font-family: courier; font-size: 12; }

.cython.tag  {  }
.cython.line { margin: 0em }
.cython.code { font-size: 9; color: #444444; display: none; margin: 0px 0px 0px 8px; border-left: 8px none; }

.cython.line .run { background-color: #B0FFB0; }
.cython.line .mis { background-color: #FFB0B0; }
.cython.code.run  { border-left: 8px solid #B0FFB0; }
.cython.code.mis  { border-left: 8px solid #FFB0B0; }

.cython.code .py_c_api  { color: red; }
.cython.code .py_macro_api  { color: #FF7000; }
.cython.code .pyx_c_api  { color: #FF3000; }
.cython.code .pyx_macro_api  { color: #FF7000; }
.cython.code .refnanny  { color: #FFA000; }
.cython.code .trace  { color: #FFA000; }
.cython.code .error_goto  { color: #FFA000; }

.cython.code .coerce  { color: #008000; border: 1px dotted #008000 }
.cython.code .py_attr { color: #FF0000; font-weight: bold; }
.cython.code .c_attr  { color: #0000FF; }
.cython.code .py_call { color: #FF0000; font-weight: bold; }
.cython.code .c_call  { color: #0000FF; }

.cython.score-0 {background-color: #FFFFff;}
.cython.score-1 {background-color: #FFFFe7;}
.cython.score-2 {background-color: #FFFFd4;}
.cython.score-3 {background-color: #FFFFc4;}
.cython.score-4 {background-color: #FFFFb6;}
.cython.score-5 {background-color: #FFFFaa;}
.cython.score-6 {background-color: #FFFF9f;}
.cython.score-7 {background-color: #FFFF96;}
.cython.score-8 {background-color: #FFFF8d;}
.cython.score-9 {background-color: #FFFF86;}
.cython.score-10 {background-color: #FFFF7f;}
.cython.score-11 {background-color: #FFFF79;}
.cython.score-12 {background-color: #FFFF73;}
.cython.score-13 {background-color: #FFFF6e;}
.cython.score-14 {background-color: #FFFF6a;}
.cython.score-15 {background-color: #FFFF66;}
.cython.score-16 {background-color: #FFFF62;}
.cython.score-17 {background-color: #FFFF5e;}
.cython.score-18 {background-color: #FFFF5b;}
.cython.score-19 {background-color: #FFFF57;}
.cython.score-20 {background-color: #FFFF55;}
.cython.score-21 {background-color: #FFFF52;}
.cython.score-22 {background-color: #FFFF4f;}
.cython.score-23 {background-color: #FFFF4d;}
.cython.score-24 {background-color: #FFFF4b;}
.cython.score-25 {background-color: #FFFF48;}
.cython.score-26 {background-color: #FFFF46;}
.cython.score-27 {background-color: #FFFF44;}
.cython.score-28 {background-color: #FFFF43;}
.cython.score-29 {background-color: #FFFF41;}
.cython.score-30 {background-color: #FFFF3f;}
.cython.score-31 {background-color: #FFFF3e;}
.cython.score-32 {background-color: #FFFF3c;}
.cython.score-33 {background-color: #FFFF3b;}
.cython.score-34 {background-color: #FFFF39;}
.cython.score-35 {background-color: #FFFF38;}
.cython.score-36 {background-color: #FFFF37;}
.cython.score-37 {background-color: #FFFF36;}
.cython.score-38 {background-color: #FFFF35;}
.cython.score-39 {background-color: #FFFF34;}
.cython.score-40 {background-color: #FFFF33;}
.cython.score-41 {background-color: #FFFF32;}
.cython.score-42 {background-color: #FFFF31;}
.cython.score-43 {background-color: #FFFF30;}
.cython.score-44 {background-color: #FFFF2f;}
.cython.score-45 {background-color: #FFFF2e;}
.cython.score-46 {background-color: #FFFF2d;}
.cython.score-47 {background-color: #FFFF2c;}
.cython.score-48 {background-color: #FFFF2b;}
.cython.score-49 {background-color: #FFFF2b;}
.cython.score-50 {background-color: #FFFF2a;}
.cython.score-51 {background-color: #FFFF29;}
.cython.score-52 {background-color: #FFFF29;}
.cython.score-53 {background-color: #FFFF28;}
.cython.score-54 {background-color: #FFFF27;}
.cython.score-55 {background-color: #FFFF27;}
.cython.score-56 {background-color: #FFFF26;}
.cython.score-57 {background-color: #FFFF26;}
.cython.score-58 {background-color: #FFFF25;}
.cython.score-59 {background-color: #FFFF24;}
.cython.score-60 {background-color: #FFFF24;}
.cython.score-61 {background-color: #FFFF23;}
.cython.score-62 {background-color: #FFFF23;}
.cython.score-63 {background-color: #FFFF22;}
.cython.score-64 {background-color: #FFFF22;}
.cython.score-65 {background-color: #FFFF22;}
.cython.score-66 {background-color: #FFFF21;}
.cython.score-67 {background-color: #FFFF21;}
.cython.score-68 {background-color: #FFFF20;}
.cython.score-69 {background-color: #FFFF20;}
.cython.score-70 {background-color: #FFFF1f;}
.cython.score-71 {background-color: #FFFF1f;}
.cython.score-72 {background-color: #FFFF1f;}
.cython.score-73 {background-color: #FFFF1e;}
.cython.score-74 {background-color: #FFFF1e;}
.cython.score-75 {background-color: #FFFF1e;}
.cython.score-76 {background-color: #FFFF1d;}
.cython.score-77 {background-color: #FFFF1d;}
.cython.score-78 {background-color: #FFFF1c;}
.cython.score-79 {background-color: #FFFF1c;}
.cython.score-80 {background-color: #FFFF1c;}
.cython.score-81 {background-color: #FFFF1c;}
.cython.score-82 {background-color: #FFFF1b;}
.cython.score-83 {background-color: #FFFF1b;}
.cython.score-84 {background-color: #FFFF1b;}
.cython.score-85 {background-color: #FFFF1a;}
.cython.score-86 {background-color: #FFFF1a;}
.cython.score-87 {background-color: #FFFF1a;}
.cython.score-88 {background-color: #FFFF1a;}
.cython.score-89 {background-color: #FFFF19;}
.cython.score-90 {background-color: #FFFF19;}
.cython.score-91 {background-color: #FFFF19;}
.cython.score-92 {background-color: #FFFF19;}
.cython.score-93 {background-color: #FFFF18;}
.cython.score-94 {background-color: #FFFF18;}
.cython.score-95 {background-color: #FFFF18;}
.cython.score-96 {background-color: #FFFF18;}
.cython.score-97 {background-color: #FFFF17;}
.cython.score-98 {background-color: #FFFF17;}
.cython.score-99 {background-color: #FFFF17;}
.cython.score-100 {background-color: #FFFF17;}
.cython.score-101 {background-color: #FFFF16;}
.cython.score-102 {background-color: #FFFF16;}
.cython.score-103 {background-color: #FFFF16;}
.cython.score-104 {background-color: #FFFF16;}
.cython.score-105 {background-color: #FFFF16;}
.cython.score-106 {background-color: #FFFF15;}
.cython.score-107 {background-color: #FFFF15;}
.cython.score-108 {background-color: #FFFF15;}
.cython.score-109 {background-color: #FFFF15;}
.cython.score-110 {background-color: #FFFF15;}
.cython.score-111 {background-color: #FFFF15;}
.cython.score-112 {background-color: #FFFF14;}
.cython.score-113 {background-color: #FFFF14;}
.cython.score-114 {background-color: #FFFF14;}
.cython.score-115 {background-color: #FFFF14;}
.cython.score-116 {background-color: #FFFF14;}
.cython.score-117 {background-color: #FFFF14;}
.cython.score-118 {background-color: #FFFF13;}
.cython.score-119 {background-color: #FFFF13;}
.cython.score-120 {background-color: #FFFF13;}
.cython.score-121 {background-color: #FFFF13;}
.cython.score-122 {background-color: #FFFF13;}
.cython.score-123 {background-color: #FFFF13;}
.cython.score-124 {background-color: #FFFF13;}
.cython.score-125 {background-color: #FFFF12;}
.cython.score-126 {background-color: #FFFF12;}
.cython.score-127 {background-color: #FFFF12;}
.cython.score-128 {background-color: #FFFF12;}
.cython.score-129 {background-color: #FFFF12;}
.cython.score-130 {background-color: #FFFF12;}
.cython.score-131 {background-color: #FFFF12;}
.cython.score-132 {background-color: #FFFF11;}
.cython.score-133 {background-color: #FFFF11;}
.cython.score-134 {background-color: #FFFF11;}
.cython.score-135 {background-color: #FFFF11;}
.cython.score-136 {background-color: #FFFF11;}
.cython.score-137 {background-color: #FFFF11;}
.cython.score-138 {background-color: #FFFF11;}
.cython.score-139 {background-color: #FFFF11;}
.cython.score-140 {background-color: #FFFF11;}
.cython.score-141 {background-color: #FFFF10;}
.cython.score-142 {background-color: #FFFF10;}
.cython.score-143 {background-color: #FFFF10;}
.cython.score-144 {background-color: #FFFF10;}
.cython.score-145 {background-color: #FFFF10;}
.cython.score-146 {background-color: #FFFF10;}
.cython.score-147 {background-color: #FFFF10;}
.cython.score-148 {background-color: #FFFF10;}
.cython.score-149 {background-color: #FFFF10;}
.cython.score-150 {background-color: #FFFF0f;}
.cython.score-151 {background-color: #FFFF0f;}
.cython.score-152 {background-color: #FFFF0f;}
.cython.score-153 {background-color: #FFFF0f;}
.cython.score-154 {background-color: #FFFF0f;}
.cython.score-155 {background-color: #FFFF0f;}
.cython.score-156 {background-color: #FFFF0f;}
.cython.score-157 {background-color: #FFFF0f;}
.cython.score-158 {background-color: #FFFF0f;}
.cython.score-159 {background-color: #FFFF0f;}
.cython.score-160 {background-color: #FFFF0f;}
.cython.score-161 {background-color: #FFFF0e;}
.cython.score-162 {background-color: #FFFF0e;}
.cython.score-163 {background-color: #FFFF0e;}
.cython.score-164 {background-color: #FFFF0e;}
.cython.score-165 {background-color: #FFFF0e;}
.cython.score-166 {background-color: #FFFF0e;}
.cython.score-167 {background-color: #FFFF0e;}
.cython.score-168 {background-color: #FFFF0e;}
.cython.score-169 {background-color: #FFFF0e;}
.cython.score-170 {background-color: #FFFF0e;}
.cython.score-171 {background-color: #FFFF0e;}
.cython.score-172 {background-color: #FFFF0e;}
.cython.score-173 {background-color: #FFFF0d;}
.cython.score-174 {background-color: #FFFF0d;}
.cython.score-175 {background-color: #FFFF0d;}
.cython.score-176 {background-color: #FFFF0d;}
.cython.score-177 {background-color: #FFFF0d;}
.cython.score-178 {background-color: #FFFF0d;}
.cython.score-179 {background-color: #FFFF0d;}
.cython.score-180 {background-color: #FFFF0d;}
.cython.score-181 {background-color: #FFFF0d;}
.cython.score-182 {background-color: #FFFF0d;}
.cython.score-183 {background-color: #FFFF0d;}
.cython.score-184 {background-color: #FFFF0d;}
.cython.score-185 {background-color: #FFFF0d;}
.cython.score-186 {background-color: #FFFF0d;}
.cython.score-187 {background-color: #FFFF0c;}
.cython.score-188 {background-color: #FFFF0c;}
.cython.score-189 {background-color: #FFFF0c;}
.cython.score-190 {background-color: #FFFF0c;}
.cython.score-191 {background-color: #FFFF0c;}
.cython.score-192 {background-color: #FFFF0c;}
.cython.score-193 {background-color: #FFFF0c;}
.cython.score-194 {background-color: #FFFF0c;}
.cython.score-195 {background-color: #FFFF0c;}
.cython.score-196 {background-color: #FFFF0c;}
.cython.score-197 {background-color: #FFFF0c;}
.cython.score-198 {background-color: #FFFF0c;}
.cython.score-199 {background-color: #FFFF0c;}
.cython.score-200 {background-color: #FFFF0c;}
.cython.score-201 {background-color: #FFFF0c;}
.cython.score-202 {background-color: #FFFF0c;}
.cython.score-203 {background-color: #FFFF0b;}
.cython.score-204 {background-color: #FFFF0b;}
.cython.score-205 {background-color: #FFFF0b;}
.cython.score-206 {background-color: #FFFF0b;}
.cython.score-207 {background-color: #FFFF0b;}
.cython.score-208 {background-color: #FFFF0b;}
.cython.score-209 {background-color: #FFFF0b;}
.cython.score-210 {background-color: #FFFF0b;}
.cython.score-211 {background-color: #FFFF0b;}
.cython.score-212 {background-color: #FFFF0b;}
.cython.score-213 {background-color: #FFFF0b;}
.cython.score-214 {background-color: #FFFF0b;}
.cython.score-215 {background-color: #FFFF0b;}
.cython.score-216 {background-color: #FFFF0b;}
.cython.score-217 {background-color: #FFFF0b;}
.cython.score-218 {background-color: #FFFF0b;}
.cython.score-219 {background-color: #FFFF0b;}
.cython.score-220 {background-color: #FFFF0b;}
.cython.score-221 {background-color: #FFFF0b;}
.cython.score-222 {background-color: #FFFF0a;}
.cython.score-223 {background-color: #FFFF0a;}
.cython.score-224 {background-color: #FFFF0a;}
.cython.score-225 {background-color: #FFFF0a;}
.cython.score-226 {background-color: #FFFF0a;}
.cython.score-227 {background-color: #FFFF0a;}
.cython.score-228 {background-color: #FFFF0a;}
.cython.score-229 {background-color: #FFFF0a;}
.cython.score-230 {background-color: #FFFF0a;}
.cython.score-231 {background-color: #FFFF0a;}
.cython.score-232 {background-color: #FFFF0a;}
.cython.score-233 {background-color: #FFFF0a;}
.cython.score-234 {background-color: #FFFF0a;}
.cython.score-235 {background-color: #FFFF0a;}
.cython.score-236 {background-color: #FFFF0a;}
.cython.score-237 {background-color: #FFFF0a;}
.cython.score-238 {background-color: #FFFF0a;}
.cython.score-239 {background-color: #FFFF0a;}
.cython.score-240 {background-color: #FFFF0a;}
.cython.score-241 {background-color: #FFFF0a;}
.cython.score-242 {background-color: #FFFF0a;}
.cython.score-243 {background-color: #FFFF0a;}
.cython.score-244 {background-color: #FFFF0a;}
.cython.score-245 {background-color: #FFFF0a;}
.cython.score-246 {background-color: #FFFF09;}
.cython.score-247 {background-color: #FFFF09;}
.cython.score-248 {background-color: #FFFF09;}
.cython.score-249 {background-color: #FFFF09;}
.cython.score-250 {background-color: #FFFF09;}
.cython.score-251 {background-color: #FFFF09;}
.cython.score-252 {background-color: #FFFF09;}
.cython.score-253 {background-color: #FFFF09;}
.cython.score-254 {background-color: #FFFF09;}
.cython .hll { background-color: #ffffcc }
.cython  { background: #f8f8f8; }
.cython .c { color: #408080; font-style: italic } /* Comment */
.cython .err { border: 1px solid #FF0000 } /* Error */
.cython .k { color: #008000; font-weight: bold } /* Keyword */
.cython .o { color: #666666 } /* Operator */
.cython .ch { color: #408080; font-style: italic } /* Comment.Hashbang */
.cython .cm { color: #408080; font-style: italic } /* Comment.Multiline */
.cython .cp { color: #BC7A00 } /* Comment.Preproc */
.cython .cpf { color: #408080; font-style: italic } /* Comment.PreprocFile */
.cython .c1 { color: #408080; font-style: italic } /* Comment.Single */
.cython .cs { color: #408080; font-style: italic } /* Comment.Special */
.cython .gd { color: #A00000 } /* Generic.Deleted */
.cython .ge { font-style: italic } /* Generic.Emph */
.cython .gr { color: #FF0000 } /* Generic.Error */
.cython .gh { color: #000080; font-weight: bold } /* Generic.Heading */
.cython .gi { color: #00A000 } /* Generic.Inserted */
.cython .go { color: #888888 } /* Generic.Output */
.cython .gp { color: #000080; font-weight: bold } /* Generic.Prompt */
.cython .gs { font-weight: bold } /* Generic.Strong */
.cython .gu { color: #800080; font-weight: bold } /* Generic.Subheading */
.cython .gt { color: #0044DD } /* Generic.Traceback */
.cython .kc { color: #008000; font-weight: bold } /* Keyword.Constant */
.cython .kd { color: #008000; font-weight: bold } /* Keyword.Declaration */
.cython .kn { color: #008000; font-weight: bold } /* Keyword.Namespace */
.cython .kp { color: #008000 } /* Keyword.Pseudo */
.cython .kr { color: #008000; font-weight: bold } /* Keyword.Reserved */
.cython .kt { color: #B00040 } /* Keyword.Type */
.cython .m { color: #666666 } /* Literal.Number */
.cython .s { color: #BA2121 } /* Literal.String */
.cython .na { color: #7D9029 } /* Name.Attribute */
.cython .nb { color: #008000 } /* Name.Builtin */
.cython .nc { color: #0000FF; font-weight: bold } /* Name.Class */
.cython .no { color: #880000 } /* Name.Constant */
.cython .nd { color: #AA22FF } /* Name.Decorator */
.cython .ni { color: #999999; font-weight: bold } /* Name.Entity */
.cython .ne { color: #D2413A; font-weight: bold } /* Name.Exception */
.cython .nf { color: #0000FF } /* Name.Function */
.cython .nl { color: #A0A000 } /* Name.Label */
.cython .nn { color: #0000FF; font-weight: bold } /* Name.Namespace */
.cython .nt { color: #008000; font-weight: bold } /* Name.Tag */
.cython .nv { color: #19177C } /* Name.Variable */
.cython .ow { color: #AA22FF; font-weight: bold } /* Operator.Word */
.cython .w { color: #bbbbbb } /* Text.Whitespace */
.cython .mb { color: #666666 } /* Literal.Number.Bin */
.cython .mf { color: #666666 } /* Literal.Number.Float */
.cython .mh { color: #666666 } /* Literal.Number.Hex */
.cython .mi { color: #666666 } /* Literal.Number.Integer */
.cython .mo { color: #666666 } /* Literal.Number.Oct */
.cython .sb { color: #BA2121 } /* Literal.String.Backtick */
.cython .sc { color: #BA2121 } /* Literal.String.Char */
.cython .sd { color: #BA2121; font-style: italic } /* Literal.String.Doc */
.cython .s2 { color: #BA2121 } /* Literal.String.Double */
.cython .se { color: #BB6622; font-weight: bold } /* Literal.String.Escape */
.cython .sh { color: #BA2121 } /* Literal.String.Heredoc */
.cython .si { color: #BB6688; font-weight: bold } /* Literal.String.Interpol */
.cython .sx { color: #008000 } /* Literal.String.Other */
.cython .sr { color: #BB6688 } /* Literal.String.Regex */
.cython .s1 { color: #BA2121 } /* Literal.String.Single */
.cython .ss { color: #19177C } /* Literal.String.Symbol */
.cython .bp { color: #008000 } /* Name.Builtin.Pseudo */
.cython .vc { color: #19177C } /* Name.Variable.Class */
.cython .vg { color: #19177C } /* Name.Variable.Global */
.cython .vi { color: #19177C } /* Name.Variable.Instance */
.cython .il { color: #666666 } /* Literal.Number.Integer.Long */
</style>
<script>
    function toggleDiv(id) {
        theDiv = id.nextElementSibling
        if (theDiv.style.display != 'block') theDiv.style.display = 'block';
        else theDiv.style.display = 'none';
    }
</script>
<div class="cython">
<p><span style="border-bottom: solid 1px grey;">Generated by Cython 0.25.2</span></p>
<p>
    <span style="background-color: #FFFF00">Yellow lines</span> hint at Python interaction.<br />
    Click on a line that starts with a "<code>+</code>" to see the C code that Cython generated for it.
</p>
<div class="cython"><pre class="cython line score-0">&#xA0;<span class="">01</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">02</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">03</span>: <span class="k">cimport</span> <span class="nn">cython</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">04</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">05</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">06</span>: <span class="nd">@cython</span><span class="o">.</span><span class="n">cdivision</span><span class="p">(</span><span class="bp">True</span><span class="p">)</span></pre>
<pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">07</span>: <span class="k">cdef</span> <span class="kr">inline</span> <span class="kt">double</span> <span class="nf">recip_square_cy2</span><span class="p">(</span><span class="nb">int</span> <span class="n">i</span><span class="p">):</span></pre>
<pre class='cython code score-0 '>static CYTHON_INLINE double __pyx_f_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_recip_square_cy2(int __pyx_v_i) {
  double __pyx_v_x;
  double __pyx_v_s;
  double __pyx_r;
  <span class='refnanny'>__Pyx_RefNannyDeclarations</span>
  <span class='refnanny'>__Pyx_RefNannySetupContext</span>("recip_square_cy2", 0);
/* … */
  /* function exit code */
  __pyx_L0:;
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return __pyx_r;
}
</pre><pre class="cython line score-0">&#xA0;<span class="">08</span>:     <span class="k">cdef</span><span class="p">:</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">09</span>:         <span class="n">double</span> <span class="n">x</span><span class="p">,</span> <span class="n">s</span></pre>
<pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">10</span>:     <span class="n">x</span> <span class="o">=</span> <span class="n">i</span><span class="o">**</span><span class="mf">2</span></pre>
<pre class='cython code score-0 '>  __pyx_v_x = __Pyx_pow_long(((long)__pyx_v_i), 2);
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">11</span>:     <span class="n">s</span> <span class="o">=</span> <span class="mf">1</span> <span class="o">/</span> <span class="n">x</span></pre>
<pre class='cython code score-0 '>  __pyx_v_s = (1.0 / __pyx_v_x);
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">12</span>:     <span class="k">return</span> <span class="n">s</span></pre>
<pre class='cython code score-0 '>  __pyx_r = __pyx_v_s;
  goto __pyx_L0;
</pre><pre class="cython line score-0">&#xA0;<span class="">13</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">14</span>: </pre>
<pre class="cython line score-23" onclick='toggleDiv(this)'>+<span class="">15</span>: <span class="k">def</span> <span class="nf">approx_pi_cy2</span><span class="p">(</span><span class="nb">int</span> <span class="n">n</span><span class="p">):</span></pre>
<pre class='cython code score-23 '>/* Python wrapper */
static PyObject *__pyx_pw_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_1approx_pi_cy2(PyObject *__pyx_self, PyObject *__pyx_arg_n); /*proto*/
static char __pyx_doc_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_approx_pi_cy2[] = "Compute an approximate value of pi.";
static PyMethodDef __pyx_mdef_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_1approx_pi_cy2 = {"approx_pi_cy2", (PyCFunction)__pyx_pw_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_1approx_pi_cy2, METH_O, __pyx_doc_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_approx_pi_cy2};
static PyObject *__pyx_pw_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_1approx_pi_cy2(PyObject *__pyx_self, PyObject *__pyx_arg_n) {
  int __pyx_v_n;
  PyObject *__pyx_r = 0;
  <span class='refnanny'>__Pyx_RefNannyDeclarations</span>
  <span class='refnanny'>__Pyx_RefNannySetupContext</span>("approx_pi_cy2 (wrapper)", 0);
  assert(__pyx_arg_n); {
    __pyx_v_n = <span class='pyx_c_api'>__Pyx_PyInt_As_int</span>(__pyx_arg_n); if (unlikely((__pyx_v_n == (int)-1) &amp;&amp; <span class='py_c_api'>PyErr_Occurred</span>())) __PYX_ERR(0, 15, __pyx_L3_error)
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L3_error:;
  <span class='pyx_c_api'>__Pyx_AddTraceback</span>("_cython_magic_c1207dcff1d0d009517ab1cb3466abf5.approx_pi_cy2", __pyx_clineno, __pyx_lineno, __pyx_filename);
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_r = __pyx_pf_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_approx_pi_cy2(__pyx_self, ((int)__pyx_v_n));

  /* function exit code */
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return __pyx_r;
}

static PyObject *__pyx_pf_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_approx_pi_cy2(CYTHON_UNUSED PyObject *__pyx_self, int __pyx_v_n) {
  int __pyx_v_k;
  double __pyx_v_val;
  double __pyx_v_x;
  double __pyx_v_pi;
  PyObject *__pyx_r = NULL;
  <span class='refnanny'>__Pyx_RefNannyDeclarations</span>
  <span class='refnanny'>__Pyx_RefNannySetupContext</span>("approx_pi_cy2", 0);
/* … */
  /* function exit code */
  __pyx_L1_error:;
  <span class='pyx_macro_api'>__Pyx_XDECREF</span>(__pyx_t_3);
  <span class='pyx_c_api'>__Pyx_AddTraceback</span>("_cython_magic_c1207dcff1d0d009517ab1cb3466abf5.approx_pi_cy2", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = NULL;
  __pyx_L0:;
  <span class='refnanny'>__Pyx_XGIVEREF</span>(__pyx_r);
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return __pyx_r;
}
/* … */
  __pyx_tuple_ = <span class='py_c_api'>PyTuple_Pack</span>(6, __pyx_n_s_n, __pyx_n_s_n, __pyx_n_s_k, __pyx_n_s_val, __pyx_n_s_x, __pyx_n_s_pi); if (unlikely(!__pyx_tuple_)) __PYX_ERR(0, 15, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_tuple_);
  <span class='refnanny'>__Pyx_GIVEREF</span>(__pyx_tuple_);
/* … */
  __pyx_t_1 = PyCFunction_NewEx(&amp;__pyx_mdef_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_1approx_pi_cy2, NULL, __pyx_n_s_cython_magic_c1207dcff1d0d00951); if (unlikely(!__pyx_t_1)) __PYX_ERR(0, 15, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_t_1);
  if (<span class='py_c_api'>PyDict_SetItem</span>(__pyx_d, __pyx_n_s_approx_pi_cy2, __pyx_t_1) &lt; 0) __PYX_ERR(0, 15, __pyx_L1_error)
  <span class='pyx_macro_api'>__Pyx_DECREF</span>(__pyx_t_1); __pyx_t_1 = 0;
</pre><pre class="cython line score-0">&#xA0;<span class="">16</span>:     <span class="sd">&quot;&quot;&quot;Compute an approximate value of pi.&quot;&quot;&quot;</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">17</span>:     <span class="k">cdef</span><span class="p">:</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">18</span>:         <span class="nb">int</span> <span class="n">k</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">19</span>:         <span class="n">double</span> <span class="n">val</span></pre>
<pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">20</span>:     <span class="n">val</span> <span class="o">=</span> <span class="mf">0</span></pre>
<pre class='cython code score-0 '>  __pyx_v_val = 0.0;
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">21</span>:     <span class="k">for</span> <span class="n">k</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="mf">1</span><span class="p">,</span> <span class="n">n</span><span class="o">+</span><span class="mf">1</span><span class="p">):</span></pre>
<pre class='cython code score-0 '>  __pyx_t_1 = (__pyx_v_n + 1);
  for (__pyx_t_2 = 1; __pyx_t_2 &lt; __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_k = __pyx_t_2;
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">22</span>:         <span class="n">x</span> <span class="o">=</span> <span class="n">recip_square_cy2</span><span class="p">(</span><span class="n">k</span><span class="p">)</span></pre>
<pre class='cython code score-0 '>    __pyx_v_x = __pyx_f_46_cython_magic_c1207dcff1d0d009517ab1cb3466abf5_recip_square_cy2(__pyx_v_k);
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">23</span>:         <span class="n">val</span> <span class="o">+=</span> <span class="n">x</span></pre>
<pre class='cython code score-0 '>    __pyx_v_val = (__pyx_v_val + __pyx_v_x);
  }
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">24</span>:     <span class="n">pi</span> <span class="o">=</span> <span class="p">(</span><span class="mf">6</span> <span class="o">*</span> <span class="n">val</span><span class="p">)</span><span class="o">**.</span><span class="mf">5</span></pre>
<pre class='cython code score-0 '>  __pyx_v_pi = pow((6.0 * __pyx_v_val), .5);
</pre><pre class="cython line score-6" onclick='toggleDiv(this)'>+<span class="">25</span>:     <span class="k">return</span> <span class="n">pi</span></pre>
<pre class='cython code score-6 '>  <span class='pyx_macro_api'>__Pyx_XDECREF</span>(__pyx_r);
  __pyx_t_3 = <span class='py_c_api'>PyFloat_FromDouble</span>(__pyx_v_pi); if (unlikely(!__pyx_t_3)) __PYX_ERR(0, 25, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_t_3);
  __pyx_r = __pyx_t_3;
  __pyx_t_3 = 0;
  goto __pyx_L0;
</pre></div></div>




{% highlight python %}
%timeit approx_pi_cy2(10000000)
{% endhighlight %}

    100 loops, best of 3: 11.6 ms per loop



{% highlight python %}
cProfile.run('approx_pi_cy2(10000000)', sort='time')
{% endhighlight %}


{% highlight python %}
l = line_profiler.LineProfiler()
# l.add_function(recip_square_cy2)
l.add_function(approx_pi_cy2)
# use a smaller value of n, otherwise line profiling takes ages
l.run('approx_pi_cy2(2000000)')
l.print_stats()
{% endhighlight %}

So the optimised Cython implementation is about 8 times faster than the NumPy implementation. You may have noticed that a couple of the optimisation steps didn't actually make much difference to performance, and so weren't really necessary.

## Cython and NumPy

Finally, here's a Cython implementation of the sum2d function, just to give an example of how to use Cython and NumPy together. As above, you may want to introduce the optimisations one at a time, and benchmark and profile after each step, to get a sense of which optimisations really make a difference.


{% highlight python %}
%%cython -a


cimport numpy as np
cimport cython


@cython.wraparound(False)
@cython.boundscheck(False)
def sum2d_cy(np.int64_t[:, :] ll):
    """Compute the sum of a 2-dimensional array of integers."""
    cdef:
        int i, j
        np.int64_t s
    s = 0
    for i in range(ll.shape[0]):
        for j in range(ll.shape[1]):
            s += ll[i, j]
    return s

{% endhighlight %}




<div class="cython">
<p><span style="border-bottom: solid 1px grey;">Generated by Cython 0.25.2</span></p>
<p>
    <span style="background-color: #FFFF00">Yellow lines</span> hint at Python interaction.<br />
    Click on a line that starts with a "<code>+</code>" to see the C code that Cython generated for it.
</p>
<div class="cython"><pre class="cython line score-0">&#xA0;<span class="">01</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">02</span>: </pre>
<pre class="cython line score-11" onclick='toggleDiv(this)'>+<span class="">03</span>: <span class="k">cimport</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span></pre>
<pre class='cython code score-11 '>  __pyx_t_1 = <span class='py_c_api'>PyDict_New</span>(); if (unlikely(!__pyx_t_1)) __PYX_ERR(0, 3, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_t_1);
  if (<span class='py_c_api'>PyDict_SetItem</span>(__pyx_d, __pyx_n_s_test, __pyx_t_1) &lt; 0) __PYX_ERR(0, 3, __pyx_L1_error)
  <span class='pyx_macro_api'>__Pyx_DECREF</span>(__pyx_t_1); __pyx_t_1 = 0;
</pre><pre class="cython line score-0">&#xA0;<span class="">04</span>: <span class="k">cimport</span> <span class="nn">cython</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">05</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">06</span>: </pre>
<pre class="cython line score-0">&#xA0;<span class="">07</span>: <span class="nd">@cython</span><span class="o">.</span><span class="n">wraparound</span><span class="p">(</span><span class="bp">False</span><span class="p">)</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">08</span>: <span class="nd">@cython</span><span class="o">.</span><span class="n">boundscheck</span><span class="p">(</span><span class="bp">False</span><span class="p">)</span></pre>
<pre class="cython line score-18" onclick='toggleDiv(this)'>+<span class="">09</span>: <span class="k">def</span> <span class="nf">sum2d_cy</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">int64_t</span><span class="p">[:,</span> <span class="p">:]</span> <span class="n">ll</span><span class="p">):</span></pre>
<pre class='cython code score-18 '>/* Python wrapper */
static PyObject *__pyx_pw_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_1sum2d_cy(PyObject *__pyx_self, PyObject *__pyx_arg_ll); /*proto*/
static char __pyx_doc_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_sum2d_cy[] = "Compute the sum of a 2-dimensional array of integers.";
static PyMethodDef __pyx_mdef_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_1sum2d_cy = {"sum2d_cy", (PyCFunction)__pyx_pw_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_1sum2d_cy, METH_O, __pyx_doc_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_sum2d_cy};
static PyObject *__pyx_pw_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_1sum2d_cy(PyObject *__pyx_self, PyObject *__pyx_arg_ll) {
  __Pyx_memviewslice __pyx_v_ll = { 0, 0, { 0 }, { 0 }, { 0 } };
  PyObject *__pyx_r = 0;
  <span class='refnanny'>__Pyx_RefNannyDeclarations</span>
  <span class='refnanny'>__Pyx_RefNannySetupContext</span>("sum2d_cy (wrapper)", 0);
  assert(__pyx_arg_ll); {
    __pyx_v_ll = __Pyx_PyObject_to_MemoryviewSlice_dsds_nn___pyx_t_5numpy_int64_t(__pyx_arg_ll); if (unlikely(!__pyx_v_ll.memview)) __PYX_ERR(0, 9, __pyx_L3_error)
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L3_error:;
  <span class='pyx_c_api'>__Pyx_AddTraceback</span>("_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24.sum2d_cy", __pyx_clineno, __pyx_lineno, __pyx_filename);
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_r = __pyx_pf_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_sum2d_cy(__pyx_self, __pyx_v_ll);

  /* function exit code */
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return __pyx_r;
}

static PyObject *__pyx_pf_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_sum2d_cy(CYTHON_UNUSED PyObject *__pyx_self, __Pyx_memviewslice __pyx_v_ll) {
  int __pyx_v_i;
  int __pyx_v_j;
  __pyx_t_5numpy_int64_t __pyx_v_s;
  PyObject *__pyx_r = NULL;
  <span class='refnanny'>__Pyx_RefNannyDeclarations</span>
  <span class='refnanny'>__Pyx_RefNannySetupContext</span>("sum2d_cy", 0);
/* … */
  /* function exit code */
  __pyx_L1_error:;
  <span class='pyx_macro_api'>__Pyx_XDECREF</span>(__pyx_t_7);
  <span class='pyx_c_api'>__Pyx_AddTraceback</span>("_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24.sum2d_cy", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = NULL;
  __pyx_L0:;
  __PYX_XDEC_MEMVIEW(&amp;__pyx_v_ll, 1);
  <span class='refnanny'>__Pyx_XGIVEREF</span>(__pyx_r);
  <span class='refnanny'>__Pyx_RefNannyFinishContext</span>();
  return __pyx_r;
}
/* … */
  __pyx_tuple__23 = <span class='py_c_api'>PyTuple_Pack</span>(5, __pyx_n_s_ll, __pyx_n_s_ll, __pyx_n_s_i, __pyx_n_s_j, __pyx_n_s_s); if (unlikely(!__pyx_tuple__23)) __PYX_ERR(0, 9, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_tuple__23);
  <span class='refnanny'>__Pyx_GIVEREF</span>(__pyx_tuple__23);
/* … */
  __pyx_t_1 = PyCFunction_NewEx(&amp;__pyx_mdef_46_cython_magic_e638c34ce8e90ea82c9e0687d99dcc24_1sum2d_cy, NULL, __pyx_n_s_cython_magic_e638c34ce8e90ea82c); if (unlikely(!__pyx_t_1)) __PYX_ERR(0, 9, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_t_1);
  if (<span class='py_c_api'>PyDict_SetItem</span>(__pyx_d, __pyx_n_s_sum2d_cy, __pyx_t_1) &lt; 0) __PYX_ERR(0, 9, __pyx_L1_error)
  <span class='pyx_macro_api'>__Pyx_DECREF</span>(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_codeobj__24 = (PyObject*)<span class='pyx_c_api'>__Pyx_PyCode_New</span>(1, 0, 5, 0, 0, __pyx_empty_bytes, __pyx_empty_tuple, __pyx_empty_tuple, __pyx_tuple__23, __pyx_empty_tuple, __pyx_empty_tuple, __pyx_kp_s_home_aliman_cache_ipython_cytho, __pyx_n_s_sum2d_cy, 9, __pyx_empty_bytes); if (unlikely(!__pyx_codeobj__24)) __PYX_ERR(0, 9, __pyx_L1_error)
</pre><pre class="cython line score-0">&#xA0;<span class="">10</span>:     <span class="sd">&quot;&quot;&quot;Compute the sum of a 2-dimensional array of integers.&quot;&quot;&quot;</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">11</span>:     <span class="k">cdef</span><span class="p">:</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">12</span>:         <span class="nb">int</span> <span class="n">i</span><span class="p">,</span> <span class="n">j</span></pre>
<pre class="cython line score-0">&#xA0;<span class="">13</span>:         <span class="n">np</span><span class="o">.</span><span class="n">int64_t</span> <span class="n">s</span></pre>
<pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">14</span>:     <span class="n">s</span> <span class="o">=</span> <span class="mf">0</span></pre>
<pre class='cython code score-0 '>  __pyx_v_s = 0;
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">15</span>:     <span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">ll</span><span class="o">.</span><span class="n">shape</span><span class="p">[</span><span class="mf">0</span><span class="p">]):</span></pre>
<pre class='cython code score-0 '>  __pyx_t_1 = (__pyx_v_ll.shape[0]);
  for (__pyx_t_2 = 0; __pyx_t_2 &lt; __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">16</span>:         <span class="k">for</span> <span class="n">j</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">ll</span><span class="o">.</span><span class="n">shape</span><span class="p">[</span><span class="mf">1</span><span class="p">]):</span></pre>
<pre class='cython code score-0 '>    __pyx_t_3 = (__pyx_v_ll.shape[1]);
    for (__pyx_t_4 = 0; __pyx_t_4 &lt; __pyx_t_3; __pyx_t_4+=1) {
      __pyx_v_j = __pyx_t_4;
</pre><pre class="cython line score-0" onclick='toggleDiv(this)'>+<span class="">17</span>:             <span class="n">s</span> <span class="o">+=</span> <span class="n">ll</span><span class="p">[</span><span class="n">i</span><span class="p">,</span> <span class="n">j</span><span class="p">]</span></pre>
<pre class='cython code score-0 '>      __pyx_t_5 = __pyx_v_i;
      __pyx_t_6 = __pyx_v_j;
      __pyx_v_s = (__pyx_v_s + (*((__pyx_t_5numpy_int64_t *) ( /* dim=1 */ (( /* dim=0 */ (__pyx_v_ll.data + __pyx_t_5 * __pyx_v_ll.strides[0]) ) + __pyx_t_6 * __pyx_v_ll.strides[1]) ))));
    }
  }
</pre><pre class="cython line score-1" onclick='toggleDiv(this)'>+<span class="">18</span>:     <span class="k">return</span> <span class="n">s</span></pre>
<pre class='cython code score-1 '>  <span class='pyx_macro_api'>__Pyx_XDECREF</span>(__pyx_r);
  __pyx_t_7 = __Pyx_PyInt_From_npy_int64(__pyx_v_s); if (unlikely(!__pyx_t_7)) __PYX_ERR(0, 18, __pyx_L1_error)
  <span class='refnanny'>__Pyx_GOTREF</span>(__pyx_t_7);
  __pyx_r = __pyx_t_7;
  __pyx_t_7 = 0;
  goto __pyx_L0;
</pre></div></div>




{% highlight python %}
%timeit sum2d_cy(big_array)
{% endhighlight %}

    10 loops, best of 3: 34.3 ms per loop



{% highlight python %}
%timeit np.sum(big_array)
{% endhighlight %}

    10 loops, best of 3: 30.4 ms per loop


In this case the Cython and NumPy implementations are nearly identical, so there is no value in using Cython.

## Further reading

There are loads of good resources on the Web on the topics covered here. Here's just a few:

* [The Python profilers](https://docs.python.org/3.6/library/profile.html)
* [A guide to analyzing Python performance](https://www.huyng.com/posts/python-performance-analysis) by Huy Nguyen
* [SnakeViz](https://jiffyclub.github.io/snakeviz/)
* [NumPy tutorial](https://docs.scipy.org/doc/numpy-dev/user/quickstart.html)
* [Cython tutorials](http://cython.readthedocs.io/en/latest/src/tutorial/)
* [The fallacy of premature optimization](http://ubiquity.acm.org/article.cfm?id=1513451) by Randall Hyde
