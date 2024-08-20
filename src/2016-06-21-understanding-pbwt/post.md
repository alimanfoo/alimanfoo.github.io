---
layout: post
title: Understanding PBWT
---


Last year a colleague pointed me at [Richard Durbin's 2014 paper on the positional Burrows Wheeler transform](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3998136/). I read it with interest but didn't really get it at the time. Recently there have been papers developing PBWT for [variation graphs](http://biorxiv.org/content/early/2016/05/02/051409.abstract) and [chromosome painting](http://biorxiv.org/content/early/2016/04/12/048280.abstract). Given the applicability of PBWT to a range of problems, I wanted to try and wrap my head around it. This notebook is an attempt to understand the basic algorithms described in Durbin (2014), I wrote it for my own benefit but thought I'd share in case it's useful to anyone and/or anyone might be able to help me fill in some of the gaps.

*2016-06-23: Updated implementation of algorithm 4 to avoid reporting 0 length matches. Also added some discussion of the value of using PBWT to improve compressibility of haplotype data versus simply transposing the haplotypes so the data are organised one variant per row.*

Let's start with some data. The variable *X* below is a list of 8 haplotypes with 6 variable sites, all biallelic. This is just some dummy data I made up myself to test out the algorithms


{% highlight python %}
X = [[0, 1, 0, 1, 0, 1],
     [1, 1, 0, 0, 0, 1],
     [1, 1, 1, 1, 1, 1],
     [0, 1, 1, 1, 1, 0],
     [0, 0, 0, 0, 0, 0],
     [1, 0, 0, 0, 1, 0],
     [1, 1, 0, 0, 0, 1],
     [0, 1, 0, 1, 1, 0]]
{% endhighlight %}

## Algorithm 1. Build prefix array

The basic idea of PBWT is to sort the haplotypes in order of reversed prefixes. This sorting then enables matches between haplotypes to be found efficiently. Algorithm 1 in Durbin (2014) builds something called the "positional prefix array", which is just the list of haplotype indices that would sort the haplotypes in reversed prefix order at some position *k*. 

Durbin uses the notation *a<sub>k</sub>* to denote the positional prefix array for position *k*, however he also uses the variable *a* for one of the intermediates within algorithm 1, so to avoid any ambiguity here I will use *ppa<sub>k</sub>* to denote the positional prefix array at position *k*.

The trick of algorithm 1 is to sweep through the data by position, building *ppa<sub>k+1</sub>* from *ppa<sub>k</sub>*. The algorithm simply sorts the haplotypes by allele at position *k*, making sure that the sort order from the previous position is retained for all haplotypes with the same allele value.

Here is a pure Python implementation of algorithm 1, purely for illustration purposes.


{% highlight python %}
def build_prefix_array(X):
    
    # M haplotypes
    M = len(X)

    # N variable sites
    N = len(X[0])
    
    # initialise positional prefix array
    ppa = list(range(M))
    
    # iterate over variants
    for k in range(N-1):

        # setup intermediates
        a = list()
        b = list()
    
        # iterate over haplotypes in reverse prefix sorted order
        for index in ppa:

            # current haplotype
            haplotype = X[index]
            
            # allele for current haplotype
            allele = haplotype[k]
            
            # update intermediates
            if allele == 0:
                a.append(index)
            else:
                b.append(index)
                
        # construct the new positional prefix array for k+1 by concatenating lists a and b
        ppa = a + b
        
    return ppa
    
{% endhighlight %}

Let's run it on the dummy data *X*.


{% highlight python %}
build_prefix_array(X)
{% endhighlight %}




    [4, 1, 6, 0, 5, 7, 3, 2]



The resulting list of integers is the list of haplotype indices that sorts the haplotypes up to the final position (*k* = *N*-1). To help visualise this, let's write a display function.


{% highlight python %}
def display_prefix_array(X):
    from IPython.display import display_html
    ppa = build_prefix_array(X)
    html = '<pre style="line-height: 100%">'
    for index in ppa:
        haplotype = X[index]
        html += str(index) + '|' + ''.join(str(allele) for allele in haplotype[:-1])
        html += '  ' + str(haplotype[-1]) + '<br/>'
    html += '</pre>'
    display_html(html, raw=True)
{% endhighlight %}


{% highlight python %}
display_prefix_array(X)
{% endhighlight %}


<pre style="line-height: 100%">4|00000  0<br/>1|11000  1<br/>6|11000  1<br/>0|01010  1<br/>5|10001  0<br/>7|01011  0<br/>3|01111  0<br/>2|11111  1<br/></pre>


The positional prefix array [4, 1, 6, ...] is the column to the left of the '``|``' showing the haplotype indices. To the right of the '``|``' are the haplotypes up to *k*-1 sorted in reversed prefix order. Then separated by a space is the next allele for each haplotype - this final column makes up what Durbin (2014) calls *y<sup>k</sup>*[*k*].

## Algorithm 2. Build prefix and divergence arrays

Sorting haplotypes in this order has a very useful property, which is that each haplotype will be adjacent to the haplotype with the longest match. If several haplotypes have equally long matches these will all be adjacent. Also, the length of match between any pair of non-adjacent haplotypes is the minimum of the lengths of the matches between all haplotypes occurring in between them in the sorted order.

Algorithm 2 shows how to find and keep track of the start position for matches between adjacent haplotypes. This can be done while sweeping through the data building the positional prefix arrays. The reason this is possible is because once a match begins between a pair of adjacent haplotypes, those two haplotypes will remain adjacent until the match breaks.

Here's a pure Python implementation of algorithm 2, again purely for illustration. Durbin (2014) uses *d<sub>k</sub>* to denote the "divergence array" at position *k*, where *d<sub>k</sub>*[*i*] is the position where a match begins between the *i*th haplotype in the sorted order and it's predecessor. Because *d* is also used for one of the intermediates, for clarity I'll use *div<sub>k</sub>* to denote the divergence array.


{% highlight python %}
def build_prefix_and_divergence_arrays(X):
    
    # M haplotypes
    M = len(X)

    # N variable sites
    N = len(X[0])
    
    # initialise positional prefix array
    ppa = list(range(M))
    
    # initialise divergence array
    div = [0] * M
    
    # iterate over variants
    for k in range(N-1):

        # setup intermediates
        a = list()
        b = list()
        d = list()
        e = list()
        p = q = k + 1
    
        # iterate over haplotypes in reverse prefix sorted order
        for index, match_start in zip(ppa, div):

            # current haplotype
            haplotype = X[index]
            
            # allele for current haplotype
            allele = haplotype[k]
            
            # update intermediates
            if match_start > p:
                p = match_start
            if match_start > q:
                q = match_start

            # update intermediates
            if allele == 0:
                a.append(index)
                d.append(p)
                p = 0
            else:
                b.append(index)
                e.append(q)
                q = 0
                
        # construct the new arrays for k+1 by concatenating intermediates
        ppa = a + b
        div = d + e
        
    return ppa, div
    
{% endhighlight %}

Let's try it out on the dummy haplotype data.


{% highlight python %}
ppa, div = build_prefix_and_divergence_arrays(X)
{% endhighlight %}


{% highlight python %}
ppa
{% endhighlight %}




    [4, 1, 6, 0, 5, 7, 3, 2]




{% highlight python %}
div
{% endhighlight %}




    [5, 2, 0, 4, 5, 4, 3, 1]



Again to help make sense of this let's write a display function.


{% highlight python %}
def display_prefix_and_divergence_arrays(X):
    from IPython.display import display_html
    ppa, div = build_prefix_and_divergence_arrays(X)
    html = '<pre style="line-height: 100%">'
    for index, match_start in zip(ppa, div):
        haplotype = X[index]
        html += str(index) + '|'
        for k, allele in enumerate(haplotype[:-1]):
            if match_start == k:
                html += '<strong><u>'
            html += str(allele)
        if match_start < len(haplotype) - 1:
            html += '</u></strong>'
        html += '  ' + str(haplotype[-1]) + '<br/>'
    html += '</pre>'
    display_html(html, raw=True)
{% endhighlight %}


{% highlight python %}
display_prefix_and_divergence_arrays(X)
{% endhighlight %}


<pre style="line-height: 100%">4|00000  0<br/>1|11<strong><u>000</u></strong>  1<br/>6|<strong><u>11000</u></strong>  1<br/>0|0101<strong><u>0</u></strong>  1<br/>5|10001  0<br/>7|0101<strong><u>1</u></strong>  0<br/>3|011<strong><u>11</u></strong>  0<br/>2|1<strong><u>1111</u></strong>  1<br/></pre>


The divergence array has been used to show underlined in bold the maximal matches between each haplotype and it's predecessor in the sorted order, as in Durbin (2014) Figure 1. 

Let's try it out on some randomly generated data, just to have another example to look at.


{% highlight python %}
import random
M = 10
N = 30
Y = [[random.randint(0, 1) for _ in range(N)] for _ in range(M)]
display_prefix_and_divergence_arrays(Y)
{% endhighlight %}


<pre style="line-height: 100%">1|11001011100011111101111000000  0<br/>9|110011011011110110010000101<strong><u>00</u></strong>  0<br/>8|11110001100100100111000001<strong><u>100</u></strong>  1<br/>4|1001111111010100011011101<strong><u>1100</u></strong>  1<br/>5|1111100110101000101000001011<strong><u>0</u></strong>  1<br/>2|00101000111100110011101010001  0<br/>6|11010001100000001<strong><u>011101010001</u></strong>  0<br/>7|10010110111010100000000011<strong><u>001</u></strong>  1<br/>3|000100011111000110110101101<strong><u>01</u></strong>  1<br/>0|1101011111011010001101001111<strong><u>1</u></strong>  0<br/></pre>


## Algorithm 3. Report long matches

Algorithm 3 shows how to find all matches between haplotypes longer than some value *L*. The trick here is to notice that all haplotypes with matches longer than L will be adjacent in the sorted order, forming a "block". So algorithm 3 iterates over the haplotypes in sorted order, collecting haplotypes with matches longer than L. When it finds a match less than *L*, it must mean a break between blocks, so it reports matches for all haplotypes in the previous block. The algorithm also collects separately for haplotypes with different alleles at position *k*, to ensure that only matches that terminate at *k* are reported. 


{% highlight python %}
def report_long_matches(X, L):
    
    # M haplotypes
    M = len(X)

    # N variable sites
    N = len(X[0])
    
    # initialise positional prefix array
    ppa = list(range(M))
    
    # initialise divergence array
    div = [0] * M
    
    # iterate over variants
    for k in range(N):

        # setup intermediates
        a = list()
        b = list()
        d = list()
        e = list()
        p = q = k + 1
        ma = list()
        mb = list()
    
        # iterate over haplotypes in reverse prefix sorted order
        for index, match_start in zip(ppa, div):
            
            # report matches
            if match_start > k - L:
                if ma and mb:
                    yield k, ma, mb
                ma = list()
                mb = list()

            # current haplotype
            haplotype = X[index]
            
            # allele for current haplotype
            allele = haplotype[k]
            
            # update intermediates
            if match_start > p:
                p = match_start
            if match_start > q:
                q = match_start

            # update intermediates
            if allele == 0:
                a.append(index)
                d.append(p)
                p = 0
                ma.append(index)
            else:
                b.append(index)
                e.append(q)
                q = 0
                mb.append(index)
                
        # report any remaining matches including final haplotype (N.B., not in Durbin 2014)
        if ma and mb:
            yield k, ma, mb
                
        # construct the new arrays for k+1
        if k < N - 1:
            ppa = a + b
            div = d + e
    
{% endhighlight %}

One minor note, I don't think algorithm 3 as presented in Durbin (2014) accounts for the case where a block of matching haplotypes extends up to the final haplotype, so I added in an extra couple of lines to account for this case.

Let's try it out on the dummy data, finding matches 3 or longer.


{% highlight python %}
for match in report_long_matches(X, 3):
    print(match)
{% endhighlight %}

    (4, [4], [5])
    (4, [0], [7])
    (5, [4], [1, 6])
    (5, [3], [2])


So 4 matches have been found. The first match terminates at position 4, and is between haplotypes 4 and 5. Again let's write a display function to help visualise.


{% highlight python %}
def display_long_matches(X, L):
    from IPython.display import display_html
    html = '<pre style="line-height: 100%">'
    for match in report_long_matches(X, L):
        k, ma, mb = match
        for i in sorted(ma):
            for j in sorted(mb):
                html += 'match ending at position %s between haplotypes %s and %s:<br/><br/>' % (k, i, j)
                h1 = X[i][k-L:k+1]
                h2 = X[j][k-L:k+1]
                html += str(i) + '|' + ''.join(str(allele) for allele in h1[:-1]) 
                html += str(h1[-1]) + '<br/>'
                html += str(j) + '|<strong><u>' + ''.join(str(allele) for allele in h2[:-1]) + '</u></strong>'
                html += str(h2[-1]) + '<br/><br/>'
    html += '</pre>'
    display_html(html, raw=True)
{% endhighlight %}


{% highlight python %}
display_long_matches(X, 3)
{% endhighlight %}


<pre style="line-height: 100%">match ending at position 4 between haplotypes 4 and 5:<br/><br/>4|0000<br/>5|<strong><u>000</u></strong>1<br/><br/>match ending at position 4 between haplotypes 0 and 7:<br/><br/>0|1010<br/>7|<strong><u>101</u></strong>1<br/><br/>match ending at position 5 between haplotypes 4 and 1:<br/><br/>4|0000<br/>1|<strong><u>000</u></strong>1<br/><br/>match ending at position 5 between haplotypes 4 and 6:<br/><br/>4|0000<br/>6|<strong><u>000</u></strong>1<br/><br/>match ending at position 5 between haplotypes 3 and 2:<br/><br/>3|1110<br/>2|<strong><u>111</u></strong>1<br/><br/></pre>



{% highlight python %}
display_long_matches(X, 4)
{% endhighlight %}


<pre style="line-height: 100%">match ending at position 4 between haplotypes 0 and 7:<br/><br/>0|01010<br/>7|<strong><u>0101</u></strong>1<br/><br/>match ending at position 5 between haplotypes 3 and 2:<br/><br/>3|11110<br/>2|<strong><u>1111</u></strong>1<br/><br/></pre>



{% highlight python %}
display_long_matches(X, 5)
{% endhighlight %}


<pre style="line-height: 100%"></pre>


You may have noticed that the dummy data *X* contains two haplotypes that are completely identical over all 6 positions, so why aren't these found? That's because the algorithm requires matches to terminate.

## Algorithm 4. Report set maximal matches

Algorithm 4 shows how to find "set maximal matches" for each haplotype. These are the longest match for each haplotype at each position *k*. Again this makes use of the fact that haplotypes with maximal matches will be adjacent in the sorted order. At each position *k*, the algorithm iterates through the haplotypes in sorted order, looking for a maximal match for each haplotype. The first thing it does is to try and find any neighbours where the match can be extended, i.e., where the alleles at position *k* are the same. If so, continue on to the next haplotype without reporting anything. Otherwise, report a match with the longest matching neighbour (taking into account the fact that several neighbours may have equally long matches). 


{% highlight python %}
def report_set_maximal_matches(X):

    # M haplotypes
    M = len(X)

    # N variable sites
    N = len(X[0])
    
    # initialise positional prefix array
    ppa = list(range(M))
    
    # initialise divergence array
    div = [0] * M
    
    # iterate over variants
    for k in range(N):
        
        # sentinel
        div.append(k+1)
        
        for i in range(M):
            
            m = i - 1
            n = i + 1
            match_continues = False
            
            if div[i] <= div[i+1]:
                # match to previous neighbour is longer, scan "down" the haplotypes (decreasing indices)
                while div[m+1] <= div[i]:
                    if X[ppa[m]][k] == X[ppa[i]][k]:
                        match_continues = True
                        break
                    m -= 1
            if match_continues:
                continue
                    
            if div[i] >= div[i+1]:
                # match to next neighbour is longer, scan "up" the haplotypes (increasing indices)
                while div[n] <= div[i+1]:
                    if X[ppa[n]][k] == X[ppa[i]][k]:
                        match_continues = True
                        break
                    n += 1
            if match_continues:
                continue
                
            for j in range(m+1, i):
                if div[i] < k:  # N.B., avoid 0 length matches, not in Durbin (2014)
                    yield ppa[i], ppa[j], div[i], k
                
            for j in range(i+1, n):
                if div[i+1] < k:  # N.B., avoid 0 length matches, not in Durbin (2014)
                    yield ppa[i], ppa[j], div[i+1], k
                    
        # build next prefix and divergence arrays
        if k < N - 1:        
                
            # setup intermediates
            a = list()
            b = list()
            d = list()
            e = list()
            p = q = k + 1

            # iterate over haplotypes in prefix sorted order
            for index, match_start in zip(ppa, div):

                # current haplotype
                haplotype = X[index]

                # allele for current haplotype
                allele = haplotype[k]

                # update intermediates
                if match_start > p:
                    p = match_start
                if match_start > q:
                    q = match_start

                # update intermediates
                if allele == 0:
                    a.append(index)
                    d.append(p)
                    p = 0
                else:
                    b.append(index)
                    e.append(q)
                    q = 0

            # construct the new arrays for k+1
            ppa = a + b
            div = d + e

{% endhighlight %}


{% highlight python %}
for match in report_set_maximal_matches(X):
    print(match)
{% endhighlight %}

    (4, 0, 0, 1)
    (4, 3, 0, 1)
    (4, 7, 0, 1)
    (5, 1, 0, 1)
    (5, 2, 0, 1)
    (5, 6, 0, 1)
    (3, 0, 0, 2)
    (3, 7, 0, 2)
    (2, 1, 0, 2)
    (2, 6, 0, 2)
    (4, 5, 1, 4)
    (5, 4, 1, 4)
    (0, 7, 0, 4)
    (7, 0, 0, 4)
    (4, 1, 2, 5)
    (4, 6, 2, 5)
    (3, 2, 1, 5)
    (2, 3, 1, 5)



{% highlight python %}
def display_set_maximal_matches(X):
    from IPython.display import display_html
    from operator import itemgetter
    from itertools import groupby
    html = '<pre style="line-height: 100%">'
    matches = sorted(report_set_maximal_matches(X), key=itemgetter(0, 2))
    for i, sub_matches in groupby(matches, key=itemgetter(0)):
        html += 'set maximal matches for haplotype %s:<br/><br/>' % i
        hi = X[i]
        html += str(i) + '|' + ''.join(map(str, hi)) + '<br/>'
        for _, j, k1, k2 in sub_matches:
            hj = X[j]
            html += str(j) + '|' + ''.join(map(str, hj[:k1]))
            html += '<strong><u>' + ''.join(map(str, hj[k1:k2])) + '</u></strong>'
            html += ''.join(map(str, hj[k2:])) + '<br/>'
        html += '<br/>'
    html += '</pre>'
    display_html(html, raw=True)
{% endhighlight %}


{% highlight python %}
display_set_maximal_matches(X)
{% endhighlight %}


<pre style="line-height: 100%">set maximal matches for haplotype 0:<br/><br/>0|010101<br/>7|<strong><u>0101</u></strong>10<br/><br/>set maximal matches for haplotype 2:<br/><br/>2|111111<br/>1|<strong><u>11</u></strong>0001<br/>6|<strong><u>11</u></strong>0001<br/>3|0<strong><u>1111</u></strong>0<br/><br/>set maximal matches for haplotype 3:<br/><br/>3|011110<br/>0|<strong><u>01</u></strong>0101<br/>7|<strong><u>01</u></strong>0110<br/>2|1<strong><u>1111</u></strong>1<br/><br/>set maximal matches for haplotype 4:<br/><br/>4|000000<br/>0|<strong><u>0</u></strong>10101<br/>3|<strong><u>0</u></strong>11110<br/>7|<strong><u>0</u></strong>10110<br/>5|1<strong><u>000</u></strong>10<br/>1|11<strong><u>000</u></strong>1<br/>6|11<strong><u>000</u></strong>1<br/><br/>set maximal matches for haplotype 5:<br/><br/>5|100010<br/>1|<strong><u>1</u></strong>10001<br/>2|<strong><u>1</u></strong>11111<br/>6|<strong><u>1</u></strong>10001<br/>4|0<strong><u>000</u></strong>00<br/><br/>set maximal matches for haplotype 7:<br/><br/>7|010110<br/>0|<strong><u>0101</u></strong>01<br/><br/></pre>


Here the bold underline indicates the regions in the haplotype representing set maximal matches to the haplotype at the top.

Again you may notice that there are no matches reported for the two haplotypes that are identical across all 6 positions. Again this is because matches are required to terminate.

Let's have a look with some random data for further illustration.


{% highlight python %}
import random
M = 3
N = 30
Y = [[random.randint(0, 1) for _ in range(N)] for _ in range(M)]
display_set_maximal_matches(Y)
{% endhighlight %}


<pre style="line-height: 100%">set maximal matches for haplotype 0:<br/><br/>0|100001010111100110111010100101<br/>2|011<strong><u>00</u></strong>0110001000111111011110011<br/>1|01101<strong><u>1</u></strong>110111001111100110000010<br/>1|0110111<strong><u>10111</u></strong>001111100110000010<br/>2|0110001100010<strong><u>0011</u></strong>1111011110011<br/>2|011000110001000111<strong><u>11101</u></strong>1110011<br/>1|0110111101110011111001<strong><u>10</u></strong>000010<br/>2|011000110001000111111011<strong><u>1</u></strong>10011<br/>1|0110111101110011111001100<strong><u>00</u></strong>010<br/><br/>set maximal matches for haplotype 1:<br/><br/>1|011011110111001111100110000010<br/>2|<strong><u>0110</u></strong>00110001000111111011110011<br/>0|10000<strong><u>1</u></strong>010111100110111010100101<br/>2|011000<strong><u>110</u></strong>001000111111011110011<br/>0|1000010<strong><u>10111</u></strong>100110111010100101<br/>2|01100011000<strong><u>100</u></strong>0111111011110011<br/>2|011000110001000<strong><u>1111</u></strong>11011110011<br/>0|1000010101111001101110<strong><u>10</u></strong>100101<br/>0|1000010101111001101110101<strong><u>00</u></strong>101<br/>2|01100011000100011111101111<strong><u>001</u></strong>1<br/><br/>set maximal matches for haplotype 2:<br/><br/>2|011000110001000111111011110011<br/>1|<strong><u>0110</u></strong>11110111001111100110000010<br/>0|100<strong><u>00</u></strong>1010111100110111010100101<br/>1|011011<strong><u>110</u></strong>111001111100110000010<br/>1|01101111011<strong><u>100</u></strong>1111100110000010<br/>0|1000010101111<strong><u>0011</u></strong>0111010100101<br/>1|011011110111001<strong><u>1111</u></strong>00110000010<br/>0|100001010111100110<strong><u>11101</u></strong>0100101<br/>0|100001010111100110111010<strong><u>1</u></strong>00101<br/>1|01101111011100111110011000<strong><u>001</u></strong>0<br/><br/></pre>


## Algorithm 5. Set maximal matches from a new sequence *z* to *X*

I haven't got my head around this yet. Any help with the basic intuition would be very welcome.

## Compact representation of *X*

I think I partly get this, at least why *PBWT* should be more compressible than the original *X*. But I don't get how you recover *X* from *PBWT*. I need to go and read up about the FM-index.

Let's at least check the assertion that the *PBWT* is more compressible than *X*, using some real haplotype data from [Ag1000G](http://www.malariagen.net/ag1000g). First, I'm going to implement algorithm 2 in a slightly different way, using NumPy and Cython to speed things up, and also constructing the *PBWT* as output.  


{% highlight python %}
%load_ext Cython
{% endhighlight %}


{% highlight python %}
%%cython


import numpy as np
cimport numpy as np


def build_pbwt(np.int8_t[:, :] H):
    cdef:
        Py_ssize_t N, M, k, i, u, v, p, q
        np.uint32_t index, match_start
        np.int8_t allele
        np.int8_t[:, :] pbwt
        np.uint32_t[:, :] ppa, div
        np.uint32_t[:] a, b, d, e
    
    # expect haplotype data transposed
    N, M = H.shape[:2]
    
    # setup pbwt
    pbwt = np.empty_like(H)

    # setup positional prefix array
    ppa = np.empty((N, M), dtype='u4')
    
    # initialise first ppa column
    for i in range(M):
        ppa[0, i] = i
    
    # setup divergence array
    div = np.zeros((N, M), dtype='u4')
    
    # setup intermediates
    a = np.zeros(M, dtype='u4')
    b = np.zeros(M, dtype='u4')
    d = np.zeros(M, dtype='u4')
    e = np.zeros(M, dtype='u4')
    
    # iterate over variants
    for k in range(N):
        
        # setup intermediates
        u = v = 0
        p = q = k + 1
    
        # iterate over haplotypes in reverse prefix sorted order
        for i in range(M):

            # index for current haplotype
            index = ppa[k, i]
            
            # match start position for current haplotype
            match_start = div[k, i]
            
            # allele for current haplotype
            allele = H[k, index]
            
            # update output
            pbwt[k, i] = allele
            
            # update intermediates
            if match_start > p:
                p = match_start
            if match_start > q:
                q = match_start

            # update intermediates
            if allele == 0:
                a[u] = index
                d[u] = p
                u += 1
                p = 0
            else:
                b[v] = index
                e[v] = q
                v += 1
                q = 0
                
        if k < N - 1:
            
            # construct the new positional prefix array for k+1
            ppa[k+1, :u] = a[:u]
            ppa[k+1, u:] = b[:v]
        
            # construct the new divergence array for k+1
            div[k+1, :u] = d[:u]
            div[k+1, u:] = e[:v]
        
    return np.asarray(pbwt), np.asarray(ppa), np.asarray(div)

{% endhighlight %}

Import some extra libraries.


{% highlight python %}
import h5py
import allel
import zarr
zarr.__version__
{% endhighlight %}




    '1.0.0'



Set up some real haplotype data.


{% highlight python %}
callset_fn = '/data/coluzzi/ag1000g/data/phase1/release/AR3/haplotypes/main/hdf5/ag1000g.phase1.ar3.haplotypes.3R.h5'
callset = h5py.File(callset_fn, mode='r')
callset
{% endhighlight %}




    <HDF5 file "ag1000g.phase1.ar3.haplotypes.3R.h5" (mode r)>




{% highlight python %}
genotypes = allel.GenotypeChunkedArray(callset['3R/calldata/genotype'])
genotypes
{% endhighlight %}




<table class='petl'>
<caption>GenotypeChunkedArray((10178802, 773, 2), int8, nbytes=14.7G, cbytes=449.7M, cratio=33.4, cname=gzip, clevel=1, shuffle=True, chunks=(678, 773, 2), data=h5py._hl.dataset.Dataset)</caption>
<thead>
<tr>
<th></th>
<th>0</th>
<th>1</th>
<th>2</th>
<th>3</th>
<th>4</th>
<th>...</th>
<th>768</th>
<th>769</th>
<th>770</th>
<th>771</th>
<th>772</th>
</tr>
</thead>
<tbody>
<tr>
<td style='font-weight: bold'>0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>...</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
</tr>
<tr>
<td style='font-weight: bold'>1</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>...</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
</tr>
<tr>
<td style='font-weight: bold'>2</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>...</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
</tr>
<tr>
<td style='font-weight: bold'>3</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>...</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
</tr>
<tr>
<td style='font-weight: bold'>4</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>...</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
<td>0/0</td>
</tr>
</tbody>
</table>
<p><strong>...</strong></p>



For comparison with results from Durbin (2014) let's take 1000 haplotypes and use the same number of segregating sites.


{% highlight python %}
H = genotypes[:500000, :500].to_haplotypes()
ac = H.count_alleles()
H = H[ac.is_segregating()]
H = H[:370264]
H
{% endhighlight %}




<table class='petl'>
<caption>HaplotypeArray((370264, 1000), dtype=int8)</caption>
<thead>
<tr>
<th></th>
<th>0</th>
<th>1</th>
<th>2</th>
<th>3</th>
<th>4</th>
<th>...</th>
<th>995</th>
<th>996</th>
<th>997</th>
<th>998</th>
<th>999</th>
</tr>
</thead>
<tbody>
<tr>
<td style='font-weight: bold'>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>...</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
<tr>
<td style='font-weight: bold'>1</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>...</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
<tr>
<td style='font-weight: bold'>2</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>...</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
<tr>
<td style='font-weight: bold'>3</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>...</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
<tr>
<td style='font-weight: bold'>4</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>...</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
</tbody>
</table>
<p><strong>...</strong></p>



Before applying the PBWT transform, check how compressible these data are already, using a standard compression algorithm (LZ4).


{% highlight python %}
H_compressed_lz4 = zarr.array(H, chunks=(679, H.shape[1]), compression='blosc', 
                              compression_opts=dict(cname='lz4', clevel=1, shuffle=2))
H_compressed_lz4
{% endhighlight %}




    zarr.core.Array((370264, 1000), int8, chunks=(679, 1000), order=C)
      compression: blosc; compression_opts: {'clevel': 1, 'cname': 'lz4', 'shuffle': 2}
      nbytes: 353.1M; nbytes_stored: 10.7M; ratio: 33.1; initialized: 546/546
      store: builtins.dict



So as they are, these haplotype data can be compressed down to 10.7M.

Now let's built the *PBWT*.


{% highlight python %}
%%time
pbwt, ppa, div = build_pbwt(H)
{% endhighlight %}

    CPU times: user 3.69 s, sys: 248 ms, total: 3.94 s
    Wall time: 3.9 s



{% highlight python %}
pbwt
{% endhighlight %}




    array([[0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           ..., 
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0]], dtype=int8)




{% highlight python %}
pbwt_compressed_lz4 = zarr.array(pbwt, chunks=(679, H.shape[1]), compression='blosc', 
                                 compression_opts=dict(cname='lz4', clevel=1, shuffle=2))
pbwt_compressed_lz4
{% endhighlight %}




    zarr.core.Array((370264, 1000), int8, chunks=(679, 1000), order=C)
      compression: blosc; compression_opts: {'clevel': 1, 'cname': 'lz4', 'shuffle': 2}
      nbytes: 353.1M; nbytes_stored: 8.2M; ratio: 43.1; initialized: 546/546
      store: builtins.dict



So in this case, *PBWT* is more compressible than the original data, but only marginally, requiring 8.2M rather than 10.7M. This seems at odds with what is reported in Durbin (2014), i.e., that PBWT is many times more compressible than raw haplotypes, so what's going on? 

I think that a lot of the benefit in terms of compressibility that is reported in Durbin (2014) is actually due to the fact that the haplotype data in the PBWT are transposed relative to the original .gz representation. I.e., If haplotypes are organised one variant per row rather than one haplotype per row, and the underlying bytes are in row-major order, then the data become more compressible as you add more haplotypes, because variants tend to be rare and hence most rows will be composed almost entirely of zeros. If you then permute each row via the PBWT you get a further benefit, because PBWT tends to bring the ones and zeros together, however the added benefit is fairly marginal and most of the gains come simply from the transposed layout.

Note that I am using LZ4 here and not the run length encoding used in Durbin (2014), I don't know how much difference that would make.

To make use of the PBWT you also need the divergence array, so how compressible is that?


{% highlight python %}
div
{% endhighlight %}




    array([[     0,      0,      0, ...,      0,      0,      0],
           [     1,      0,      0, ...,      0,      0,      1],
           [     2,      0,      0, ...,      0,      1,      2],
           ..., 
           [370261, 369928, 369973, ..., 370259, 370260, 370261],
           [370262, 369928, 369973, ..., 370260, 370261, 370262],
           [370263, 369928, 369973, ..., 370261, 370262, 370263]], dtype=uint32)




{% highlight python %}
div_compressed_lz4 = zarr.array(div, chunks=(679, H.shape[1]), compression='blosc', 
                                compression_opts=dict(cname='lz4', clevel=1, shuffle=1))
div_compressed_lz4
{% endhighlight %}




    zarr.core.Array((370264, 1000), uint32, chunks=(679, 1000), order=C)
      compression: blosc; compression_opts: {'clevel': 1, 'cname': 'lz4', 'shuffle': 1}
      nbytes: 1.4G; nbytes_stored: 249.8M; ratio: 5.7; initialized: 546/546
      store: builtins.dict



The divergence array is not so compressible with this compression configuration. However Durbin (2014) suggests using a delta filter. Let's try this via the Python LZMA library.


{% highlight python %}
%%time
import lzma
filters = [dict(id=lzma.FILTER_DELTA, dist=4),
           dict(id=lzma.FILTER_LZMA2, preset=1)]
div_compressed_lzma = zarr.array(div, chunks=(679, H.shape[1]), compression='lzma', 
                                 compression_opts=dict(filters=filters))
{% endhighlight %}

    CPU times: user 29.3 s, sys: 12 ms, total: 29.3 s
    Wall time: 29.3 s



{% highlight python %}
div_compressed_lzma
{% endhighlight %}




    zarr.core.Array((370264, 1000), uint32, chunks=(679, 1000), order=C)
      compression: lzma; compression_opts: {'check': 0, 'filters': [{'id': 3, 'dist': 4}, {'id': 33, 'preset': 1}], 'format': 1, 'preset': None}
      nbytes: 1.4G; nbytes_stored: 16.0M; ratio: 88.5; initialized: 546/546
      store: builtins.dict



Unfortunately LZMA is painfully slow, however it does at least demonstrate that the divergence arrays are also highly compressible in principle, here taking 16.0M.

The raw .ipynb for this post is [here](https://github.com/alimanfoo/alimanfoo.github.io/tree/master/_posts).


