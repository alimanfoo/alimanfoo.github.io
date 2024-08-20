---
layout: post
title: Introducing scikit-allel
---


## Genome variation

My work involves the analysis of genetic variation from next-generation sequencing experiments. After running the sequence data through some variant calling pipeline like [GATK](https://www.broadinstitute.org/gatk/) I end up with a data set containing genotype calls for millions variants (SNPs, INDELS, etc.) on hundreds (soon to be thousands) of individuals. I study the malaria-carrying mosquito *Anopheles gambiae* but the principles are similar no matter what organism you're working on.

Once these data are in hand, there are many analyses that can be performed. These typically begin with quality-control (looking for poor-quality features in the data), but from there the possibilities are endless. Because DNA mutates and recombines as it passes from one generation to the next, every genome captures information about its own history. When you have data on many genomes from different populations, it's pretty amazing what you can begin to infer, from recent selection on insecticide resistance mutations to ancient population expansions and crashes.

## Exploratory analysis

Dealing with a data set of this richness, complexity and scale, you have to explore the data. When you're doing exploratory analyses, you don't start out with a concrete analysis plan, rather the plan emerges and evolves as you discover different features of the data. It's a very fluid process, and there are several things that can get in the way. 

First, things need to run fast. Typically this means seconds, minutes at a push. That means you can come up with an idea, run it, visualise the result, realise you did something stupid, re-run it a different way, visualise again, etc. If you're always having to go off for a cup of tea while something runs, it breaks the flow and you don't get very far.

Second, you want to avoid context-switching. That means you want to stay within the same computing environment, working with the same programming language, as much as possible. If you are always having to switch from one tool to another, writing out and parsing different text files, life becomes painful. You spend your time debugging weird errors and waiting for I/O.

## Scientific Python

I used to be a software engineer, and I did a lot of work in Python. When I started doing data analysis, naturally I stuck with what I knew. I was very happy to discover the rich ecosystem of general-purpose scientific software libraries available for Python.

At the heart of scientific computing in Python is [NumPy](http://www.numpy.org). NumPy is a library providing array-based numerical computing. Basically, this means you can do fast numerical computing whilst writing concise, readable code. 

As I started to work with large-scale genome variation data, I also started building some tools to work with data, based on NumPy. Partly this was out of necessity, partly frustration, partly as a learning exercise (you never really understand something until you've tried to code it yourself). Recently I decided to follow the example of many others developing NumPy-based programming libraries for specific scientific domains and package my tools up as a [scikit](https://scikits.appspot.com/scikits).

## scikit-allel

[`scikit-allel`](http://scikit-allel.readthedocs.org/en/latest/) is a Python library for exploratory analysis of large-scale genome variation data. It's still in an early stage of development, but there are some useful features, and I hope it will make life easier for others doing similar work to me.

Here's a quick illustration of a couple of features.


{% highlight python %}
import numpy as np
import allel; print('scikit-allel', allel.__version__)
import matplotlib.pyplot as plt
import zarr
import urllib
import shutil
from pathlib import Path
import zipfile
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
%matplotlib inline
{% endhighlight %}

    scikit-allel 1.2.0


I'm going to use data from the [Ag1000G phase 1 AR3 release](http://www.malariagen.net/data/ag1000g-phase1-AR3). Let's download some data from the [Ag1000G public FTP site](ftp://ngs.sanger.ac.uk/production/ag1000g/). The files we're downloading are ~200-400 Mb so this may take a little while, depending on your internet connection.


{% highlight python %}
def download(source_url, dest_path):
    """Helper function to download a remote file to the local file system."""
    with urllib.request.urlopen(source_url) as source:
        with open(dest_path, mode='wb') as dest:
            shutil.copyfileobj(source, dest)

            
# base URL for files to be downloaded
base_url = 'ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/zarr/'
# download files to a local directory called "data"
dest_dir = Path('data')
dest_dir.mkdir(exist_ok=True)

# files to be downloaded
files = [
    'ag1000g.phase1.ar3.pass.biallelic.metadata.zip',
    'ag1000g.phase1.ar3.pass.biallelic.2L.variants.zip',
    'ag1000g.phase1.ar3.pass.biallelic.2L.calldata.GT.zip',
]

for f in files:
    print('downloading', f)
    dest_path = dest_dir / f
    if not dest_path.exists():
        source_url = base_url + f
        download(source_url, dest_path)
    print('extracting', f)
    with zipfile.ZipFile(dest_path, mode='r') as z:
        z.extractall(dest_dir)

{% endhighlight %}

    downloading ag1000g.phase1.ar3.pass.biallelic.metadata.zip
    extracting ag1000g.phase1.ar3.pass.biallelic.metadata.zip
    downloading ag1000g.phase1.ar3.pass.biallelic.2L.variants.zip
    extracting ag1000g.phase1.ar3.pass.biallelic.2L.variants.zip
    downloading ag1000g.phase1.ar3.pass.biallelic.2L.calldata.GT.zip
    extracting ag1000g.phase1.ar3.pass.biallelic.2L.calldata.GT.zip


The data are stored in [Zarr files](https://zarr.readthedocs.io). Let's take a look.


{% highlight python %}
callset = zarr.open_consolidated('data/ag1000g.phase1.ar3.pass.biallelic')
callset
{% endhighlight %}




    <zarr.hierarchy.Group '/'>




{% highlight python %}
# pick a chromosome to work with
chrom = '2L'
{% endhighlight %}


{% highlight python %}
# access variants
variants = callset[chrom]['variants']
variants
{% endhighlight %}




    <zarr.hierarchy.Group '/2L/variants'>



Every genetic variant (in this case they are all SNPs) has a position on the genome. Working with these positions is a very common task.


{% highlight python %}
pos = allel.SortedIndex(variants['POS'])
pos
{% endhighlight %}




<div class="allel allel-DisplayAs1D"><span>&lt;SortedIndex shape=(8296600,) dtype=int32&gt;</span><table><thead><tr><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">8296595</th><th style="text-align: center">8296596</th><th style="text-align: center">8296597</th><th style="text-align: center">8296598</th><th style="text-align: center">8296599</th></tr></thead><tbody><tr><td style="text-align: center">44688</td><td style="text-align: center">44691</td><td style="text-align: center">44732</td><td style="text-align: center">44736</td><td style="text-align: center">44756</td><td style="text-align: center">...</td><td style="text-align: center">49356424</td><td style="text-align: center">49356425</td><td style="text-align: center">49356426</td><td style="text-align: center">49356429</td><td style="text-align: center">49356435</td></tr></tbody></table></div>



Let's plot the density of these variants over the chromosome.


{% highlight python %}
bin_width = 100000
bins = np.arange(0, pos.max(), bin_width)
# set X coordinate as bin midpoints
x = (bins[1:] + bins[:-1])/2
# compute variant density
h, _ = np.histogram(pos, bins=bins)
y = h / bin_width
# plot
fig, ax = plt.subplots(figsize=(9, 2))
sns.despine(ax=ax, offset=5)
ax.plot(x, y)
ax.set_xlabel('Position (bp)')
ax.set_ylabel('Density (bp$^{-1}$)')
ax.set_title('Variant density');
{% endhighlight %}


![png](/assets/2015-09-15-introducing-scikit-allel_files/2015-09-15-introducing-scikit-allel_11_0.png)


Let's say I have a gene of interest. I know what position it starts and ends, and I want to find variants within the gene.


{% highlight python %}
start, stop = 2358158, 2431617
loc = pos.locate_range(start, stop)
loc
{% endhighlight %}




    slice(24471, 26181, None)



I can use this slice to load genotype data for the region of interest.


{% highlight python %}
calldata = callset[chrom]['calldata']
g = allel.GenotypeArray(calldata['GT'][loc])
g
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(1710, 765, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">760</th><th style="text-align: center">761</th><th style="text-align: center">762</th><th style="text-align: center">763</th><th style="text-align: center">764</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1707</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1708</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1709</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></tbody></table></div>



There are various manipulations that can be done on a genotype array, e.g., convert to the number of alternate alleles per call.


{% highlight python %}
gn = g.to_n_alt()
gn
{% endhighlight %}




    array([[0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           ...,
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0]], dtype=int8)



From there we could compute genetic distance between each pair of individuals...


{% highlight python %}
dist = allel.pairwise_distance(gn, metric='euclidean')
dist
{% endhighlight %}




    array([12.4498996 ,  1.41421356,  2.        , ...,  2.        ,
            2.        ,  0.        ])



...which we could quickly visualise...


{% highlight python %}
allel.plot_pairwise_distance(dist);
{% endhighlight %}


![png](/assets/2015-09-15-introducing-scikit-allel_files/2015-09-15-introducing-scikit-allel_21_0.png)


## Further reading

I will leave it there for now, but check out the [scikit-allel docs](http://scikit-allel.readthedocs.org/en/latest/index.html) for more information. There is a section on [data structures](http://scikit-allel.readthedocs.org/en/latest/model.html), which includes both contiguous in-memory and compressed data structures for dealing with large arrays. There is also a [statistics](http://scikit-allel.readthedocs.org/en/latest/stats.html) section with various functions for computing diversity, Fst, LD, running PCA, and doing admixture tests, as well as a few useful plotting functions.

It is just a beginning, but hopefully a step in a good direction.

<hr/>


{% highlight python %}
import datetime
print(datetime.datetime.now().isoformat())
{% endhighlight %}

    2019-06-13T22:59:13.953334

