---
layout: post
title: Fast PCA
---


Principal components analysis (PCA) is a mainstay of population genetics, providing a model-free method for exploring patterns of relatedness within a collection of individuals. PCA was introduced as a tool for genetic genetic analysis by [Patterson, Price & Reich (2006)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190). Subsequently [Gil McVean (2009)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000686) provided an analytical framework for understanding PCA in terms of genetic ancestry. However, although PCA is widely used and the analytical details are worked out, there are a number of practical issues that come up when trying to run PCA on large SNP datasets from next-generation sequencing experiments. For example, small changes in how you prepare the input data can make a big difference to the outputs. The [Ag1000G phase 1 data](http://www.malariagen.net/data/ag1000g-phase1-ar3) provide a concrete illustration of some of these issues, so I thought I'd try to bring together some experiences here.

Also, while PCA is fairly quick to run on smaller datasets, it can become slow and memory-intensive with larger data. A few months ago I discovered that [scikit-learn](http://scikit-learn.org/stable/) includes a [randomized SVD](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.RandomizedPCA.html) implementation, which is a potentially faster and more scalable method for approximating the top N components than using a conventional singular value decomposition. To evaluate randomized PCA I implemented some functions in [scikit-allel](http://scikit-allel.readthedocs.org) which provide a convenience layer between underlying SVD implementations in NumPy and scikit-learn and the typical data structures I used to store genotype data. I know others have also started working with randomized PCA for genotype data ([Galinsky et al. 2015](http://biorxiv.org/content/early/2015/04/16/018143)) so I thought it would be interesting to apply both conventional and randomized SVD implementations to a non-human dataset and report some performance data.

## Setup


{% highlight python %}
import random
random.seed(42)
import time
import numpy as np
np.random.seed(42)
import h5py
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas
import allel; print('scikit-allel', allel.__version__)
%reload_ext memory_profiler
{% endhighlight %}

    scikit-allel 1.0.3


I have a copy of the [Ag1000G phase 1 AR3 data release](http://www.malariagen.net/data/ag1000g-phase1-ar3) on a local drive. The SNP genotype data is available in an HDF5 file.


{% highlight python %}
callset_fn = 'data/2015-09-28/ag1000g.phase1.ar3.pass.h5'
callset = h5py.File(callset_fn, mode='r')
callset
{% endhighlight %}




    <HDF5 file "ag1000g.phase1.ar3.pass.h5" (mode r)>



Let's work with chromosome arm 3L.


{% highlight python %}
chrom = '3L'
{% endhighlight %}

Setup the genotype data.


{% highlight python %}
g = allel.GenotypeChunkedArray(callset[chrom]['calldata']['genotype'])
g
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeChunkedArray shape=(9643193, 765, 2) dtype=int8 chunks=(6553, 10, 2)
   nbytes=13.7G cbytes=548.0M cratio=25.7
   compression=gzip compression_opts=3
   values=h5py._hl.dataset.Dataset&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">760</th><th style="text-align: center">761</th><th style="text-align: center">762</th><th style="text-align: center">763</th><th style="text-align: center">764</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">9643190</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">9643191</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">9643192</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></table></div>



Count alleles at each variant. This takes a minute or so.


{% highlight python %}
# the '[:]' syntax pulls the data from compressed storage into a numpy array
ac = g.count_alleles()[:]
ac
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;AlleleCountsArray shape=(9643193, 4) dtype=int32&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">1527</td><td style="text-align: center">   3</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">1529</td><td style="text-align: center">   1</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">1528</td><td style="text-align: center">   2</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="5">...</td></tr><tr><th style="text-align: center">9643190</th><td style="text-align: center">1512</td><td style="text-align: center">  16</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center">9643191</th><td style="text-align: center">1527</td><td style="text-align: center">   1</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center">9643192</th><td style="text-align: center">1507</td><td style="text-align: center">  18</td><td style="text-align: center">   1</td><td style="text-align: center">   0</td></tr></table></div>



Before going any further, I'm going to remove singletons and multiallelic SNPs. Singletons are not informative for PCA, and the analysis is simpler if we restrict to biallelic SNPs.

For interest, how many multiallelic SNPs are there?


{% highlight python %}
np.count_nonzero(ac.max_allele() > 1)
{% endhighlight %}




    2193707



How many biallelic singletons?


{% highlight python %}
np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1))
{% endhighlight %}




    2622060



Apply the filtering.


{% highlight python %}
flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
gf = g.compress(flt, axis=0)
gf
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeChunkedArray shape=(4825329, 765, 2) dtype=int8 chunks=(590, 765, 2)
   nbytes=6.9G cbytes=442.8M cratio=15.9
   compression=blosc compression_opts={'shuffle': 1, 'cname': 'lz4', 'clevel': 5}
   values=zarr.core.Array&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">760</th><th style="text-align: center">761</th><th style="text-align: center">762</th><th style="text-align: center">763</th><th style="text-align: center">764</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">4825326</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">4825327</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">4825328</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></table></div>



Finally, transform the genotype data into a 2-dimensional matrix where each cell has the number of non-reference alleles per call. This is what we'll use as the input to PCA.


{% highlight python %}
gn = gf.to_n_alt()
gn
{% endhighlight %}




    <ChunkedArrayWrapper shape=(4825329, 765) dtype=int8 chunks=(1179, 765)
       nbytes=3.4G cbytes=331.4M cratio=10.6
       compression=blosc compression_opts={'shuffle': 1, 'cname': 'lz4', 'clevel': 5}
       values=zarr.core.Array>



Note that we are still working with reasonably large amounts of data here, and so we are using chunked compressed arrays for storage. The `gn` variable above uses a [Zarr](http://zarr.readthedocs.io) array to store the data, one of several chunked storage containers that can be used with scikit-allel.

## Removing correlated features (LD pruning)

As I understand it, PCA works best when the features you provide as input are independent from each other. Here each SNP is a feature, however, because DNA is transmitted from one generation to the next with some recombination between parents, genotypes at nearby SNPs tend to be correlated, with the correlation (linkage disequlibrium) decaying as you increase the separation between SNPs.

We can get a sense of that correlation structure by visualising pairwise linkage disequilibrium in the first 1000 SNPs.


{% highlight python %}
def plot_ld(gn, title):
    m = allel.stats.rogers_huff_r(gn) ** 2
    ax = allel.plot.pairwise_ld(m)
    ax.set_title(title)
{% endhighlight %}


{% highlight python %}
plot_ld(gn[:1000], 'Figure 1. Pairwise LD.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_23_0.png)


The darker regions in the plot above indicate pairs of SNPs where genotypes are correlated.

Before I deal with this correlation directly, I'm going to thin down the data a bit. There are 4,825,329 SNPs left after the initial filtering steps above, and analysing this many features would be slow. Here we are more concerned with running an exploratory analysis, so I'm going to randomly choose a subset of these SNPs to work with. This should still reveal the main signals in the data, while making runtime faster. 


{% highlight python %}
n = 100000  # number of SNPs to choose randomly
vidx = np.random.choice(gn.shape[0], n, replace=False)
vidx.sort()
gnr = gn.take(vidx, axis=0)
gnr
{% endhighlight %}




    <ChunkedArrayWrapper shape=(100000, 765) dtype=int8 chunks=(391, 765)
       nbytes=73.0M cbytes=7.1M cratio=10.3
       compression=blosc compression_opts={'shuffle': 1, 'cname': 'lz4', 'clevel': 5}
       values=zarr.core.Array>



By randomly downsampling SNPs, this should have dealt with much of the correlation between nearby features. Let's take a look at the first 1000.


{% highlight python %}
plot_ld(gnr[:1000], 'Figure 2. Pairwise LD after random downsampling.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_27_0.png)


You can see that much of the correlation is gone. However, depending how dusty your screen is, you may be able to see some speckling, indicating that there are still some correlated SNPs in the dataset.

To remove this remaining correlation, I'm going to explicitly locate SNPs that are not correlated with each other, using the `locate_unlinked` function from scikit-allel. This is known as LD pruning, and works by sliding a window along the data, computing pairwise LD between all SNPs within each window, then removing one SNP from each correlated pair.

Conventionally, LD pruning is run just once, however I'm going to run several iterations. In some cases this may make a difference to the results, in others it may not, probably depending on how much long-range LD is present in your samples. Running multiple iterations does slow things down a bit, but it's interesting to demonstrate and see what the effect is.


{% highlight python %}
def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn
{% endhighlight %}


{% highlight python %}
gnu = ld_prune(gnr, size=500, step=200, threshold=.1, n_iter=5)
{% endhighlight %}

    iteration 1 retaining 56651 removing 43349 variants
    iteration 2 retaining 47465 removing 9186 variants
    iteration 3 retaining 44583 removing 2882 variants
    iteration 4 retaining 43266 removing 1317 variants
    iteration 5 retaining 42509 removing 757 variants


5 iterations is probably more than necessary for this dataset, as you can see not many SNPs are removed after the first few iterations.

I've used a sliding window size of 500 SNPs here, which is larger than others typically use. Out of interest, how many SNPs would be removed if we used a smaller window and just one iteration?


{% highlight python %}
ld_prune(gnr, size=100, step=20, threshold=.1, n_iter=1);
{% endhighlight %}

    iteration 1 retaining 75762 removing 24238 variants


So with this dataset, using a larger window and multiple iterations finds and removes a lot more correlated SNPs. This is probably related to the fact that there are a lot of rare variants in the data, and so a larger window is required to find variants in linkage.

Let's take a look at how much LD is left after LD pruning.


{% highlight python %}
plot_ld(gnu[:1000], 'Figure 3. Pairwise LD after LD pruning.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_34_0.png)


The data are relatively small now after downsampling and LD-pruning, so we can bring the data out of chunked storage and into memory uncompressed, which is necessary for PCA.


{% highlight python %}
gnu = gnu[:]
gnu
{% endhighlight %}




    array([[0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0],
           ..., 
           [0, 0, 0, ..., 0, 0, 0],
           [0, 0, 1, ..., 0, 0, 0],
           [0, 0, 0, ..., 0, 0, 0]], dtype=int8)



## PCA via conventional SVD

Let's run a conventional PCA analysis of the LD-pruned genotype data.


{% highlight python %}
coords1, model1 = allel.pca(gnu, n_components=10, scaler='patterson')
{% endhighlight %}

To help visualise the results, I need to pull in some metadata about which population each individual mosquito belongs to.


{% highlight python %}
df_samples = pandas.read_csv('data/2015-09-28/samples.meta.txt', delimiter='\t', index_col='index')
df_samples.head()
{% endhighlight %}




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ox_code</th>
      <th>src_code</th>
      <th>sra_sample_accession</th>
      <th>population</th>
      <th>country</th>
      <th>region</th>
      <th>contributor</th>
      <th>contact</th>
      <th>year</th>
      <th>m_s</th>
      <th>sex</th>
      <th>n_sequences</th>
      <th>mean_coverage</th>
      <th>latitude</th>
      <th>longitude</th>
    </tr>
    <tr>
      <th>index</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AB0085-C</td>
      <td>BF2-4</td>
      <td>ERS223996</td>
      <td>BFS</td>
      <td>Burkina Faso</td>
      <td>Pala</td>
      <td>Austin Burt</td>
      <td>Sam O'Loughlin</td>
      <td>2012</td>
      <td>S</td>
      <td>F</td>
      <td>89905852</td>
      <td>28.01</td>
      <td>11.150</td>
      <td>-4.235</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AB0087-C</td>
      <td>BF3-3</td>
      <td>ERS224013</td>
      <td>BFM</td>
      <td>Burkina Faso</td>
      <td>Bana</td>
      <td>Austin Burt</td>
      <td>Sam O'Loughlin</td>
      <td>2012</td>
      <td>M</td>
      <td>F</td>
      <td>116706234</td>
      <td>36.76</td>
      <td>11.233</td>
      <td>-4.472</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AB0088-C</td>
      <td>BF3-5</td>
      <td>ERS223991</td>
      <td>BFM</td>
      <td>Burkina Faso</td>
      <td>Bana</td>
      <td>Austin Burt</td>
      <td>Sam O'Loughlin</td>
      <td>2012</td>
      <td>M</td>
      <td>F</td>
      <td>112090460</td>
      <td>23.30</td>
      <td>11.233</td>
      <td>-4.472</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AB0089-C</td>
      <td>BF3-8</td>
      <td>ERS224031</td>
      <td>BFM</td>
      <td>Burkina Faso</td>
      <td>Bana</td>
      <td>Austin Burt</td>
      <td>Sam O'Loughlin</td>
      <td>2012</td>
      <td>M</td>
      <td>F</td>
      <td>145350454</td>
      <td>41.36</td>
      <td>11.233</td>
      <td>-4.472</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AB0090-C</td>
      <td>BF3-10</td>
      <td>ERS223936</td>
      <td>BFM</td>
      <td>Burkina Faso</td>
      <td>Bana</td>
      <td>Austin Burt</td>
      <td>Sam O'Loughlin</td>
      <td>2012</td>
      <td>M</td>
      <td>F</td>
      <td>105012254</td>
      <td>34.64</td>
      <td>11.233</td>
      <td>-4.472</td>
    </tr>
  </tbody>
</table>
</div>




{% highlight python %}
populations = df_samples.population.unique()
populations
{% endhighlight %}




    array(['BFS', 'BFM', 'UGS', 'GWA', 'KES', 'CMS', 'AOM', 'GAS', 'GNS'], dtype=object)




{% highlight python %}
pop_colours = {
    'BFM': '#FF0000',
    'GAS': '#008000',
    'GNS': '#00FFFF',
    'UGS': '#90EE90',
    'GWA': '#FFA500',
    'AOM': '#8B0000',
    'BFS': '#1E90FF',
    'KES': '#808080',
    'CMS': '#0000FF',
}
{% endhighlight %}


{% highlight python %}
def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))
    

def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = df_samples.population.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    
{% endhighlight %}


{% highlight python %}
fig_pca(coords1, model1, 'Figure 4. Conventional PCA.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_45_0.png)


Looking at the left-hand plot of PC1 versus PC2, there is a clear separation of individuals into 6 different clusters. This indicates there are at least 6 genetically distinct populations represented by the mosquitoes we've sequenced. The plot of PC3 vs PC4 gives us additional evidence that certain populations (GAS and KES) are genetically distinct from each other and the rest, but doesn't reveal any new clusters.

## Effect of LD pruning

What would happen if we ran PCA on the data **without** removing correlated SNPs?


{% highlight python %}
coords2, model2 = allel.pca(gnr, n_components=10, scaler='patterson')
fig_pca(coords2, model2, 'Figure 5. Conventional PCA without LD pruning.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_49_0.png)


Although all of the same population sub-divisions are visible in the first four components, they are resolved in a very different way. The first two components are now driven strongly by two populations, Angola (AOM) and Kenya (KES), and further population structure is not clearly resolved until PC3 and PC4. 

It is interesting to note that the Kenyan and Angolan populations are the two populations with the lowest heterozygosity. In particular, almost all Kenyan samples have very long runs of homozygosity, suggesting a recent population crash. I would hazard a guess that, in particular for Kenya, there is long-range LD which is affecting the PCA. When we used the aggressively LD-pruned data in Figure 4 above, this effect is reduced.

## Effect of scaling

Patterson et al. (2006) proposed scaling the data to unit variance at each SNP, assuming that the alleles are approximately binomially distributed. McVean (2009) remarks that scaling the data in this way should have little effect, although it will upweight rare variants (i.e., SNPs where the minor allele is at low frequency in the dataset). Let's return to using the LD pruned data, and see what happens if we **don't** use Patterson's scaling method but instead just centre the data.


{% highlight python %}
coords3, model3 = allel.pca(gnu, n_components=10, scaler=None)
fig_pca(coords3, model3, 'Figure 6. Conventional PCA without variance scaling.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_53_0.png)


Here again the same clusters are visible but are resolved in a different way. Also, note more of the total variance is explained by the first four components than when using the Patterson scaler. As McVean (2009) suggests, I would guess that these effects are both due to the weighting of rare variants. When rare variants are upweighted, this resolves more clearly any subtle population structure in the data. However, there a lot of rare variants in this dataset, and so the total amount of variance explained by the first few components goes down.

## Effect of unequal sample sizes

McVean (2009) provides a very elegant demonstration of what can happen if different populations are not equally represented in your dataset. If there are many more samples from one particular population, this has the effect of warping the principal components. 

In Ag1000G phase one there are a lot more samples from Cameroon (CMS) than any of the other locations.


{% highlight python %}
df_samples.groupby('population').population.count()
{% endhighlight %}




    population
    AOM     60
    BFM     69
    BFS     81
    CMS    275
    GAS     56
    GNS     31
    GWA     46
    KES     44
    UGS    103
    Name: population, dtype: int64



What would happen if we randomly pick a subset of CMS samples, to achieve a more even representation? 


{% highlight python %}
sidx_cms = df_samples[df_samples.population == 'CMS'].index
sidx_other = df_samples[df_samples.population != 'CMS'].index
sidx = sorted(list(sidx_other) + random.sample(list(sidx_cms), 50))
len(sidx)
{% endhighlight %}




    540




{% highlight python %}
df_samples.take(sidx).groupby('population').population.count()
{% endhighlight %}




    population
    AOM     60
    BFM     69
    BFS     81
    CMS     50
    GAS     56
    GNS     31
    GWA     46
    KES     44
    UGS    103
    Name: population, dtype: int64




{% highlight python %}
gnus = gnu.take(sidx, axis=1)
# also remove any non-segregating variants after removing samples
gnus = gnus.compress(np.any(gnus > 0, axis=1) & np.any(gnus < 2, axis=1), axis=0)
{% endhighlight %}


{% highlight python %}
coords4, model4 = allel.pca(gnus, n_components=10, scaler='patterson')
fig_pca(coords4, model4, 'Figure 7. Conventional PCA, CMS downsampled to 50.', 
        sample_population=df_samples.take(sidx).population.values)
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_62_0.png)


Now the results are similar to the original PCA we plotted in Figure 4, however PC2 appears to be more balanced. So sample size clearly does matter. However, there is a chicken-and-egg problem here. If you are using PCA to discover population structure in some collection of individuals, you won't know *a priori* if any particular population is overrepresented. Perhaps in that situation, an initial round of PCA to discover population structure can be followed up with a second round, downsampling any populations within which you observe no differentiation.

## Randomized PCA

Randomized PCA is an alternative to conventional PCA. I don't claim to understand the details, but apparently it uses an approximation to estimate the top principal components only, rather than evaluating all principal components as in a conventional SVD. So it should be faster and use less memory. 

Let's run a randomized PCA on the Ag1000G data and compare the results with the conventional PCA.


{% highlight python %}
coords5, model5 = allel.randomized_pca(gnu, n_components=10, scaler='patterson')
fig_pca(coords5, model5, 'Figure 8. Randomized PCA.')
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_66_0.png)


For the first four components at least, the results are indistinguishable from the conventional PCA.

Let's compare performance.


{% highlight python %}
n_variants = np.arange(5000, gnu.shape[0], 5000)

pca_time_v = []
for n in n_variants:
    gx = gnu[:n]
    t1 = time.time()
    allel.stats.pca(gx, n_components=10, scaler='patterson')
    t2 = time.time()
    dur = t2-t1
    pca_time_v.append(dur)

random_pca_time_v = []
for n in n_variants:
    gx = gnu[:n]
    t1 = time.time()
    allel.stats.randomized_pca(gx, n_components=10, scaler='patterson')
    t2 = time.time()
    dur = t2-t1
    random_pca_time_v.append(dur)
    
{% endhighlight %}


{% highlight python %}
fig, ax = plt.subplots()
sns.despine(ax=ax, offset=5)
ax.plot(n_variants, pca_time_v, label='Conventional PCA', marker='o')
ax.plot(n_variants, random_pca_time_v, label='Randomized PCA', marker='o')
ax.set_xlabel('no. variants')
ax.set_ylabel('runtime (seconds)')
ax.legend(bbox_to_anchor=[1, 1], loc='upper left')
ax.set_title('Figure 9. PCA performance with number of variants');
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_69_0.png)



{% highlight python %}
n_samples = np.arange(50, gnu.shape[1], 50)

pca_time_s = []
for n in n_samples:
    gx = gnu[:, :n]
    gx = gx.compress(np.any(gx > 0, axis=1) & np.any(gx < 2, axis=1), axis=0)
    t1 = time.time()
    allel.pca(gx, n_components=10, scaler='patterson')
    t2 = time.time()
    dur = t2-t1
    pca_time_s.append(dur)

random_pca_time_s = []
for n in n_samples:
    gx = gnu[:, :n]
    gx = gx.compress(np.any(gx > 0, axis=1) & np.any(gx < 2, axis=1), axis=0)
    t1 = time.time()
    allel.randomized_pca(gx, n_components=10, scaler='patterson')
    t2 = time.time()
    dur = t2-t1
    random_pca_time_s.append(dur)
    
{% endhighlight %}


{% highlight python %}
fig, ax = plt.subplots()
sns.despine(ax=ax, offset=5)
ax.plot(n_samples, pca_time_s, label='Conventional PCA', marker='o')
ax.plot(n_samples, random_pca_time_s, label='Randomized PCA', marker='o')
ax.set_xlabel('no. samples')
ax.set_ylabel('runtime (seconds)')
ax.legend(bbox_to_anchor=[1, 1], loc='upper left')
ax.set_title('Figure 10. PCA performance with number of samples');
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_71_0.png)


What about memory usage?


{% highlight python %}
%memit allel.pca(gnu, n_components=10, scaler='patterson');
{% endhighlight %}

    peak memory: 3397.09 MiB, increment: 589.41 MiB



{% highlight python %}
%memit allel.randomized_pca(gnu, n_components=10, scaler='patterson');
{% endhighlight %}

    peak memory: 3148.99 MiB, increment: 310.09 MiB


So the randomized PCA is faster, scales better with more samples, and uses around half of the memory required by conventional SVD.

## Conclusions

* LD pruning makes a difference. For NGS data, a larger window size and/or multiple rounds of pruning may be required to deal with regions of long-range LD. LD pruning may also impact different populations in different ways, if populations have different levels of LD.
* Scaling input data to unit variance using the method of Patterson et al. (2006) makes a small but noticeable difference, increasing the ability to resolve distinct populations within higher principal components.
* Unequal sample sizes warps the principal components, as predicted by McVean (2009).
* Randomized PCA produces results that are almost indistinguishable from conventional PCA, while running faster and using less memory. However, preparing the data (LD pruning) can also take a long time, so it would be good to find a way to optimise that step too.

## Further reading

* Patterson, N., Price, A. L., & Reich, D. (2006). [Population structure and eigenanalysis](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190). PLoS Genetics, 2(12), 2074–2093.
* McVean, G. (2009). [A genealogical interpretation of principal components analysis](http://doi.org/10.1371/journal.pgen.1000686). PLoS Genetics, 5(10), e1000686. 
* Galinsky, K. J., Bhatia, G., Loh, P.-R., Georgiev, S., Mukherjee, S., Patterson, N. J., & Price, A. L. (2015). [Fast principal components analysis reveals independent evolution of ADH1B gene in Europe and East Asia](http://biorxiv.org/content/early/2015/08/24/018143.abstract). bioRxiv.
* [`scikit-allel` PCA functions](http://scikit-allel.readthedocs.org/en/latest/stats/decomposition.html)

<hr/>

## Post-script: Randomized PCA and lower components

When I originally wrote this post, I only looked at the first four components. These account for most of the variance and so capture the major signals of population subdivision. 


{% highlight python %}
fig, ax = plt.subplots()
sns.despine(ax=ax)
x = np.arange(10)
y = model1.explained_variance_ratio_ * 100
ax.bar(x+.6, y, width=.8)
ax.set_xticks(x+1)
ax.set_xlim(0, 11)
ax.set_xlabel('component')
ax.set_ylabel('% variance explained');
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_82_0.png)


For the first four components, conventional and randomized PCA are basically the same. However, recently I looked into the lower components, where there are some interesting signals of population structure, but less variance is captured. For these lower components, results from conventional and randomized PCA are not so similar.


{% highlight python %}
fig = plt.figure(figsize=(8, 4))

ax = fig.add_subplot(1, 2, 1)
plot_pca_coords(coords1, model1, 4, 5, ax=ax, sample_population=df_samples.population.values)
ax.set_title('Conventional PCA')

ax = fig.add_subplot(1, 2, 2)
plot_pca_coords(coords5, model5, 4, 5, ax=ax, sample_population=df_samples.population.values)
ax.set_title('Randomized PCA')

fig.tight_layout();
{% endhighlight %}


![png](/assets/2015-09-28-fast-pca_files/2015-09-28-fast-pca_84_0.png)


So an important caveat when using randomized PCA is that lower components may not be resolved very well, compared with conventional PCA.

## Post-script: OpenBLAS

The original version of this post was run with NumPy built without any linear algebra optimisations. Following a comment from Andreas Noack I have re-run this notebook with NumPy built with OpenBLAS. This significantly speeds up the runtime for both conventional and randomized PCA. With the data size used in this notebook, the bottleneck is now more data preparation (LD pruning) than running the PCA itself.


{% highlight python %}
np.__config__.show()
{% endhighlight %}

    openblas_lapack_info:
        libraries = ['openblas', 'openblas']
        language = c
        define_macros = [('HAVE_CBLAS', None)]
        library_dirs = ['/home/aliman/miniconda3/envs/biipy240/lib']
    blas_mkl_info:
      NOT AVAILABLE
    blas_opt_info:
        libraries = ['openblas', 'openblas']
        language = c
        define_macros = [('HAVE_CBLAS', None)]
        library_dirs = ['/home/aliman/miniconda3/envs/biipy240/lib']
    openblas_info:
        libraries = ['openblas', 'openblas']
        language = c
        define_macros = [('HAVE_CBLAS', None)]
        library_dirs = ['/home/aliman/miniconda3/envs/biipy240/lib']
    lapack_opt_info:
        libraries = ['openblas', 'openblas']
        language = c
        define_macros = [('HAVE_CBLAS', None)]
        library_dirs = ['/home/aliman/miniconda3/envs/biipy240/lib']



{% highlight python %}
!lscpu
{% endhighlight %}

    Architecture:          x86_64
    CPU op-mode(s):        32-bit, 64-bit
    Byte Order:            Little Endian
    CPU(s):                8
    On-line CPU(s) list:   0-7
    Thread(s) per core:    2
    Core(s) per socket:    4
    Socket(s):             1
    NUMA node(s):          1
    Vendor ID:             GenuineIntel
    CPU family:            6
    Model:                 94
    Model name:            Intel(R) Xeon(R) CPU E3-1505M v5 @ 2.80GHz
    Stepping:              3
    CPU MHz:               3317.453
    CPU max MHz:           3700.0000
    CPU min MHz:           800.0000
    BogoMIPS:              5615.30
    Virtualisation:        VT-x
    L1d cache:             32K
    L1i cache:             32K
    L2 cache:              256K
    L3 cache:              8192K
    NUMA node0 CPU(s):     0-7
    Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc art arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc aperfmperf eagerfpu pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch epb intel_pt tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm mpx rdseed adx smap clflushopt xsaveopt xsavec xgetbv1 dtherm ida arat pln pts hwp hwp_notify hwp_act_window hwp_epp



{% highlight python %}
import datetime
print(datetime.datetime.now().isoformat())
{% endhighlight %}

    2016-11-01T19:57:14.214815

