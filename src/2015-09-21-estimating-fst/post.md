---
layout: post
title: Estimating F<sub>ST</sub>
---


In phase 1 of the [Ag1000G project](http://www.malariagen.net/ag1000g) we have whole genome sequence data for mosquitoes from 9 African countries. As part of our analysis of population structure, I recently needed to calculate average F<sub>ST</sub> between each pair of populations. I also needed to calculate F<sub>ST</sub> in windows over the genome, to look for genome regions that are particularly differentiated between certain populations.

F<sub>ST</sub> is a statistic which seems simple at first yet quickly becomes very technical when you start reading the literature. I asked around my lab for advice and [George](http://www.well.ox.ac.uk/george-busby) pointed me to [Bhatia et al. (2013)](http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3759727&tool=pmcentrez&rendertype=abstract) which provides some clear advice on how to estimate F<sub>ST</sub>. However, Bhatia et al. were working with relatively well-studied human populations, and mosquito population genetics can get pretty extreme by comparison, so I didn't want to take anything for granted.

To help explore the impact of different F<sub>ST</sub> estimators and SNP ascertainment schemes, I implemented both the Weir and Cockerham estimator and the Hudson estimator in  [`scikit-allel`](http://scikit-allel.readthedocs.org/en/latest/stats/fst.html). This post gives some examples of using these functions with large scale SNP data, and some practical experiences from applying them to mosquito populations.

## Setup


{% highlight python %}
import numpy as np
import h5py
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas
import allel; print('scikit-allel', allel.__version__)
{% endhighlight %}

    scikit-allel 1.0.3


I have a copy of the [Ag1000G phase 1 AR3 data release](http://www.malariagen.net/data/ag1000g-phase1-ar3) on a local drive. The SNP genotype data is available in an HDF5 file.


{% highlight python %}
callset_fn = 'data/2015-09-21/ag1000g.phase1.ar3.pass.h5'
callset = h5py.File(callset_fn, mode='r')
callset
{% endhighlight %}




    <HDF5 file "ag1000g.phase1.ar3.pass.h5" (mode r)>



Let's work with chromosome arm 3L.


{% highlight python %}
chrom = '3L'
# load all variant positions
pos_all = allel.SortedIndex(callset[chrom]['variants']['POS'])
pos_all
{% endhighlight %}




<div class="allel allel-DisplayAs1D"><span>&lt;SortedIndex shape=(9643193,) dtype=int32&gt;</span><table><tr><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">9643188</th><th style="text-align: center">9643189</th><th style="text-align: center">9643190</th><th style="text-align: center">9643191</th><th style="text-align: center">9643192</th></tr><tr><td style="text-align: center">9790</td><td style="text-align: center">9798</td><td style="text-align: center">9812</td><td style="text-align: center">9818</td><td style="text-align: center">9829</td><td style="text-align: center">...</td><td style="text-align: center">41956532</td><td style="text-align: center">41956537</td><td style="text-align: center">41956541</td><td style="text-align: center">41956551</td><td style="text-align: center">41956556</td></tr></table></div>



There are 9,643,193 SNPs on this chromosome, genotyped in 765 individuals, so this is a relatively big dataset, too big to work with in memory uncompressed. To cope with larger datasets, scikit-allel can run computations directly against the data on disk, and can also make use of compression to reduce data size.


{% highlight python %}
genotype_all = allel.GenotypeChunkedArray(callset[chrom]['calldata']['genotype'])
genotype_all
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeChunkedArray shape=(9643193, 765, 2) dtype=int8 chunks=(6553, 10, 2)
   nbytes=13.7G cbytes=548.0M cratio=25.7
   compression=gzip compression_opts=3
   values=h5py._hl.dataset.Dataset&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">760</th><th style="text-align: center">761</th><th style="text-align: center">762</th><th style="text-align: center">763</th><th style="text-align: center">764</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">9643190</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">9643191</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">9643192</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></table></div>



This array would be 13.7G uncompressed, but [genotype data compress very well](http://alimanfoo.github.io/2016/09/21/genotype-compression-benchmark.html), so the actual size on disk is only 548M.

There is also a table of sample metadata which we'll need because it tells us which mosquito comes from which population.


{% highlight python %}
df_samples = pandas.read_csv('data/2015-09-21/samples.meta.txt', sep='\t', index_col='index')
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



The 'index' column in this table corresponds to the order of columns in the genotype array.

Let's pick two populations to work with.


{% highlight python %}
pop1 = 'BFM'
pop2 = 'AOM'
n_samples_pop1 = np.count_nonzero(df_samples.population == pop1)
n_samples_pop2 = np.count_nonzero(df_samples.population == pop2)
print(pop1, n_samples_pop1, pop2, n_samples_pop2)
{% endhighlight %}

    BFM 69 AOM 60


I've chosen BFM (*Anopheles coluzzii* from Burkina Faso) and AOM (*Anopheles coluzzii* from Angola) because F<sub>ST</sub> is reasonably high between these two populations.

Now compute allele counts in each population. This takes a minute or so.


{% highlight python %}
# dictionary mapping population names to sample indices
subpops = {
    pop1: df_samples[df_samples.population == pop1].index,
    pop2: df_samples[df_samples.population == pop2].index,
}
# allele counts
acs = genotype_all.count_alleles_subpops(subpops)
acs
{% endhighlight %}




<div class="allel allel-DisplayAsTable"><span>&lt;AlleleCountsChunkedTable shape=(9643193,) dtype=[('AOM', '&lt;i4', (4,)), ('BFM', '&lt;i4', (4,))]
   nbytes=294.3M cbytes=23.5M cratio=12.5
   values=allel.chunked.storage_zarr.ZarrTable&gt;</span><table><tr><th></th><th style="text-align: center">AOM</th><th style="text-align: center">BFM</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="3">...</td></tr><tr><th style="text-align: center">9643190</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr><tr><th style="text-align: center">9643191</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr><tr><th style="text-align: center">9643192</th><td style="text-align: center">[120   0   0   0]</td><td style="text-align: center">[138   0   0   0]</td></tr></table></div>



Finally, we can filter out variants that aren't segregating in the union of our two populations. Let's also filter out multiallelic variants for simplicity.


{% highlight python %}
acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
flt = acu.is_segregating() & (acu.max_allele() == 1)
print('retaining', np.count_nonzero(flt), 'SNPs')
{% endhighlight %}

    retaining 3177369 SNPs



{% highlight python %}
pos = pos_all.compress(flt)
ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
genotype = genotype_all.compress(flt, axis=0)
genotype
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeChunkedArray shape=(3177369, 765, 2) dtype=int8 chunks=(776, 765, 2)
   nbytes=4.5G cbytes=412.7M cratio=11.2
   compression=blosc compression_opts={'cname': 'lz4', 'clevel': 5, 'shuffle': 1}
   values=zarr.core.Array&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">760</th><th style="text-align: center">761</th><th style="text-align: center">762</th><th style="text-align: center">763</th><th style="text-align: center">764</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">3177366</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">3177367</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">3177368</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></table></div>



Note that the new `genotype` array we have created is in-memory, however we are making use of [compressed arrays](http://zarr.readthedocs.io/) to store these intermediate data, so the amount of memory required is small. Btw don't get confused with the name of the 'compress' method, this comes from [numpy](http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.compress.html) and means making a selection along a particular axis, it has nothing to do with data compression which is handled automatically by `scikit-allel`.

## Comparing F<sub>ST</sub> estimators

### Per-SNP estimates

Let's first compute the per-SNP F<sub>ST</sub> value from each of the two estimators. The Weir & Cockerham estimator takes a little longer because it has to revisit the genotype data. The Hudson estimator is faster because it only needs the allele counts, which we've already computed. 


{% highlight python %}
# sample indices for population 1
pop1_idx = subpops[pop1]
# sample indices for population 2
pop2_idx = subpops[pop2]
a, b, c = allel.stats.weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], max_allele=1)
snp_fst_wc = (a / (a + b + c))[:, 0]
snp_fst_wc
{% endhighlight %}




    array([-0.00102459,  0.01277143,  0.08574083, ..., -0.00102459,
            0.0554781 ,  0.08247673])




{% highlight python %}
num, den = allel.stats.hudson_fst(ac1, ac2)
snp_fst_hudson = num / den
snp_fst_hudson
{% endhighlight %}




    array([ 0.        ,  0.01459854,  0.08403361, ...,  0.        ,
            0.05042017,  0.07563025])




{% highlight python %}
fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_wc, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (%s) vs %s (%s), SNP $F_{ST}$' % (pop1, n_samples_pop1, pop2, n_samples_pop2));
{% endhighlight %}


![png](/assets/2015-09-21-estimating-fst_files/2015-09-21-estimating-fst_29_0.png)


With a couple of exceptions, the two estimators are virtually identical for all SNPs. However, one thing that Bhatia et al. warn is that the Weir & Cockerham estimator can give different results if sample sizes are unequal. We've chosen two populations with similar sample sizes, but what happens if we fake one of the populations to have a much smaller sample size?


{% highlight python %}
# keep only 20 samples from first population
pop1_idx_ds = subpops[pop1][:20]
a, b, c = allel.stats.weir_cockerham_fst(genotype, subpops=[pop1_idx_ds, pop2_idx], max_allele=1)
# there may be some non-segregating variants after down-sampling, suppress errors about zero division
with np.errstate(divide='ignore', invalid='ignore'):
    snp_fst_wc_ds = (a / (a + b + c))[:, 0]
snp_fst_wc_ds
{% endhighlight %}




    array([        nan,  0.03306823,  0.03605869, ...,         nan,
            0.02435467,  0.04264706])




{% highlight python %}
np.count_nonzero(np.isnan(snp_fst_wc_ds))
{% endhighlight %}




    1105405




{% highlight python %}
# recompute allele counts for downsampled population
ac1_ds = genotype.count_alleles(subpop=pop1_idx_ds, max_allele=1)
num, den = allel.stats.hudson_fst(ac1_ds, ac2)
# there may be some non-segregating variants after down-sampling, suppress errors about zero division
with np.errstate(divide='ignore', invalid='ignore'):
    snp_fst_hudson_ds = num / den
snp_fst_hudson_ds
{% endhighlight %}




    array([        nan,  0.        ,  0.08403361, ...,         nan,
            0.05042017,  0.07563025])




{% highlight python %}
fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson_ds, snp_fst_wc_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (20) vs %s (%s), SNP $F_{ST}$' % (pop1, pop2, n_samples_pop2));
{% endhighlight %}


![png](/assets/2015-09-21-estimating-fst_files/2015-09-21-estimating-fst_34_0.png)



{% highlight python %}
fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_hudson_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Hudson $F_{ST}$ (one population downsampled)')
ax.set_title('%s vs %s, SNP $F_{ST}$' % (pop1, pop2));
{% endhighlight %}


![png](/assets/2015-09-21-estimating-fst_files/2015-09-21-estimating-fst_35_0.png)



{% highlight python %}
fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_wc, snp_fst_wc_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Weir & Cockerham $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$ (one population downsampled)')
ax.set_title('%s vs %s, SNP $F_{ST}$' % (pop1, pop2));
{% endhighlight %}


![png](/assets/2015-09-21-estimating-fst_files/2015-09-21-estimating-fst_36_0.png)


When the sample sizes are unequal, the correspondance between the two estimators is clearly much less. Also, the Weir & Cockerham estimator appears to be systematically different with and without one population down-sampled.

### Chromosome average estimates

Now compute chromosome-wide average F<sub>ST</sub> with standard errors approximated via a block-jackknife.


{% highlight python %}
fst_wc, se_wc, vb_wc, _ = allel.stats.blockwise_weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], 
                                                                   blen=10000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc, se_wc))
{% endhighlight %}

    0.1057 +/- 0.0030 (Weir & Cockerham)



{% highlight python %}
fst_hudson, se_hudson, vb_hudson, _ = allel.stats.blockwise_hudson_fst(ac1, ac2, blen=10000)
print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))
{% endhighlight %}

    0.1065 +/- 0.0030 (Hudson)


The two estimates are very close, well within one standard error.

How about with one population downsampled to 20?


{% highlight python %}
fst_wc_ds, se_wc_ds, _, _ = allel.stats.blockwise_weir_cockerham_fst(genotype, subpops=[pop1_idx_ds, pop2_idx], 
                                                                     blen=10000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc_ds, se_wc_ds))
{% endhighlight %}

    0.1120 +/- 0.0030 (Weir & Cockerham)



{% highlight python %}
fst_hudson_ds, se_hudson_ds, _, _ = allel.stats.blockwise_hudson_fst(ac1_ds, ac2, blen=10000)
print('%.04f +/- %.04f (Hudson)' % (fst_hudson_ds, se_hudson_ds))
{% endhighlight %}

    0.1057 +/- 0.0030 (Hudson)


The two estimates are now separated by about two standard errors, with the Weir & Cockerham estimator inflated relative to the estimate with full samples from both populations

## SNP ascertainment

Another issue that Bhatia et al. discuss is SNP ascertainment. Basically, how you choose which SNPs to use when estimating F<sub>ST</sub> can make a difference. As I understand it, when computing average F<sub>ST</sub> you want to use a set of SNPs which segregated in the ancestral population, because changes in allele frequency at these SNPs will tell you something about genetic drift. 

Bhatia et al. recommend ascertaining SNPs by choosing SNPs that are segregating in a third "outgroup" population. However, we don't really have an obvious outgroup population in Ag1000G. So we then have four choices: (1) choose SNPs segregating in the first population; (2) choose SNPs segregating in the second population; (3) choose SNPs segregating in either population; (4) choose SNPs segregating in both populations.

Let's explore the impact of different ascertainment schemes, using the Hudson estimator.


{% highlight python %}
def compute_fst(ac1, ac2, scheme):
    
    if scheme == 'first':
        loc_asc = ac1.is_segregating()
    elif scheme == 'second':
        loc_asc = ac2.is_segregating()
    elif scheme == 'either':
        loc_asc = ac1.is_segregating() | ac2.is_segregating()
    elif scheme == 'both':
        loc_asc = ac1.is_segregating() & ac2.is_segregating()    
    n_snps = np.count_nonzero(loc_asc)
    
    ac1 = ac1.compress(loc_asc, axis=0)
    ac2 = ac2.compress(loc_asc, axis=0)
    
    fst, se, _, _ = allel.stats.blockwise_hudson_fst(ac1, ac2, blen=10000)
    
    print('%.04f +/- %.04f (using %s SNPs segregating in %s population)' % (fst, se, n_snps, scheme))
{% endhighlight %}


{% highlight python %}
for scheme in 'first', 'second', 'either', 'both':
    compute_fst(ac1, ac2, scheme)
{% endhighlight %}

    0.1027 +/- 0.0031 (using 2797341 SNPs segregating in first population)
    0.1140 +/- 0.0040 (using 1086498 SNPs segregating in second population)
    0.1065 +/- 0.0030 (using 3177369 SNPs segregating in either population)
    0.1101 +/- 0.0045 (using 706470 SNPs segregating in both population)


The spread of values here is more than three standard errors, so clearly ascertainment makes a difference. Here I'd be inclined to use SNPs segregating in both populations as it is a stricter criterion, however comments very welcome.

## Genome plot

Finally, let's plot F<sub>ST</sub> over the chromosome, to see if any regions are particularly differentiated.


{% highlight python %}
def plot_fst(ac1, ac2, pos, blen=2000):
    
    fst, se, vb, _ = allel.stats.blockwise_hudson_fst(ac1, ac2, blen=blen)
    
    # use the per-block average Fst as the Y coordinate
    y = vb
    
    # use the block centres as the X coordinate
    x = allel.stats.moving_statistic(pos, statistic=lambda v: (v[0] + v[-1]) / 2, size=blen)
    
    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    
{% endhighlight %}


{% highlight python %}
plot_fst(ac1, ac2, pos)
{% endhighlight %}


![png](/assets/2015-09-21-estimating-fst_files/2015-09-21-estimating-fst_55_0.png)


This plot suggests some genome regions where F<sub>ST</sub> is higher than the chromosome-wide average, which are interesting to follow up.

## Conclusions

Hudson's F<sub>ST</sub> estimator is more robust to unequal sample sizes, and faster to compute because it only requires allele counts as input. 

SNP ascertainment also makes a difference. It's probably a good idea to try different ascertainment schemes to see what impact they have on the results. 

## Further reading

* Bhatia, G., Patterson, N., Sankararaman, S., & Price, A. L. (2013). [Estimating and interpreting FST: the impact of rare variants](http://doi.org/10.1101/gr.154831.113). Genome Research, 23(9), 1514–21. 
* [`scikit-allel` F<sub>ST</sub> functions](http://scikit-allel.readthedocs.org/en/latest/stats/fst.html)

<hr/>


{% highlight python %}
import datetime
print(datetime.datetime.now().isoformat())
{% endhighlight %}

    2016-11-01T19:32:06.817362

