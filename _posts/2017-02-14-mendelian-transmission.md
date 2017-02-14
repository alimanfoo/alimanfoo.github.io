---
layout: post
title: Mendelian transmission
---


*This post has some examples of analysing a genetic cross, using [scikit-allel](@@) and standard scientific Python libraries ([NumPy](@@), [matplotlib](@@), etc.). As usual, if you spot any errors or have any suggestions, please drop a comment below.*

## Setup


{% highlight python %}
import numpy as np
import pandas
import h5py
import allel
{% endhighlight %}

I'm going to use data from the [Ag1000G](@@) project [phase 1 data releases](@@), which includes genotype calls for four genetic crosses. Each cross involves two parents (a mother and a father) and up to 20 offspring (progeny). These are mosquito crosses, but mosquitoes are diploid (like us), so the genetics are the same as if analysing a cross or family of any other diploid species.

Here's some information about the crosses.


{% highlight python %}
samples = pandas.read_csv('data/phase1.AR3.1/samples/cross.samples.meta.txt', sep='\t')
samples.head()
{% endhighlight %}




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ox_code</th>
      <th>cross</th>
      <th>role</th>
      <th>n_reads</th>
      <th>median_cov</th>
      <th>mean_cov</th>
      <th>sex</th>
      <th>colony_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>AD0231-C</td>
      <td>29-2</td>
      <td>parent</td>
      <td>451.762</td>
      <td>20.0</td>
      <td>19.375</td>
      <td>F</td>
      <td>ghana</td>
    </tr>
    <tr>
      <th>1</th>
      <td>AD0232-C</td>
      <td>29-2</td>
      <td>parent</td>
      <td>572.326</td>
      <td>25.0</td>
      <td>24.370</td>
      <td>M</td>
      <td>kisumu</td>
    </tr>
    <tr>
      <th>2</th>
      <td>AD0234-C</td>
      <td>29-2</td>
      <td>progeny</td>
      <td>489.057</td>
      <td>16.0</td>
      <td>15.742</td>
      <td>F</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>AD0235-C</td>
      <td>29-2</td>
      <td>progeny</td>
      <td>539.649</td>
      <td>17.0</td>
      <td>17.364</td>
      <td>F</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>AD0236-C</td>
      <td>29-2</td>
      <td>progeny</td>
      <td>537.237</td>
      <td>17.0</td>
      <td>17.284</td>
      <td>F</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




{% highlight python %}
samples.cross.value_counts()
{% endhighlight %}




    29-2    22
    46-9    22
    36-9    20
    42-4    16
    Name: cross, dtype: int64



So there are four crosses. The two largest (crosses '29-2' and '46-9') each have 22 individuals (2 parents, 20 progeny), and the smallest ('42-4') has 16 individuals (2 parents, 14 progeny).

All individuals in all crosses have been sequenced on Illumina HiSeq machines, and then have had genotypes called at variant sites discovered in a cohort of wild specimens. The genotype data were originally in [VCF format](@@), however for ease of analysis we've converted the data to [HDF5 format](@@).

Open the file containing genotype data for chromosome arm 3R.


{% highlight python %}
callset = h5py.File('data/phase1.AR3/variation/crosses/ar3/hdf5/ag1000g.crosses.phase1.ar3sites.3R.h5',
                    mode='r')
callset
{% endhighlight %}




    <HDF5 file "ag1000g.crosses.phase1.ar3sites.3R.h5" (mode r)>



To analyse your own data using the examples shown below, you would need to convert the genotype data to either NumPy or HDF5 format. If you have the data in VCF format then you can use the [vcfnp](@@) utility to perform the conversion. There is some documentation in the [vcfnp README](@@) but please feel free to email me if you run into any difficulty. This data conversion is the painful step, if you can get over it then the rest should be relatively plain sailing.

Here I am going to start from unphased genotype data. If you have already phased the data that's fine, convert to NumPy or HDF5 as you would for unphased data.

In total I have genotype calls in 80 individuals at 22,632,425 SNPs on chromosome 3R. However, I'm only going to analyse data for a single cross, '29-2', between a mother from the 'Ghana' colony and a father from the 'Kisumu' colony. I can subset out the genotype data for just this cross, and keep only SNPs that are segregating in this cross. I'm also only going to keep SNPs that passed all quality filters.


{% highlight python %}
genotypes = allel.GenotypeChunkedArray(callset['3R/calldata/genotype'])
genotypes
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeChunkedArray shape=(22632425, 80, 2) dtype=int8 chunks=(6553, 10, 2)
   nbytes=3.4G cbytes=131.6M cratio=26.2
   compression=gzip compression_opts=1
   values=h5py._hl.dataset.Dataset&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">75</th><th style="text-align: center">76</th><th style="text-align: center">77</th><th style="text-align: center">78</th><th style="text-align: center">79</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">1/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">...</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">22632422</th><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">22632423</th><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">22632424</th><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">...</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">./.</td><td style="text-align: center">1/1</td></tr></table></div>




{% highlight python %}
# locate the indices of the samples within the callset
sample_indices = samples[samples.cross == '29-2'].index.values.tolist()
sample_indices
{% endhighlight %}




    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]




{% highlight python %}
# do an allele count to find segregating variants
ac = genotypes.count_alleles(max_allele=3, subpop=sample_indices)[:]
ac
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;AlleleCountsArray shape=(22632425, 4) dtype=int32&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">19</td><td style="text-align: center">25</td><td style="text-align: center"> 0</td><td style="text-align: center"> 0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">44</td><td style="text-align: center"> 0</td><td style="text-align: center"> 0</td><td style="text-align: center"> 0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">44</td><td style="text-align: center"> 0</td><td style="text-align: center"> 0</td><td style="text-align: center"> 0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="5">...</td></tr><tr><th style="text-align: center">22632422</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">22632423</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">22632424</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr></table></div>




{% highlight python %}
# how many SNPs are segregating within the cross?
loc_seg = ac.is_segregating()
np.count_nonzero(loc_seg)
{% endhighlight %}




    2142258




{% highlight python %}
# locate SNPs that passed all quality filters
loc_pass = callset['3R/variants/FILTER_PASS'][:]
{% endhighlight %}


{% highlight python %}
# perform the subset and load the results into memory uncompressed
genotypes_cross = genotypes.subset(loc_seg & loc_pass, sample_indices)[:]
genotypes_cross
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(709399, 22, 2) dtype=int8&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">17</th><th style="text-align: center">18</th><th style="text-align: center">19</th><th style="text-align: center">20</th><th style="text-align: center">21</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">./.</td><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">./.</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">709396</th><td style="text-align: center">1/1</td><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">...</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center">709397</th><td style="text-align: center">0/0</td><td style="text-align: center">1/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">...</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center">709398</th><td style="text-align: center">0/0</td><td style="text-align: center">1/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">...</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td><td style="text-align: center">0/1</td></tr></table></div>



Now I have an array of genotype calls at 709,399 segregating SNPs in 22 individuals. The mother is the first column, the father is the second column, and the progeny are the remaining columns. You'll notice that the mother's genotype call is missing at the first two SNPs: we could remove these, but we'll leave them in, just to check that the analyses are robust to some missing data.

## Phasing by transmission

I'm starting from unphased genotype calls, so the first thing to do is phase the calls to generate haplotypes. There are several options for phasing a cross. Here I'm going to use the [`phase_by_transmission()`](http://scikit-allel.readthedocs.io/en/latest/stats/mendel.html#allel.stats.mendel.phase_by_transmission) function from [scikit-allel](@@), because it's convenient and fast (couple of seconds). We've found this function works well for crosses with relatively large numbers of progeny. However, if you have a smaller family with only a couple of progeny, and/or you have a more complicated pedigree with multiple generations, try phasing with [SHAPEIT2 + DuoHMM](@@).


{% highlight python %}
genotypes_cross_phased = allel.phase_by_transmission(genotypes_cross, window_size=100)
genotypes_cross_phased
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(709399, 22, 2) dtype=int8&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">17</th><th style="text-align: center">18</th><th style="text-align: center">19</th><th style="text-align: center">20</th><th style="text-align: center">21</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">./.</td><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">./.</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0|1</td><td style="text-align: center">0|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">0|0</td><td style="text-align: center">...</td><td style="text-align: center">0|0</td><td style="text-align: center">0|0</td><td style="text-align: center">0|0</td><td style="text-align: center">0|0</td><td style="text-align: center">0|0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">709396</th><td style="text-align: center">1|1</td><td style="text-align: center">0|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">...</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td><td style="text-align: center">1|0</td></tr><tr><th style="text-align: center">709397</th><td style="text-align: center">0|0</td><td style="text-align: center">1|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">...</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td></tr><tr><th style="text-align: center">709398</th><td style="text-align: center">0|0</td><td style="text-align: center">1|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">...</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td><td style="text-align: center">0|1</td></tr></table></div>



Notice that most of the genotype calls in the snippet shown above now have a pipe character ('`|`') as the allele separator, indicating the call is phased. The first two SNPs are not phased, however, because the genotype call for one of the parents is missing.

## Visualising transmission

Now we have phased data, let's plot the transmission of alleles from parents to progeny. This is a useful diagnostic for assessing the quality of the phasing, and also gives an indication of how much recombination has occurred.

Before plotting, I'm going to separate out the data into maternal and paternal haplotypes. The maternal haplotypes are the two haplotypes carried by the mother, and the haplotype in each of the progeny inherited from the mother. The paternal haplotypes are the same but for the father.


{% highlight python %}
# pull out mother's genotypes from the first column
genotypes_mother = genotypes_cross_phased[:, 0]
# convert to haplotype array
haplotypes_mother = genotypes_mother.to_haplotypes()
# pull out maternal haplotypes from the progeny
haplotypes_progeny_maternal = allel.HaplotypeArray(genotypes_cross_phased[:, 2:, 0])
# stack mother's haplotypes alongside haplotypes she transmitted to her progeny
haplotypes_maternal = haplotypes_mother.concatenate(haplotypes_progeny_maternal, axis=1)
haplotypes_maternal
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;HaplotypeArray shape=(709399, 22) dtype=int8&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">17</th><th style="text-align: center">18</th><th style="text-align: center">19</th><th style="text-align: center">20</th><th style="text-align: center">21</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">.</td><td style="text-align: center">.</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">.</td><td style="text-align: center">.</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">709396</th><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">...</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td></tr><tr><th style="text-align: center">709397</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">709398</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr></table></div>



Let's fix on the mother for a moment. The mother has two haplotypes. Because recombination occurs during gamete formation, each haplotype the mother passes on to her progeny is a unique mosaic of her own two haplotypes. For any SNP where the two maternal haplotypes carry a different allele, we can "paint" the maternal haplotypes within the progeny according to which of the mother's two alleles were inherited, using the [`paint_transmission()`](@@) function.


{% highlight python %}
painting_maternal = allel.paint_transmission(haplotypes_mother, haplotypes_progeny_maternal)
painting_maternal
{% endhighlight %}




    array([[6, 6, 6, ..., 6, 6, 6],
           [6, 6, 6, ..., 6, 6, 6],
           [2, 2, 1, ..., 1, 1, 1],
           ..., 
           [4, 4, 4, ..., 4, 4, 4],
           [3, 3, 3, ..., 3, 3, 3],
           [3, 3, 3, ..., 3, 3, 3]], dtype=uint8)



This new "painting" array is an array of integer codes, where each number means something:


{% highlight python %}
help(allel.paint_transmission)
{% endhighlight %}

    Help on function paint_transmission in module allel.stats.mendel:
    
    paint_transmission(parent_haplotypes, progeny_haplotypes)
        Paint haplotypes inherited from a single diploid parent according to
        their allelic inheritance.
        
        Parameters
        ----------
        parent_haplotypes : array_like, int, shape (n_variants, 2)
            Both haplotypes from a single diploid parent.
        progeny_haplotypes : array_like, int, shape (n_variants, n_progeny)
            Haplotypes found in progeny of the given parent, inherited from the
            given parent. I.e., haplotypes from gametes of the given parent.
        
        Returns
        -------
        painting : ndarray, uint8, shape (n_variants, n_progeny)
            An array of integers coded as follows: 1 = allele inherited from
            first parental haplotype; 2 = allele inherited from second parental
            haplotype; 3 = reference allele, also carried by both parental
            haplotypes; 4 = non-reference allele, also carried by both parental
            haplotypes; 5 = non-parental allele; 6 = either or both parental
            alleles missing; 7 = missing allele; 0 = undetermined.
        
        Examples
        --------
        >>> import allel
        >>> haplotypes = allel.HaplotypeArray([
        ...     [0, 0, 0, 1, 2, -1],
        ...     [0, 1, 0, 1, 2, -1],
        ...     [1, 0, 0, 1, 2, -1],
        ...     [1, 1, 0, 1, 2, -1],
        ...     [0, 2, 0, 1, 2, -1],
        ...     [0, -1, 0, 1, 2, -1],
        ...     [-1, 1, 0, 1, 2, -1],
        ...     [-1, -1, 0, 1, 2, -1],
        ... ], dtype='i1')
        >>> painting = allel.stats.paint_transmission(haplotypes[:, :2],
        ...                                           haplotypes[:, 2:])
        >>> painting
        array([[3, 5, 5, 7],
               [1, 2, 5, 7],
               [2, 1, 5, 7],
               [5, 4, 5, 7],
               [1, 5, 2, 7],
               [6, 6, 6, 7],
               [6, 6, 6, 7],
               [6, 6, 6, 7]], dtype=uint8)
    


We are particularly interested in plotting the "1" and "2" values, because these occur where the mother's haplotypes carried two different alleles, and so we have information about which allele has been transmitted.

We can do the same for the father.


{% highlight python %}
genotypes_father = genotypes_cross_phased[:, 1]
haplotypes_father = genotypes_father.to_haplotypes()
haplotypes_progeny_paternal = allel.HaplotypeArray(genotypes_cross_phased[:, 2:, 1])
haplotypes_paternal = haplotypes_father.concatenate(haplotypes_progeny_paternal, axis=1)
haplotypes_paternal
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;HaplotypeArray shape=(709399, 22) dtype=int8&gt;</span><table><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">17</th><th style="text-align: center">18</th><th style="text-align: center">19</th><th style="text-align: center">20</th><th style="text-align: center">21</th></tr><tr><th style="text-align: center">0</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">1</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">1</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">2</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center">709396</th><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">...</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center">709397</th><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">...</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td></tr><tr><th style="text-align: center">709398</th><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">...</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">1</td></tr></table></div>




{% highlight python %}
painting_paternal = allel.paint_transmission(haplotypes_father, haplotypes_progeny_paternal)
painting_paternal
{% endhighlight %}




    array([[5, 3, 3, ..., 3, 3, 3],
           [3, 3, 3, ..., 3, 3, 3],
           [3, 3, 3, ..., 3, 3, 3],
           ..., 
           [3, 3, 3, ..., 3, 3, 3],
           [4, 4, 4, ..., 4, 4, 4],
           [4, 4, 4, ..., 4, 4, 4]], dtype=uint8)



Notice the "5" code in this snippet. This code indicates an allele that is found on a progeny haplotype but not present on either of the parent's haplotypes. This is not possible according to the rules of Mendelian transmission, and hence is a form of "Mendelian error".

Now we have these "paintings", it is fairly straightforward to plot them.


{% highlight python %}

{% endhighlight %}
