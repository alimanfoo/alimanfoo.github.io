# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# In phase 1 of the [Ag1000G project](http://www.malariagen.net/ag1000g) we have whole genome sequence data for mosquitoes from 9 African countries. As part of our analysis of population structure, I recently needed to calculate average F<sub>ST</sub> between each pair of populations. I also needed to calculate F<sub>ST</sub> in windows over the genome, to look for genome regions that are particularly differentiated between certain populations.

# %% [markdown]
# F<sub>ST</sub> is a statistic which seems simple at first yet quickly becomes very technical when you start reading the literature. I asked around my lab for advice and [George](http://www.well.ox.ac.uk/george-busby) pointed me to [Bhatia et al. (2013)](http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3759727&tool=pmcentrez&rendertype=abstract) which provides some clear advice on how to estimate F<sub>ST</sub>. However, Bhatia et al. were working with relatively well-studied human populations, and mosquito population genetics can get pretty extreme by comparison, so I didn't want to take anything for granted.

# %% [markdown]
# To help explore the impact of different F<sub>ST</sub> estimators and SNP ascertainment schemes, I implemented both the Weir and Cockerham estimator and the Hudson estimator in  [`scikit-allel`](http://scikit-allel.readthedocs.org/en/latest/stats/fst.html). This post gives some examples of using these functions with large scale SNP data, and some practical experiences from applying them to mosquito populations.

# %% [markdown]
# ## Setup

# %%
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
import zarr
import urllib
import shutil
from pathlib import Path
import zipfile
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas
import allel; print('scikit-allel', allel.__version__)


# %% [markdown]
# I'm going to use data from the [Ag1000G phase 1 AR3 release](http://www.malariagen.net/data/ag1000g-phase1-AR3). Let's download some data from the [Ag1000G public FTP site](ftp://ngs.sanger.ac.uk/production/ag1000g/). The files we're downloading are ~200-400 Mb so this may take a little while, depending on your internet connection.

# %%
def download(source_url, dest_path):
    """Helper function to download a remote file to the local file system."""
    with urllib.request.urlopen(source_url) as source:
        with open(dest_path, mode='wb') as dest:
            shutil.copyfileobj(source, dest)

            
# base URL for variation data files to be downloaded
variation_base_url = 'ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/zarr/'
# download files to a local directory called "data"
dest_dir = Path('data')
dest_dir.mkdir(exist_ok=True)

# files to be downloaded
variation_files = [
    'ag1000g.phase1.ar3.pass.biallelic.metadata.zip',
    'ag1000g.phase1.ar3.pass.biallelic.3L.variants.zip',
    'ag1000g.phase1.ar3.pass.biallelic.3L.calldata.GT.zip',
]

for f in variation_files:
    print('downloading', f)
    dest_path = dest_dir / f
    if not dest_path.exists():
        source_url = variation_base_url + f
        download(source_url, dest_path)
    print('extracting', f)
    with zipfile.ZipFile(dest_path, mode='r') as z:
        z.extractall(dest_dir)

# also grab some sample metadata
samples_url = 'ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/samples/samples.meta.txt'
download(samples_url, dest_dir / 'ag1000g.phase1.ar3.samples.meta.txt')

# %%
callset = zarr.open_consolidated('data/ag1000g.phase1.ar3.pass.biallelic')
callset

# %% [markdown]
# Let's work with chromosome arm 3L.

# %%
chrom = '3L'
# load all variant positions
pos_all = allel.SortedIndex(callset[chrom]['variants']['POS'])
pos_all

# %% [markdown]
# There are 7,449,486 SNPs on this chromosome, genotyped in 765 individuals, so this is a relatively big dataset, too big to work with in memory uncompressed. To cope with larger datasets, scikit-allel can run computations directly against the data on disk, and can also make use of compression to reduce data size.

# %%
gt_all = allel.GenotypeDaskArray(callset[chrom]['calldata']['GT'])
gt_all

# %% [markdown]
# This array would be 11.4G uncompressed, but [genotype data compress very well](http://alimanfoo.github.io/2016/09/21/genotype-compression-benchmark.html), so the actual size on disk is only 182M.

# %% [markdown]
# There is also a table of sample metadata which we'll need because it tells us which mosquito comes from which population.

# %%
df_samples = pandas.read_csv('data/ag1000g.phase1.ar3.samples.meta.txt', sep='\t', index_col='index')
df_samples.head()

# %% [markdown]
# The 'index' column in this table corresponds to the order of columns in the genotype array.

# %% [markdown]
# Let's pick two populations to work with.

# %%
pop1 = 'BFM'
pop2 = 'AOM'
n_samples_pop1 = np.count_nonzero(df_samples.population == pop1)
n_samples_pop2 = np.count_nonzero(df_samples.population == pop2)
print(pop1, n_samples_pop1, pop2, n_samples_pop2)

# %% [markdown]
# I've chosen BFM (*Anopheles coluzzii* from Burkina Faso) and AOM (*Anopheles coluzzii* from Angola) because F<sub>ST</sub> is reasonably high between these two populations.

# %% [markdown]
# Now compute allele counts in each population. This may take a few seconds.

# %%
import dask.array as da
from dask.diagnostics import ProgressBar

# %%
# sample indices for each population
pop1_idx = df_samples[df_samples.population == pop1].index
pop2_idx = df_samples[df_samples.population == pop2].index

# compute allele counts
with ProgressBar():
    ac1_all = gt_all.count_alleles(max_allele=1, subpop=pop1_idx).compute()
    ac2_all = gt_all.count_alleles(max_allele=1, subpop=pop2_idx).compute()


# %% [markdown]
# Finally, we can filter out variants that aren't segregating in the union of our two populations.

# %%
acu = ac1_all + ac2_all
acu

# %%
flt = acu.is_segregating()
print('retaining', np.count_nonzero(flt), 'SNPs')

# %%
pos = pos_all.compress(flt)
ac1 = ac1_all.compress(flt, axis=0)
ac2 = ac2_all.compress(flt, axis=0)
gt = gt_all.compress(flt, axis=0)
gt

# %% [markdown]
# ## Comparing F<sub>ST</sub> estimators

# %% [markdown]
# ### Per-SNP estimates

# %% [markdown]
# Let's first compute the per-SNP F<sub>ST</sub> value from each of the two estimators. The Weir & Cockerham estimator takes a little longer because it has to revisit the genotype data. The Hudson estimator is faster because it only needs the allele counts, which we've already computed. 

# %%
# TODO fully daskify implementation of weir_cockerham_fst
a, b, c = allel.weir_cockerham_fst(gt, subpops=[pop1_idx, pop2_idx], max_allele=1, blen=200000)
snp_fst_wc = (a / (a + b + c))[:, 0]
snp_fst_wc

# %%
num, den = allel.hudson_fst(ac1, ac2)
snp_fst_hudson = num / den
snp_fst_hudson

# %%
fig, ax = plt.subplots(figsize=(7, 7))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_wc, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (%s) vs %s (%s), SNP $F_{ST}$' % (pop1, n_samples_pop1, pop2, n_samples_pop2));

# %% [markdown]
# With a couple of exceptions, the two estimators are virtually identical for all SNPs. However, one thing that Bhatia et al. warn is that the Weir & Cockerham estimator can give different results if sample sizes are unequal. We've chosen two populations with similar sample sizes, but what happens if we fake one of the populations to have a much smaller sample size?

# %%
# keep only 20 samples from first population
pop1_idx_ds = subpops[pop1][:20]
a, b, c = allel.weir_cockerham_fst(gt, subpops=[pop1_idx_ds, pop2_idx], max_allele=1, blen=200000)
# there may be some non-segregating variants after down-sampling, suppress errors about zero division
with np.errstate(divide='ignore', invalid='ignore'):
    snp_fst_wc_ds = (a / (a + b + c))[:, 0]
snp_fst_wc_ds

# %%
np.count_nonzero(np.isnan(snp_fst_wc_ds))

# %%
# recompute allele counts for downsampled population
ac1_ds = gt.count_alleles(subpop=pop1_idx_ds, max_allele=1).compute()
num, den = allel.hudson_fst(ac1_ds, ac2)
# there may be some non-segregating variants after down-sampling, suppress errors about zero division
with np.errstate(divide='ignore', invalid='ignore'):
    snp_fst_hudson_ds = num / den
snp_fst_hudson_ds

# %%
fig, ax = plt.subplots(figsize=(7, 7))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson_ds, snp_fst_wc_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (20) vs %s (%s), SNP $F_{ST}$' % (pop1, pop2, n_samples_pop2));

# %%
fig, ax = plt.subplots(figsize=(7, 7))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_hudson_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Hudson $F_{ST}$ (one population downsampled)')
ax.set_title('%s vs %s, SNP $F_{ST}$' % (pop1, pop2));

# %%
fig, ax = plt.subplots(figsize=(7, 7))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_wc, snp_fst_wc_ds, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Weir & Cockerham $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$ (one population downsampled)')
ax.set_title('%s vs %s, SNP $F_{ST}$' % (pop1, pop2));

# %% [markdown]
# When the sample sizes are unequal, the correspondance between the two estimators is clearly much less. Also, the Weir & Cockerham estimator appears to be systematically different with and without one population down-sampled.

# %% [markdown]
# ### Chromosome average estimates

# %% [markdown]
# Now compute chromosome-wide average F<sub>ST</sub> with standard errors approximated via a block-jackknife.

# %%
# TODO fix this!
fst_wc, se_wc, vb_wc, _ = allel.average_weir_cockerham_fst(gt, subpops=[pop1_idx, pop2_idx], 
                                                           blen=10000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc, se_wc))

# %%
fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))

# %% [markdown]
# The two estimates are very close, well within one standard error.

# %% [markdown]
# How about with one population downsampled to 20?

# %%
# TODO fix this!
fst_wc_ds, se_wc_ds, _, _ = allel.average_weir_cockerham_fst(gt, subpops=[pop1_idx_ds, pop2_idx], 
                                                             blen=10000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc_ds, se_wc_ds))

# %%
fst_hudson_ds, se_hudson_ds, _, _ = allel.blockwise_hudson_fst(ac1_ds, ac2, blen=10000)
print('%.04f +/- %.04f (Hudson)' % (fst_hudson_ds, se_hudson_ds))


# %% [markdown]
# The two estimates are now separated by about two standard errors, with the Weir & Cockerham estimator inflated relative to the estimate with full samples from both populations

# %% [markdown]
# ## SNP ascertainment

# %% [markdown]
# Another issue that Bhatia et al. discuss is SNP ascertainment. Basically, how you choose which SNPs to use when estimating F<sub>ST</sub> can make a difference. As I understand it, when computing average F<sub>ST</sub> you want to use a set of SNPs which segregated in the ancestral population, because changes in allele frequency at these SNPs will tell you something about genetic drift. 
#
# Bhatia et al. recommend ascertaining SNPs by choosing SNPs that are segregating in a third "outgroup" population. However, we don't really have an obvious outgroup population in Ag1000G. So we then have four choices: (1) choose SNPs segregating in the first population; (2) choose SNPs segregating in the second population; (3) choose SNPs segregating in either population; (4) choose SNPs segregating in both populations.
#
# Let's explore the impact of different ascertainment schemes, using the Hudson estimator.

# %%
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
    
    fst, se, _, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
    
    print('%.04f +/- %.04f (using %s SNPs segregating in %s population)' % (fst, se, n_snps, scheme))


# %%
for scheme in 'first', 'second', 'either', 'both':
    compute_fst(ac1, ac2, scheme)


# %% [markdown]
# The spread of values here is more than three standard errors, so clearly ascertainment makes a difference. Here I'd be inclined to use SNPs segregating in both populations as it is a stricter criterion, however comments very welcome.

# %% [markdown]
# ## Genome plot

# %% [markdown]
# Finally, let's plot F<sub>ST</sub> over the chromosome, to see if any regions are particularly differentiated.

# %%
def plot_fst(ac1, ac2, pos, blen=2000):
    
    fst, se, vb, _ = allel.average_hudson_fst(ac1, ac2, blen=blen)
    
    # use the per-block average Fst as the Y coordinate
    y = vb
    
    # use the block centres as the X coordinate
    x = allel.moving_statistic(pos, statistic=lambda v: (v[0] + v[-1]) / 2, size=blen)
    
    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    


# %%
plot_fst(ac1, ac2, pos)

# %% [markdown]
# This plot suggests some genome regions where F<sub>ST</sub> is higher than the chromosome-wide average, which are interesting to follow up.

# %% [markdown]
# ## Conclusions

# %% [markdown]
# Hudson's F<sub>ST</sub> estimator is more robust to unequal sample sizes, and faster to compute because it only requires allele counts as input. 
#
# SNP ascertainment also makes a difference. It's probably a good idea to try different ascertainment schemes to see what impact they have on the results. 

# %% [markdown]
# ## Further reading

# %% [markdown]
# * Bhatia, G., Patterson, N., Sankararaman, S., & Price, A. L. (2013). [Estimating and interpreting FST: the impact of rare variants](http://doi.org/10.1101/gr.154831.113). Genome Research, 23(9), 1514–21. 
# * [`scikit-allel` F<sub>ST</sub> functions](http://scikit-allel.readthedocs.org/en/latest/stats/fst.html)
#
# <hr/>

# %%
import datetime
print(datetime.datetime.now().isoformat())

# %%
