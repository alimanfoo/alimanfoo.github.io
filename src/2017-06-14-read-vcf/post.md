---
layout: post
title: Extracting data from VCF files
---


*This post gives an introduction to functions for extracting data from [Variant Call Format (VCF)](https://samtools.github.io/hts-specs/VCFv4.3.pdf) files and loading into [NumPy](http://www.numpy.org/) arrays, [pandas](http://pandas.pydata.org/) data frames, [HDF5](https://support.hdfgroup.org/HDF5/) files or [Zarr](http://zarr.readthedocs.io) arrays for ease of analysis. These functions are available in [scikit-allel](http://scikit-allel.readthedocs.io/en/latest/) version 1.1 or later. Any feedback or bug reports welcome.*

*Update 2018-10-12: This post has been updated to use scikit-allel 1.1.10 and zarr 2.2, and adds examples of how to store data grouped by chromosome.*

## Introduction

### Variant Call Format (VCF)

VCF is a widely-used file format for genetic variation data. Here is an example of a small VCF file, based on the example given in the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf):


{% highlight python %}
with open('example.vcf', mode='r') as vcf:
    print(vcf.read())
{% endhighlight %}

    ##fileformat=VCFv4.3
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
    ##FILTER=<ID=q10,Description="Quality below 10">
    ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
    20	14370	rs6054257	G	A	29	PASS	DP=14;AF=0.5;DB	GT:DP	0/0:1	0/1:8	1/1:5
    20	17330	.	T	A	3	q10	DP=11;AF=0.017	GT:DP	0/0:3	0/1:5	0/0:41
    20	1110696	rs6040355	A	G,T	67	PASS	DP=10;AF=0.333,0.667;DB	GT:DP	0/2:6	1/2:0	2/2:4
    20	1230237	.	T	.	47	PASS	DP=13	GT:DP	0/0:7	0/0:4	./.:.
    20	1234567	microsat1	GTC	G,GTCT	50	PASS	DP=9	GT:DP	0/1:4	0/2:2	1/1:3


A VCF file begins with a number of meta-information lines, which start with two hash ('##') characters. Then there is a single header line beginning with a single hash ('#') character. After the header line there are data lines, with each data line describing a genetic variant at a particular position relative to the reference genome of whichever species you are studying. Each data line is divided into fields separated by tab characters. There are 9 fixed fields, labelled "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" and "FORMAT". Following these are fields containing data about samples, which usually contain a genotype call for each sample plus some associated data. 

For example, the first data line in the file above describes a variant on chromosome 20 at position 14370 relative to the B36 assembly of the human genome. The reference allele is 'G' and the alternate allele is 'A', so this is a single nucleotide polymorphism (SNP). In this file there are three fields with data about samples labelled 'NA00001', 'NA00002' and 'NA00003'. The genotype call in the first sample is '0/0', meaning that individual 'NA0001' is homozygous for the reference allele at this position. The genotype call for the second sample is '0/1' (you may need to scroll across to see this), which means that individual 'NA00002' is heterozygous for the reference and alternate alleles at this position.

### NumPy, pandas, HDF5, Zarr, ...

There are a number of software tools that can read VCF files and perform various analyses. However, if your dataset is large and/or you need to do some bespoke analysis, then it can be faster and more convenient to first extract the necessary data from the VCF file and load into a more efficient storage container.

For analysis and plotting of numerical data in Python, it is very convenient to load data into [NumPy arrays](https://docs.scipy.org/doc/numpy-dev/user/quickstart.html). A NumPy array is an in-memory data structure that provides support for fast arithmetic and data manipulation. For analysing tables of data, [pandas DataFrames](https://pandas.pydata.org/pandas-docs/stable/dsintro.html) provide useful features such as querying, aggregation and joins. When data are too large to fit into main memory, [HDF5 files](http://docs.h5py.org/en/latest/quick.html) and [Zarr arrays](http://zarr.readthedocs.io/en/latest/tutorial.html) can provide fast on-disk storage and retrieval of numerical arrays. 

[scikit-allel](http://scikit-allel.readthedocs.io/en/latest/) is a Python package intended to enable exploratory analysis of large-scale genetic variation data. Version 1.1.0 of scikit-allel adds some new functions for extracting data from VCF files and loading the data into NumPy arrays, pandas DataFrames or HDF5 files. Once you have extracted these data, there are many analyses that can be run interactively on a commodity laptop or desktop computer, even with large-scale datasets from population resequencing studies. To give a flavour of what analyses can be done, there are a few previous articles on my blog, touching on topics including [variant and sample QC](http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html), [allele frequency differentiation](http://alimanfoo.github.io/2015/09/21/estimating-fst.html), [population structure](http://alimanfoo.github.io/2015/09/28/fast-pca.html), and [genetic crosses](http://alimanfoo.github.io/2017/02/14/mendelian-transmission.html).

Until now, getting data out of VCF files and into NumPy etc. has been a bit of a pain point. Hopefully the new scikit-allel functions will make this a bit less of a hurdle. Let's take a look at the new functions...

## [`read_vcf()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.read_vcf)

Let's start with the scikit-allel function [`read_vcf()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.read_vcf). First, some imports:


{% highlight python %}
# import scikit-allel
import allel
# check which version is installed
print(allel.__version__)
{% endhighlight %}

    1.1.10


Read the example VCF file shown above, using default parameters:


{% highlight python %}
callset = allel.read_vcf('example.vcf')
{% endhighlight %}

The `callset` object returned by `read_vcf()` is a Python dictionary (`dict`). It contains several NumPy arrays, each of which can be accessed via a key. Here are the available keys: 


{% highlight python %}
sorted(callset.keys())
{% endhighlight %}




    ['calldata/GT',
     'samples',
     'variants/ALT',
     'variants/CHROM',
     'variants/FILTER_PASS',
     'variants/ID',
     'variants/POS',
     'variants/QUAL',
     'variants/REF']



The 'samples' array contains sample identifiers extracted from the header line in the VCF file.


{% highlight python %}
callset['samples']
{% endhighlight %}




    array(['NA00001', 'NA00002', 'NA00003'], dtype=object)



All arrays with keys beginning 'variants/' come from the fixed fields in the VCF file. For example, here is the data from the 'CHROM' field:


{% highlight python %}
callset['variants/CHROM']
{% endhighlight %}




    array(['20', '20', '20', '20', '20'], dtype=object)



Here is the data from the 'POS' field:


{% highlight python %}
callset['variants/POS']
{% endhighlight %}




    array([  14370,   17330, 1110696, 1230237, 1234567], dtype=int32)



Here is the data from the 'QUAL' field:


{% highlight python %}
callset['variants/QUAL']
{% endhighlight %}




    array([ 29.,   3.,  67.,  47.,  50.], dtype=float32)



All arrays with keys beginning 'calldata/' come from the sample fields in the VCF file. For example, here are the actual genotype calls from the 'GT' field:


{% highlight python %}
callset['calldata/GT']
{% endhighlight %}




    array([[[ 0,  0],
            [ 0,  1],
            [ 1,  1]],
    
           [[ 0,  0],
            [ 0,  1],
            [ 0,  0]],
    
           [[ 0,  2],
            [ 1,  2],
            [ 2,  2]],
    
           [[ 0,  0],
            [ 0,  0],
            [-1, -1]],
    
           [[ 0,  1],
            [ 0,  2],
            [ 1,  1]]], dtype=int8)



Note the -1 values for one of the genotype calls. By default scikit-allel uses -1 to indicate a missing value for any array with a signed integer data type (although you can change this if you want).

### Aside: genotype arrays

Because working with genotype calls is a very common task, scikit-allel has a [`GenotypeArray`](http://scikit-allel.readthedocs.io/en/latest/model/ndarray.html#allel.model.ndarray.GenotypeArray) class which adds some convenient functionality to an array of genotype calls. To use this class, pass the raw NumPy array into the `GenotypeArray` class constructor, e.g.:


{% highlight python %}
gt = allel.GenotypeArray(callset['calldata/GT'])
gt
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 3, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">1/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td><td style="text-align: center">2/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">./.</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">0/2</td><td style="text-align: center">1/1</td></tr></tbody></table></div>



One of the things that the `GenotypeArray` class does is provide a slightly more visually-appealing representation when used in a Jupyter notebook, as can be seen above. There are also methods for making various computations over the genotype calls. For example, the `is_het()` method locates all heterozygous genotype calls:


{% highlight python %}
gt.is_het()
{% endhighlight %}




    array([[False,  True, False],
           [False,  True, False],
           [ True,  True, False],
           [False, False, False],
           [ True,  True, False]], dtype=bool)



To give another example, the `count_het()` method will count heterozygous calls, summing over variants (axis=0) or samples (axis=1) if requested. E.g., to count the number of het calls per variant:


{% highlight python %}
gt.count_het(axis=1)
{% endhighlight %}




    array([1, 1, 2, 0, 2])



One more example, here is how to perform an allele count, i.e., count the number times each allele (0=reference, 1=first alternate, 2=second alternate, etc.) is observed for each variant:


{% highlight python %}
ac = gt.count_alleles()
ac
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;AlleleCountsArray shape=(5, 3) dtype=int32&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">3</td><td style="text-align: center">3</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">5</td><td style="text-align: center">1</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">1</td><td style="text-align: center">1</td><td style="text-align: center">4</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">4</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">2</td><td style="text-align: center">3</td><td style="text-align: center">1</td></tr></tbody></table></div>



### Fields

VCF files can often contain many fields of data, and you may only need to extract some of them to perform a particular analysis. You can select which fields to extract by passing a list of strings as the **fields** parameter. For example, let's extract the 'DP' field from within the 'INFO' field, and let's also extract the 'DP' field from the genotype call data:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['variants/DP', 'calldata/DP'])
sorted(callset.keys())
{% endhighlight %}




    ['calldata/DP', 'variants/DP']



Here is the data that we just extracted:


{% highlight python %}
callset['variants/DP']
{% endhighlight %}




    array([14, 11, 10, 13,  9], dtype=int32)




{% highlight python %}
callset['calldata/DP']
{% endhighlight %}




    array([[ 1,  8,  5],
           [ 3,  5, 41],
           [ 6,  0,  4],
           [ 7,  4, -1],
           [ 4,  2,  3]], dtype=int16)



I chose these two fields to illustrate the point that sometimes the same field name (e.g., 'DP') can be used both within the INFO field of a VCF and also within the genotype call data. When selecting fields, to make sure there is no ambiguity, you can include a prefix which is either 'variants/' or 'calldata/'. For example, if you provide 'variants/DP', then the `read_vcf()` function will look for an INFO field named 'DP'. If you provide 'calldata/DP' then `read_vcf()` will look for a FORMAT field named 'DP' within the call data. 

If you are feeling lazy, you can drop the 'variants/' and 'calldata/' prefixes, in which case `read_vcf()` will assume you mean 'variants/' if there is any ambiguity. E.g.:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['DP', 'GT'])
sorted(callset.keys())
{% endhighlight %}




    ['calldata/GT', 'variants/DP']



If you want to extract absolutely everything from a VCF file, then you can provide a special value ``'*'`` as the **fields** parameter: 


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields='*')
sorted(callset.keys())
{% endhighlight %}




    ['calldata/DP',
     'calldata/GT',
     'samples',
     'variants/AF',
     'variants/ALT',
     'variants/CHROM',
     'variants/DB',
     'variants/DP',
     'variants/FILTER_PASS',
     'variants/FILTER_q10',
     'variants/FILTER_s50',
     'variants/ID',
     'variants/POS',
     'variants/QUAL',
     'variants/REF',
     'variants/is_snp',
     'variants/numalt',
     'variants/svlen']



You can also provide the special value ``'variants/*'`` to request all variants fields (including all INFO), and the special value ``'calldata/*'`` to request all call data fields. 

If you don't specify the **fields** parameter, scikit-allel will default to extracting data from the fixed variants fields (but no INFO) and the GT genotype field if present (but no other call data).

### Types

NumPy arrays can have various [data types](https://docs.scipy.org/doc/numpy-1.12.0/user/basics.types.html), including signed integers ('int8', 'int16', 'int32', 'int64'), unsigned integers ('uint8', 'uint16', 'uint32', 'uint64'), floating point numbers ('float32', 'float64'), variable length strings ('object') and fixed length strings (e.g., 'S4' for a 4-character ASCII string). scikit-allel will try to choose a sensible default data type for the fields you want to extract, based on the meta-information in the VCF file, but you can override these via the **types** parameter. 

For example, by default the 'DP' INFO field is loaded into a 32-bit integer array:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['DP'])
callset['variants/DP']
{% endhighlight %}




    array([14, 11, 10, 13,  9], dtype=int32)



If you know the maximum value this field is going to contain, to save some memory you could choose a 16-bit integer array instead:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['DP'], types={'DP': 'int16'})
callset['variants/DP']
{% endhighlight %}




    array([14, 11, 10, 13,  9], dtype=int16)



You can also choose a floating-point data type, even for fields that are declared as type 'Integer' in the VCF meta-information, and vice versa. E.g.:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['DP'], types={'DP': 'float32'})
callset['variants/DP']
{% endhighlight %}




    array([ 14.,  11.,  10.,  13.,   9.], dtype=float32)



For fields containing textual data, there are two choices for data type. By default, scikit-allel will use an 'object' data type, which means that values are stored as an array of Python strings. E.g.:


{% highlight python %}
callset = allel.read_vcf('example.vcf')
callset['variants/REF']
{% endhighlight %}




    array(['G', 'T', 'A', 'T', 'GTC'], dtype=object)



The advantage of using 'object' dtype is that strings can be of any length. Alternatively, you can use a fixed-length string dtype, e.g.:


{% highlight python %}
callset = allel.read_vcf('example.vcf', types={'REF': 'S3'})
callset['variants/REF']
{% endhighlight %}




    array([b'G', b'T', b'A', b'T', b'GTC'],
          dtype='|S3')



Note that fixed-length string dtypes will cause any string values longer than the requested number of characters to be truncated. I.e., there can be some data loss. E.g., if using a single-character string for the 'REF' field, the correct value of 'GTC' for the final variant will get truncated to 'G':


{% highlight python %}
callset = allel.read_vcf('example.vcf', types={'REF': 'S1'})
callset['variants/REF']
{% endhighlight %}




    array([b'G', b'T', b'A', b'T', b'G'],
          dtype='|S1')



### Numbers

Some fields like 'ALT', 'AC' and 'AF' can have a variable number of values. I.e., each variant may have a different number of data values for these fields. One trade-off you have to make when loading data into NumPy arrays is that you cannot have arrays with a variable number of items per row. Rather, you have to fix the maximum number of possible items. While you lose some flexibility, you gain speed of access.

For fields like 'ALT', scikit-allel will choose a default number of expected values, which is set at 3. E.g., here is what you get by default:


{% highlight python %}
callset = allel.read_vcf('example.vcf')
callset['variants/ALT']
{% endhighlight %}




    array([['A', '', ''],
           ['A', '', ''],
           ['G', 'T', ''],
           ['', '', ''],
           ['G', 'GTCT', '']], dtype=object)



In this case, 3 is more that we need, because no variant has more than 2 ALT values. However, some VCF files (especially those including INDELs) may have more than 3 ALT values. 

If you need to increase or decrease the expected number of values for any field, you can do this via the **numbers** parameter. E.g., increase the number of ALT values to 5:


{% highlight python %}
callset = allel.read_vcf('example.vcf', numbers={'ALT': 5})
callset['variants/ALT']
{% endhighlight %}




    array([['A', '', '', '', ''],
           ['A', '', '', '', ''],
           ['G', 'T', '', '', ''],
           ['', '', '', '', ''],
           ['G', 'GTCT', '', '', '']], dtype=object)



Some care is needed here, because if you choose a value that is lower than the maximum number of values in the VCF file, then any extra values will not get extracted. E.g., the following would be fine if all variants were biallelic, but would lose information for any multi-allelic variants:


{% highlight python %}
callset = allel.read_vcf('example.vcf', numbers={'ALT': 1})
callset['variants/ALT']
{% endhighlight %}




    array(['A', 'A', 'G', '', 'G'], dtype=object)



#### Number of alternate alleles

Often there will be several fields within a VCF that all have a number of values that depends on the number of alternate alleles (declared with number 'A' or 'R' in the VCF meta-information). You can set the expected number of values simultaneously for all such fields via the **alt_number** parameter:


{% highlight python %}
callset = allel.read_vcf('example.vcf', fields=['ALT', 'AF'], alt_number=2)
{% endhighlight %}


{% highlight python %}
callset['variants/ALT']
{% endhighlight %}




    array([['A', ''],
           ['A', ''],
           ['G', 'T'],
           ['', ''],
           ['G', 'GTCT']], dtype=object)




{% highlight python %}
callset['variants/AF']
{% endhighlight %}




    array([[ 0.5  ,    nan],
           [ 0.017,    nan],
           [ 0.333,  0.667],
           [   nan,    nan],
           [   nan,    nan]], dtype=float32)



#### Genotype ploidy

By default, scikit-allel assumes you are working with a diploid organism, and so expects to parse out 2 alleles for each genotype call. If you are working with an organism with some other ploidy, you can change the expected number of alleles via the **numbers** parameter.

For example, here is an example VCF with tetraploid genotype calls:


{% highlight python %}
with open('example_polyploid.vcf', mode='r') as f:
    print(f.read())
{% endhighlight %}

    ##fileformat=VCFv4.3
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
    20	14370	.	G	A	.	.	.	GT	0/0/0/0	0/0/0/1	0/0/1/1
    20	17330	.	T	A,C,G	.	.	.	GT	1/1/2/2	0/1/2/3	3/3/3/3
    


Here is how to indicate the ploidy:


{% highlight python %}
callset = allel.read_vcf('example_polyploid.vcf', numbers={'GT': 4})
callset['calldata/GT']
{% endhighlight %}




    array([[[0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 1]],
    
           [[1, 1, 2, 2],
            [0, 1, 2, 3],
            [3, 3, 3, 3]]], dtype=int8)



As shown earlier for diploid calls, the `GenotypeArray` class can provide some extra functionality over a plain NumPy array, e.g.:


{% highlight python %}
gt = allel.GenotypeArray(callset['calldata/GT'])
gt
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(2, 3, 4) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0/0/0</td><td style="text-align: center">0/0/0/1</td><td style="text-align: center">0/0/1/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">1/1/2/2</td><td style="text-align: center">0/1/2/3</td><td style="text-align: center">3/3/3/3</td></tr></tbody></table></div>




{% highlight python %}
gt.is_het()
{% endhighlight %}




    array([[False,  True,  True],
           [ True,  True, False]], dtype=bool)




{% highlight python %}
ac = gt.count_alleles()
ac
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;AlleleCountsArray shape=(2, 4) dtype=int32&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">9</td><td style="text-align: center">3</td><td style="text-align: center">0</td><td style="text-align: center">0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">1</td><td style="text-align: center">3</td><td style="text-align: center">3</td><td style="text-align: center">5</td></tr></tbody></table></div>



### Region

You can extract data for only a specific chromosome or genome region via the **region** parameter. The value of the parameter should be a region string of the format '{chromosome}:{begin}-{end}', just like you would give to tabix or samtools. E.g.:


{% highlight python %}
callset = allel.read_vcf('example.vcf', region='20:1000000-1231000')
callset['variants/POS']
{% endhighlight %}




    array([1110696, 1230237], dtype=int32)



By default, scikit-allel will try to use tabix to extract the data from the requested region. If tabix is not available on your system or the VCF file has not been tabix indexed, scikit-allel will fall back to scanning through the VCF file from the beginning, which may be slow. If you have tabix installed but it is in a non-standard location, you can specify the path to the tabix executable via the **tabix** parameter.

Note that if you are using tabix, then tabix needs to be at least version 0.2.5. Some older versions of tabix do not support the "-h" option to output the headers from the VCF, which scikit-allel needs to get the meta-information for parsing. If you get an error message like "RuntimeError: VCF file is missing mandatory header line ("#CHROM...")" then check your tabix version and upgrade if necessary. If you have conda installed, a recent version tabix can be installed via the following command: ``conda install -c bioconda htslib``.

### Samples

You can extract data for only specific samples via the **samples** parameter. E.g., extract data for samples 'NA00001' and 'NA00003':


{% highlight python %}
callset = allel.read_vcf('example.vcf', samples=['NA00001', 'NA00003'])
callset['samples']
{% endhighlight %}




    array(['NA00001', 'NA00003'], dtype=object)




{% highlight python %}
allel.GenotypeArray(callset['calldata/GT'])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 2, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">1/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">2/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">./.</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">1/1</td></tr></tbody></table></div>



Note that the genotype array now only has two columns, corresponding to the two samples requested.

## [`vcf_to_npz()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_npz)

NumPy arrays are stored in main memory (a.k.a., RAM), which means that as soon as you end your Python session or restart your Jupyter notebook kernel, any data stored in a NumPy array will be lost. 

If your VCF file is not too big, you can extract data from the file into NumPy arrays then save those arrays to disk via the [`vcf_to_npz()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_npz) function. This function has most of the same parameters as the `read_vcf()` function, except that you also specify an output path, which is the name of the file you want to save the extracted data to. For example:


{% highlight python %}
allel.vcf_to_npz('example.vcf', 'example.npz', fields='*', overwrite=True)
{% endhighlight %}

This extraction only needs to be done once, then you can load the data any time directly into NumPy arrays via the NumPy [`load()`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.load.html) function, e.g.:


{% highlight python %}
import numpy as np
callset = np.load('example.npz')
callset
{% endhighlight %}




    <numpy.lib.npyio.NpzFile at 0x7f829b07d6d8>




{% highlight python %}
sorted(callset.keys())
{% endhighlight %}




    ['calldata/DP',
     'calldata/GT',
     'samples',
     'variants/AF',
     'variants/ALT',
     'variants/CHROM',
     'variants/DB',
     'variants/DP',
     'variants/FILTER_PASS',
     'variants/FILTER_q10',
     'variants/FILTER_s50',
     'variants/ID',
     'variants/POS',
     'variants/QUAL',
     'variants/REF',
     'variants/is_snp',
     'variants/numalt',
     'variants/svlen']




{% highlight python %}
callset['variants/POS']
{% endhighlight %}




    array([  14370,   17330, 1110696, 1230237, 1234567], dtype=int32)




{% highlight python %}
callset['calldata/GT']
{% endhighlight %}




    array([[[ 0,  0],
            [ 0,  1],
            [ 1,  1]],
    
           [[ 0,  0],
            [ 0,  1],
            [ 0,  0]],
    
           [[ 0,  2],
            [ 1,  2],
            [ 2,  2]],
    
           [[ 0,  0],
            [ 0,  0],
            [-1, -1]],
    
           [[ 0,  1],
            [ 0,  2],
            [ 1,  1]]], dtype=int8)



Note that although the data have been saved to disk, all of the data must be loaded into main memory first during the extraction process, so this function is not suitable if you have a dataset that is too large to fit into main memory.

## [`vcf_to_hdf5()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_hdf5)

For large datasets, the [`vcf_to_hdf5()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_hdf5) function is available. This function again takes similar parameters to `read_vcf()`, but will store extracted data into an HDF5 file stored on disk. The extraction process works through the VCF file in chunks, and so the entire dataset is never loaded entirely into main memory. A bit further below I give worked examples with a large dataset, but for now here is a simple example:


{% highlight python %}
allel.vcf_to_hdf5('example.vcf', 'example.h5', fields='*', overwrite=True)
{% endhighlight %}

The saved data can be accessed via the [`h5py`](http://www.h5py.org/) library, e.g.:


{% highlight python %}
import h5py
callset = h5py.File('example.h5', mode='r')
callset
{% endhighlight %}




    <HDF5 file "example.h5" (mode r)>



The one difference to be aware of here is that accessing data via a key like 'variants/POS' does not return a NumPy array, instead you get an HDF5 dataset object. 


{% highlight python %}
chrom = callset['variants/CHROM']
chrom
{% endhighlight %}




    <HDF5 dataset "CHROM": shape (5,), type "|O">




{% highlight python %}
pos = callset['variants/POS']
pos
{% endhighlight %}




    <HDF5 dataset "POS": shape (5,), type "<i4">




{% highlight python %}
gt = callset['calldata/GT']
gt
{% endhighlight %}




    <HDF5 dataset "GT": shape (5, 3, 2), type "|i1">



This dataset object is useful because you can then load all or only part of the underlying data into main memory via slicing. E.g.:


{% highlight python %}
# load second to fourth items into NumPy array
chrom[1:3]
{% endhighlight %}




    array(['20', '20'], dtype=object)




{% highlight python %}
# load second to fourth items into NumPy array
pos[1:3]
{% endhighlight %}




    array([  17330, 1110696], dtype=int32)




{% highlight python %}
# load all items into NumPy array
chrom[:]
{% endhighlight %}




    array(['20', '20', '20', '20', '20'], dtype=object)




{% highlight python %}
# load all items into NumPy array
pos[:]
{% endhighlight %}




    array([  14370,   17330, 1110696, 1230237, 1234567], dtype=int32)



This is particularly useful for the large data coming from the sample fields, e.g., the genotype calls. With these data you can use slicing to pull out particular rows or columns without having to load all data into memory. E.g.:


{% highlight python %}
# load genotype calls into memory for second to fourth variants, all samples
allel.GenotypeArray(gt[1:3, :])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(2, 3, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td><td style="text-align: center">2/2</td></tr></tbody></table></div>




{% highlight python %}
# load genotype calls into memory for all variants, first and second samples
allel.GenotypeArray(gt[:, 0:2])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 2, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">0/2</td></tr></tbody></table></div>




{% highlight python %}
# load all genotype calls into memory
allel.GenotypeArray(gt[:, :])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 3, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">1/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td><td style="text-align: center">2/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">./.</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">0/2</td><td style="text-align: center">1/1</td></tr></tbody></table></div>



### Storing strings

By default, scikit-allel will store data from string fields like CHROM, ID, REF and ALT in the HDF5 file as variable length strings. Unfortunately HDF5 uses a lot of space to store variable length strings, and so for larger datasets this can lead to very large HDF5 file sizes. You can disable this option by setting the **vlen** parameter to ``False``, which will force all data from string fields to be stored as fixed length strings. For example:


{% highlight python %}
allel.vcf_to_hdf5('example.vcf', 'example_novlen.h5', fields='*', overwrite=True, vlen=False)
{% endhighlight %}


{% highlight python %}
callset = h5py.File('example_novlen.h5', mode='r')
{% endhighlight %}


{% highlight python %}
callset['variants/ID']
{% endhighlight %}




    <HDF5 dataset "ID": shape (5,), type "|S9">




{% highlight python %}
callset['variants/ID'][:]
{% endhighlight %}




    array([b'rs6054257', b'.', b'rs6040355', b'.', b'microsat1'],
          dtype='|S9')



Note that an appropriate length to use for the fixed length string data type for each string field will be guessed from the longest string found for each field in the first chunk of the VCF. In some cases there may be longer strings found later in a VCF, in which case string values may get truncated. If you already know something about the longest value in a given field, you can specify **types** manually to avoid truncation, e.g.:


{% highlight python %}
allel.vcf_to_hdf5('example.vcf', 'example_novlen_types.h5', 
                  fields=['CHROM', 'ID'],
                  types={'CHROM': 'S2', 'ID': 'S10'},
                  overwrite=True, vlen=False)
{% endhighlight %}


{% highlight python %}
callset = h5py.File('example_novlen_types.h5', mode='r')
{% endhighlight %}


{% highlight python %}
callset['variants/ID']
{% endhighlight %}




    <HDF5 dataset "ID": shape (5,), type "|S10">



## [`vcf_to_zarr()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_zarr)

An alternative to HDF5 is a newer storage library called [Zarr](http://zarr.readthedocs.io/en/latest/). Currently there is only a Python implementation of Zarr, whereas you can access HDF5 files using a variety of different programming languages, so if you need portability then HDF5 is a better option. However, if you know you're only going to be using Python, then Zarr is a good option. In general it is a bit faster than HDF5, provides more storage and compression options, and also plays better with parellel computing libraries like [Dask](http://dask.pydata.org/en/latest/). 

To install Zarr, from the command line, run: `conda install -c conda-forge zarr`

To extract data to Zarr, for example:


{% highlight python %}
allel.vcf_to_zarr('example.vcf', 'example.zarr', fields='*', overwrite=True)
{% endhighlight %}

Then, to access saved data:


{% highlight python %}
import zarr
import numcodecs
print('zarr', zarr.__version__, 'numcodecs', numcodecs.__version__)
{% endhighlight %}

    zarr 2.2.0 numcodecs 0.5.5



{% highlight python %}
callset = zarr.open_group('example.zarr', mode='r')
callset
{% endhighlight %}




    <zarr.hierarchy.Group '/' read-only>



Zarr group objects provide a `tree()` method which can be useful for exploring the hierarchy of groups and arrays, e.g.:


{% highlight python %}
callset.tree(expand=True)
{% endhighlight %}




<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/themes/default/style.min.css"/><div id="af2fed46-a4fc-43fa-bbcf-0bc02ebd5d5b" class="zarr-tree"><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>/</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>calldata</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DP (5, 3) int16</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>GT (5, 3, 2) int8</span></li></ul></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>samples (3,) object</span></li><li data-jstree='{"type": "Group"}' class='jstree-open'><span>variants</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AF (5, 3) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ALT (5, 3) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CHROM (5,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DB (5,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DP (5,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_PASS (5,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_q10 (5,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_s50 (5,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ID (5,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>POS (5,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>QUAL (5,) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>REF (5,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>is_snp (5,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>numalt (5,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>svlen (5, 3) int32</span></li></ul></li></ul></li></ul></div>
<script>
    if (!require.defined('jquery')) {
        require.config({
            paths: {
                jquery: '//cdnjs.cloudflare.com/ajax/libs/jquery/1.12.1/jquery.min'
            },
        });
    }
    if (!require.defined('jstree')) {
        require.config({
            paths: {
                jstree: '//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/jstree.min'
            },
        });
    }
    require(['jstree'], function() {
        $('#af2fed46-a4fc-43fa-bbcf-0bc02ebd5d5b').jstree({
            types: {
                Group: {
                    icon: "fa fa-folder"
                },
                Array: {
                    icon: "fa fa-table"
                }
            },
            plugins: ["types"]
        });
    });
</script>




As with h5py, arrays can be accessed via a path-like notation, e.g.:


{% highlight python %}
chrom = callset['variants/CHROM']
chrom
{% endhighlight %}




    <zarr.core.Array '/variants/CHROM' (5,) object read-only>




{% highlight python %}
pos = callset['variants/POS']
pos
{% endhighlight %}




    <zarr.core.Array '/variants/POS' (5,) int32 read-only>




{% highlight python %}
gt = callset['calldata/GT']
gt
{% endhighlight %}




    <zarr.core.Array '/calldata/GT' (5, 3, 2) int8 read-only>



Again as with h5py, data can be loaded into memory as numpy arrays via slice notation:


{% highlight python %}
chrom[1:3]
{% endhighlight %}




    array(['20', '20'], dtype=object)




{% highlight python %}
pos[1:3]
{% endhighlight %}




    array([  17330, 1110696], dtype=int32)




{% highlight python %}
chrom[:]
{% endhighlight %}




    array(['20', '20', '20', '20', '20'], dtype=object)




{% highlight python %}
pos[:]
{% endhighlight %}




    array([  14370,   17330, 1110696, 1230237, 1234567], dtype=int32)




{% highlight python %}
allel.GenotypeArray(gt[1:3, :])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(2, 3, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td><td style="text-align: center">2/2</td></tr></tbody></table></div>




{% highlight python %}
allel.GenotypeArray(gt[:, 0:2])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 2, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">0/2</td></tr></tbody></table></div>




{% highlight python %}
allel.GenotypeArray(gt[:])
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeArray shape=(5, 3, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">1/1</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/1</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/2</td><td style="text-align: center">1/2</td><td style="text-align: center">2/2</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">3</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">./.</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">4</th><td style="text-align: center">0/1</td><td style="text-align: center">0/2</td><td style="text-align: center">1/1</td></tr></tbody></table></div>



## [`vcf_to_dataframe()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_dataframe)

For some analyses it can be useful to think of the data in a VCF file as a table or data frame, especially if you are only analysing data from the fixed fields and don't need the genotype calls or any other call data. The [`vcf_to_dataframe()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_dataframe) function extracts data from a VCF and loads into a [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html). E.g.:


{% highlight python %}
df = allel.vcf_to_dataframe('example.vcf')
df
{% endhighlight %}




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>ALT_2</th>
      <th>ALT_3</th>
      <th>QUAL</th>
      <th>FILTER_PASS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>20</td>
      <td>14370</td>
      <td>rs6054257</td>
      <td>G</td>
      <td>A</td>
      <td></td>
      <td></td>
      <td>29.0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>20</td>
      <td>17330</td>
      <td>.</td>
      <td>T</td>
      <td>A</td>
      <td></td>
      <td></td>
      <td>3.0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>20</td>
      <td>1110696</td>
      <td>rs6040355</td>
      <td>A</td>
      <td>G</td>
      <td>T</td>
      <td></td>
      <td>67.0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>20</td>
      <td>1230237</td>
      <td>.</td>
      <td>T</td>
      <td></td>
      <td></td>
      <td></td>
      <td>47.0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>20</td>
      <td>1234567</td>
      <td>microsat1</td>
      <td>GTC</td>
      <td>G</td>
      <td>GTCT</td>
      <td></td>
      <td>50.0</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



Note that the 'ALT' field has been broken into three separate columns, labelled 'ALT_1', 'ALT_2' and 'ALT_3'. When loading data into a data frame, any field with multiple values will be broken into multiple columns in this way.

Let's extract all fields, and also reduce number of values per field given that we know there are at most 2 alternate alleles in our example VCF file:


{% highlight python %}
df = allel.vcf_to_dataframe('example.vcf', fields='*', alt_number=2)
df
{% endhighlight %}




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>ALT_2</th>
      <th>QUAL</th>
      <th>DP</th>
      <th>AF_1</th>
      <th>AF_2</th>
      <th>DB</th>
      <th>FILTER_PASS</th>
      <th>FILTER_q10</th>
      <th>FILTER_s50</th>
      <th>numalt</th>
      <th>svlen_1</th>
      <th>svlen_2</th>
      <th>is_snp</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>20</td>
      <td>14370</td>
      <td>rs6054257</td>
      <td>G</td>
      <td>A</td>
      <td></td>
      <td>29.0</td>
      <td>14</td>
      <td>0.500</td>
      <td>NaN</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>20</td>
      <td>17330</td>
      <td>.</td>
      <td>T</td>
      <td>A</td>
      <td></td>
      <td>3.0</td>
      <td>11</td>
      <td>0.017</td>
      <td>NaN</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>20</td>
      <td>1110696</td>
      <td>rs6040355</td>
      <td>A</td>
      <td>G</td>
      <td>T</td>
      <td>67.0</td>
      <td>10</td>
      <td>0.333</td>
      <td>0.667</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>20</td>
      <td>1230237</td>
      <td>.</td>
      <td>T</td>
      <td></td>
      <td></td>
      <td>47.0</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>20</td>
      <td>1234567</td>
      <td>microsat1</td>
      <td>GTC</td>
      <td>G</td>
      <td>GTCT</td>
      <td>50.0</td>
      <td>9</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>2</td>
      <td>-2</td>
      <td>1</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



In case you were wondering, the 'numalt', 'svlen' and 'is_snp' fields are computed by scikit-allel, they are not present in the original VCF.

Pandas DataFrames have many useful features. For example, you can query:


{% highlight python %}
df.query('DP > 10 and QUAL > 20')
{% endhighlight %}




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>ALT_2</th>
      <th>QUAL</th>
      <th>DP</th>
      <th>AF_1</th>
      <th>AF_2</th>
      <th>DB</th>
      <th>FILTER_PASS</th>
      <th>FILTER_q10</th>
      <th>FILTER_s50</th>
      <th>numalt</th>
      <th>svlen_1</th>
      <th>svlen_2</th>
      <th>is_snp</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>20</td>
      <td>14370</td>
      <td>rs6054257</td>
      <td>G</td>
      <td>A</td>
      <td></td>
      <td>29.0</td>
      <td>14</td>
      <td>0.5</td>
      <td>NaN</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>20</td>
      <td>1230237</td>
      <td>.</td>
      <td>T</td>
      <td></td>
      <td></td>
      <td>47.0</td>
      <td>13</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



## [`vcf_to_recarray()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_recarray)

If you prefer to work with NumPy structured arrays rather than pandas DataFrames, try the [`vcf_to_recarray()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_recarray) function, e.g.: 


{% highlight python %}
ra = allel.vcf_to_recarray('example.vcf')
ra
{% endhighlight %}




    array([('20',   14370, 'rs6054257', 'G', 'A', '', '',  29.,  True),
           ('20',   17330, '.', 'T', 'A', '', '',   3., False),
           ('20', 1110696, 'rs6040355', 'A', 'G', 'T', '',  67.,  True),
           ('20', 1230237, '.', 'T', '', '', '',  47.,  True),
           ('20', 1234567, 'microsat1', 'GTC', 'G', 'GTCT', '',  50.,  True)],
          dtype=(numpy.record, [('CHROM', 'O'), ('POS', '<i4'), ('ID', 'O'), ('REF', 'O'), ('ALT_1', 'O'), ('ALT_2', 'O'), ('ALT_3', 'O'), ('QUAL', '<f4'), ('FILTER_PASS', '?')]))



## [`vcf_to_csv()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_csv)

Finally, the [`vcf_to_csv()`](http://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_csv) function is available if you need to dump data from a VCF file out to a generic CSV file (e.g., to load into a database). E.g.:


{% highlight python %}
allel.vcf_to_csv('example.vcf', 'example.csv', fields=['CHROM', 'POS', 'DP'])
{% endhighlight %}


{% highlight python %}
with open('example.csv', mode='r') as f:
    print(f.read())
{% endhighlight %}

    CHROM,POS,DP
    20,14370,14
    20,17330,11
    20,1110696,10
    20,1230237,13
    20,1234567,9
    


This function uses the pandas [`to_csv()`](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_csv.html) function under the hood to write the CSV, so you can control various output parameters (e.g., field separator) by passing through keyword arguments.

## Worked example: human 1000 genomes phase 3

I've downloaded a [VCF file](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz) with genotype data for Chromosome 22 from the 1000 genomes project phase 3 [FTP site](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 


{% highlight python %}
vcf_path = 'data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
{% endhighlight %}

Before we start processing, let's see how big the file is (the '!' is special Jupyter notebook syntax for running a command via the operating system shell):


{% highlight python %}
!ls -lh {vcf_path}
{% endhighlight %}

    -rw-r--r-- 1 aliman aliman 205M Jun 20  2017 data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


...and how many lines: 


{% highlight python %}
!zcat {vcf_path} | wc -l
{% endhighlight %}

    1103800


So there are more than a million variants.

### Preparation

When processing larger VCF files it's useful to get some feedback on how fast things are going. Let's import the `sys` module so we can log to standard output:


{% highlight python %}
import sys
{% endhighlight %}

Ultimately I am going to extract all the data from this VCF file into a Zarr store. However, before I do that, I'm going to check how many alternate alleles I should expect. I'm going to do that by extracting just the 'numalt' field, which scikit-allel will compute from the number of values in the 'ALT' field:


{% highlight python %}
callset = allel.read_vcf(vcf_path, fields=['numalt'], log=sys.stdout)
{% endhighlight %}

    [read_vcf] 65536 rows in 8.19s; chunk in 8.19s (8001 rows/s); 22 :18539397
    [read_vcf] 131072 rows in 16.05s; chunk in 7.86s (8341 rows/s); 22 :21016127
    [read_vcf] 196608 rows in 22.73s; chunk in 6.69s (9800 rows/s); 22 :23236362
    [read_vcf] 262144 rows in 27.52s; chunk in 4.78s (13696 rows/s); 22 :25227844
    [read_vcf] 327680 rows in 32.98s; chunk in 5.46s (12001 rows/s); 22 :27285434
    [read_vcf] 393216 rows in 38.07s; chunk in 5.09s (12879 rows/s); 22 :29572822
    [read_vcf] 458752 rows in 43.13s; chunk in 5.06s (12940 rows/s); 22 :31900536
    [read_vcf] 524288 rows in 47.94s; chunk in 4.80s (13640 rows/s); 22 :34069864
    [read_vcf] 589824 rows in 53.03s; chunk in 5.09s (12869 rows/s); 22 :36053392
    [read_vcf] 655360 rows in 59.10s; chunk in 6.07s (10804 rows/s); 22 :38088395
    [read_vcf] 720896 rows in 66.34s; chunk in 7.25s (9044 rows/s); 22 :40216200
    [read_vcf] 786432 rows in 72.77s; chunk in 6.43s (10192 rows/s); 22 :42597446
    [read_vcf] 851968 rows in 78.85s; chunk in 6.08s (10775 rows/s); 22 :44564263
    [read_vcf] 917504 rows in 84.40s; chunk in 5.54s (11823 rows/s); 22 :46390672
    [read_vcf] 983040 rows in 89.51s; chunk in 5.12s (12807 rows/s); 22 :48116697
    [read_vcf] 1048576 rows in 95.96s; chunk in 6.45s (10162 rows/s); 22 :49713436
    [read_vcf] 1103547 rows in 102.48s; chunk in 6.52s (8436 rows/s)
    [read_vcf] all done (10768 rows/s)


Let's see what the largest number of alternate alleles is:


{% highlight python %}
numalt = callset['variants/numalt']
np.max(numalt)
{% endhighlight %}




    8



Out of interest, how many variants are multi-allelic?


{% highlight python %}
count_numalt = np.bincount(numalt)
count_numalt
{% endhighlight %}




    array([      0, 1097199,    6073,     224,      38,       9,       3,
                 0,       1])




{% highlight python %}
n_multiallelic = np.sum(count_numalt[2:])
n_multiallelic
{% endhighlight %}




    6348



So there are only a very small number of multi-allelic variants (6,348), the vast majority (1,097,199) have just one alternate allele.

### Extract to Zarr

Now we know how many alternate alleles to expect, let's go ahead and extract everything out into a Zarr on-disk store. Time for a cup of tea:


{% highlight python %}
zarr_path = 'data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.zarr'
{% endhighlight %}


{% highlight python %}
allel.vcf_to_zarr(vcf_path, zarr_path, group='22',
                  fields='*', alt_number=8, log=sys.stdout,
                  compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))
{% endhighlight %}

    [vcf_to_zarr] 65536 rows in 15.14s; chunk in 15.14s (4328 rows/s); 22 :18539397
    [vcf_to_zarr] 131072 rows in 30.19s; chunk in 15.05s (4354 rows/s); 22 :21016127
    [vcf_to_zarr] 196608 rows in 42.07s; chunk in 11.88s (5515 rows/s); 22 :23236362
    [vcf_to_zarr] 262144 rows in 53.31s; chunk in 11.23s (5833 rows/s); 22 :25227844
    [vcf_to_zarr] 327680 rows in 64.61s; chunk in 11.30s (5800 rows/s); 22 :27285434
    [vcf_to_zarr] 393216 rows in 78.52s; chunk in 13.91s (4709 rows/s); 22 :29572822
    [vcf_to_zarr] 458752 rows in 92.61s; chunk in 14.09s (4650 rows/s); 22 :31900536
    [vcf_to_zarr] 524288 rows in 104.47s; chunk in 11.86s (5526 rows/s); 22 :34069864
    [vcf_to_zarr] 589824 rows in 117.29s; chunk in 12.82s (5111 rows/s); 22 :36053392
    [vcf_to_zarr] 655360 rows in 128.32s; chunk in 11.03s (5942 rows/s); 22 :38088395
    [vcf_to_zarr] 720896 rows in 139.75s; chunk in 11.43s (5735 rows/s); 22 :40216200
    [vcf_to_zarr] 786432 rows in 153.38s; chunk in 13.63s (4808 rows/s); 22 :42597446
    [vcf_to_zarr] 851968 rows in 165.66s; chunk in 12.28s (5336 rows/s); 22 :44564263
    [vcf_to_zarr] 917504 rows in 178.74s; chunk in 13.08s (5008 rows/s); 22 :46390672
    [vcf_to_zarr] 983040 rows in 190.38s; chunk in 11.64s (5630 rows/s); 22 :48116697
    [vcf_to_zarr] 1048576 rows in 201.35s; chunk in 10.96s (5977 rows/s); 22 :49713436
    [vcf_to_zarr] 1103547 rows in 211.84s; chunk in 10.49s (5238 rows/s)
    [vcf_to_zarr] all done (5202 rows/s)


By default Zarr compresses data using LZ4 via the Blosc compression library, but here I've decided to use Blosc with Zstandard instead of LZ4. LZ4 is faster but has slightly lower compression ratio, so is a good option if disk space is not a major issue and you are working from a storage device with high bandwidth (e.g., an SSD). Zstd is slightly slower but compresses a bit better.

Zarr stores the data in multiple files within a directory hierarchy. How much storage is required in total?


{% highlight python %}
!du -hs {zarr_path}
{% endhighlight %}

    140M	data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.zarr


To check things have worked as expected, let's run a few quick diagnostics:


{% highlight python %}
callset_h1k = zarr.open_group(zarr_path, mode='r')
callset_h1k
{% endhighlight %}




    <zarr.hierarchy.Group '/' read-only>




{% highlight python %}
callset_h1k.tree(expand=True)
{% endhighlight %}




<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/themes/default/style.min.css"/><div id="8a8162d9-3aa5-4110-acf7-80fb2ac0355f" class="zarr-tree"><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>/</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>22</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>calldata</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>GT (1103547, 2504, 2) int8</span></li></ul></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>samples (2504,) object</span></li><li data-jstree='{"type": "Group"}' class='jstree-open'><span>variants</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AA (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AC (1103547, 8) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AFR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ALT (1103547, 8) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AMR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CHROM (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIEND (1103547, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIPOS (1103547, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CS (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DP (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EAS_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>END (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EUR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EX_TARGET (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_PASS (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ID (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>IMPRECISE (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MC (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEINFO (1103547, 4) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEND (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MLEN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MSTART (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MULTI_ALLELIC (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>NS (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>POS (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>QUAL (1103547,) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>REF (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SAS_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVLEN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVTYPE (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>TSD (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>VT (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>is_snp (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>numalt (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>svlen (1103547, 8) int32</span></li></ul></li></ul></li></ul></li></ul></div>
<script>
    if (!require.defined('jquery')) {
        require.config({
            paths: {
                jquery: '//cdnjs.cloudflare.com/ajax/libs/jquery/1.12.1/jquery.min'
            },
        });
    }
    if (!require.defined('jstree')) {
        require.config({
            paths: {
                jstree: '//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/jstree.min'
            },
        });
    }
    require(['jstree'], function() {
        $('#8a8162d9-3aa5-4110-acf7-80fb2ac0355f').jstree({
            types: {
                Group: {
                    icon: "fa fa-folder"
                },
                Array: {
                    icon: "fa fa-table"
                }
            },
            plugins: ["types"]
        });
    });
</script>




Note here that I have loaded the data into a group named '22', i.e., the data are grouped by chromosome. This is typically how we organise data in HDF5 or Zarr files in our own work with mosquitoes.


{% highlight python %}
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
{% endhighlight %}


{% highlight python %}
pos = allel.SortedIndex(callset_h1k['22/variants/POS'])
pos
{% endhighlight %}




<div class="allel allel-DisplayAs1D"><span>&lt;SortedIndex shape=(1103547,) dtype=int32&gt;</span><table><thead><tr><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">1103542</th><th style="text-align: center">1103543</th><th style="text-align: center">1103544</th><th style="text-align: center">1103545</th><th style="text-align: center">1103546</th></tr></thead><tbody><tr><td style="text-align: center">16050075</td><td style="text-align: center">16050115</td><td style="text-align: center">16050213</td><td style="text-align: center">16050319</td><td style="text-align: center">16050527</td><td style="text-align: center">...</td><td style="text-align: center">51241342</td><td style="text-align: center">51241386</td><td style="text-align: center">51244163</td><td style="text-align: center">51244205</td><td style="text-align: center">51244237</td></tr></tbody></table></div>




{% highlight python %}
def plot_windowed_variant_density(pos, window_size, title=None):
    
    # setup windows 
    bins = np.arange(0, pos.max(), window_size)
    
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size
    
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
{% endhighlight %}


{% highlight python %}
plot_windowed_variant_density(pos, window_size=100000, title='Variant density')
{% endhighlight %}


![png](/assets/2017-06-14-read-vcf_files/2017-06-14-read-vcf_169_0.png)


When working with large genotype arrays, scikit-allel has a [`GenotypeDaskArray`](@@TODO) class that is like the [`GenotypeArray`](@@TODO) class we met earlier but can handle data stored on-disk in HDF5 or Zarr files and compute over them without loading all data into memory.   


{% highlight python %}
gt = allel.GenotypeDaskArray(callset_h1k['22/calldata/GT'])
gt
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;GenotypeDaskArray shape=(1103547, 2504, 2) dtype=int8&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">...</th><th style="text-align: center">2499</th><th style="text-align: center">2500</th><th style="text-align: center">2501</th><th style="text-align: center">2502</th><th style="text-align: center">2503</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">...</th><td style="text-align: center" colspan="12">...</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103544</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103545</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103546</th><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">...</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td><td style="text-align: center">0/0</td></tr></tbody></table></div>



The main difference with the Dask-backed arrays is that you need to explicitly call `compute()` to run a computation, e.g.:


{% highlight python %}
%%time
ac = gt.count_alleles(max_allele=8).compute()
{% endhighlight %}

    CPU times: user 28 s, sys: 651 ms, total: 28.7 s
    Wall time: 4.23 s



{% highlight python %}
ac
{% endhighlight %}




<div class="allel allel-DisplayAs2D"><span>&lt;AlleleCountsArray shape=(1103547, 9) dtype=int64&gt;</span><table><thead><tr><th></th><th style="text-align: center">0</th><th style="text-align: center">1</th><th style="text-align: center">2</th><th style="text-align: center">3</th><th style="text-align: center">4</th><th style="text-align: center">5</th><th style="text-align: center">6</th><th style="text-align: center">7</th><th style="text-align: center">8</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">5007</td><td style="text-align: center">   1</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">4976</td><td style="text-align: center">  32</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">4970</td><td style="text-align: center">  38</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">...</th><td style="text-align: center" colspan="10">...</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103544</th><td style="text-align: center">4969</td><td style="text-align: center">  39</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103545</th><td style="text-align: center">5007</td><td style="text-align: center">   1</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1103546</th><td style="text-align: center">4989</td><td style="text-align: center">  19</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td><td style="text-align: center">   0</td></tr></tbody></table></div>



## Further reading

If you have any questions, please send an email to the [scikit-allel mailing list](https://groups.google.com/forum/#!forum/scikit-allel). If you're wondering where to go next, some of these articles may be of interest:

* [Installing Python for data analysis](http://alimanfoo.github.io/2017/05/18/installing-python.html)
* [A tour of scikit-allel](http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html)
* [Fast PCA](http://alimanfoo.github.io/2015/09/28/fast-pca.html)
* [Estimating F<sub>ST</sub>](http://alimanfoo.github.io/2015/09/21/estimating-fst.html)
* [Mendelian transmission](http://alimanfoo.github.io/2017/02/14/mendelian-transmission.html)

Further documentation about scikit-allel is available from the [scikit-allel API docs](http://scikit-allel.readthedocs.io/en/latest/).

## Post-script: grouping by chromosome

In the worked example above, I had a VCF file with data from a single chromosome, and I used the "group" argument to extract the data into the Zarr store under a group named after the chromosome ("22"). This wasn't strictly necessary for the worked example, but if I then wanted to extract data for other chromosomes, grouping the data by chromosome provides a convenient way to keep the data organised.

For example, if I also wanted to extract data for Chromosome 21, I would do the following:


{% highlight python %}
vcf_path = 'data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

allel.vcf_to_zarr(vcf_path, zarr_path, group='21',
                  fields='*', alt_number=8, log=sys.stdout,
                  compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))
{% endhighlight %}

    [vcf_to_zarr] 65536 rows in 12.83s; chunk in 12.83s (5109 rows/s); 21 :15597733
    [vcf_to_zarr] 131072 rows in 23.83s; chunk in 11.01s (5953 rows/s); 21 :17688253
    [vcf_to_zarr] 196608 rows in 34.98s; chunk in 11.15s (5877 rows/s); 21 :19697865
    [vcf_to_zarr] 262144 rows in 46.24s; chunk in 11.26s (5821 rows/s); 21 :21745371
    [vcf_to_zarr] 327680 rows in 57.22s; chunk in 10.98s (5969 rows/s); 21 :23694335
    [vcf_to_zarr] 393216 rows in 68.53s; chunk in 11.31s (5795 rows/s); 21 :25609367
    [vcf_to_zarr] 458752 rows in 81.13s; chunk in 12.60s (5200 rows/s); 21 :27761870
    [vcf_to_zarr] 524288 rows in 93.15s; chunk in 12.02s (5450 rows/s); 21 :29885097
    [vcf_to_zarr] 589824 rows in 104.42s; chunk in 11.27s (5816 rows/s); 21 :32050808
    [vcf_to_zarr] 655360 rows in 115.55s; chunk in 11.13s (5888 rows/s); 21 :34251524
    [vcf_to_zarr] 720896 rows in 127.02s; chunk in 11.47s (5714 rows/s); 21 :36455104
    [vcf_to_zarr] 786432 rows in 142.15s; chunk in 15.13s (4330 rows/s); 21 :38563725
    [vcf_to_zarr] 851968 rows in 153.81s; chunk in 11.65s (5623 rows/s); 21 :40757642
    [vcf_to_zarr] 917504 rows in 165.23s; chunk in 11.42s (5737 rows/s); 21 :42768192
    [vcf_to_zarr] 983040 rows in 179.12s; chunk in 13.89s (4717 rows/s); 21 :44738104
    [vcf_to_zarr] 1048576 rows in 191.46s; chunk in 12.34s (5309 rows/s); 21 :46520145
    [vcf_to_zarr] 1105538 rows in 201.12s; chunk in 9.66s (5898 rows/s)
    [vcf_to_zarr] all done (5488 rows/s)


Note that the `vcf_path` has changed because I'm extracting from a different VCF file, but I'm using the same `zarr_path` as before because the data will be grouped by chromosome within the Zarr store.

Here's the new data hierarchy after extracting data for Chromosome 21:


{% highlight python %}
callset_h1k.tree(expand=True)
{% endhighlight %}




<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/themes/default/style.min.css"/><div id="3f621278-8095-4c16-95d8-344b69f1a19e" class="zarr-tree"><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>/</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>21</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>calldata</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>GT (1105538, 2504, 2) int8</span></li></ul></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>samples (2504,) object</span></li><li data-jstree='{"type": "Group"}' class='jstree-open'><span>variants</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AA (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AC (1105538, 8) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AFR_AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ALT (1105538, 8) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AMR_AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AN (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CHROM (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIEND (1105538, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIPOS (1105538, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CS (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DP (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EAS_AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>END (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EUR_AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EX_TARGET (1105538,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_PASS (1105538,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ID (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>IMPRECISE (1105538,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MC (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEINFO (1105538, 4) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEND (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MLEN (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MSTART (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MULTI_ALLELIC (1105538,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>NS (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>POS (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>QUAL (1105538,) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>REF (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SAS_AF (1105538, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVLEN (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVTYPE (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>TSD (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>VT (1105538,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>is_snp (1105538,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>numalt (1105538,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>svlen (1105538, 8) int32</span></li></ul></li></ul></li><li data-jstree='{"type": "Group"}' class='jstree-open'><span>22</span><ul><li data-jstree='{"type": "Group"}' class='jstree-open'><span>calldata</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>GT (1103547, 2504, 2) int8</span></li></ul></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>samples (2504,) object</span></li><li data-jstree='{"type": "Group"}' class='jstree-open'><span>variants</span><ul><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AA (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AC (1103547, 8) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AFR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ALT (1103547, 8) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AMR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>AN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CHROM (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIEND (1103547, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CIPOS (1103547, 2) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>CS (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>DP (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EAS_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>END (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EUR_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>EX_TARGET (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>FILTER_PASS (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>ID (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>IMPRECISE (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MC (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEINFO (1103547, 4) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MEND (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MLEN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MSTART (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>MULTI_ALLELIC (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>NS (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>POS (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>QUAL (1103547,) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>REF (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SAS_AF (1103547, 8) float32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVLEN (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>SVTYPE (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>TSD (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>VT (1103547,) object</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>is_snp (1103547,) bool</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>numalt (1103547,) int32</span></li><li data-jstree='{"type": "Array"}' class='jstree-open'><span>svlen (1103547, 8) int32</span></li></ul></li></ul></li></ul></li></ul></div>
<script>
    if (!require.defined('jquery')) {
        require.config({
            paths: {
                jquery: '//cdnjs.cloudflare.com/ajax/libs/jquery/1.12.1/jquery.min'
            },
        });
    }
    if (!require.defined('jstree')) {
        require.config({
            paths: {
                jstree: '//cdnjs.cloudflare.com/ajax/libs/jstree/3.3.3/jstree.min'
            },
        });
    }
    require(['jstree'], function() {
        $('#3f621278-8095-4c16-95d8-344b69f1a19e').jstree({
            types: {
                Group: {
                    icon: "fa fa-folder"
                },
                Array: {
                    icon: "fa fa-table"
                }
            },
            plugins: ["types"]
        });
    });
</script>




If you want to save time, you could run the extraction for each chromosome in parallel, it is fine to write to a Zarr store from multiple processes.

The 1000 genomes project has built a separate VCF file for each chromosome, but you may have a single VCF file with data for multiple chromosomes. In this case you can still group the data by chromosome in the Zarr output, but you need to use the `region` argument when doing the extraction, and the VCF file needs to be tabix indexed. E.g.:


{% highlight python %}
vcf_path = '/path/to/input.vcf.gz'  # tabix-indexed VCF with data for multiple chromosomes
zarr_path = '/path/to/output.zarr'
chromosomes = [str(i) for i in range(1, 23)]
for chrom in chromosomes:
    allel.vcf_to_zarr(vcf_path, zarr_path, group=chrom, region=chrom, ...)    
{% endhighlight %}

## Post-script: SNPEFF annotations

If you have a VCF with annotations from SNPEFF, you can post-process these to extract out data from within the ANN annotations. For example, here's an example VCF with SNPEFF annotations:


{% highlight python %}
with open('example_snpeff.vcf', mode='r') as f:
    print(f.read())
{% endhighlight %}

    ##fileformat=VCFv4.1
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##SnpEffVersion="4.1b (build 2015-02-13), by Pablo Cingolani"
    ##SnpEffCmd="SnpEff  agam4.2 -noStats -lof - "
    ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    2L	103	.	C	T	53.36	.	AC=4;ANN=T|intergenic_region|MODIFIER|AGAP004677|AGAP004677|intergenic_region|AGAP004677||||||||3000|
    2L	192	.	G	A	181.16	.	AC=8;
    2L	13513722	.	A	T,G	284.06	PASS	AC=1;ANN=T|missense_variant|MODERATE|AGAP005273|AGAP005273|transcript|AGAP005273-RA|VectorBase|1/4|n.17A>T|p.Asp6Val|17/4788|17/4788|6/1596||,G|synonymous_variant|LOW|AGAP005273|AGAP005273|transcript|AGAP005273-RA|VectorBase|1/4|n.17A>G|p.Asp6Asp|12/4768|12/4768|4/1592||
    


The main thing to note is that the values of the ANN field contain multiple sub-fields. E.g., in the first variant, the ANN value is "``T|intergenic_region|MODIFIER|AGAP004677|AGAP004677|intergenic_region|AGAP004677||||||||3000|``". The sub-fields are delimited by the pipe character ("``|``"). The first sub-field gives the allele to which the effect annotation applies ("T"); the second sub-field gives the effect annotation ("intergenic_region"); etc.

By default, scikit-allel parses the ANN field like any other string field, e.g.:


{% highlight python %}
callset = allel.read_vcf('example_snpeff.vcf', fields='ANN')
ann = callset['variants/ANN']
ann
{% endhighlight %}




    array([ 'T|intergenic_region|MODIFIER|AGAP004677|AGAP004677|intergenic_region|AGAP004677||||||||3000|',
           '',
           'T|missense_variant|MODERATE|AGAP005273|AGAP005273|transcript|AGAP005273-RA|VectorBase|1/4|n.17A>T|p.Asp6Val|17/4788|17/4788|6/1596||'], dtype=object)



However, like this the data aren't much use. If you want to access the sub-fields, scikit-allel provides a **transformers** parameter, which can be used like this:


{% highlight python %}
callset = allel.read_vcf('example_snpeff.vcf', fields='ANN', transformers=allel.ANNTransformer())
list(callset)
{% endhighlight %}




    ['variants/ANN_AA_length',
     'variants/ANN_AA_pos',
     'variants/ANN_Allele',
     'variants/ANN_Annotation',
     'variants/ANN_Annotation_Impact',
     'variants/ANN_CDS_length',
     'variants/ANN_CDS_pos',
     'variants/ANN_Distance',
     'variants/ANN_Feature_ID',
     'variants/ANN_Feature_Type',
     'variants/ANN_Gene_ID',
     'variants/ANN_Gene_Name',
     'variants/ANN_HGVS_c',
     'variants/ANN_HGVS_p',
     'variants/ANN_Rank',
     'variants/ANN_Transcript_BioType',
     'variants/ANN_cDNA_length',
     'variants/ANN_cDNA_pos']



Note that even though I only requested the "ANN" field, the resulting callset has a number of fields, each of which represents a sub-field from the ANN annotation. E.g.:


{% highlight python %}
callset['variants/ANN_Allele']
{% endhighlight %}




    array(['T', '', 'T'], dtype=object)




{% highlight python %}
callset['variants/ANN_Annotation']
{% endhighlight %}




    array(['intergenic_region', '', 'missense_variant'], dtype=object)



You can also use the **transformers** parameter with other functions, e.g.:


{% highlight python %}
df = allel.vcf_to_dataframe('example_snpeff.vcf', fields='ANN', transformers=allel.ANNTransformer())
df
{% endhighlight %}




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ANN_Allele</th>
      <th>ANN_Annotation</th>
      <th>ANN_Annotation_Impact</th>
      <th>ANN_Gene_Name</th>
      <th>ANN_Gene_ID</th>
      <th>ANN_Feature_Type</th>
      <th>ANN_Feature_ID</th>
      <th>ANN_Transcript_BioType</th>
      <th>ANN_Rank</th>
      <th>ANN_HGVS_c</th>
      <th>ANN_HGVS_p</th>
      <th>ANN_cDNA_pos</th>
      <th>ANN_cDNA_length</th>
      <th>ANN_CDS_pos</th>
      <th>ANN_CDS_length</th>
      <th>ANN_AA_pos</th>
      <th>ANN_AA_length</th>
      <th>ANN_Distance</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>T</td>
      <td>intergenic_region</td>
      <td>MODIFIER</td>
      <td>AGAP004677</td>
      <td>AGAP004677</td>
      <td>intergenic_region</td>
      <td>AGAP004677</td>
      <td></td>
      <td>-1</td>
      <td></td>
      <td></td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>3000</td>
    </tr>
    <tr>
      <th>1</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td>-1</td>
      <td></td>
      <td></td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
      <td>-1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>T</td>
      <td>missense_variant</td>
      <td>MODERATE</td>
      <td>AGAP005273</td>
      <td>AGAP005273</td>
      <td>transcript</td>
      <td>AGAP005273-RA</td>
      <td>VectorBase</td>
      <td>1</td>
      <td>17A&gt;T</td>
      <td>Asp6Val</td>
      <td>17</td>
      <td>4788</td>
      <td>17</td>
      <td>4788</td>
      <td>6</td>
      <td>1596</td>
      <td>-1</td>
    </tr>
  </tbody>
</table>
</div>



One thing to watch out for, some variants can have multiple SNPEFF annotations, i.e., there can be multiple values for the ANN annotation. This can happen if a variant has multiple alleles, or if a variant is within a gene that has multiple transcripts, or if a variant is within two overlapping genes (e.g., on opposite strands, or one gene within an intron of another), or if genes are closely spaced and so you get an upstream or downstream annotation for one gene and another annotation for the gene the variant is in, or some combination of any or all of these.

E.g., in the example above, the third variant has two alternate alleles, and therefore has two ANN values, and by default only the first value has been parsed. You can increase the number of values via the **numbers** parameter as usual, e.g.:


{% highlight python %}
callset = allel.read_vcf('example_snpeff.vcf', fields='ANN', numbers={'ANN': 2},
                         transformers=allel.ANNTransformer())
{% endhighlight %}


{% highlight python %}
callset['variants/ANN_Allele']
{% endhighlight %}




    array([['T', ''],
           ['', ''],
           ['T', 'G']], dtype=object)




{% highlight python %}
callset['variants/ANN_Annotation']
{% endhighlight %}




    array([['intergenic_region', ''],
           ['', ''],
           ['missense_variant', 'synonymous_variant']], dtype=object)



Dealing with multiple SNPEFF annotations can be tricky. FWIW my suggestion is, when generating the SNPEFF annotations, use the option to only call annotations against the canonical (longest) transcript for each gene, and disable upstream and downstream annotations. This will at least reduce some of the complexity, although you'll still need to handle multiallelic variants and the odd case of overlapping genes.

## Post-script: changes from `vcfnp`

The new functions available in `scikit-allel` supercede a package I previously wrote for extracting data from VCF files called [`vcfnp`](@@TODO). I rewrote this functionality from the ground up and ported the functionality to `scikit-allel` for two main reasons. Firstly, `vcfnp` was slow and so you needed a cluster to parse big VCF files, which is obviously a pain. The new functions in `scikit-allel` should be up to ~40 times faster. Secondly, the `vcfnp` API was somewhat complicated, requiring three separate steps to get data from VCF into an HDF5 file or Zarr store. The new functions in `scikit-allel` hopefully simplify this process, enabling data to be extracted from VCF and loaded into any of a variety of storage containers via a single function call.

If you previously used `vcfnp` here are a few notes on some of the things that have changed.

* No need for separate function calls to extract data from variants and calldata fields, both can be extracted via a single call to `read_vcf()` or any of the `vcf_to_...()` functions described above.
* Data can be extracted from VCF and loaded into HDF5 with a single function call to `vcf_to_hdf5()`; i.e., no need to first extract parts of the data out to .npy files then load into HDF5.
* No need to use a cluster or do any parallelisation, it should be possible to run `vcf_to_hdf5()` or `vcf_to_zarr()` on a whole VCF on a half-decent desktop or laptop computer, although big VCF files might take a couple of hours and require a reasonably large hard disk.
* The default NumPy data type for string fields has changed to use 'object' dtype, which means that strings of any length will be stored automatically (i.e., no need to configure separate dtypes for each string field) and there will be no truncation of long strings.
* Previously in `vcfnp` the genotype calls were extracted into a special field called 'genotype' separate from the 'GT' calldata field if requested. In `scikit-allel` the default behaviour is to parse the 'GT' field as a 3-dimensional integer array and return simply as 'calldata/GT'. If you really want to process the 'GT' field as a string then you can override this by setting the type for 'calldata/GT' to 'S3' or 'object'.
* The "arity" argument in `vcfnp` is instead called "numbers" in scikit-allel, to match better with the terminology used in VCF meta-information headers.
* Genotype ploidy is now specified via the "numbers" argument, there is no special "ploidy" argument.

