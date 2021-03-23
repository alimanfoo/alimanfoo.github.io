---
layout: post
title: Querying multidimensional data with xarray
---

<a href="https://colab.research.google.com/github/alimanfoo/alimanfoo.github.io/blob/master/_posts/2021-03-16-xarray-query.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

It felt great getting my [first pull request to the xarray project](https://github.com/pydata/xarray/pull/4984) merged this week. The PR adds a `query()` method to the `Dataset` class, which provides the ability to apply a selection to all arrays within a dataset, based on evaluating a query expression against one or more of those arrays. This is analogous to the `query()` method that the pandas `DataFrame` class provides, but for multidimensional data. In fact it leverages a lot on functionality already available in pandas and xarray, the [actual implementation](https://github.com/pydata/xarray/blob/fba449b135fa69477a48967d66cdf3a2de58c633/xarray/core/dataset.py#L7004) was straightforward.

What's the use case? In our mosquito population genomics work we have datasets comprising genotype calls at ~100 million variants in ~10,000 mosquito samples. However, not all of those variants are reliable, and we very often begin an analysis by selecting a subset of high quality variants, or variants above a given allele frequency, or some combination of multiple criteria. We also often want to analyse a subset of mosquito samples, e.g., samples coming from a particular country, collected in a given year. Here's an illustration with a dummy dataset:


{% highlight python %}
!pip install -q -U git+git://github.com/pydata/xarray.git
{% endhighlight %}


{% highlight python %}
import numpy as np
import xarray as xr
{% endhighlight %}


{% highlight python %}
# dimension names
DIM_VARIANTS = "variants"
DIM_SAMPLES = "samples"
DIM_PLOIDY = "ploidy"
DIM_ALLELES = "alleles"

# coordinate arrays
coords = dict(

    # genomic positions of the variants
    variant_position = ([DIM_VARIANTS], 
                        np.array([1, 3, 7, 12, 25])),

    # sample identifiers
    sample_id = ([DIM_SAMPLES], np.array(['foo', 'bar', 'baz', 'qux'])),

)

# data variables
data_vars = dict(

    # reference and alternate alleles (these are biallelic SNPs)
    variant_alleles = ([DIM_VARIANTS, DIM_ALLELES],
                       np.array([['A', 'T'], 
                                 ['C', 'A'], 
                                 ['G', 'T'], 
                                 ['A', 'G'], 
                                 ['C', 'T']])),
                 
    # variant average mapping quality
    variant_MQ = ([DIM_VARIANTS], 
                  np.array([45, 34, 12, 50, 55])),

    # variant quality by depth
    variant_QD = ([DIM_VARIANTS],
                  np.array([23.4, 3.2, 34.9, 7.6, 15.7])),

    # country of collection
    sample_country = ([DIM_SAMPLES], 
                      np.array(['Burkina Faso', 'Burkina Faso', 'Cameroon', 'Angola'])),
    
    # year of collection
    sample_year = ([DIM_SAMPLES],
                   np.array([2019, 2018, 2018, 2021])),

    # simulate some genotype calls
    call_genotype = ([DIM_VARIANTS, DIM_SAMPLES, DIM_PLOIDY],
                     np.random.randint(0, 2, size=(5, 4, 2))),
)

ds = xr.Dataset(data_vars=data_vars, coords=coords)
ds
{% endhighlight %}




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:           (alleles: 2, ploidy: 2, samples: 4, variants: 5)
Coordinates:
    variant_position  (variants) int64 1 3 7 12 25
    sample_id         (samples) &lt;U3 &#x27;foo&#x27; &#x27;bar&#x27; &#x27;baz&#x27; &#x27;qux&#x27;
Dimensions without coordinates: alleles, ploidy, samples, variants
Data variables:
    variant_alleles   (variants, alleles) &lt;U1 &#x27;A&#x27; &#x27;T&#x27; &#x27;C&#x27; &#x27;A&#x27; ... &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;
    variant_MQ        (variants) int64 45 34 12 50 55
    variant_QD        (variants) float64 23.4 3.2 34.9 7.6 15.7
    sample_country    (samples) &lt;U12 &#x27;Burkina Faso&#x27; &#x27;Burkina Faso&#x27; ... &#x27;Angola&#x27;
    sample_year       (samples) int64 2019 2018 2018 2021
    call_genotype     (variants, samples, ploidy) int64 0 0 0 0 1 ... 0 0 1 1 1</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-8587856c-745e-4c46-990f-d3b3d26eb7b0' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8587856c-745e-4c46-990f-d3b3d26eb7b0' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span>alleles</span>: 2</li><li><span>ploidy</span>: 2</li><li><span>samples</span>: 4</li><li><span>variants</span>: 5</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-a087e066-1cf4-430e-9889-c7a42f1140d4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-a087e066-1cf4-430e-9889-c7a42f1140d4' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_position</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 3 7 12 25</div><input id='attrs-493f3fb7-e81b-413e-8351-26063c733595' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-493f3fb7-e81b-413e-8351-26063c733595' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-da98a922-39fb-40c9-852e-bc0a17cd83e8' class='xr-var-data-in' type='checkbox'><label for='data-da98a922-39fb-40c9-852e-bc0a17cd83e8' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 1,  3,  7, 12, 25])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_id</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U3</div><div class='xr-var-preview xr-preview'>&#x27;foo&#x27; &#x27;bar&#x27; &#x27;baz&#x27; &#x27;qux&#x27;</div><input id='attrs-39c717d8-0322-4f84-b3d4-73a3e5470d2d' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-39c717d8-0322-4f84-b3d4-73a3e5470d2d' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b6214731-dbd4-4a70-84de-de87842dfee4' class='xr-var-data-in' type='checkbox'><label for='data-b6214731-dbd4-4a70-84de-de87842dfee4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;foo&#x27;, &#x27;bar&#x27;, &#x27;baz&#x27;, &#x27;qux&#x27;], dtype=&#x27;&lt;U3&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-8be012d9-f6cb-414e-ab8e-a08ef10c3a12' class='xr-section-summary-in' type='checkbox'  checked><label for='section-8be012d9-f6cb-414e-ab8e-a08ef10c3a12' class='xr-section-summary' >Data variables: <span>(6)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_alleles</span></div><div class='xr-var-dims'>(variants, alleles)</div><div class='xr-var-dtype'>&lt;U1</div><div class='xr-var-preview xr-preview'>&#x27;A&#x27; &#x27;T&#x27; &#x27;C&#x27; &#x27;A&#x27; ... &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;</div><input id='attrs-2839212c-e378-47aa-8d48-7eabc8547296' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-2839212c-e378-47aa-8d48-7eabc8547296' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-348ec4c2-52a5-4a7b-ad23-794412b5a1fe' class='xr-var-data-in' type='checkbox'><label for='data-348ec4c2-52a5-4a7b-ad23-794412b5a1fe' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[&#x27;A&#x27;, &#x27;T&#x27;],
       [&#x27;C&#x27;, &#x27;A&#x27;],
       [&#x27;G&#x27;, &#x27;T&#x27;],
       [&#x27;A&#x27;, &#x27;G&#x27;],
       [&#x27;C&#x27;, &#x27;T&#x27;]], dtype=&#x27;&lt;U1&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_MQ</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>45 34 12 50 55</div><input id='attrs-390b5bcd-e35b-4a39-8ed4-2109c9b18871' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-390b5bcd-e35b-4a39-8ed4-2109c9b18871' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-8222416e-896c-4d65-a2ba-8f7d0c0b11d0' class='xr-var-data-in' type='checkbox'><label for='data-8222416e-896c-4d65-a2ba-8f7d0c0b11d0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([45, 34, 12, 50, 55])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_QD</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>23.4 3.2 34.9 7.6 15.7</div><input id='attrs-e14cc0e4-51bd-4698-b9d0-a17ff28c224f' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-e14cc0e4-51bd-4698-b9d0-a17ff28c224f' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c556dfa9-6a07-4958-8f66-60dbb5d5c9df' class='xr-var-data-in' type='checkbox'><label for='data-c556dfa9-6a07-4958-8f66-60dbb5d5c9df' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([23.4,  3.2, 34.9,  7.6, 15.7])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_country</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U12</div><div class='xr-var-preview xr-preview'>&#x27;Burkina Faso&#x27; ... &#x27;Angola&#x27;</div><input id='attrs-57669cb9-d518-4bac-a7a5-02293a0e7dc8' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-57669cb9-d518-4bac-a7a5-02293a0e7dc8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-20372246-f485-426a-8cfe-539bbf1d59cb' class='xr-var-data-in' type='checkbox'><label for='data-20372246-f485-426a-8cfe-539bbf1d59cb' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;Burkina Faso&#x27;, &#x27;Burkina Faso&#x27;, &#x27;Cameroon&#x27;, &#x27;Angola&#x27;], dtype=&#x27;&lt;U12&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_year</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>2019 2018 2018 2021</div><input id='attrs-0007ceb6-d488-4dba-993c-9126ff2360e2' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-0007ceb6-d488-4dba-993c-9126ff2360e2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-331ae348-c6af-493c-ab03-e87636d89e69' class='xr-var-data-in' type='checkbox'><label for='data-331ae348-c6af-493c-ab03-e87636d89e69' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([2019, 2018, 2018, 2021])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>call_genotype</span></div><div class='xr-var-dims'>(variants, samples, ploidy)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>0 0 0 0 1 0 1 1 ... 1 0 1 0 0 1 1 1</div><input id='attrs-21dd61b1-a9c7-4439-8a17-d74de6c6a6aa' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-21dd61b1-a9c7-4439-8a17-d74de6c6a6aa' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-7cc95603-c98c-41cf-90a0-0e4c6a25b7a5' class='xr-var-data-in' type='checkbox'><label for='data-7cc95603-c98c-41cf-90a0-0e4c6a25b7a5' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[0, 0],
        [0, 0],
        [1, 0],
        [1, 1]],

       [[1, 1],
        [1, 1],
        [1, 0],
        [0, 0]],

       [[0, 1],
        [0, 0],
        [1, 1],
        [0, 0]],

       [[0, 0],
        [1, 0],
        [0, 0],
        [0, 1]],

       [[1, 0],
        [1, 0],
        [0, 1],
        [1, 1]]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-757672ad-6fa1-415c-b92f-a67379b4e35a' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-757672ad-6fa1-415c-b92f-a67379b4e35a' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



Now we can make a query to select variants. E.g., we could select variants where MQ and QD are above some threshold.


{% highlight python %}
ds.query(variants="variant_MQ > 30 and variant_QD > 5")
{% endhighlight %}




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:           (alleles: 2, ploidy: 2, samples: 4, variants: 3)
Coordinates:
    variant_position  (variants) int64 1 12 25
    sample_id         (samples) &lt;U3 &#x27;foo&#x27; &#x27;bar&#x27; &#x27;baz&#x27; &#x27;qux&#x27;
Dimensions without coordinates: alleles, ploidy, samples, variants
Data variables:
    variant_alleles   (variants, alleles) &lt;U1 &#x27;A&#x27; &#x27;T&#x27; &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;
    variant_MQ        (variants) int64 45 50 55
    variant_QD        (variants) float64 23.4 7.6 15.7
    sample_country    (samples) &lt;U12 &#x27;Burkina Faso&#x27; &#x27;Burkina Faso&#x27; ... &#x27;Angola&#x27;
    sample_year       (samples) int64 2019 2018 2018 2021
    call_genotype     (variants, samples, ploidy) int64 0 0 0 0 1 ... 0 0 1 1 1</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-37e19d06-d3a8-4f16-a229-b0ca39ae5241' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-37e19d06-d3a8-4f16-a229-b0ca39ae5241' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span>alleles</span>: 2</li><li><span>ploidy</span>: 2</li><li><span>samples</span>: 4</li><li><span>variants</span>: 3</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-2900d8f1-8c2e-40e6-a8fa-4d57dffedf75' class='xr-section-summary-in' type='checkbox'  checked><label for='section-2900d8f1-8c2e-40e6-a8fa-4d57dffedf75' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_position</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 12 25</div><input id='attrs-0461fea0-9b26-4e3e-9cc0-29914ae12f98' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-0461fea0-9b26-4e3e-9cc0-29914ae12f98' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-166aac85-4d8b-4235-9b22-29ad6609e839' class='xr-var-data-in' type='checkbox'><label for='data-166aac85-4d8b-4235-9b22-29ad6609e839' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 1, 12, 25])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_id</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U3</div><div class='xr-var-preview xr-preview'>&#x27;foo&#x27; &#x27;bar&#x27; &#x27;baz&#x27; &#x27;qux&#x27;</div><input id='attrs-08fba647-290e-44bc-9591-6cc43868081b' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-08fba647-290e-44bc-9591-6cc43868081b' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-9c13cc4e-b813-4c6c-b41f-34cc8bf41272' class='xr-var-data-in' type='checkbox'><label for='data-9c13cc4e-b813-4c6c-b41f-34cc8bf41272' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;foo&#x27;, &#x27;bar&#x27;, &#x27;baz&#x27;, &#x27;qux&#x27;], dtype=&#x27;&lt;U3&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-706ff6e6-adef-4ab1-8f5c-3ba305b833fa' class='xr-section-summary-in' type='checkbox'  checked><label for='section-706ff6e6-adef-4ab1-8f5c-3ba305b833fa' class='xr-section-summary' >Data variables: <span>(6)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_alleles</span></div><div class='xr-var-dims'>(variants, alleles)</div><div class='xr-var-dtype'>&lt;U1</div><div class='xr-var-preview xr-preview'>&#x27;A&#x27; &#x27;T&#x27; &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;</div><input id='attrs-6fca0e5c-2485-4f95-8eeb-ec37ece7e832' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-6fca0e5c-2485-4f95-8eeb-ec37ece7e832' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-49584f2b-8d34-463d-89a4-f1c22015fc7d' class='xr-var-data-in' type='checkbox'><label for='data-49584f2b-8d34-463d-89a4-f1c22015fc7d' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[&#x27;A&#x27;, &#x27;T&#x27;],
       [&#x27;A&#x27;, &#x27;G&#x27;],
       [&#x27;C&#x27;, &#x27;T&#x27;]], dtype=&#x27;&lt;U1&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_MQ</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>45 50 55</div><input id='attrs-f7757ad2-b5c2-45be-9dc0-bb291e545f2e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-f7757ad2-b5c2-45be-9dc0-bb291e545f2e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-0855d4cd-3348-4bb5-b05d-6306245a895a' class='xr-var-data-in' type='checkbox'><label for='data-0855d4cd-3348-4bb5-b05d-6306245a895a' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([45, 50, 55])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_QD</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>23.4 7.6 15.7</div><input id='attrs-a0d792d6-cff0-4acb-a129-b4d61639f219' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-a0d792d6-cff0-4acb-a129-b4d61639f219' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2a7a2fbc-5b3d-497b-b629-fdb111e0376b' class='xr-var-data-in' type='checkbox'><label for='data-2a7a2fbc-5b3d-497b-b629-fdb111e0376b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([23.4,  7.6, 15.7])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_country</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U12</div><div class='xr-var-preview xr-preview'>&#x27;Burkina Faso&#x27; ... &#x27;Angola&#x27;</div><input id='attrs-d3f90293-c6eb-4547-a036-2368e553b737' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-d3f90293-c6eb-4547-a036-2368e553b737' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-e3d0baf0-049d-471f-a660-38b620782f52' class='xr-var-data-in' type='checkbox'><label for='data-e3d0baf0-049d-471f-a660-38b620782f52' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;Burkina Faso&#x27;, &#x27;Burkina Faso&#x27;, &#x27;Cameroon&#x27;, &#x27;Angola&#x27;], dtype=&#x27;&lt;U12&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_year</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>2019 2018 2018 2021</div><input id='attrs-0e79198a-4c90-49b1-959e-4bba3709d129' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-0e79198a-4c90-49b1-959e-4bba3709d129' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-6d97f28d-fe41-435a-a245-42da82a8ea70' class='xr-var-data-in' type='checkbox'><label for='data-6d97f28d-fe41-435a-a245-42da82a8ea70' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([2019, 2018, 2018, 2021])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>call_genotype</span></div><div class='xr-var-dims'>(variants, samples, ploidy)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>0 0 0 0 1 0 1 1 ... 1 0 1 0 0 1 1 1</div><input id='attrs-64187f05-02ea-4b23-a2ca-b7836f915be8' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-64187f05-02ea-4b23-a2ca-b7836f915be8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2cc42444-21bb-46ae-af10-7a4b4d359a94' class='xr-var-data-in' type='checkbox'><label for='data-2cc42444-21bb-46ae-af10-7a4b4d359a94' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[0, 0],
        [0, 0],
        [1, 0],
        [1, 1]],

       [[0, 0],
        [1, 0],
        [0, 0],
        [0, 1]],

       [[1, 0],
        [1, 0],
        [0, 1],
        [1, 1]]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-b9e4d840-3b00-4221-8bbd-7d372937af0c' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-b9e4d840-3b00-4221-8bbd-7d372937af0c' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



We could also select samples, e.g., collected in Burkina Faso since 2018.


{% highlight python %}
ds.query(samples="sample_country == 'Burkina Faso' and sample_year > 2018")
{% endhighlight %}




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:           (alleles: 2, ploidy: 2, samples: 1, variants: 5)
Coordinates:
    variant_position  (variants) int64 1 3 7 12 25
    sample_id         (samples) &lt;U3 &#x27;foo&#x27;
Dimensions without coordinates: alleles, ploidy, samples, variants
Data variables:
    variant_alleles   (variants, alleles) &lt;U1 &#x27;A&#x27; &#x27;T&#x27; &#x27;C&#x27; &#x27;A&#x27; ... &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;
    variant_MQ        (variants) int64 45 34 12 50 55
    variant_QD        (variants) float64 23.4 3.2 34.9 7.6 15.7
    sample_country    (samples) &lt;U12 &#x27;Burkina Faso&#x27;
    sample_year       (samples) int64 2019
    call_genotype     (variants, samples, ploidy) int64 0 0 1 1 0 1 0 0 1 0</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-b3a8e1c2-2102-4f1d-9829-3a7ef5a07c72' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-b3a8e1c2-2102-4f1d-9829-3a7ef5a07c72' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span>alleles</span>: 2</li><li><span>ploidy</span>: 2</li><li><span>samples</span>: 1</li><li><span>variants</span>: 5</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-c76fa396-597e-4a6b-81c9-49b926ab97b6' class='xr-section-summary-in' type='checkbox'  checked><label for='section-c76fa396-597e-4a6b-81c9-49b926ab97b6' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_position</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 3 7 12 25</div><input id='attrs-769dd983-3724-46b3-a3e2-7b7bcb499ba6' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-769dd983-3724-46b3-a3e2-7b7bcb499ba6' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f141270d-3fbb-436a-9abe-a1daca02e4a8' class='xr-var-data-in' type='checkbox'><label for='data-f141270d-3fbb-436a-9abe-a1daca02e4a8' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 1,  3,  7, 12, 25])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_id</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U3</div><div class='xr-var-preview xr-preview'>&#x27;foo&#x27;</div><input id='attrs-57945b86-b8f3-4ab6-8a77-f900c87d0630' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-57945b86-b8f3-4ab6-8a77-f900c87d0630' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-dc6455c9-150d-4c09-84ba-d4e7088b38be' class='xr-var-data-in' type='checkbox'><label for='data-dc6455c9-150d-4c09-84ba-d4e7088b38be' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;foo&#x27;], dtype=&#x27;&lt;U3&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-41955a2d-9714-488c-88ec-66fc50284ffa' class='xr-section-summary-in' type='checkbox'  checked><label for='section-41955a2d-9714-488c-88ec-66fc50284ffa' class='xr-section-summary' >Data variables: <span>(6)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_alleles</span></div><div class='xr-var-dims'>(variants, alleles)</div><div class='xr-var-dtype'>&lt;U1</div><div class='xr-var-preview xr-preview'>&#x27;A&#x27; &#x27;T&#x27; &#x27;C&#x27; &#x27;A&#x27; ... &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;</div><input id='attrs-cd9db239-c104-420d-aaea-5f82157ff1b1' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-cd9db239-c104-420d-aaea-5f82157ff1b1' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-cecb9349-2041-44dc-91cc-b4f97f56626a' class='xr-var-data-in' type='checkbox'><label for='data-cecb9349-2041-44dc-91cc-b4f97f56626a' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[&#x27;A&#x27;, &#x27;T&#x27;],
       [&#x27;C&#x27;, &#x27;A&#x27;],
       [&#x27;G&#x27;, &#x27;T&#x27;],
       [&#x27;A&#x27;, &#x27;G&#x27;],
       [&#x27;C&#x27;, &#x27;T&#x27;]], dtype=&#x27;&lt;U1&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_MQ</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>45 34 12 50 55</div><input id='attrs-4780977e-a144-4089-a2f9-70432a0d1734' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-4780977e-a144-4089-a2f9-70432a0d1734' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-6dbcab33-696b-4c96-9f12-95e2a9deb5bf' class='xr-var-data-in' type='checkbox'><label for='data-6dbcab33-696b-4c96-9f12-95e2a9deb5bf' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([45, 34, 12, 50, 55])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_QD</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>23.4 3.2 34.9 7.6 15.7</div><input id='attrs-275b9521-f711-4683-a010-ef6a8555edc8' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-275b9521-f711-4683-a010-ef6a8555edc8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2e73da11-28b5-480d-aa9e-a2638a5ebc77' class='xr-var-data-in' type='checkbox'><label for='data-2e73da11-28b5-480d-aa9e-a2638a5ebc77' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([23.4,  3.2, 34.9,  7.6, 15.7])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_country</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U12</div><div class='xr-var-preview xr-preview'>&#x27;Burkina Faso&#x27;</div><input id='attrs-8cded44a-c3ee-46ed-8a8e-9397e2bb7215' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-8cded44a-c3ee-46ed-8a8e-9397e2bb7215' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f4c0e82c-54fa-4b68-9d14-8584d32ffb93' class='xr-var-data-in' type='checkbox'><label for='data-f4c0e82c-54fa-4b68-9d14-8584d32ffb93' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;Burkina Faso&#x27;], dtype=&#x27;&lt;U12&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_year</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>2019</div><input id='attrs-eec8fca2-de47-4f2c-a29a-8fd72c6b3cc6' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-eec8fca2-de47-4f2c-a29a-8fd72c6b3cc6' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-d7bbf966-c468-4255-84da-7ea1153ca1bf' class='xr-var-data-in' type='checkbox'><label for='data-d7bbf966-c468-4255-84da-7ea1153ca1bf' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([2019])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>call_genotype</span></div><div class='xr-var-dims'>(variants, samples, ploidy)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>0 0 1 1 0 1 0 0 1 0</div><input id='attrs-2192cbf8-362a-4ff3-b681-d049981e1bff' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-2192cbf8-362a-4ff3-b681-d049981e1bff' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-152facbc-09aa-4d60-a36c-7d01cf12722d' class='xr-var-data-in' type='checkbox'><label for='data-152facbc-09aa-4d60-a36c-7d01cf12722d' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[0, 0]],

       [[1, 1]],

       [[0, 1]],

       [[0, 0]],

       [[1, 0]]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-1580c1d1-e5b0-4bfa-ae11-467252235ec7' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-1580c1d1-e5b0-4bfa-ae11-467252235ec7' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



Also, both of these queries could be applied at the same time.


{% highlight python %}
ds.query(
    variants="variant_MQ > 30 and variant_QD > 5",
    samples="sample_country == 'Burkina Faso' and sample_year > 2018"
)
{% endhighlight %}




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:           (alleles: 2, ploidy: 2, samples: 1, variants: 3)
Coordinates:
    variant_position  (variants) int64 1 12 25
    sample_id         (samples) &lt;U3 &#x27;foo&#x27;
Dimensions without coordinates: alleles, ploidy, samples, variants
Data variables:
    variant_alleles   (variants, alleles) &lt;U1 &#x27;A&#x27; &#x27;T&#x27; &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;
    variant_MQ        (variants) int64 45 50 55
    variant_QD        (variants) float64 23.4 7.6 15.7
    sample_country    (samples) &lt;U12 &#x27;Burkina Faso&#x27;
    sample_year       (samples) int64 2019
    call_genotype     (variants, samples, ploidy) int64 0 0 0 0 1 0</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-6aeb696e-679f-40a3-9ee5-d27cf162c708' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-6aeb696e-679f-40a3-9ee5-d27cf162c708' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span>alleles</span>: 2</li><li><span>ploidy</span>: 2</li><li><span>samples</span>: 1</li><li><span>variants</span>: 3</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-541316a6-fe0f-40d9-a215-4cc844eea884' class='xr-section-summary-in' type='checkbox'  checked><label for='section-541316a6-fe0f-40d9-a215-4cc844eea884' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_position</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 12 25</div><input id='attrs-5a73d1a5-ce4e-475c-83b0-cc3b2e1ba512' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-5a73d1a5-ce4e-475c-83b0-cc3b2e1ba512' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-8dd50e4c-3066-4b90-acc6-40ca451f1183' class='xr-var-data-in' type='checkbox'><label for='data-8dd50e4c-3066-4b90-acc6-40ca451f1183' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 1, 12, 25])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_id</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U3</div><div class='xr-var-preview xr-preview'>&#x27;foo&#x27;</div><input id='attrs-04609f2b-3f5c-4611-bb89-5d2587694b50' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-04609f2b-3f5c-4611-bb89-5d2587694b50' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-19fe7100-a1f2-4afd-882a-9f886bed10d0' class='xr-var-data-in' type='checkbox'><label for='data-19fe7100-a1f2-4afd-882a-9f886bed10d0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;foo&#x27;], dtype=&#x27;&lt;U3&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-e2054450-5937-4f9b-a449-da9ba8133ca1' class='xr-section-summary-in' type='checkbox'  checked><label for='section-e2054450-5937-4f9b-a449-da9ba8133ca1' class='xr-section-summary' >Data variables: <span>(6)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>variant_alleles</span></div><div class='xr-var-dims'>(variants, alleles)</div><div class='xr-var-dtype'>&lt;U1</div><div class='xr-var-preview xr-preview'>&#x27;A&#x27; &#x27;T&#x27; &#x27;A&#x27; &#x27;G&#x27; &#x27;C&#x27; &#x27;T&#x27;</div><input id='attrs-dfbb9508-7c90-4cd9-812b-5a2641a32e92' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-dfbb9508-7c90-4cd9-812b-5a2641a32e92' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-ab30b182-7af2-44e0-9a13-e57ccdd548f0' class='xr-var-data-in' type='checkbox'><label for='data-ab30b182-7af2-44e0-9a13-e57ccdd548f0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[&#x27;A&#x27;, &#x27;T&#x27;],
       [&#x27;A&#x27;, &#x27;G&#x27;],
       [&#x27;C&#x27;, &#x27;T&#x27;]], dtype=&#x27;&lt;U1&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_MQ</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>45 50 55</div><input id='attrs-134c6afb-3ced-4dd1-bedc-0e42a3276b65' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-134c6afb-3ced-4dd1-bedc-0e42a3276b65' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-6eb9548b-9118-4fe0-b6b0-01d37ac0d8d4' class='xr-var-data-in' type='checkbox'><label for='data-6eb9548b-9118-4fe0-b6b0-01d37ac0d8d4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([45, 50, 55])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>variant_QD</span></div><div class='xr-var-dims'>(variants)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>23.4 7.6 15.7</div><input id='attrs-b63b357f-c855-491a-9614-f8c62318cab3' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-b63b357f-c855-491a-9614-f8c62318cab3' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-82d83628-49a8-4dbd-84f0-e8e42fe2009b' class='xr-var-data-in' type='checkbox'><label for='data-82d83628-49a8-4dbd-84f0-e8e42fe2009b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([23.4,  7.6, 15.7])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_country</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>&lt;U12</div><div class='xr-var-preview xr-preview'>&#x27;Burkina Faso&#x27;</div><input id='attrs-4b968d36-f5b5-4201-bbdc-49376a67586f' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-4b968d36-f5b5-4201-bbdc-49376a67586f' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-cbca9fed-3881-4cf8-a15f-68e49911f2e0' class='xr-var-data-in' type='checkbox'><label for='data-cbca9fed-3881-4cf8-a15f-68e49911f2e0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;Burkina Faso&#x27;], dtype=&#x27;&lt;U12&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>sample_year</span></div><div class='xr-var-dims'>(samples)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>2019</div><input id='attrs-54bc7614-2268-4afa-974b-10bf8071ed69' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-54bc7614-2268-4afa-974b-10bf8071ed69' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-cc3e4df3-3bbc-46de-8206-c6b7fb877cbb' class='xr-var-data-in' type='checkbox'><label for='data-cc3e4df3-3bbc-46de-8206-c6b7fb877cbb' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([2019])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>call_genotype</span></div><div class='xr-var-dims'>(variants, samples, ploidy)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>0 0 0 0 1 0</div><input id='attrs-71194e5d-3c32-4269-ad90-11b04d5a32f4' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-71194e5d-3c32-4269-ad90-11b04d5a32f4' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-30f3a5b4-9d7c-4286-a3b6-5c908c4a4cd4' class='xr-var-data-in' type='checkbox'><label for='data-30f3a5b4-9d7c-4286-a3b6-5c908c4a4cd4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[0, 0]],

       [[0, 0]],

       [[1, 0]]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-1c6011ee-8feb-4a36-8a5e-82034429f148' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-1c6011ee-8feb-4a36-8a5e-82034429f148' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



Of course the same thing can be achieved by other means using the existing `isel()` method, but I've found having this query functionality can be useful for several reasons. It is concise and relatively easy to explain to newer users. It is also convenient to be able to have the query as a string, as this can be stored in metadata files if you ever need to keep a record of queries applied for different analyses.

I'm looking forward to using this new feature, and making more use of xarray generally, it's a great package for managing scientific data of all different shapes and sizes. [The xarray docs are here](http://xarray.pydata.org/en/stable/) if you'd like to learn more about it.


{% highlight python %}

{% endhighlight %}
