---
layout: post
title: CPU blues
---

<script src="https://cdn.pydata.org/bokeh/release/bokeh-0.11.1.min.js"></script>
<script src="https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.11.1.min.js"></script>
<script src="https://cdn.pydata.org/bokeh/release/bokeh-compiler-0.11.1.min.js"></script>
<link href="https://cdn.pydata.org/bokeh/release/bokeh-0.11.1.min.css" rel="stylesheet" type="text/css">
<link href="https://cdn.pydata.org/bokeh/release/bokeh-widgets-0.11.1.min.css" rel="stylesheet" type="text/css">

*Notes on reducing I/O bottlenecks in parallel computations with [Dask](http://dask.pydata.org/en/latest/) and [Zarr](http://zarr.readthedocs.io/en/latest/).*

When I'm analysing data I tend to keep one eye on the system monitor at the top of my screen. The height of the green bar tells me how much RAM I'm using and the height of the blue bar tells me how much CPU...

<p><image src='/assets/Selection_058.png' alt='status bar'/></p>

I'll kick off a computation then look up to see what's happening. Most of the time the blue bar chugs away at 1/8 height, which means that 1 of the 8 logical cores on my computer is fully utilised and the others are idle. I spend a fair amount of time waiting for computations to finish and it started to bug me that 7/8 CPU capacity was going spare. Many of the computations I run are simple and could be parallelized, so I started looking into ways of adapting my analysis code to make better use of multiple CPUs.

## Dask + HDF5

One solution I really like is [Dask](http://dask.pydata.org). Using [``dask.array``](http://dask.pydata.org/en/latest/array.html) it's simple to parallelize a wide range of operations over numerical arrays, using either multiple threads or multiple processes. To evaluate Dask I wrote an alternative [Dask-backed implementation](http://scikit-allel.readthedocs.io/en/latest/model/dask.html) of some of the basic genetic data transformations I use every day. I then ran these on some [data from the Ag1000G project](https://www.malariagen.net/data/ag1000g-phase1-ar3) using Dask's multi-threaded scheduler, hoping to see my idle CPU fire up to maximum capacity...


{% highlight python %}
import h5py
import multiprocessing
import dask; print('dask', dask.__version__)
import dask.array as da
from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
from dask.diagnostics.profile_visualize import visualize
from cachey import nbytes
import allel; print('scikit-allel', allel.__version__)
import bokeh
from bokeh.io import output_notebook
output_notebook()
from functools import reduce
import operator
{% endhighlight %}

    dask 0.11.1
    scikit-allel 1.0.3



{% highlight python %}
# Setup input data, using data on local SSD
callset = h5py.File('data/2016-05-16/ag1000g.phase1.ar3.pass.3R.h5', mode='r')
genotype = callset['3R/calldata/genotype']
genotype
{% endhighlight %}




    <HDF5 dataset "genotype": shape (13167162, 765, 2), type "|i1">




{% highlight python %}
# check how the data were compressed
genotype.compression, genotype.compression_opts
{% endhighlight %}




    ('gzip', 3)




{% highlight python %}
# how many cores on this computer?
multiprocessing.cpu_count()
{% endhighlight %}




    8




{% highlight python %}
# when I made the HDF5 file I set the chunks too small, so let's operate on bigger chunks
chunks = (genotype.chunks[0], genotype.chunks[1] * 20, genotype.chunks[2])
chunks
{% endhighlight %}




    (6553, 200, 2)




{% highlight python %}
import humanize
import numpy as np
print('chunk size:', humanize.naturalsize(np.product(chunks)))
{% endhighlight %}

    chunk size: 2.6 MB



{% highlight python %}
!cat /usr/local/bin/drop_caches
{% endhighlight %}

    #!/bin/bash
    # This must be run as sudo, so to avoid passwords this 
    # script should be set in the sudoers file with NOPASSWD.
    echo 1 > /proc/sys/vm/drop_caches
    



{% highlight python %}
# ensure OS page cache is cleared 
!sudo drop_caches
{% endhighlight %}


{% highlight python %}
# run the allele count computation via Dask
gd = allel.GenotypeDaskArray(genotype, chunks=chunks)
ac = gd.count_alleles(max_allele=3)
with ResourceProfiler(dt=1) as rprof:
    ac.compute()
visualize(rprof);
{% endhighlight %}




<div class="plotdiv" id="19b1b3ab-cc84-4e20-8ced-a836c3f430cc"></div>
<script type="text/javascript">
  
  (function(global) {
    function now() {
      return new Date();
    }
  
    if (typeof (window._bokeh_onload_callbacks) === "undefined") {
      window._bokeh_onload_callbacks = [];
    }
  
    function run_callbacks() {
      window._bokeh_onload_callbacks.forEach(function(callback) { callback() });
      delete window._bokeh_onload_callbacks
      console.info("Bokeh: all callbacks have finished");
    }
  
    function load_libs(js_urls, callback) {
      window._bokeh_onload_callbacks.push(callback);
      if (window._bokeh_is_loading > 0) {
        console.log("Bokeh: BokehJS is being loaded, scheduling callback at", now());
        return null;
      }
      if (js_urls == null || js_urls.length === 0) {
        run_callbacks();
        return null;
      }
      console.log("Bokeh: BokehJS not loaded, scheduling load and callback at", now());
      window._bokeh_is_loading = js_urls.length;
      for (var i = 0; i < js_urls.length; i++) {
        var url = js_urls[i];
        var s = document.createElement('script');
        s.src = url;
        s.async = false;
        s.onreadystatechange = s.onload = function() {
          window._bokeh_is_loading--;
          if (window._bokeh_is_loading === 0) {
            console.log("Bokeh: all BokehJS libraries loaded");
            run_callbacks()
          }
        };
        s.onerror = function() {
          console.warn("failed to load library " + url);
        };
        console.log("Bokeh: injecting script tag for BokehJS library: ", url);
        document.getElementsByTagName("head")[0].appendChild(s);
      }
    };var element = document.getElementById("19b1b3ab-cc84-4e20-8ced-a836c3f430cc");
    if (element == null) {
      console.log("Bokeh: ERROR: autoload.js configured with elementid '19b1b3ab-cc84-4e20-8ced-a836c3f430cc' but no matching script tag was found. ")
      return false;
    }
  
    var js_urls = [];
  
    var inline_js = [
      function(Bokeh) {
        Bokeh.$(function() {
            var docs_json = {"5fa95652-5e47-49ac-a58d-a27ba9d352aa":{"roots":{"references":[{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[0.0,103.6,133.8,133.8,132.9,132.8,132.8,132.8,133.8,133.8,132.8,131.8,133.8,131.9,133.8,133.8,131.8,133.9,133.9,132.8,132.8,132.9,132.9,132.8,131.8,132.9,133.8,131.9,133.8,131.9,131.9,133.8,132.8,132.9,132.8,132.8,132.8,132.8,133.8,132.8,133.9,131.8,133.8,131.8,134.8,132.8,131.8,132.8,132.9,132.8,132.8,133.8,132.0,131.8,132.8,132.8,133.8,131.8,131.8,132.8,131.9,131.8,131.8,131.9,130.8,132.8,132.8,130.8,134.8,131.8,129.8,132.8,134.8,132.8,132.8,131.8,131.8,133.8,131.8,131.9,132.8,130.8,134.8,131.8,132.8,133.8,132.8,131.9,133.8,132.8,132.8,132.8,133.8,132.9,131.8,131.8,133.8,132.8,133.8,132.8,131.8,133.0,133.8,131.8,134.8,131.8,133.8,132.8,132.8,131.8,131.8,133.8,129.8,132.8,133.8,131.8,132.8,132.9,134.8,131.8,130.8,132.8,130.8,131.9,133.8]}},"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"},{"attributes":{"callback":null,"end":124.16280556499999},"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},{"attributes":{"dimension":1,"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"attributes":{},"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},{"attributes":{},"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},"y_range_name":"memory"},"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"},{"attributes":{},"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},{"attributes":{},"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"attributes":{},"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"},{"attributes":{"callback":null,"end":134.8},"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[268.312576,304.861184,317.472768,319.01696,322.7648,325.931008,332.898304,337.629184,338.86208,345.014272,348.28288,351.514624,359.612416,357.261312,369.62304,370.520064,373.76,379.14624,384.540672,387.198976,390.43072,397.90592,403.357696,410.107904,415.506432,414.834688,417.472512,419.90144,424.751104,427.184128,434.208768,440.46336,444.100608,445.976576,451.371008,454.344704,458.514432,463.306752,468.578304,468.672512,474.341376,479.465472,485.945344,489.730048,491.58144,492.658688,494.923776,503.296,507.61728,511.668224,511.807488,517.484544,521.805824,524.98432,534.163456,536.993792,540.50816,546.537472,550.85056,550.191104,554.491904,557.73184,561.78688,568.799232,570.88,573.8496,577.089536,583.573504,589.791232,592.224256,593.563648,598.58944,606.154752,606.961664,614.52288,616.415232,622.559232,625.803264,626.880512,626.733056,633.761792,639.971328,643.096576,644.169728,647.954432,651.190272,656.044032,663.605248,672.251904,676.573184,677.920768,682.31168,685.551616,690.688,692.842496,695.812096,696.89344,700.858368,710.311936,712.761344,718.70464,720.32256,726.540288,726.540288,730.3168,735.911936,739.69664,743.940096,749.461504,753.512448,753.713152,761.757696,769.048576,770.932736,773.029888,774.815744,781.139968,784.654336,785.637376,789.159936,790.171648,789.970944,799.0272,803.59424,809.971712]}},"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},{"attributes":{"below":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"}},"left":[{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"},{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"},{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"right":[{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},"tools":[{"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"}],"x_range":{"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},"y_range":{"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"}},"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"attributes":{"legends":[["% CPU",[{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"}]],["Memory",[{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"}]]],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"attributes":{"callback":null,"end":809.971712,"start":268.312576},"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"},{"attributes":{},"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},{"attributes":{"data_source":{"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},"glyph":{"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"attributes":{"data_source":{"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},"glyph":{"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},"selection_glyph":null},"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"}],"root_ids":["4317cd29-f8a6-4e35-a040-49d7fe673b23"]},"title":"Bokeh Application","version":"0.11.1"}};
            var render_items = [{"docid":"5fa95652-5e47-49ac-a58d-a27ba9d352aa","elementid":"19b1b3ab-cc84-4e20-8ced-a836c3f430cc","modelid":"4317cd29-f8a6-4e35-a040-49d7fe673b23","notebook_comms_target":"342ba418-8982-4ae8-9b7d-9cd17e02ddef"}];
            
            Bokeh.embed.embed_items(docs_json, render_items);
        });
      },
      function(Bokeh) {
      }
    ];
  
    function run_inline_js() {
      for (var i = 0; i < inline_js.length; i++) {
        inline_js[i](window.Bokeh);
      }
    }
  
    if (window._bokeh_is_loading === 0) {
      console.log("Bokeh: BokehJS loaded, going straight to plotting");
      run_inline_js();
    } else {
      load_libs(js_urls, function() {
        console.log("Bokeh: BokehJS plotting callback run at", now());
        run_inline_js();
      });
    }
  }(this));
</script>


As you can see from the plot above, this computation uses just over 1 core (~130% CPU). The limiting factor is related to [h5py](h5py.org) which I've used to pull input data out of an [HDF5 file](https://www.hdfgroup.org/HDF5/). The h5py library is a totally awesome piece of software that I use every day, but HDF5 is not designed to support multi-threaded data access. Also, h5py doesn't release the [GIL](https://wiki.python.org/moin/GlobalInterpreterLock), a Python technicality which means other threads cannot run while h5py is doing anything, even if the other threads want to do something unrelated to HDF5 I/O.

## Dask + Zarr

Recently I've been working on [Zarr](zarr.readthedocs.io/), a new Python library for chunked, compressed, N-dimensional data. Previously I [introduced Zarr](http://alimanfoo.github.io/2016/04/14/to-hdf5-and-beyond.html) and showed how it can be used to get fast access into large multi-dimensional arrays. The other thing Zarr can do is let you read or write to an array from multiple threads or processes in parallel. Also, Zarr releases the GIL during compression and decompression, so other threads can carry on working. Here's the allele count example again, but this time using a Zarr array as the input data source:


{% highlight python %}
import zarr
zarr.__version__
{% endhighlight %}




    '2.1.3'




{% highlight python %}
# Setup a Zarr array, copying in genotype data from the HDF5 file.
# N.B., let's use the similar compression options as the HDF5 file for a fairer
# comparison, although other compressors might be faster.
# Let's also use SSD, same as where HDF5 was stored above.
genotype_zarr = zarr.open_like(genotype, path='data/2016-05-16/genotype.zarr', mode='w', 
                               chunks=chunks, compression='blosc',
                               compression_opts=dict(cname='zlib', clevel=3, shuffle=0))
genotype_zarr[:] = genotype
genotype_zarr
{% endhighlight %}




    Array((13167162, 765, 2), int8, chunks=(6553, 200, 2), order=C)
      nbytes: 18.8G; nbytes_stored: 640.7M; ratio: 30.0; initialized: 8040/8040
      compressor: Blosc(cname='zlib', clevel=3, shuffle=0)
      store: DirectoryStore




{% highlight python %}
# ensure OS pagecache is cleared 
!sudo drop_caches
{% endhighlight %}


{% highlight python %}
# run allele count computation via dask
gdz = allel.model.dask.GenotypeDaskArray(genotype_zarr)
acz = gdz.count_alleles(max_allele=3)
with ResourceProfiler(dt=1) as rprof:
    acz.compute()
visualize(rprof);
{% endhighlight %}




<div class="plotdiv" id="a37327e0-4d10-4f2e-be85-7a5fdc8198b9"></div>
<script type="text/javascript">
  
  (function(global) {
    function now() {
      return new Date();
    }
  
    if (typeof (window._bokeh_onload_callbacks) === "undefined") {
      window._bokeh_onload_callbacks = [];
    }
  
    function run_callbacks() {
      window._bokeh_onload_callbacks.forEach(function(callback) { callback() });
      delete window._bokeh_onload_callbacks
      console.info("Bokeh: all callbacks have finished");
    }
  
    function load_libs(js_urls, callback) {
      window._bokeh_onload_callbacks.push(callback);
      if (window._bokeh_is_loading > 0) {
        console.log("Bokeh: BokehJS is being loaded, scheduling callback at", now());
        return null;
      }
      if (js_urls == null || js_urls.length === 0) {
        run_callbacks();
        return null;
      }
      console.log("Bokeh: BokehJS not loaded, scheduling load and callback at", now());
      window._bokeh_is_loading = js_urls.length;
      for (var i = 0; i < js_urls.length; i++) {
        var url = js_urls[i];
        var s = document.createElement('script');
        s.src = url;
        s.async = false;
        s.onreadystatechange = s.onload = function() {
          window._bokeh_is_loading--;
          if (window._bokeh_is_loading === 0) {
            console.log("Bokeh: all BokehJS libraries loaded");
            run_callbacks()
          }
        };
        s.onerror = function() {
          console.warn("failed to load library " + url);
        };
        console.log("Bokeh: injecting script tag for BokehJS library: ", url);
        document.getElementsByTagName("head")[0].appendChild(s);
      }
    };var element = document.getElementById("a37327e0-4d10-4f2e-be85-7a5fdc8198b9");
    if (element == null) {
      console.log("Bokeh: ERROR: autoload.js configured with elementid 'a37327e0-4d10-4f2e-be85-7a5fdc8198b9' but no matching script tag was found. ")
      return false;
    }
  
    var js_urls = [];
  
    var inline_js = [
      function(Bokeh) {
        Bokeh.$(function() {
            var docs_json = {"20dc5c9d-0501-44c0-aa96-fa3d058c67f8":{"roots":{"references":[{"attributes":{"below":[{"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"2e732485-d13c-4557-bbea-69da07c1c371","type":"Range1d"}},"left":[{"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"},{"id":"f8a5709e-a198-4788-8bf7-4883dd26ebc0","type":"Grid"},{"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"},{"id":"5c5089e1-08cb-4817-aba8-7a894ca2f004","type":"Grid"},{"id":"7ab2d2a2-a95c-493c-8a95-65e493130694","type":"Legend"},{"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"},{"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"},{"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"}],"right":[{"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"2d04983b-c4d9-4414-a9d9-0269221a98ea","type":"ToolEvents"},"tools":[{"id":"d1547512-70fa-4bd7-a362-b29ce3aa912e","type":"PreviewSaveTool"},{"id":"0f24786c-df78-4a5e-b1c7-a75286cabaa4","type":"ResetTool"},{"id":"5befe7f4-ba41-4f0c-b520-7c8e5f618abb","type":"ResizeTool"},{"id":"8f63143d-3e06-4282-97f3-859d19700f85","type":"WheelZoomTool"},{"id":"87eb214a-1d74-4960-8341-1fa8277a4d8c","type":"PanTool"}],"x_range":{"id":"a6a65e67-e09d-47ca-960f-2693f4b2878f","type":"Range1d"},"y_range":{"id":"75272091-efa7-41a9-8047-6b001b21f75d","type":"Range1d"}},"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},{"attributes":{"data_source":{"id":"81ac2a2c-64e7-4dc7-9fbe-92ebd667c571","type":"ColumnDataSource"},"glyph":{"id":"e37126ab-e7e2-4c51-8c8a-a969618f2115","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"565ffb4f-8413-4fa6-b5c6-6eee6963841d","type":"Line"},"selection_glyph":null},"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"d003b7e3-1f54-41dd-b7e8-f6e0ddb8a8ab","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"}},"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"},{"attributes":{},"id":"26abe704-8c64-4ca3-92e1-b8be4617c179","type":"BasicTickFormatter"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0190051049999056,2.0202676689999635,3.021584463999943,4.02281131999996,5.024128939999969,6.02546096399999,7.026767620999976,8.028078234999953,9.029396844999951,10.030752998999901,11.032060126999909,12.03339271699997,13.034689307999997],"y":[1052.823552,853.66784,861.872128,867.831808,876.904448,901.791744,931.49184,960.098304,992.649216,1024.729088,1058.029568,1091.784704,1126.567936,1159.774208]}},"id":"3ef02d94-35f4-4ae0-9c38-2b2678f0347b","type":"ColumnDataSource"},{"attributes":{"dimension":1,"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"}},"id":"5c5089e1-08cb-4817-aba8-7a894ca2f004","type":"Grid"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},{"attributes":{},"id":"8281868a-449f-48e5-ae10-d85864434cff","type":"BasicTicker"},{"attributes":{},"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"},{"attributes":{"callback":null,"end":763.0},"id":"75272091-efa7-41a9-8047-6b001b21f75d","type":"Range1d"},{"attributes":{},"id":"2d04983b-c4d9-4414-a9d9-0269221a98ea","type":"ToolEvents"},{"attributes":{"dimensions":["width"],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"8f63143d-3e06-4282-97f3-859d19700f85","type":"WheelZoomTool"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"}},"id":"f8a5709e-a198-4788-8bf7-4883dd26ebc0","type":"Grid"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"d1547512-70fa-4bd7-a362-b29ce3aa912e","type":"PreviewSaveTool"},{"attributes":{"legends":[["% CPU",[{"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"}]],["Memory",[{"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"}]]],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"7ab2d2a2-a95c-493c-8a95-65e493130694","type":"Legend"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"aa09d41c-e052-4eb6-8ba9-7b1c291a5321","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"}},"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"attributes":{"callback":null,"end":134.8},"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"},{"attributes":{"callback":null,"end":13.034689307999997},"id":"a6a65e67-e09d-47ca-960f-2693f4b2878f","type":"Range1d"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"0f24786c-df78-4a5e-b1c7-a75286cabaa4","type":"ResetTool"},{"attributes":{},"id":"aa09d41c-e052-4eb6-8ba9-7b1c291a5321","type":"BasicTickFormatter"},{"attributes":{"callback":null,"end":1159.774208,"start":853.66784},"id":"2e732485-d13c-4557-bbea-69da07c1c371","type":"Range1d"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"2f20233c-f8ab-4aab-b72c-fe34fa6c03d3","type":"Line"},{"attributes":{"below":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"}},"left":[{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"},{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"},{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"right":[{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},"tools":[{"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"}],"x_range":{"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},"y_range":{"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"}},"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},{"attributes":{"data_source":{"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},"glyph":{"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"attributes":{"data_source":{"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},"glyph":{"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},"selection_glyph":null},"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[0.0,103.6,133.8,133.8,132.9,132.8,132.8,132.8,133.8,133.8,132.8,131.8,133.8,131.9,133.8,133.8,131.8,133.9,133.9,132.8,132.8,132.9,132.9,132.8,131.8,132.9,133.8,131.9,133.8,131.9,131.9,133.8,132.8,132.9,132.8,132.8,132.8,132.8,133.8,132.8,133.9,131.8,133.8,131.8,134.8,132.8,131.8,132.8,132.9,132.8,132.8,133.8,132.0,131.8,132.8,132.8,133.8,131.8,131.8,132.8,131.9,131.8,131.8,131.9,130.8,132.8,132.8,130.8,134.8,131.8,129.8,132.8,134.8,132.8,132.8,131.8,131.8,133.8,131.8,131.9,132.8,130.8,134.8,131.8,132.8,133.8,132.8,131.9,133.8,132.8,132.8,132.8,133.8,132.9,131.8,131.8,133.8,132.8,133.8,132.8,131.8,133.0,133.8,131.8,134.8,131.8,133.8,132.8,132.8,131.8,131.8,133.8,129.8,132.8,133.8,131.8,132.8,132.9,134.8,131.8,130.8,132.8,130.8,131.9,133.8]}},"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"565ffb4f-8413-4fa6-b5c6-6eee6963841d","type":"Line"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"e37126ab-e7e2-4c51-8c8a-a969618f2115","type":"Line"},{"attributes":{"callback":null,"end":124.16280556499999},"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},{"attributes":{},"id":"d003b7e3-1f54-41dd-b7e8-f6e0ddb8a8ab","type":"BasicTickFormatter"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0190051049999056,2.0202676689999635,3.021584463999943,4.02281131999996,5.024128939999969,6.02546096399999,7.026767620999976,8.028078234999953,9.029396844999951,10.030752998999901,11.032060126999909,12.03339271699997,13.034689307999997],"y":[0.0,132.5,746.0,736.0,743.1,752.0,763.0,760.0,753.0,744.0,738.0,747.0,747.0,723.1]}},"id":"81ac2a2c-64e7-4dc7-9fbe-92ebd667c571","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"3ef02d94-35f4-4ae0-9c38-2b2678f0347b","type":"ColumnDataSource"},"glyph":{"id":"2f20233c-f8ab-4aab-b72c-fe34fa6c03d3","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"ec565e49-f33a-49a5-b973-98c000f5773b","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"},{"attributes":{},"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"ec565e49-f33a-49a5-b973-98c000f5773b","type":"Line"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"attributes":{"legends":[["% CPU",[{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"}]],["Memory",[{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"}]]],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"attributes":{"callback":null,"end":809.971712,"start":268.312576},"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"},{"attributes":{},"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},"y_range_name":"memory"},"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"},{"attributes":{},"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"5befe7f4-ba41-4f0c-b520-7c8e5f618abb","type":"ResizeTool"},{"attributes":{},"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},{"attributes":{},"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"},{"attributes":{"dimension":1,"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[268.312576,304.861184,317.472768,319.01696,322.7648,325.931008,332.898304,337.629184,338.86208,345.014272,348.28288,351.514624,359.612416,357.261312,369.62304,370.520064,373.76,379.14624,384.540672,387.198976,390.43072,397.90592,403.357696,410.107904,415.506432,414.834688,417.472512,419.90144,424.751104,427.184128,434.208768,440.46336,444.100608,445.976576,451.371008,454.344704,458.514432,463.306752,468.578304,468.672512,474.341376,479.465472,485.945344,489.730048,491.58144,492.658688,494.923776,503.296,507.61728,511.668224,511.807488,517.484544,521.805824,524.98432,534.163456,536.993792,540.50816,546.537472,550.85056,550.191104,554.491904,557.73184,561.78688,568.799232,570.88,573.8496,577.089536,583.573504,589.791232,592.224256,593.563648,598.58944,606.154752,606.961664,614.52288,616.415232,622.559232,625.803264,626.880512,626.733056,633.761792,639.971328,643.096576,644.169728,647.954432,651.190272,656.044032,663.605248,672.251904,676.573184,677.920768,682.31168,685.551616,690.688,692.842496,695.812096,696.89344,700.858368,710.311936,712.761344,718.70464,720.32256,726.540288,726.540288,730.3168,735.911936,739.69664,743.940096,749.461504,753.512448,753.713152,761.757696,769.048576,770.932736,773.029888,774.815744,781.139968,784.654336,785.637376,789.159936,790.171648,789.970944,799.0272,803.59424,809.971712]}},"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"26abe704-8c64-4ca3-92e1-b8be4617c179","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"8281868a-449f-48e5-ae10-d85864434cff","type":"BasicTicker"},"y_range_name":"memory"},"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"},{"attributes":{},"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"},{"attributes":{},"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},{"attributes":{},"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},{"attributes":{"dimensions":["width"],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"87eb214a-1d74-4960-8341-1fa8277a4d8c","type":"PanTool"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"}],"root_ids":["4317cd29-f8a6-4e35-a040-49d7fe673b23","2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3"]},"title":"Bokeh Application","version":"0.11.1"}};
            var render_items = [{"docid":"20dc5c9d-0501-44c0-aa96-fa3d058c67f8","elementid":"a37327e0-4d10-4f2e-be85-7a5fdc8198b9","modelid":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","notebook_comms_target":"e9a51277-7f6a-4ea1-9ef9-7faafbc6662e"}];
            
            Bokeh.embed.embed_items(docs_json, render_items);
        });
      },
      function(Bokeh) {
      }
    ];
  
    function run_inline_js() {
      for (var i = 0; i < inline_js.length; i++) {
        inline_js[i](window.Bokeh);
      }
    }
  
    if (window._bokeh_is_loading === 0) {
      console.log("Bokeh: BokehJS loaded, going straight to plotting");
      run_inline_js();
    } else {
      load_libs(js_urls, function() {
        console.log("Bokeh: BokehJS plotting callback run at", now());
        run_inline_js();
      });
    }
  }(this));
</script>


This time I get over 700% CPU usage. Also the computation is about 8 times faster, which is about what you'd expect given the higher CPU usage. 

Zlib is a fairly slow compressor, what happens if we use something faster like LZ4?


{% highlight python %}
genotype_zarr_lz4 = zarr.open_like(genotype, path='data/2016-05-16/genotype.zarr.lz4', mode='w', 
                                   chunks=chunks, compression='blosc',
                                   compression_opts=dict(cname='lz4', clevel=5, shuffle=0))
genotype_zarr_lz4[:] = genotype_zarr
genotype_zarr_lz4
{% endhighlight %}




    Array((13167162, 765, 2), int8, chunks=(6553, 200, 2), order=C)
      nbytes: 18.8G; nbytes_stored: 1.0G; ratio: 18.2; initialized: 8040/8040
      compressor: Blosc(cname='lz4', clevel=5, shuffle=0)
      store: DirectoryStore




{% highlight python %}
# ensure OS pagecache is cleared 
!sudo drop_caches
{% endhighlight %}


{% highlight python %}
# run allele count computation via dask
gdz = allel.model.dask.GenotypeDaskArray(genotype_zarr_lz4)
acz = gdz.count_alleles(max_allele=3)
with ResourceProfiler(dt=1) as rprof:
    acz.compute()
visualize(rprof);
{% endhighlight %}




<div class="plotdiv" id="47baa5ba-affd-4b88-8fa0-f3156f6817f0"></div>
<script type="text/javascript">
  
  (function(global) {
    function now() {
      return new Date();
    }
  
    if (typeof (window._bokeh_onload_callbacks) === "undefined") {
      window._bokeh_onload_callbacks = [];
    }
  
    function run_callbacks() {
      window._bokeh_onload_callbacks.forEach(function(callback) { callback() });
      delete window._bokeh_onload_callbacks
      console.info("Bokeh: all callbacks have finished");
    }
  
    function load_libs(js_urls, callback) {
      window._bokeh_onload_callbacks.push(callback);
      if (window._bokeh_is_loading > 0) {
        console.log("Bokeh: BokehJS is being loaded, scheduling callback at", now());
        return null;
      }
      if (js_urls == null || js_urls.length === 0) {
        run_callbacks();
        return null;
      }
      console.log("Bokeh: BokehJS not loaded, scheduling load and callback at", now());
      window._bokeh_is_loading = js_urls.length;
      for (var i = 0; i < js_urls.length; i++) {
        var url = js_urls[i];
        var s = document.createElement('script');
        s.src = url;
        s.async = false;
        s.onreadystatechange = s.onload = function() {
          window._bokeh_is_loading--;
          if (window._bokeh_is_loading === 0) {
            console.log("Bokeh: all BokehJS libraries loaded");
            run_callbacks()
          }
        };
        s.onerror = function() {
          console.warn("failed to load library " + url);
        };
        console.log("Bokeh: injecting script tag for BokehJS library: ", url);
        document.getElementsByTagName("head")[0].appendChild(s);
      }
    };var element = document.getElementById("47baa5ba-affd-4b88-8fa0-f3156f6817f0");
    if (element == null) {
      console.log("Bokeh: ERROR: autoload.js configured with elementid '47baa5ba-affd-4b88-8fa0-f3156f6817f0' but no matching script tag was found. ")
      return false;
    }
  
    var js_urls = [];
  
    var inline_js = [
      function(Bokeh) {
        Bokeh.$(function() {
            var docs_json = {"b01bd468-82e5-4abe-bf59-f4b1fa201fd5":{"roots":{"references":[{"attributes":{"data_source":{"id":"81ac2a2c-64e7-4dc7-9fbe-92ebd667c571","type":"ColumnDataSource"},"glyph":{"id":"e37126ab-e7e2-4c51-8c8a-a969618f2115","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"565ffb4f-8413-4fa6-b5c6-6eee6963841d","type":"Line"},"selection_glyph":null},"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"9605ba79-6cf8-47ca-b55c-c6afef05a3c6","type":"BasicTickFormatter"},"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},"ticker":{"id":"fe116f05-a0d6-4b2c-85f9-b16adc071201","type":"BasicTicker"},"y_range_name":"memory"},"id":"1177c305-d1f0-4c5b-9808-c46e3a0cd7d1","type":"LinearAxis"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"29349351-e02f-4e8c-8ce9-a2bebf17fc4d","type":"BasicTickFormatter"},"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},"ticker":{"id":"7cef8189-8216-4c44-b278-5d7fef27c88a","type":"BasicTicker"}},"id":"2393f418-32de-4255-bc1a-7e703cf137bd","type":"LinearAxis"},{"attributes":{"dimension":1,"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"}},"id":"5c5089e1-08cb-4817-aba8-7a894ca2f004","type":"Grid"},{"attributes":{"data_source":{"id":"3ef02d94-35f4-4ae0-9c38-2b2678f0347b","type":"ColumnDataSource"},"glyph":{"id":"2f20233c-f8ab-4aab-b72c-fe34fa6c03d3","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"ec565e49-f33a-49a5-b973-98c000f5773b","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"},{"attributes":{},"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"},{"attributes":{"dimensions":["width"],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"8f63143d-3e06-4282-97f3-859d19700f85","type":"WheelZoomTool"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"565ffb4f-8413-4fa6-b5c6-6eee6963841d","type":"Line"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"d1547512-70fa-4bd7-a362-b29ce3aa912e","type":"PreviewSaveTool"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},{"attributes":{"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"43421aa2-a938-4c3d-b498-1d9991f3b0ad","type":"ResizeTool"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"attributes":{"dimensions":["width"],"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"8c2dca8e-bdfa-47e3-ae3d-9f9389c3658f","type":"WheelZoomTool"},{"attributes":{"callback":null,"end":134.8},"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},{"attributes":{"data_source":{"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},"glyph":{"id":"681c69f3-79e6-484b-9e2b-f2c66e73c51b","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},"selection_glyph":null},"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"ed4a578e-6252-477d-b4fb-24c403c3976b","type":"Line"},{"attributes":{},"id":"9e7b78b2-16e8-4710-a3ea-4bc6f6aaa855","type":"BasicTicker"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"0f24786c-df78-4a5e-b1c7-a75286cabaa4","type":"ResetTool"},{"attributes":{},"id":"dce5cafc-8e87-41c7-a4ba-c829e358979b","type":"BasicTickFormatter"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"a2057f1d-eb49-4b41-86d3-bec65a3bd611","type":"BasicTickFormatter"},"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},"ticker":{"id":"9e7b78b2-16e8-4710-a3ea-4bc6f6aaa855","type":"BasicTicker"}},"id":"0133a549-47b4-4b1f-be9d-dce808ffb16c","type":"LinearAxis"},{"attributes":{},"id":"d003b7e3-1f54-41dd-b7e8-f6e0ddb8a8ab","type":"BasicTickFormatter"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0190051049999056,2.0202676689999635,3.021584463999943,4.02281131999996,5.024128939999969,6.02546096399999,7.026767620999976,8.028078234999953,9.029396844999951,10.030752998999901,11.032060126999909,12.03339271699997,13.034689307999997],"y":[0.0,132.5,746.0,736.0,743.1,752.0,763.0,760.0,753.0,744.0,738.0,747.0,747.0,723.1]}},"id":"81ac2a2c-64e7-4dc7-9fbe-92ebd667c571","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"3cf987e0-6cc5-4539-80cd-a8d8b1bb3935","type":"Line"},{"attributes":{"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"5ec3890c-bb0e-4850-8d8b-2c93696a4e49","type":"ResetTool"},{"attributes":{},"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"ec565e49-f33a-49a5-b973-98c000f5773b","type":"Line"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"attributes":{"dimension":1,"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},"ticker":{"id":"7cef8189-8216-4c44-b278-5d7fef27c88a","type":"BasicTicker"}},"id":"6b754aac-a226-49ea-a1a3-2788adf3b7b3","type":"Grid"},{"attributes":{"callback":null,"end":809.971712,"start":268.312576},"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0198318770001151,2.021112637999977,3.0224466010001834,4.023769853999966,5.025097741000081,6.025634717000003,7.026700150000124,8.028008294000074,9.029321900000014,10.030628143000058],"y":[0.0,117.7,739.0,734.0,733.0,733.0,732.6,728.2,737.0,736.0,734.0]}},"id":"afc91503-dbf8-4bb9-b26b-0cf4b56dea97","type":"ColumnDataSource"},{"attributes":{},"id":"7cef8189-8216-4c44-b278-5d7fef27c88a","type":"BasicTicker"},{"attributes":{"callback":null,"end":124.16280556499999},"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"2f20233c-f8ab-4aab-b72c-fe34fa6c03d3","type":"Line"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"6d46f92a-f11e-46be-897b-9e952b8097d2","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"attributes":{"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"9aa0b661-87e5-4d3c-b27e-5137edfc3713","type":"PreviewSaveTool"},{"attributes":{"data_source":{"id":"d127548c-19c7-4816-acb7-4f71a2ebc063","type":"ColumnDataSource"},"glyph":{"id":"0c7a474f-4234-4108-986a-1250da910d9c","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"bfcf7253-6be0-4c7c-b62e-f71cbd5b6437","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"4e26634f-845f-4f70-a233-b8c3fad13643","type":"GlyphRenderer"},{"attributes":{},"id":"8281868a-449f-48e5-ae10-d85864434cff","type":"BasicTicker"},{"attributes":{},"id":"fe116f05-a0d6-4b2c-85f9-b16adc071201","type":"BasicTicker"},{"attributes":{"callback":null,"end":1120.681984,"start":769.257472},"id":"0eb7b1b3-e7f2-43d2-b835-a59aa0c7eb0e","type":"Range1d"},{"attributes":{},"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"4d6ecfe9-6538-4050-a5fb-23e42b151078","type":"Line"},{"attributes":{"axis_label":"Time (s)","formatter":{"id":"d003b7e3-1f54-41dd-b7e8-f6e0ddb8a8ab","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"}},"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0198318770001151,2.021112637999977,3.0224466010001834,4.023769853999966,5.025097741000081,6.025634717000003,7.026700150000124,8.028008294000074,9.029321900000014,10.030628143000058],"y":[897.937408,769.257472,780.701696,796.643328,836.653056,882.46272,925.724672,974.426112,1024.90112,1075.376128,1120.681984]}},"id":"d127548c-19c7-4816-acb7-4f71a2ebc063","type":"ColumnDataSource"},{"attributes":{"below":[{"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"2e732485-d13c-4557-bbea-69da07c1c371","type":"Range1d"}},"left":[{"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"397acca2-39a4-4016-9d45-3c7a97b5e309","type":"LinearAxis"},{"id":"f8a5709e-a198-4788-8bf7-4883dd26ebc0","type":"Grid"},{"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"},{"id":"5c5089e1-08cb-4817-aba8-7a894ca2f004","type":"Grid"},{"id":"7ab2d2a2-a95c-493c-8a95-65e493130694","type":"Legend"},{"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"},{"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"},{"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"}],"right":[{"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"2d04983b-c4d9-4414-a9d9-0269221a98ea","type":"ToolEvents"},"tools":[{"id":"d1547512-70fa-4bd7-a362-b29ce3aa912e","type":"PreviewSaveTool"},{"id":"0f24786c-df78-4a5e-b1c7-a75286cabaa4","type":"ResetTool"},{"id":"5befe7f4-ba41-4f0c-b520-7c8e5f618abb","type":"ResizeTool"},{"id":"8f63143d-3e06-4282-97f3-859d19700f85","type":"WheelZoomTool"},{"id":"87eb214a-1d74-4960-8341-1fa8277a4d8c","type":"PanTool"}],"x_range":{"id":"a6a65e67-e09d-47ca-960f-2693f4b2878f","type":"Range1d"},"y_range":{"id":"75272091-efa7-41a9-8047-6b001b21f75d","type":"Range1d"}},"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"26abe704-8c64-4ca3-92e1-b8be4617c179","type":"BasicTickFormatter"},{"attributes":{"callback":null,"end":10.030628143000058},"id":"2967b4b0-fdb3-4cb3-8de3-a58eeea83549","type":"Range1d"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.0190051049999056,2.0202676689999635,3.021584463999943,4.02281131999996,5.024128939999969,6.02546096399999,7.026767620999976,8.028078234999953,9.029396844999951,10.030752998999901,11.032060126999909,12.03339271699997,13.034689307999997],"y":[1052.823552,853.66784,861.872128,867.831808,876.904448,901.791744,931.49184,960.098304,992.649216,1024.729088,1058.029568,1091.784704,1126.567936,1159.774208]}},"id":"3ef02d94-35f4-4ae0-9c38-2b2678f0347b","type":"ColumnDataSource"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"attributes":{"callback":null,"end":1159.774208,"start":853.66784},"id":"2e732485-d13c-4557-bbea-69da07c1c371","type":"Range1d"},{"attributes":{"callback":null,"end":763.0},"id":"75272091-efa7-41a9-8047-6b001b21f75d","type":"Range1d"},{"attributes":{},"id":"2d04983b-c4d9-4414-a9d9-0269221a98ea","type":"ToolEvents"},{"attributes":{"dimension":1,"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"}},"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"attributes":{"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},"ticker":{"id":"9e7b78b2-16e8-4710-a3ea-4bc6f6aaa855","type":"BasicTicker"}},"id":"7bf69ce3-9cde-4dab-a8d4-03510cde2040","type":"Grid"},{"attributes":{},"id":"9605ba79-6cf8-47ca-b55c-c6afef05a3c6","type":"BasicTickFormatter"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"}},"id":"f8a5709e-a198-4788-8bf7-4883dd26ebc0","type":"Grid"},{"attributes":{"legends":[["% CPU",[{"id":"1515f122-38e5-4447-8a98-5a1840922398","type":"GlyphRenderer"}]],["Memory",[{"id":"18cc2340-fb7b-4231-8673-83872926939f","type":"GlyphRenderer"}]]],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"7ab2d2a2-a95c-493c-8a95-65e493130694","type":"Legend"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"0c7a474f-4234-4108-986a-1250da910d9c","type":"Line"},{"attributes":{"axis_label":"% CPU","formatter":{"id":"aa09d41c-e052-4eb6-8ba9-7b1c291a5321","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"11810a31-7149-4901-958c-67279cbe3f44","type":"BasicTicker"}},"id":"128d6cc8-a043-42be-897e-e1c7b31a4a4b","type":"LinearAxis"},{"attributes":{"callback":null,"end":13.034689307999997},"id":"a6a65e67-e09d-47ca-960f-2693f4b2878f","type":"Range1d"},{"attributes":{"callback":null,"end":739.0},"id":"911dd33b-4395-478b-845f-8bdc18966cae","type":"Range1d"},{"attributes":{},"id":"aa09d41c-e052-4eb6-8ba9-7b1c291a5321","type":"BasicTickFormatter"},{"attributes":{"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"4afe7a47-54bc-4016-8fb2-ff0bab13b21e","type":"BasicTicker"}},"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"},{"attributes":{},"id":"a2057f1d-eb49-4b41-86d3-bec65a3bd611","type":"BasicTickFormatter"},{"attributes":{"legends":[["% CPU",[{"id":"50989b1e-c571-493e-a02d-fc79ae05ff08","type":"GlyphRenderer"}]],["Memory",[{"id":"4e26634f-845f-4f70-a233-b8c3fad13643","type":"GlyphRenderer"}]]],"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"1d275b6b-f1a2-4ae1-9585-1588481d72d6","type":"Legend"},{"attributes":{"below":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"f913c37f-f3af-4446-b6cb-7d2b6493e4aa","type":"Range1d"}},"left":[{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"7a5e5a63-1d26-404e-ad84-1b5ce49dfdc4","type":"LinearAxis"},{"id":"32a5fcfa-ce96-45e9-b098-04976fd59bce","type":"Grid"},{"id":"ce88ebcf-fece-4097-b115-b3d3a06e511b","type":"LinearAxis"},{"id":"92e6bc32-1660-4da6-a949-82832b9b820b","type":"Grid"},{"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"},{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"right":[{"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},"tools":[{"id":"def710c2-ffe1-4931-ba4e-8a32d94d8589","type":"PreviewSaveTool"},{"id":"7649a88b-7466-431a-aa63-0e4a8d51dee6","type":"ResetTool"},{"id":"d64d2edf-afeb-45ef-847d-8fc29f635007","type":"ResizeTool"},{"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"}],"x_range":{"id":"08f28539-3242-47dd-8779-a92ca0fabbf1","type":"Range1d"},"y_range":{"id":"a54a935a-a8bd-4f4a-bb36-f9feded0bc02","type":"Range1d"}},"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},{"attributes":{"data_source":{"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},"glyph":{"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"0a01ff01-3ed8-46a3-986c-90feecaad3a2","type":"Line"},"selection_glyph":null,"y_range_name":"memory"},"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"fb36bf78-e46f-4307-8e41-e6ba86363175","type":"WheelZoomTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[0.0,103.6,133.8,133.8,132.9,132.8,132.8,132.8,133.8,133.8,132.8,131.8,133.8,131.9,133.8,133.8,131.8,133.9,133.9,132.8,132.8,132.9,132.9,132.8,131.8,132.9,133.8,131.9,133.8,131.9,131.9,133.8,132.8,132.9,132.8,132.8,132.8,132.8,133.8,132.8,133.9,131.8,133.8,131.8,134.8,132.8,131.8,132.8,132.9,132.8,132.8,133.8,132.0,131.8,132.8,132.8,133.8,131.8,131.8,132.8,131.9,131.8,131.8,131.9,130.8,132.8,132.8,130.8,134.8,131.8,129.8,132.8,134.8,132.8,132.8,131.8,131.8,133.8,131.8,131.9,132.8,130.8,134.8,131.8,132.8,133.8,132.8,131.9,133.8,132.8,132.8,132.8,133.8,132.9,131.8,131.8,133.8,132.8,133.8,132.8,131.8,133.0,133.8,131.8,134.8,131.8,133.8,132.8,132.8,131.8,131.8,133.8,129.8,132.8,133.8,131.8,132.8,132.9,134.8,131.8,130.8,132.8,130.8,131.9,133.8]}},"id":"191ca85e-352a-4d0d-ba42-a95b4d8f3362","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"value":"#1f77b4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"bfcf7253-6be0-4c7c-b62e-f71cbd5b6437","type":"Line"},{"attributes":{"line_color":{"value":"#253494"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"e37126ab-e7e2-4c51-8c8a-a969618f2115","type":"Line"},{"attributes":{},"id":"28b12531-34bc-479b-bc1e-c5f3be1d5144","type":"ToolEvents"},{"attributes":{},"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},{"attributes":{"data_source":{"id":"afc91503-dbf8-4bb9-b26b-0cf4b56dea97","type":"ColumnDataSource"},"glyph":{"id":"4d6ecfe9-6538-4050-a5fb-23e42b151078","type":"Line"},"hover_glyph":null,"nonselection_glyph":{"id":"ed4a578e-6252-477d-b4fb-24c403c3976b","type":"Line"},"selection_glyph":null},"id":"50989b1e-c571-493e-a02d-fc79ae05ff08","type":"GlyphRenderer"},{"attributes":{"line_color":{"value":"#41b6c4"},"line_width":{"value":4},"x":{"field":"x"},"y":{"field":"y"}},"id":"b50fdf47-3a55-4f7e-87e0-a8c119eccc56","type":"Line"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"},"ticker":{"id":"ca7af755-9f0c-440c-8760-a7967bc3c40e","type":"BasicTicker"},"y_range_name":"memory"},"id":"660cef3e-978e-4dc0-a92c-8210d5dc03ad","type":"LinearAxis"},{"attributes":{},"id":"e500e5c9-44ba-49f4-b143-a23fe544d46e","type":"ToolEvents"},{"attributes":{"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"5befe7f4-ba41-4f0c-b520-7c8e5f618abb","type":"ResizeTool"},{"attributes":{},"id":"704388ec-331d-4de4-99b3-424ba24021ee","type":"BasicTicker"},{"attributes":{"dimensions":["width"],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"a1d56c59-4fdb-4fca-80c7-f8ad41c63970","type":"PanTool"},{"attributes":{"callback":null,"column_names":["y","x"],"data":{"x":[0.0,1.014188417000014,2.0154294809999556,3.016738250000003,4.017277289000049,5.018627856999956,6.019928776000029,7.021277900999962,8.02261046000001,9.023941553999975,10.025255754,11.026850345999947,12.028178293999986,13.028941281000016,14.03023930400002,15.031555677000028,16.03284427999995,17.033500260999972,18.03394214900004,19.035234228000036,20.03654364700003,21.037096778999967,22.037506318000055,23.03882286700002,24.040171299000008,25.041163642000015,26.042451005999965,27.042876620000015,28.04416660200002,29.04523222900002,30.046199636999972,31.047680549999995,32.04898403899995,33.04937198499999,34.050693902000035,35.052048399,36.05344595500003,37.054809216999956,38.056097987000044,39.05741041900001,40.05825603599999,41.059548631999974,42.06083237400003,43.06219031399996,44.06364561600003,45.064962741000045,46.06625671799998,47.06755891199998,48.06803835699998,49.069360279999955,50.07064495199995,51.07200168999998,52.07237172199996,53.07350930200005,54.074809078000044,55.07611321599995,56.077436678000026,57.078722063999976,58.08001983500003,59.08131795899999,60.08183769100003,61.083213459000035,62.08449821700003,63.08554264300005,64.08680745100003,65.08816960299998,66.089493189,67.09077945399997,68.09209643199995,69.09342901100001,70.09485882900003,71.09617007600002,72.09749028600004,73.09880496899996,74.1002555,75.10157527900003,76.10288713299997,77.10436809700002,78.10565651700006,79.106761119,80.10809843699997,81.10931925800003,82.110654748,83.11199392000003,84.11328420899997,85.11463560000004,86.11592398100004,87.11651629000005,88.11792202100003,89.11927366199996,90.12055926699998,91.12185950000003,92.12314320999997,93.123908594,94.12522140099998,95.12657526999999,96.12792411299995,97.12922770199998,98.13057869500005,99.13189241400005,100.13318787000003,101.13350977300001,102.13479593600005,103.13612645600006,104.13745701200003,105.13876582700004,106.14007211600006,107.14135435699995,108.14290023800004,109.14433436399997,110.14567750200001,111.14700049099997,112.148388899,113.14969484599999,114.15099289199998,115.15231488899997,116.15349427299998,117.15423889199997,118.155556822,119.15685521199998,120.158146203,121.15944299499995,122.160782143,123.16149741599997,124.16280556499999],"y":[268.312576,304.861184,317.472768,319.01696,322.7648,325.931008,332.898304,337.629184,338.86208,345.014272,348.28288,351.514624,359.612416,357.261312,369.62304,370.520064,373.76,379.14624,384.540672,387.198976,390.43072,397.90592,403.357696,410.107904,415.506432,414.834688,417.472512,419.90144,424.751104,427.184128,434.208768,440.46336,444.100608,445.976576,451.371008,454.344704,458.514432,463.306752,468.578304,468.672512,474.341376,479.465472,485.945344,489.730048,491.58144,492.658688,494.923776,503.296,507.61728,511.668224,511.807488,517.484544,521.805824,524.98432,534.163456,536.993792,540.50816,546.537472,550.85056,550.191104,554.491904,557.73184,561.78688,568.799232,570.88,573.8496,577.089536,583.573504,589.791232,592.224256,593.563648,598.58944,606.154752,606.961664,614.52288,616.415232,622.559232,625.803264,626.880512,626.733056,633.761792,639.971328,643.096576,644.169728,647.954432,651.190272,656.044032,663.605248,672.251904,676.573184,677.920768,682.31168,685.551616,690.688,692.842496,695.812096,696.89344,700.858368,710.311936,712.761344,718.70464,720.32256,726.540288,726.540288,730.3168,735.911936,739.69664,743.940096,749.461504,753.512448,753.713152,761.757696,769.048576,770.932736,773.029888,774.815744,781.139968,784.654336,785.637376,789.159936,790.171648,789.970944,799.0272,803.59424,809.971712]}},"id":"f00e361d-a3f4-4c78-89de-480cb38c2cdf","type":"ColumnDataSource"},{"attributes":{"axis_label":"Memory (MB)","formatter":{"id":"26abe704-8c64-4ca3-92e1-b8be4617c179","type":"BasicTickFormatter"},"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"},"ticker":{"id":"8281868a-449f-48e5-ae10-d85864434cff","type":"BasicTicker"},"y_range_name":"memory"},"id":"25223ada-82a8-4f7a-a696-5fa6cac4fe24","type":"LinearAxis"},{"attributes":{"dimensions":["width"],"plot":{"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"}},"id":"2f427b29-15b0-484a-af57-efc989d350a2","type":"PanTool"},{"attributes":{},"id":"0785866c-8b43-4627-ad01-62e2fddbe82f","type":"BasicTicker"},{"attributes":{"below":[{"id":"0133a549-47b4-4b1f-be9d-dce808ffb16c","type":"LinearAxis"}],"extra_y_ranges":{"memory":{"id":"0eb7b1b3-e7f2-43d2-b835-a59aa0c7eb0e","type":"Range1d"}},"left":[{"id":"2393f418-32de-4255-bc1a-7e703cf137bd","type":"LinearAxis"}],"plot_height":300,"plot_width":800,"renderers":[{"id":"0133a549-47b4-4b1f-be9d-dce808ffb16c","type":"LinearAxis"},{"id":"7bf69ce3-9cde-4dab-a8d4-03510cde2040","type":"Grid"},{"id":"2393f418-32de-4255-bc1a-7e703cf137bd","type":"LinearAxis"},{"id":"6b754aac-a226-49ea-a1a3-2788adf3b7b3","type":"Grid"},{"id":"1d275b6b-f1a2-4ae1-9585-1588481d72d6","type":"Legend"},{"id":"50989b1e-c571-493e-a02d-fc79ae05ff08","type":"GlyphRenderer"},{"id":"4e26634f-845f-4f70-a233-b8c3fad13643","type":"GlyphRenderer"},{"id":"1177c305-d1f0-4c5b-9808-c46e3a0cd7d1","type":"LinearAxis"}],"right":[{"id":"1177c305-d1f0-4c5b-9808-c46e3a0cd7d1","type":"LinearAxis"}],"title":"Profile Results","tool_events":{"id":"28b12531-34bc-479b-bc1e-c5f3be1d5144","type":"ToolEvents"},"tools":[{"id":"9aa0b661-87e5-4d3c-b27e-5137edfc3713","type":"PreviewSaveTool"},{"id":"5ec3890c-bb0e-4850-8d8b-2c93696a4e49","type":"ResetTool"},{"id":"43421aa2-a938-4c3d-b498-1d9991f3b0ad","type":"ResizeTool"},{"id":"8c2dca8e-bdfa-47e3-ae3d-9f9389c3658f","type":"WheelZoomTool"},{"id":"2f427b29-15b0-484a-af57-efc989d350a2","type":"PanTool"}],"x_range":{"id":"2967b4b0-fdb3-4cb3-8de3-a58eeea83549","type":"Range1d"},"y_range":{"id":"911dd33b-4395-478b-845f-8bdc18966cae","type":"Range1d"}},"id":"07a4a937-14c6-45e1-944b-a3affb94c046","subtype":"Figure","type":"Plot"},{"attributes":{"legends":[["% CPU",[{"id":"88307054-b40f-45bf-8fe7-f44f94992e28","type":"GlyphRenderer"}]],["Memory",[{"id":"9f20068d-63c5-47d7-bd67-6e13f3b3008d","type":"GlyphRenderer"}]]],"plot":{"id":"4317cd29-f8a6-4e35-a040-49d7fe673b23","subtype":"Figure","type":"Plot"}},"id":"02cae5a5-b90a-4b6a-bae4-c2c359657db8","type":"Legend"},{"attributes":{},"id":"6e6050c1-e849-4668-90f3-2a0a7f4da6e0","type":"BasicTickFormatter"},{"attributes":{},"id":"29349351-e02f-4e8c-8ce9-a2bebf17fc4d","type":"BasicTickFormatter"},{"attributes":{"dimensions":["width"],"plot":{"id":"2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","subtype":"Figure","type":"Plot"}},"id":"87eb214a-1d74-4960-8341-1fa8277a4d8c","type":"PanTool"}],"root_ids":["4317cd29-f8a6-4e35-a040-49d7fe673b23","2bdac8a7-a88b-4edc-9d51-9ea46abdb1f3","07a4a937-14c6-45e1-944b-a3affb94c046"]},"title":"Bokeh Application","version":"0.11.1"}};
            var render_items = [{"docid":"b01bd468-82e5-4abe-bf59-f4b1fa201fd5","elementid":"47baa5ba-affd-4b88-8fa0-f3156f6817f0","modelid":"07a4a937-14c6-45e1-944b-a3affb94c046","notebook_comms_target":"4dec070f-f668-4278-a474-00b8c71a86b5"}];
            
            Bokeh.embed.embed_items(docs_json, render_items);
        });
      },
      function(Bokeh) {
      }
    ];
  
    function run_inline_js() {
      for (var i = 0; i < inline_js.length; i++) {
        inline_js[i](window.Bokeh);
      }
    }
  
    if (window._bokeh_is_loading === 0) {
      console.log("Bokeh: BokehJS loaded, going straight to plotting");
      run_inline_js();
    } else {
      load_libs(js_urls, function() {
        console.log("Bokeh: BokehJS plotting callback run at", now());
        run_inline_js();
      });
    }
  }(this));
</script>


This goes even faster, and I'm still getting nearly full CPU utilisation, so probably I could push my SSD harder.

## Distributed Dask + Zarr

I'm currently focused on making better use of multi-core processors, but others like [Matt Rocklin](https://github.com/mrocklin) are working on frameworks for large-scale distributed computing. After I released Zarr v0.4.0 in April, Matt got in touch to suggest a reorganization of the code so that Zarr arrays could be stored in distributed systems like S3 or HDFS. Earlier this week I [released Zarr v1.0.0](http://zarr.readthedocs.io/en/latest/release.html#release-1-0-0) which includes a new storage architecture to support this. Here is Matt on using the new version of Zarr with Dask and S3 on a 20 node (160 core) cluster...

<iframe width="560" height="315" src="https://www.youtube.com/embed/8WtaYvqhxHc" frameborder="0" allowfullscreen></iframe>

Dask’s distributed scheduler looks seriously cool. It’s exciting to think that computations I’m currently coding to run in parallel via Dask on my multi-core desktop could in future be scaled up to a large compute cluster without any extra work at all.

## Further reading

To go with the new Zarr release there is some [new documentation](http://zarr.readthedocs.io/), including a [tutorial](http://zarr.readthedocs.io/en/latest/tutorial.html), [API reference](http://zarr.readthedocs.io/en/latest/api.html) and [storage specification](http://zarr.readthedocs.io/en/latest/spec/v1.html). Please bear in mind that Zarr is still a young project, so if you do take it for a spin, any [feedback is appreciated](https://github.com/alimanfoo/zarr/issues).

See also:

* [Dask](http://dask.pydata.org)
* [h5py](http://www.h5py.org/)
* [h5py parallel](http://docs.h5py.org/en/latest/mpi.html)


{% highlight python %}
import datetime
print(datetime.datetime.now().isoformat())
{% endhighlight %}

    2016-11-01T19:38:06.960872
