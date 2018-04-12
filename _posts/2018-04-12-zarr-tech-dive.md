---
layout: post
title: Zarr 2.2 tech dive
---

I recently released version 2.2.0 of [Zarr](http://zarr.readthedocs.io), a
Python package for numerical data storage. One of the new features available
with Zarr is the ability to store and retrieve data directly from cloud object
storage systems like Amazon S3 or Google Cloud Storage (GCS). One of the coolest
things about working on Zarr is that although I work on genomics, I get to meet
people from a number of different scientific disciplines, with similar issues
and interests in being able to run interactive analyses over large numerical
datasets. In particular, there is some really exciting work going on within the
geoscience community around the
[Pangeo project](https://github.com/pangeo-data/pangeo), which is working to
enable better use of cloud infrastructure for ocean, atmosphere and climate data
science. The Pangeo project has put together some very nice demos of using Zarr
in combination with other packages like
[Dask](https://dask.pydata.org/en/latest/) and
[Xarray](https://xarray.pydata.org/en/stable/) to run interactive analyses on
large datasets via Google cloud. Off the back of that work, I was invited to
give a webinar as part of the
[ESIP tech dive](http://wiki.esipfed.org/index.php/Interoperability_and_Technology/Tech_Dive_Webinar_Series)
series. In the webinar I tried to give an overview the main architectural
elements of Zarr, with some worked examples. Hopefully it's useful, here it is
on YouTube:

<iframe width="740" height="430"
src="https://www.youtube.com/embed/np_p4JBAIYI?rel=0" frameborder="0"
allow="autoplay; encrypted-media" allowfullscreen></iframe>


