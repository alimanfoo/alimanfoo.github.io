---
layout: post
title: Web-based exploration of genome variation data
---

Last week our paper on [Panoptes](https://github.com/cggh/panoptes)
was
[published in Bioinformatics](https://doi.org/10.1093/bioinformatics/btx410). Panoptes
is the software we use to power several web applications on the
[MalariaGEN](http://www.malariagen.net) site, including the
[Ag1000G browser](https://www.malariagen.net/apps/ag1000g/). Although
it's taken us a while to get this paper in print, these web
applications have been in production for a couple of years now. To
celebrate, I thought I'd post links to some screencasts featuring the
Ag1000G browser which we made around the time it was first launched
back in 2014/2015. We've made some updates to the Ag1000G browser
since these videos were made, but they're still relevant, and
hopefully provide a useful introduction to some of the capabilities of
the Panoptes framework.

## 5 minute tour

Below is a screencast that
[Steve Pritchard](https://twitter.com/topcat3005) made, giving a brief
tour of the main features in the Ag1000G browser:

<iframe width="560" height="315" src="https://www.youtube.com/embed/LWCbi8t9Zug" frameborder="0" allowfullscreen></iframe>

## Extended tour

Below are 4 videos I made while giving a longer tour of the Ag1000G
browser. These are much more rough-and-ready (apologies) but go into
more depth on various features. 

### Part 1 - tables, queries, plots

The first video looks at the samples table and illustrates some basic
features around querying the table and making interactive plots:

<iframe width="560" height="315" src="https://www.youtube.com/embed/z4yxsE0hWXo?list=PLKbXDtRY2ZfWkgBwCBk-Hbp8dbPQwKUwL" frameborder="0" allowfullscreen></iframe>

### Part 2 - larger tables

The second video looks at the variants table and illustrates some
features available for working with larger tables:

<iframe width="560" height="315" src="https://www.youtube.com/embed/yHtmmqQDciQ?list=PLKbXDtRY2ZfWkgBwCBk-Hbp8dbPQwKUwL" frameborder="0" allowfullscreen></iframe>

### Part 3 - genome browser

The third video looks at the genome browser functionality, which has
been fantastic for it's ability to zoom in and out and summarise data
at any genomic scale:

<iframe width="560" height="315" src="https://www.youtube.com/embed/jLWMQoN3oJY?list=PLKbXDtRY2ZfWkgBwCBk-Hbp8dbPQwKUwL" frameborder="0" allowfullscreen></iframe>

### Part 4 - visual analytics

The fourth video looks at some of the visual analytics features,
including how to select data items and see how items are related
across different representations of the data:

<iframe width="560" height="315" src="https://www.youtube.com/embed/XUmmI8iewvY?list=PLKbXDtRY2ZfWkgBwCBk-Hbp8dbPQwKUwL" frameborder="0" allowfullscreen></iframe>

## Future work

We are working on a new (2.0) version of the Panoptes framework, which
addresses some of the scalability limitations in the current (1.6)
production release, as well as improving the data import process and
providing new capabilities in the user interface. The new version is
still in beta, but please feel free to get in touch if you are
interested in trying it out.

I should add that although Panoptes has been deployed in production
for a while, I still don't think of it as a finished product by any
stretch. Rather it is an exploration of some ideas about how to
visualise, explore and interact with large genome variation
datasets. Working with large and complex data is a major challenge,
and I'm very conscious that there are many questions people might like
to ask of the Ag1000G data that are very difficult or impossible to
answer with the current Ag1000G browser. All of us involved in ongoing
development of Panoptes would love to make connections with others
working on similar problems, and to hear from anyone who has ideas,
comments or frustrations about the software as it stands.
