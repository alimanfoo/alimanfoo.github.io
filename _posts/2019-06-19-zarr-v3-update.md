---
layout: post
title:  "Zarr protocol v3 design update"
date:   2019-06-19
categories: zarr specs
---

*This article was originally posted on the [Zarr
blog](https://zarr-developers.github.io/zarr/specs/2019/06/19/zarr-v3-update.html).*

Today I put together some [slides summarising the current state of
exploratory work on the Zarr v3 protocol
spec](https://zarr-developers.github.io/slides/v3-update-20190619.html). The
purpose of this blog post is to share those slides more widely, and to
provide some context explaining why work has started on a v3 spec.

## Why work on a v3 spec?

The [current (v2) Zarr
spec](https://zarr.readthedocs.io/en/stable/spec/v2.html) is
implemented in a number of software libraries, and is a stable and
robust protocol that is used in production in a number of different
scientific communities. If you need to store and compute in parallel
against large array-like data, it's a good solution. So why start
thinking about a new protocol version?

### Language-agnostic

One reason is that the v2 protocol is somewhat Python-centric, and
includes some features which are not straightforward to implement in
other languages. This has meant that implementations do not all
support the same feature set. It would be good to have a minimal v3
protocol spec that could be fully implemented in any language, so all
implementations have parity around a core feature set.

### Unifying Zarr and N5

Another reason is that we would like to merge development efforts
between the Zarr and N5 communities, and so a goal for the v3 spec is
to unify the two approaches and provide a common implementation
target.

### Extensibility

A third reason is that a number of different groups have started
experimenting and extending the Zarr protocol in interesting ways, but
it's not always clear how to extend the v2 protocol to support new
features. It would be good if the v3 spec provided a variety of clear
extension points and extension mechanisms.

### Cloud storage

Finally, while the v2 spec can be used very effectively with
distributed storage systems like Amazon S3 or Google Cloud Storage,
there is room for improvement, particularly regarding how metadata is
stored and organised.

## Zarr v3 design update

I you are interested in knowing more about the current status of work
on the v3 spec, please take a look at the [v3 design update
slides](https://zarr-developers.github.io/slides/v3-update-20190619.html). The
slides use reveal.js and have both horizontal and vertical
navigation - if you haven't seen that before, then navigate downwards
first wherever you can, before navigating to the right.

As I mention in the slides, the current v3 spec is just a straw man,
meant to illustrate some ideas and potential solutions, but everything
is up for discussion. So if you have any comments or ideas, please do
get in touch, anyone is welcome to participate.
