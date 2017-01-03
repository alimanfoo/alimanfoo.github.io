---
layout: page
title: Reading
permalink: /reading/
---

This page is a reading diary, inspired by [#365papers](https://twitter.com/hashtag/365papers?src=hash) (or 
[#230papers](https://twitter.com/hashtag/230papers?src=hash) realistically).

<div id="papers">

<div style="margin: 0 0 1em 0">

  <p>
    <input class="search" placeholder="Search" />

    <button class="sort" data-sort="author">
      Sort by first author
    </button>

    <button class="sort" data-sort="year">
      Sort by year
    </button>

    <button class="sort" data-sort="read">
      Sort by date read
    </button>

    <input id="abstracts" type="checkbox" name="abstracts" onclick="abstracts();">
    Show abstracts
  </p>

</div>

  <ul class="list">

    <li>
      <span class="citation">
        <span class="author">Keightley et al.</span>
        (<span class="year">2015</span>)
      </span>
      <a href="https://www.ncbi.nlm.nih.gov/pubmed/25371432">
        <span class="title">Estimation of the spontaneous mutation rate in 
Heliconius melpomene</span>
      </a>
      [<span class="read">2017-01-03</span>]
      <blockquote class="abstract" style="display: none">We estimated the spontaneous mutation rate in Heliconius 
melpomene by genome sequencing of a pair of parents and 30 of their offspring, 
based on the ratio of number of de novo heterozygotes to the number of 
callable site-individuals. We detected nine new mutations, each one affecting 
a single site in a single offspring. This yields an estimated mutation rate of 
2.9 × 10(-9) (95% confidence interval, 1.3 × 10(-9)-5.5 × 10(-9)), which is 
similar to recent estimates in Drosophila melanogaster, the only other insect 
species in which the mutation rate has been directly estimated. We infer that 
recent effective population size of H. melpomene is about 2 million, a 
substantially lower value than its census size, suggesting a role for natural 
selection reducing diversity. We estimate that H. melpomene diverged from its 
Müllerian comimic H. erato about 6 Ma, a somewhat later date than estimates 
based on a local molecular clock.</blockquote>
    </li>

  </ul>

</div>

<style type="text/css">
.citation {
  font-weight: bold;
}
</style>

<script type="text/javascript" src="/assets/list.min.js"></script>

<script type="text/javascript">

var options = {
  valueNames: [ 'author', 'year', 'title', 'abstract', 'read' ]
};

var userList = new List('papers', options);

function abstracts() {
  if (document.getElementById('abstracts').checked) {
    var display = "block";
  } else {
    var display = "none";
  }
  var nodes = document.querySelectorAll(".abstract");
  for (var i=0; i < nodes.length; i++) {
    nodes[i].style.display = display;
  }  
}
</script>
