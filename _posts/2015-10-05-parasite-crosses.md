---
layout: post
title: Parasite crosses
---

The malaria parasite *Plasmodium falciparum* is a tiny, single-celled
organism, uniquely adapted to living inside humans and mosquitoes. It
is, however, a eukaryote, which means it shares more in common with
animals and plants than it does with bacteria or viruses. For example,
malaria parasites reproduce sexually. When they are taken up by a
mosquito bite, they differentiate into male and female gametes, which
fuse to form a zygote. As a zygote, the parasite is diploid, meaning
it has two copies of each of the chromosomes that make up its genome,
like the cells in our bodies. However, it doesn't stay like that for
long, undergoing meiosis shortly afterwards. During meiosis the cell
divides, transmitting only a single genome copy to each daughter. The
rest of the parasite's life cycle is spent as a haploid. If you've
ever heard the joke about how man is but a brief stage in the life
cycle of a sperm, for malaria parasites that is not so far from the
truth, at least from a genomic point of view.

Because malaria parasites have sex, they can be crossed. A parasite
cross is similar in principle to a cross with any other
organism. Choose two "parents" that differ in some trait of interest
(e.g., resistance to an anti-malarial drug), give them the opportunity
to mate, collect their "children", and study how the trait has been
inherited. In practice it's not quite so simple, because you cannot
actually cross two individual parasites. Rather you have to cross two
parasite *clones*, meaning that each "parent" is actually many
thousands of parasites which have been cultured to be genetically
identical. However, Walliker et al. (1987) showed that it can be done,
performing the first cross between clones 3D7 and HB3. Wellems et
al. (1991) performed a second cross between clones HB3 and Dd2, which
they used to discover the gene responsible for chloroquine
resistance. Hayton et al. (200X) performed a third cross between
clones 7G8 and GB4, discovering a gene involved in parasite invasion
of human red blood cells, now a promising vaccine candidate.

The MalariaGEN *P. falciparum* genetic crosses was set up in 201X with
the goal of deep sequencing the parents and progeny of these three
crosses. The project has given us some unique insights into the
unusual world of *P. falciparum* genome biology, and it's been a great
privilege to take part in the analysis. We recently made a public data
release from the project, which includes SNP and INDEL calls from two
different methods of variant discovery, as well as the underlying
sequence data itself. We also created a web application to help
explore the data, and posted a preprint on biorxiv describing the data
and what it can tell us about how parasites evolve to face new
challenges, such as pressure from anti-malarial drugs. This post gives
an informal tour of the data and web application, and brings together
a few experiences that might be useful to others studying genome
variation in *Plasmodium* and other species.

## Variant discovery

The crosses are a great testing ground for variant discovery
methods. Each cross comprises two parent clones and up to 35 progeny
clones. We sequence at the haploid merozoite life cycle stage, so
there is not heterozygosity to worry about, calling a genotype means
calling a single allele. And for the progeny, we know their genotype
must be inherited from one parent or the other. True *de novo*
mutations are relatively rare, so any allele called in a progeny clone
but not called in either parent is probably a "Mendelian"
error. Furthermore, for some progeny there are multiple biological
replicates, which means we sequenced the same culture several
times. These replicates should be identical, so any disagreement
between them is also probably an error.

To call SNPs and small INDELs, we used two methods. The first was
based on alignment of sequence reads to the 3D7 reference genome using
BWA, followed by GATK best practice (BQSR, INDEL realignment,
UnifiedGenotyper, VQSR). The second used Cortex via the independent
workflow, which involves a De Bruijn graph assembly for each clone,
followed by mapping of assembled variants back to the reference genome
to obtain coordinates.

@@IMAGE







@@TODO


## Sandbox

The fact that parasites have sex is a blessing and a curse. It is a
curse because sex means recombination, and recombination generates
genetic novelty, 

You could say that genetics began with a cross. Two pea plants, one
tall, one short. Their children, all tall. But their children's
children, some short again. From patterns of inheritance like these,
Mendel inferred laws of sexual reproduction. The first law could be
put this way: you have two copies of your genetic material (DNA), one
inherited from your father, the other from your mother. Hold that
thought.

Mendel also found that different traits, like height and seed colour,
were inherited independently. But Morgan and Sturtevant, studying
fruit flies, found some traits that were often inherited
together. @@EXAMPLE. To explain their data, they proposed that genetic
material was carried by chromosomes. Before being transmitted to the
next generation, the two copies of each chromosome recombine, forming
a unique patchwork of the originals.

@@IMAGE

In 1987 Walliker et al. showed that what goes for plants and flies
also goes for Plasmodium falciparum. Malaria parasites have sex,
recombining chromosomes with each cycle of transmission from one host
to the next. To demonstrate this, Walliker et al. used a cross. They
selected two "parent" parasite clones, gave them the opportunity to
mate, then collected their "progeny". However, crossing parasites
isn't as easy as crossing plants. Malaria parasites have a complex
life cycle, and only have sex during the brief time they spend inside
a mosquito's midgut. 

This was fundamental biology, but for
malaria there was also a medical application. In 1990 Wellems et al
used a second Plasmodium falciparum cross to discover the gene
responsible for resistance to chloroquine, which at the time was in
widespread use for treatment of malaria but whose efficacy was
dropping rapidly. In 200X Hayton et al. performed a third cross, using
it to discover a gene involved in red blood cell invasion, which is
now a promising vaccine candidate.

Walliker et al (1987), Wellems et al. (1990) and Hayton et al. (200X)
each performed a single parasite cross. The concept is the
fundamentally the same as for crossing plants or flies: pick two
parents with different phenotypes, induce them to mate, then study how
the different parental traits are transmitted to the progeny. However,
for microscopic unicellular parasites the practicalities are very
different. 