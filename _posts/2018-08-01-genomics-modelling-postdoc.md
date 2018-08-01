---
layout: post
title: Senior post-doc in malaria vector genomics and modelling
---

*TL;DR we have a job opportunity working at the interface of genomics
and disease modelling for malaria vector control. Anyone with a strong
background in science and/or computing is welcome to apply. See the
[job
description](https://www.ndm.ox.ac.uk/current-job-vacancies/vacancy/136286-Senior-Postdoctoral-Research-Scientist-in-Malaria-Vector-Genomics-and-Modelling)
for details of how to apply. Read on for more background and
information on the post.*

## Current challenges in malaria control

Malaria is still a massive global health problem, particularly in
Africa. The best way to control malaria in Africa is to target the
mosquitoes (a.k.a. vectors) that transmit the disease from one person
to another. You can protect a household from malaria by giving
everyone an insecticide-treated bed net to sleep under, or by spraying
the walls with insecticides. Both of these measures have been very
effective, but there are challenges.

The first challenge is cost. To have an impact, literally millions of
bed nets have to be produced and distributed, and have to reach
everyone who needs them, even in remote regions. Nets last up to three
years at best, so new nets have to be distributed regularly. Spraying
is also a major logistical challenge. Hundreds of thousands of
households need to be sprayed, and spraying needs to be repeated every
6 months or so. Teams of people need to be trained and provided with
equipment and insecticide. A lot of money is spent on malaria control
each year, but it is a some way short of what would be required to
eliminate malaria within any reasonable time frame, and sustaining
funding over a long period is a big ask.

The second challenge is resistance. Mosquitoes evolve insecticide
resistance, which means that genetic changes occur that alter the
mosquitoes’ physiology and/or internal chemistry, rendering
insecticides ineffective. Those genetic changes are spreading
throughout mosquito populations across the African continent. There is
a broad consensus that something needs to be done, but there are no
easy answers, and the global health community is struggling to mount a
coordinated response.

## Advances in malaria vector control

Broadly speaking, there are two ways of responding to these
challenges. One is to use existing insecticide-based mosquito control
tools more wisely, to avoid or delay the spread of insecticide
resistance. The other is to develop new tools for mosquito control.

Insecticide resistance does not appear in a mosquito population
instantaneously. It is a process, whereby resistance mutations occur,
increase in frequency over a number of generations, recombine with
other mutations, and spread via the movement and mating of mosquitoes
from different locations. If a single insecticide is used continuously
year after year over a broad geographical area, then resistance will
emerge and spread quickly. But if more than one insecticide is
available, there are options for how these can be used. One option is
to rotate insecticides, switching after a period of one or two
years. Another is to use different insecticides in different areas, in
a mosaic pattern. A third option is to use more than one insecticide
at the same time. The good news is that there are a number of
different insecticides now available for spraying, and new nets are
also available that combine an insecticide with a “synergist” that
counteracts resistance. So these are practical options, and serious
efforts are underway to plan and implement a coordinated strategy for
insecticide resistance management.

With careful management we might be able to delay the spread of
insecticide resistance, but that doesn’t change the fact that bed nets
and spraying campaigns are expensive and logistically challenging. For
these reasons, serious investments are being made in the development
of alternative technologies for mosquito control. In particular, new
tools based on genome editing are in development that can be used to
introduce a genetic modification into a mosquito population. The
modification can be designed to crash the population, or to change its
ability to transmit malaria parasites. The trick is to engineer the
modifications so they are able to selfishly propagate themselves
through a mosquito population, a process known as gene drive. This is
not science fiction, the technology exists and has been proven in the
lab. There are still many hurdles to overcome before this technology
could be deployed in the field, but the potential value is enormous,
because these modifications are self-spreading, meaning that they are
carried from one place to another by mosquitoes themselves. The
logistics required to deploy them are therefore simpler, and they
could reach remote areas which are difficult to access for net
distribution or spraying.

## Why modelling?

If we want to manage insecticide resistance, how do we know if
rotations, mosaics or mixtures are best? Which will maximise the time
for resistance to emerge and spread? On what time frame should
rotations be done? On what geographical scale should mosaics be
implemented? Is it better to do both rotations and mosaics at the same
time? Which insecticide and/or synergist combinations are best? Should
we combine nets and spraying? If we want to deploy a gene drive
construct, how many mosquitoes do we need to release? Where do we need
to release them? How often do we need to release them? How far will
the gene drive spread? By how much will the mosquito population size
be reduced? What impact will this have on the rate of malaria?

These are not simple questions, and cannot be answered by simple
analysis. Instead, we have to build computer models of mosquito
populations, which capture all relevant aspects of climate, geography,
ecology, behaviour, epidemiology and evolution. These models are a way
to synthesize everything we know about mosquitoes, and allow us to try
out different vector control strategies in the virtual world of a
computer simulation, investigate the likely outcomes, and then make
informed decisions about how to proceed in the real world. For
example, models are now available that can be used to explore and
compare the impact of rotations, mosaics and mixtures of
insecticide-based interventions. Models have also been developed that
help us to understand how effective gene drive constructs have to be
before they will spread throughout a mosquito population in a
self-sustaining way. Models like these are a vital bridge between
theory and practice, and need to be as complete and accurate as we can
make them.

## Known unknowns

If you’re modelling a mosquito population, there are certain things
you need to know. For example, there are some basic questions about
mosquito biology, such as how long does a mosquito live, how many
times will a female lay eggs, how many eggs will she lay each time,
how many offspring will survive, etc. There are also questions about
the environment in which mosquitoes live, like where are the breeding
sites, where are people to feed on, are breeding sites there all year
round or only during the wet season, if so when is the wet season, and
so on.

Many of these questions can be answered with careful field and/or lab
work. However, other parameters are harder to determine directly. For
example, how big are mosquito populations? How far do mosquitoes fly
on average? Do some mosquitoes use aestivation (a kind of hibernation)
to sit out the dry season? Previous wisdom was that mosquitoes don’t
move very far, up to around 5 km in their lifetime. But recent studies
have shown that many insect species engage in purposeful,
wind-assisted migration, travelling hundreds of kilometres. There have
also been studies suggesting that different mosquito species use
different strategies to cope with the dry season, with some
aestivating and others migrating. These are important parameters when
it comes to modelling insecticide resistance or gene drive, because
they have a big impact on how fast and far these could spread.

## Mosquito genome sequencing

An alternate way to resolve some of these unknowns is via genome
sequencing. If you're not a population geneticist then this may seem
like a knight’s move, but stick with me. You are probably familiar
with the idea that DNA contains genes which encode proteins, and thus
a genome (the total DNA sequence of an individual organism) provides a
blueprint for growth and function. You may also have heard that a
genome contains a lot of “junk” DNA, which doesn't encode any
functional proteins and is silently passed down from one generation to
the next, mostly hidden from the buffeting forces of natural
selection. In fact this junk DNA is a treasure trove, because the
processes of random mutation and recombination that occur with each
cycle of mating and reproduction create a living record within this
DNA of each organism’s ancestry. By comparing these DNA sequences
between lots of individual mosquitoes, it is possible to reconstruct a
surprisingly detailed account of the major events that occurred to the
ancestors of those individuals, stretching from the previous
generation back an astonishingly long way into the past. This is
because events and processes like an increase or decrease in
population size, or migration between populations, affect the patterns
of mating between individuals, which in turn affect the patterns of
genetic diversity and relatedness observed in DNA sequences from
present day individuals.

This process of inferring the demographic history of a population from
genome sequence data is a major field of statistical and computational
research. There are a wealth of papers, methods and software packages
available (e.g., dadi, Stairway Plot, MSMC, ABC, to name a few). The
major drivers behind this research have largely been human history,
where genomic analysis has been used to reconstruct major migration
events or infer interbreeding between modern humans and neanderthals;
and conservation biology, where genetic estimates of population size
are used to determine whether a species is endangered, and if so, what
action is needed. However, this is still a young field, wrestling to
make best use of the data deluge becoming available thanks to recent
developments in high-throughout, whole-genome sequencing
technology. Another complication is that mosquitoes are quite unlike
many of the species that have been studied in this way, because of
their large population sizes and high levels of genetic
diversity. This means that there are no well-calibrated, off-the-shelf
solutions for inferring mosquito population size or migration rates
from genome sequence data. Substantial work is needed to understand,
adapt, improve and tune methods for demographic inference from genomic
data before we can get the answers we need.

## The role

We are looking for someone to work at the intersection between genome
sequence data analysis, statistical methods for population demographic
inference, and computational modelling and simulation of mosquito
populations and control. That probably sounds like an intense place to
be, requiring a lot of specialist knowledge, but that's not the
case. If you do have experience in one or more of these areas then
that will undoubtedly help. But this is fundamentally a
cross-disciplinary role, and so the most important requirement is a
willingness to step into an unfamiliar field, and to dive in without
losing sight of the bigger picture or the connections between
different areas. That's why we're open to applications from anyone
with a strong track record in any scientific or computational area.

Our team’s strength is in the analysis of genome sequence data, and so
a major part of the role will involve exploring, visualizing,
summarizing and understanding the genomic data being generated by
MalariaGEN (more on that below). But an important component will be to
work with and improve statistical methods and tools for inferring
demographic features of mosquito populations. And although we are not
a modelling group, some time spent gaining an intuition for infectious
disease modelling, and in particular, current approaches to modelling
of vector-borne diseases and vector populations, will be beneficial.

## The data

The genomic data you’ll be working with come from the Malaria Genomic
Epidemiology Network (MalariaGEN), which is a multi-centre
collaboration generating the world’s largest resources of genome
sequence data on malaria parasites and mosquitoes. Last year we
published the largest ever genetic study of mosquitoes (or in fact any
arthropod), in which we sequenced the genomes of 765 mosquitoes
collected from 8 African countries. That paper is just the beginning
of a much larger programme of mosquito genome sequencing. Data from
the second project phase are already curated and available online,
comprising genomes of 1,142 mosquitoes from 13 African countries. We
have already sequenced all ~3,000 mosquitoes from the third project
phase and are curating those data now. And we have launched a
follow-on project which will aim to sequence ~10,000 mosquito genomes
per year. All of this is made possible by our close partnership with
the Wellcome Sanger Institute, a world-leading centre for genomic
research.

## The code

This is both a scientific and a computational role, so hopefully you
will already have some programming ability, and will be keen to
develop new skills. I expect you'll have to work with existing
software written in a variety of languages, so some flexibility and
open-mindedness will be needed. People in our group code in a number
of different languages so there is a diverse set of experiences to
draw on, but we are converging on Python and making heavy use of the
scientific Python / PyData ecosystem of packages for our day-to-day
genomic data analysis. We are big fans of open source software,
everything we create goes up on GitHub, and we try to align our code
with other packages and best practices in the scientific Python
ecosystem as much as possible. In recent years we’ve created a couple
of open source Python packages, including scikit-allel and Zarr, and
it’s been a great pleasure to interact and collaborate with others
working in the scientific Python community. The culture of
professionalism, collaboration and high quality software engineering
is so strong now in the open source community, and we are benefiting a
lot from bringing that culture into our own work as a team. I
encourage everyone in our team to get involved in open source
software, it is an invaluable learning experience.

## The team

Professor Dominic Kwiatkowski is the head of our group at the Big Data
Institute (BDI) in Oxford where this role is based, and also heads our
sister group at the Sanger Institute. The work in our group spans
human, parasite and mosquito genomics, and I lead the mosquito
(vector) programme. This role will also involve working in close
collaboration with the Target Malaria project led by Austin Burt at
Imperial College, London. Target Malaria is leading the development of
gene drive technologies for malaria vector control, and is making
great strides.  Our groups in Oxford and Sanger also include teams
working on the genomics of malaria parasites, and there are many
parallels, so we keep in close touch. We are also very fortunate to be
supported by amazing lab, informatics, data, admin and communications
teams, and there is a great ethos and spirit of cooperation across the
piece. They are genuinely a great bunch of people to work with.

## The building

The job is based at the newly built Big Data Institute (BDI), which is
part of the Li Ka Shing Centre for Health Information and Discovery at
the University of Oxford. The director of the BDI is professor Gil
McVean, one of the most respected people in the field of statistical
genetics. Our group at the BDI is part of the infectious disease
programme, which also includes the Malaria Atlas Project led by Pete
Gething, and Christophe Fraser's group who work on HIV, as well as
Gil's research group who work on a variety of problems in statistical
genetics. The building is only a year old and still not full, but it
is a great environment to work in. There is a weekly infections
seminar where you can hear anything from me talking about the history
of insecticide resistance, to Gil or Christophe giving a whiteboard
talk on coalescent theory, to the head of the HIV programme at the
Gates Foundation talking about the plan for HIV control in Africa, to
how Bayesian geostatistics can be used to model the spread of
terrorism.

And the social areas and all of the seminar rooms have been painted
with whiteboard paint, so you can write on the walls. That's probably
all you needed to know.

## Work/life balance

This is an academic position, which means that we cannot compete with
the commercial sector on salary. There has to be some compensation,
and apart from all of the above, there is a lot of flexibility around
working hours. We all make an effort to be at the office between 10 am
and 4pm at least 4 days a week, but otherwise as long as the work is
getting done you are free to work where and when suits best. I have
three children aged 7, 2 and 9 months, and I usually stay at home in
the morning long enough to help get everybody ready and out the door
in time for school/pre-school. At the other end of the day, I always
leave in time to get home for bath and bedtime. I do sometimes work in
the evening, but the weekends are sacred. I work from home one day a
week, and I also travel over to the Sanger Institute one day a week,
so there is plenty of variety and opportunity to get focused work
done, as well as catching up with colleagues and collaborators,
without missing out on family time.

## How to apply

See the [official job
description](https://www.ndm.ox.ac.uk/current-job-vacancies/vacancy/136286-Senior-Postdoctoral-Research-Scientist-in-Malaria-Vector-Genomics-and-Modelling)
for details of how to apply. Don't forget to include a tailored cover
letter and a CV tailored to address the requirements listed in the job
description. If you can, try to back up the claims you make in your CV
with some evidence. E.g., if you say you are a fluent Python
programmer, give a brief description of a non-trivial program you have
implemented, or even better, a link to a GitHub repository with code
you wrote.

If you have read this far, thank you for taking an interest. If you
have any further questions, please feel free to get in touch.