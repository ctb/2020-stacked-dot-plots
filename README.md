# Stacked dot plots

## an exploratory repository

Run [mashmap](https://github.com/marbl/MashMap) or
[nucmer](https://github.com/mummer4/mummer) to find alignments between
a query genome and one or more target genomes; plot, dot-plot style.

See [the source in `alignplot.py`](./alignplot.py) and
[the `dotplot` notebook](./dotplot.ipynb) for more info.

The first two examples below are of "non-stacked" dotplots, as in,
they are just pairwise dotplots.

The third example below "stacks" the dotplots on the x axis, so that you can
see two dotplots against the same query genome in the query genome
coordinates, but with the target genomes split out.

## installation and running

For now, if you want to try this out for yourself --

* you'll need python3 installed
* you'll need mummer4 and/or mashmap installed, along with matplotlib.
* you'll need to put genomes somewhere (e.g. the `genomes/` directory)
* you can add names or replace names in info.csv
* see [the source in `alignplot.py`](./alignplot.py) and
[the `dotplot` notebook](./dotplot.ipynb) for more info on actually running.

(More work to be done here -- will binderize eventually :)

## examples

### Shared sequences (contamination?) between two Genbank genomes.

Here, the Rokubacterium (red, x axis) is probably contaminated by sequences
from the Acidobacterium (y axis coordinates).

![](images/example1.png)

### Shared sequences between a TARA MAG and a Genbank genome.

The TARA binned genome/MAG is on the y axis, the Genbank genome is on x in
red.

![](images/example2.png)

### Shared sequences between a TARA MAG and two Genbank genomes.

The TARA binned genome/MAG is on the y axis, the first Genbank genome
is on x in red, the second is on x in blue.

![](images/example3.png)

## contact

Titus Brown, titus@idyll.org

CTB Nov 2020
