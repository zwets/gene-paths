# gene-paths

_Determine gene order and orientation in assemblies._


## Installation

#### Requirements

* A C++ compiler

#### Building

* `make`


## Quick Usage

* `gene-paths --help`
* `gene-paths -g examples/assembly.gfa
              -f examples assembly.fna
              -q examples/query1.gfa`


## Background

#### Motivation

These tools  were developed to help answer the question "where are genes
relative to each other?"  For instance, when typing _Staphyloccos aureus_
SCCmec cassettes, it is relevant whether IS431 is up- or downstream of the
_mecA_ gene, and whether it is inverted.

This question is simple to answer if the features are on a single contig in
the genome assembly.  Unfortunately this is often not the case, especially
when repetitions and mobile elements are involved.

It would seem that if we find features A and B on separate contigs, we are
out of luck.  Not knowing what is between the contigs, or even their order,
we can't tell which of A or B is upstream, whether one is inverted relative
to the other, and we certainly can't know the distance between them.  Or so
you would think.

You may be surprised to learn that quite the opposite is true.

Contigs do not generally end because of a lack of data, i.e. because what
comes next is not covered by reads.  In most cases their left and right ends
are covered by reads, but these imply that there are multiple continuations.

In the assembly graph this could look like this:

                       ___ contig 2 ___                
                      /                \                
     --- contig 1 ---+                  +--- contig 4 ---
                      \___ contig 3 ___/

Here contig 1 cannot be extended further because in one place on the genome
it is followed by the sequence captured in contig 2, whereas elsewhere it is
followed by the (different) sequence in contig 3.

However, if we know that our genes of interest are on contigs 1 and 2, then,
with the knowledge of the assembly graph, we can actually find their genomic
distance, relative orientation, and order.  In fact, we know the exact
nucleotide sequence connecting the two genes!

Even if the features were located on contigs 1 and 4, we could still put a
bound on their genomic distance, and (with a bit of luck) know their relative
orientation and order with certainty.

It is remarkable how much information we discard by routinely working with
assembled contigs rather than assembly graphs.

Clearly, as the number of edges between the features of interest increases,
the number of possible paths connecting them rapidly explodes.  It is already
impossible to predict from this continuation of the graph:

     __ contig 2 ___                    ___ contig 5 __
                    \                  /
                     +--- contig 4 ---+
     __ contig 3 ___/                  \___ contig 6 __

which of the sequences `2-4-5`, `2-4-6`, `3-4-5`, `3-4-6` are present with
certainty on the genome (though we know that at least two must be).

Nevertheless, for the analysis of localised features such as arrangement of
gene cassettes, insertion sites, promoter regions, etc. the assembly graph
provides indispensable information.

This is the motivation for writing `gene-paths` (and yes, the work is still
in progress).

#### How are assembly graphs stored?

The GFA format is the de facto standard for assembly graphs.  Its spec is
maintained at <https://github.com/GFA-spec/GFA-spec>.
[GFA2](https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md) is the
current recommended format, and can be easily upgraded to from GFA v1.

#### How do I get an assembly graph?

* [SPAdes](http://cab.spbu.ru/software/spades/) by default writes the assembly
  graph to `contigs.gfa`, but note the points made about overlap removal in
  the [Unicycler documentation](https://github.com/rrwick/Unicycler#background).
* [Unicycler](https://github.com/rrwick/Unicycler) and its successor,
  [Trycycler](https://github.com/rrwick/Trycycler/wiki) optimise SPAdes and
  can perform hybrid short and long read assembly.
* [SKESA](https://github.com/ncbi/SKESA.git) can generate assembly graphs for
  assemblies and reads using its `gfa_connector` utility.  SKESA's assemblies
  tend to be more fragmented than SPAdes's, but there is a
  [trade-off with speed and the number of misassemblies](https://cab.spbu.ru/benchmarking-tools-for-de-novo-microbial-assembly/).
* [SAUTE](https://github.com/ncbi/SKESA#saute---sequence-assembly-using-target-enrichment)
  was recently added to SKESA, and combines _de novo_ assembly with alignment
  to a reference, producing all paths through the assembly graph consistent
  with the reference.

#### Can I look at assembly graphs?

Yes, and it's very insightful.  Get the fabulous
[Bandage](https://github.com/rrwick/Bandage) and enjoy!

#### What other tools work with GFA?

* [gfatools](https://github.com/lh3/gfatools) is a GFA parser with a succinct
  C implementation of data structures to hold assembly graphs.  It has tools
  to convert GFA to FASTA, BED, and SQL, and to perform miniasm-like
  transformations on the graph.
* [gfakluge](https://github.com/edawson/gfakluge) is a GFA parser and C++
  class library that maps one-to-one to the entities in GFA.  Its `gfak` tool
  converts between GFA formats, and can sort, subset, merge, and trim graphs.
* The [page on modularising assemblers](https://github.com/GFA-spec/assembler-components)
  will have more references of interest.
* [vg](https://github.com/vgteam/vg) is the nuclear option.  Whereas the tools
  above target GFA specifically, and (thus far) provide primarily parsers and
  basic data structures, VG is the tool of the trade for working with [genome
  variation graphs]().


---
#### Licence

gene-paths - determine gene order and orientation in assemblies  
Copyright (C) 2021  Marco van Zwetselaar <io@zwets.it>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

