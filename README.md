# gene-paths

_Determine gene order and orientation in assemblies._


## Motivation

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
the number of possible paths between them rapidly explodes.  It is already
impossible to predict from this continuation of the graph:

     ___ contig 2 ___                    ___ contig 5 ___
                     \                  /
                      +--- contig 4 ---+
     ___ contig 3 ___/                  \___ contig 6 ___

which of the sequences `2-4-5`, `2-4-6`, `3-4-5`, `3-4-6` are present with
certainty on the genome (though we know that at least two must be).

Nevertheless, for the analysis of localised features such as arrangement of
gene cassettes, insertion sites, promoter regions, etc. the assembly graph
provides indispensable information.

This is the motivation for writing `gene-paths` (and yes, the work is still
in progress).


## Usage

#### Where to get assembly graphs?

SPAdes by default writes the assembly graph to  `contigs.gfa`.
[Unicycler](https://github.com/rrwick/Unicycler) can be used to further
polish assemblies.  SKESA can generate the assembly graph post-hoc using
`gfa-connector`.  It's recent addition `saute` generates GFA by default.

#### Can I look at assembly graphs?

Yes, and it's very insightful.  Get the fabulous
[Bandage](https://github.com/rrwick/Bandage) and enjoy!


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

