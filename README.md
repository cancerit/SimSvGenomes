SimSvGenomes
============

Exhaustive enumeration or derivative genome structures from simple
rearrangement types. The simple rearrangement types are the following.

* Deletion
* Tandem duplication
* Unbalanced translocation (can happen intra-chromosomally and between sister
  chromatids)
* Balanced translocation (can happen intra-chromosomally and between sister
  chromatics)
* Telomeric break
* Telomeric break and fusion (i.e. one BFB cycle)
* Whole chromosome duplication
* Whole chromosome deletion
* Whole genome duplication

This program supports haploid or diploid wild type genomes with one or many
distinct chromosomes. 

SimSvGenomes is an unofficial placeholder name for this program.


Compiling
=========

This program requires glib. 

I am not an expert in compiling C code, but here are commands that work for me on Ubuntu Linux. 

    # Compiling with maximum compiler optimisation (O3)
    gcc -O3 rg_enumerator.multi_chr.main.c -lglib-2.0 -I/nfs/users/nfs_y/yl3/programs/glib-2.38.2/glib -I/nfs/users/nfs_y/yl3/programs/glib-2.38.2/ -o rg_enumerator.multi_chr.O3

    # Compiling with debugging information - I use Valgrind
    gcc -g rg_enumerator.multi_chr.main.c -lglib-2.0 -I/nfs/users/nfs_y/yl3/programs/glib-2.38.2/glib -I/nfs/users/nfs_y/yl3/programs/glib-2.38.2/ -o rg_enumerator.multi_chr


Usage
=====

    ./rg_enumerator.multi_chr.O3 <n_chrs> <diploid> <max_dup_depth> <max_overall_depth>

    Options:

    n_chrs - integer, number of different chromosomes in the wild type genome. 

    diploid - 0 or 1, whether the wild type genome is haploid or diploid.

    max_dup_depth - integer, maximum number of duplication-type rearrangements
        to enumerate.

    max_overall_depth - integer, maximum overall number of SVs to exhausively
        enumerate.

Duplicating rearrangements are tandem duplication, BFB-induced fold-back
inversion and whole-chromosome duplication. 

**Whole-genome duplication is hardcoded to occur at most once in the history of
a simulated derivative genome.**

As an example, feel free to try the following command.

    ./rg_enumerator.multi_chr.O3 2 1 2 2


Output format
=============

The output of the program is a space-delimited file with the following columns.

1. Detailed history of a simulated genome
2. Simple history of a simulated genome
3. Allelic copy number portion of the normalised rearrangement pattern string
4. Rearrangement portion of the rearrangement pattern string
5. Genome string describing the karyotype of the simulated genome

### 1. Detailed history of a simulated genome
This is the simple history, but each individual rearrangement is appended by an
index indicating the order in which it was generated. The two breakpoints of a
rearrangement are enumerated in a given somatic genome in a defined order,
so detailed history of a simulated genome enable the precise reconstruction of
the breakpoint positions of each rearrangement to generate a given simulated
genome. 

### 2. Simple history of a simulated genome
Same as detailed history, but without the indices. 

### 3. Allelic copy number portion of the normalised rearrangement pattern string
Chromosomes are semicolon-delimited. Segments within chromosomes are delimited
by forward slashes ('/'). When run in diploid mode, maternal and paternal copy
number of each segment are comma-separated. 

### 4. Rearrangement portion of the rearrangement pattern string
Segments above are identified using a zero-based index in the order in which
they appear. Rearrangements are forward slash-separated and the low and high
end of each rearrangement is separated by a comma. 

### 5. Genome string describing the karyotype of the simulated genome
The normalised genome string of the current derivative chromosome. Every
individual *haploid* chromosome is described within a pair of curly brackets.
That is, there will be two chromosomes in a wild type genome containing a
single diploid chromosome. Each chromosome contains a semicolon-separated list
of segments, identified by three comma-separated integers. The meaning of the
integers are `<segment ID>,<segment is paternal>,<segment is inverted>`. Note
that segment IDs in the normalised genome string are distinct from segment IDs
in the rearrangement pattern string.

When this column contains a detailed history notation instead of a genome,
string, it means that the detailed history specified in the first column
produces a derivative genome with an identical structure to a previously
encountered derivative genome that was generated using the detailed history in
this column. 


Notes
=====

There are no chromosomal identities beyond rearrangements patterns present (or
simulated) in each chromosome. For the purposes of this program "chr1" and
"chr2" are the same thing and produce the same genome string as long as they
have the same rearrangement patterns. 

Balanced translocations, and direct inversions and fold-back rearrangements are
generated with balanced breakpoints. For example, the rearrangement pattern
string of '0+,1+/1-,2-'.

When annotating real rearrangement breakpoints in the PCAWG project, we had no
way to distinguish balanced from unbalanced breakpoints (do breakpoints have to
touch? or have a <10bp distance? or <100bp?). Therefore no breakpoints are
annotated as balanced in PCAWG analysis, and therefore a direct inversion in
actual PCAWG data would get the rearrangement pattern string '0+,2+/2-,4-',
with deleted segments 1 and 3.

There are two ways to consolidate this difference. First way is to ignore it,
because a direct inversion '0+,1+/1-,2-' will become '0+,2+/2-,4-' when a
subsequent deletion overlaps with one of its breakpoints. So with a
sufficient simulation depth, the "unbalanced versions" of balanced breakpoints
will be guaranteed to be covered. Secondly, one could simply parse the
rearrangement simulator output to convert all balanced breakpoint into their
"unbalanced versions".

The maximum depth I have managed to reach with this program is on a single
haploid chromosome down to five rearrangements, of which at most four are
duplicating rearrangements. 

Different depth thresholds are provided for duplicative and non-duplicative
rearrangements because duplicative events duplicate existing segments and
thus significantly increase the space in which subsequent breakpoints can
occur. Whole-genome duplication events are limited to one per history
because of the same reason. 


LICENCE
========
Copyright (c) 2014-2017 Genome Research Ltd.

Author: Yilong Li <yl3@sanger.ac.uk>, <liyilong623@gmail.com>

SimSvGenomes is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
