# Index
 * [What is pufferfish?](#whatis)
 * [Using pufferfish?](#using)
 
# What is pufferfish? <a name="whatis"></a>

**short answer** : Pufferfish is a new time and memory-efficient data structure for indexing a compacted, colored de Bruijn graph (ccdBG).  You can read more about pufferfish in our [pre-print on bioRxiv](https://www.biorxiv.org/content/early/2017/09/21/191874).

**long answer** : 
Though the de Bruijn Graph (dBG) has enjoyed tremendous popularity as an assembly and sequence comparison data structure, it has only relatively recently begun to see use as an index of the reference sequences (e.g. [deBGA](https://github.com/HongzheGuo/deBGA), [kallisto](https://github.com/pachterlab/kallisto)). Particularly, these tools index the compacted dBG (cdBG), in which all non-branching paths are collapsed into individual nodes and labeled with the string they spell out. This data structure is particularly well-suited for representing repetitive reference sequences, since a single contig in the cdBG represents all occurrences of the repeated sequence. The original positions in the reference can be recovered with the help of an auxiliary "contig table" that maps each contig to the reference sequence, position, and orientation where it appears as a substring. The [deBGA paper](https://academic.oup.com/bioinformatics/article-abstract/32/21/3224/2415064/deBGA-read-alignment-with-de-Bruijn-graph-based?redirectedFrom=fulltext) has a nice description how this kind of index looks (they call it a unipath index, because the contigs we index are unitigs in the cdBG), and how all the pieces fit together to be able to resolve the queries we care about.  Moreover, the cdBG can be built on multiple reference sequences (transcripts, chromosomes, genomes), where each reference is given a distinct color (or colour, if you're of the British persuasion). The resulting structure, which also encodes the relationships between the cdBGs of the underlying reference sequences, is called the compacted, colored de Bruijn graph (ccdBG).  This is not, of course, the only variant of the dBG that has proven useful from an indexing perspective. The (pruned) dBG has also proven useful as a graph upon which to build a path index of arbitrary variation / sequence graphs, which has enabled very interesting and clever indexing schemes like that adopted in [GCSA2](https://github.com/jltsiren/gcsa2) (which we won't discuss further here, but which I hope to cover in more detail in a future post).  Also, thinking about sequence search in terms of the dBG has led to interesting representations for variation-aware sequence search backed by indexes like the vBWT (implemented in the excellent [gramtools](https://github.com/iqbal-lab-org/gramtools) package).

While existing hash-based indices based on the cdBG (and ccdBG) are very efficient for search, they typically occupy a large amount of space in memory (both during construction and even when built). As a result, to make use of such data structures on large reference sequences (e.g., the human genome) or collections of reference sequences (e.g., in a metagenomic context), one typically requires a very large memory machine — if the structures can be built at all. Pufferfish implements a new and much more compact data structure for indexing the ccdBG. While maintaining very efficient queries, this allows Pufferfish to index reference sequences while reducing the memory requirements considerably (by an order-of-magnitude or more). This greatly reduces the memory burden for indexing reference sequences and makes it possible to build hash-based indexes of sequences of size that were not previously feasible.

**about pufferfish development:**
Currently, Pufferfish is the software implementing this efficient ccdBG index, and allowing point (i.e., k-mer) queries.  Pufferfish is under active development, but we want to be as open (and as useful to as many people) as possible early on. However, we are also in the process of building higher-level tools (e.g., read mappers and aligners) around this index, so stay tuned!

**Dependency:** 
Pufferfish depends on sdsl-lite, to install from the pufferfish root directory execute following, 
>git clone https://github.com/simongog/sdsl-lite.git

>cd sdsl-lite

>./install.sh ../


**branches:**
The **master** branch of pufferfish is _not_ necessarily stable, but it should, at any given time contain a working version of the index.  That is, breaking changes should not be pushed to master.  The **develop** branch of pufferfish is guaranteed to be neither stable nor working at any given point, but a best-faith effort will be made to not commit broken code to this branch.  For feature branches, all bets are off.

For more details about pufferfish, please check out our [pre-print on bioRxiv](https://www.biorxiv.org/content/early/2017/09/21/191874), as well as the (updating) blog post [here](http://robpatro.com/blog/?p=494).


# Using Pufferfish <a name="using"></a>

## Using Pufferfish with BCALM2

You can use pufferfish with the unitig file provided by [BCALM2](https://github.com/GATB/bcalm).  Once you have downloaded and built bcalm, you can run it on your reference sequence file to produce a list of compacted unitigs like so:

```
bcalm -abundance-min 1 -max-memory <mem_in_MB> -nb-cores <cores_to_use> -in reference.fa  -out out_prefix -k <k>
```

This will generate, in addition to some other files, a file called `out_prefix.unitigs.fa`.  These are, unfortunately, not quite in the format required by pufferfish yet (unitigs can span reference boundaries).  To fix this, we must `pufferize` the file.  We can do that as such:

```
bcalm_pufferizer -k <k> -r reference.fa -u out_prefix.unitigs.fa
```

This will create a file called `out_prefix.unitigs.fa.pufferized.gfa`, on which you can then build the pufferfish index.
