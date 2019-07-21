# Index
 * [What is pufferfish?](#whatis)
 * [Using pufferfish?](#using)
 
# What is pufferfish? <a name="whatis"></a>

**short answer** : Pufferfish is a new time and memory-efficient data structure for indexing a compacted, colored de Bruijn graph (ccdBG).  You can read more about pufferfish in the [paper](https://academic.oup.com/bioinformatics/article/34/13/i169/5045749), which appeared at ISMB 2018.

**long answer** : 
Though the de Bruijn Graph (dBG) has enjoyed tremendous popularity as an assembly and sequence comparison data structure, it has only relatively recently begun to see use as an index of the reference sequences (e.g. [deBGA](https://github.com/HongzheGuo/deBGA), [kallisto](https://github.com/pachterlab/kallisto)). Particularly, these tools index the compacted dBG (cdBG), in which all non-branching paths are collapsed into individual nodes and labeled with the string they spell out. This data structure is particularly well-suited for representing repetitive reference sequences, since a single contig in the cdBG represents all occurrences of the repeated sequence. The original positions in the reference can be recovered with the help of an auxiliary "contig table" that maps each contig to the reference sequence, position, and orientation where it appears as a substring. The [deBGA paper](https://academic.oup.com/bioinformatics/article-abstract/32/21/3224/2415064/deBGA-read-alignment-with-de-Bruijn-graph-based?redirectedFrom=fulltext) has a nice description how this kind of index looks (they call it a unipath index, because the contigs we index are unitigs in the cdBG), and how all the pieces fit together to be able to resolve the queries we care about.  Moreover, the cdBG can be built on multiple reference sequences (transcripts, chromosomes, genomes), where each reference is given a distinct color (or colour, if you're of the British persuasion). The resulting structure, which also encodes the relationships between the cdBGs of the underlying reference sequences, is called the compacted, colored de Bruijn graph (ccdBG).  This is not, of course, the only variant of the dBG that has proven useful from an indexing perspective. The (pruned) dBG has also proven useful as a graph upon which to build a path index of arbitrary variation / sequence graphs, which has enabled very interesting and clever indexing schemes like that adopted in [GCSA2](https://github.com/jltsiren/gcsa2) (which we won't discuss further here, but which I hope to cover in more detail in a future post).  Also, thinking about sequence search in terms of the dBG has led to interesting representations for variation-aware sequence search backed by indexes like the vBWT (implemented in the excellent [gramtools](https://github.com/iqbal-lab-org/gramtools) package).

While existing hash-based indices based on the cdBG (and ccdBG) are very efficient for search, they typically occupy a large amount of space in memory (both during construction and even when built). As a result, to make use of such data structures on large reference sequences (e.g., the human genome) or collections of reference sequences (e.g., in a metagenomic context), one typically requires a very large memory machine — if the structures can be built at all. Pufferfish implements a new and much more compact data structure for indexing the ccdBG. While maintaining very efficient queries, this allows Pufferfish to index reference sequences while reducing the memory requirements considerably (by an order-of-magnitude or more). This greatly reduces the memory burden for indexing reference sequences and makes it possible to build hash-based indexes of sequences of size that were not previously feasible.

**about pufferfish development:**
Currently, Pufferfish is the software implementing this efficient ccdBG index, and allowing point (i.e., k-mer) queries.  Pufferfish is under active development, but we want to be as open (and as useful to as many people) as possible early on. However, we are also in the process of building higher-level tools (e.g., read mappers and aligners) around this index, so stay tuned!


**branches:**
The **master** branch of pufferfish is _not_ necessarily stable, but it should, at any given time contain a working version of the index.  That is, breaking changes should not be pushed to master.  The **develop** branch of pufferfish is guaranteed to be neither stable nor working at any given point, but a best-faith effort will be made to not commit broken code to this branch.  For feature branches, all bets are off.

For more details about pufferfish, please check out our [paper](https://academic.oup.com/bioinformatics/article/34/13/i169/5045749), as well as the blog post [here](http://robpatro.com/blog/?p=494).

# Building Pufferfish <a name="building"></a>
To build the pufferfish do the following,

```
>git clone git@github.com:COMBINE-lab/pufferfish.git
>cd pufferfish
> mkdir build
> cd build
> cmake ../
> make
```

# Using Pufferfish <a name="using"></a>

**Programs used within pufferfish:**

Building a pufferfish index requires first having a compacted de Bruijn graph, for which we use a modified version of [TwoPaCo](https://github.com/medvedevgroup/TwoPaCo).  However, some modification of the TwoPaCo output is required for pufferfish to properly index the graph (e.g. a k-mer must appear at most once in the graph and palindromic contigs output by TwoPaCo must be removed).  Thus we rely on a modified version of TwoPaCo which we bundle with pufferfish in the `external` directory.

To choose an appropriate filter size to pass to TwoPaCo to build the compacted dBG, we make use the the hyper-log-log implementation of [ntCard](https://github.com/bcgsc/ntCard).  Because we use this as a library instead of an executable, and to avoid an external dependency to simply call one function, we bundle a modified version of that code with pufferfish and also include it in the `external` directory.

We are also dependent on [SeqLib](https://github.com/walaj/SeqLib) and hence all the libraries that it is dependent on such as `bz2`, `lzma`, and `z` for mapping part. So it is required to install these libraries on the system as well and also update the CMakeLists.txt file in `src` directory and change the line for setting variable `SEQLIBDIR` statistically.

## Core Operations

**Building a pufferfish index**

To build a pufferfish inded, you can use the `index` command.  It is used like so:

```
pufferfish index -k <ksize> -r <fasta_file_to_index> -o <pufferfish index directory>
```

There are also optional parameters including `-s` (the ability to build a sparser and smaller index), `-p` (control the number of threads used during construction), and `-f` (to provide an explicit filter size for TwoPaCo dBG construction).

# Mapping using Pufferfish
We can generate different types of output including sam. If you have samtools installed on your system you can run the first command in the puff_align.bash file to generate a bam file for mapping a set of reads to the pufferfish index. The rest of the commands do quantification on the bam file using Salmon and then validation of the results on reference set and taxonomic tree using two python scripts in scripts section.

For the address of the reads, index, truth files, etc. you need to have a config file like `microbiome_config.json`.

