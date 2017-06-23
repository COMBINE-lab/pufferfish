# What is pufferfish?

**short answer** : Pufferfish is a new time and memory-efficient data structure for indexing the compacted, colored de Bruijn graph (ccdBG).  

**long answer** : Though the de Bruijn Graph (dBG) has enjoyed tremendous popularity as an assembly and comparision data structure, it has only relatively recently begun to see use as an index of the reference sequences (e.g. [deBGA](https://github.com/HongzheGuo/deBGA), [kallisto](https://github.com/pachterlab/kallisto)).  Particularly, these tools index the _compacted_ dBG (cdBG), in which all non-branching paths are collapsed into individual nodes and labeled with the string they spell out.  This data structure is particularly well-suited for representing repetitive reference sequences, since a single contig in the cdBG represents all occurences of the repeated sequence.  The original positions in the reference can be recovered with the help of an auxiliary "contig table", that maps each contig to the reference sequence, position (and orientation) where it appears as a substring.  Moreoever, the cdBG can be built on multiple reference sequences (transcripts, chromosomes, genomes), where each reference is given a distinct color.  The resulting structure, which also encodes the relationships between the cdBGs of the underlying refernece sequences is called the compacted, colored de Bruijn graph.

While existing indices based on the cdBG (and ccdBG) are very efficient for search, they typically occupy a large amount of space in memory (both during construction and even when built).  As a result, to make use of such data structures on large reference sequences (e.g., the human genome) or collections of reference sequences (e.g., in a metagenomic context), one typically requires a very large memory machine --- if the structures can be built at all.  Pufferfish implements a new and much more compact data structure for indexing the ccdBG.  While maintaining very efficient queries, this allows Pufferfish to index reference sequences while reducing the memory requirements considerably (by an order-of-magnitude or more).  This considerably reduces the memory burden for indexing reference sequences and makes it possible to index sequences of size that were not previously feasible.

**about pufferfish development:**
Currently, Pufferfish is the software implementing this efficient ccdBG index, and allowing point (i.e. k-mer) queries.  Pufferfish is under active development, but we want to be as open (and as useful to as many people) as possible early on. However, we are also in the process of building higher-level tools (e.g. read mappers and aligners) around this index, so stay tuned!

**branches:**
The **master** branch of pufferfish is _not_ necessarily stable, but it should, at any given time contain a working version of the index.  That is, breaking changes should not be pushed to master.  The **develop** branch of pufferfish is guaranteed to be neither stable nor working at any given point, but a best-faith effort will be made to not commit broken code to this branch.  For feature branches, all bets are off.

# Some details about pufferfish

Pufferfish supports (i.e., implements) two types of indices — dense and sparse.  Both of these make use of an efficient [minimum perfect hash implementation](https://github.com/rizkg/BBHash) to vastly reduce the memory requirements compared to other hashing-based dBG, colored-dBG and compacted-colored-dBG indices.  The basic (dense) pufferfish index can be visualized as the image below:


