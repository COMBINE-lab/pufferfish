#!/usr/bin/env python
import sys, os
if len(sys.argv) < 4:
    print("alters BCALM's unitigs so that they fit pufferfish input. more specifically:")
    print("  the script considers the set B and E of all k-mers that are extremities of the reference genomes, respectively beginning and end of genomes.")
    print("  output is: modified unitigs such that each k-mer in B should be the beginning of an unitig, and each kmer in E should be end of an unitig.")
    print("  in order words, unitigs are split at kmers that are extremities of the reference sequences")
    exit("arguments: references.fa unitigs.fa k")

references=sys.argv[1]
unitigs=sys.argv[2]
k=int(sys.argv[3])

# https://www.biostars.org/p/710/#1412
from itertools import groupby
def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def normalize(kmer):
    return kmer if kmer < revcomp(kmer) else revcomp(kmer)


#ref_kmers=set()
# parse references
#for header, ref in fasta_iter(references):
#  for kmer in [ref[:k], ref[-k:]]:
#        ref_kmers.add(normalize(kmer))


ref_kmers={}
# parse references
# information we keep for each ref_kmer : start/end[0/1], unitigID, unitigLen
for header, ref in fasta_iter(references):
    ref_kmers[normalize(ref[:k])] = [0, -1, -1]
    ref_kmers[normalize(ref[-k:])] = [1, -1, -1]


# parse unitigs and split if necessary
# ASSUMPTION: we will never have to split a unitig twice, because really, what are the chances that a unitig spans an entire reference sequence?
# NOTE: unitigs are renumbered consecutively
# NOTE: former unitigs links are discarded
#output = open(unitigs+".pufferized.fa","w")
output = open(unitigs+".pufferized.gfa", "w")
nb_unitigs=0
def save_unitig(header,seq):
    global output, p, nb_unitigs
    #if "L:" in header:
    #    header = header.split("L:")[0] # strip links
    #assert(header.split(' ')[0].isdigit())
    #new_header = ">%d " % nb_unitigs + ' '.join(header.split(' ')[1:]) # renumber unitigs
    #output.write("S\t%s\t%s\n" % (new_header, seq))
    output.write("S\t%s\t%s\n" % (nb_unitigs, seq))
    nb_unitigs += 1
    return(nb_unitigs)

print("Started parsing and spliting unitigs .. ")
for header, unitig in fasta_iter(unitigs):
    seq = unitig
    already_split = False
    prev = 0
    for i in range(1,len(unitig)-k+1):
        kmer = normalize(unitig[i:i+k])
        if kmer in ref_kmers:
            #if already_split:
            #    exit("Error: this dataset violates the assumption that a unitig only needs to be split once.")
            #if i > 0 and i < len(unitig)-k+1:
            unitig_id=save_unitig(header,unitig[prev:i+k-1])
            if (ref_kmers[kmer][1]) != -1:
                exit("Error: a kmer cannot belong to more than one unitig")
            ref_kmers[kmer][1] = unitig_id
            ref_kmers[kmer][2] = i+k-1-prev
            #save_unitig(header,unitig[i:])
            #print("split",unitig,"at",i)
            #already_split = True
            prev = i            
    #if not already_split:
    save_unitig(header,unitig[prev:])


print("Started reconstructing the path .. ")
for ref in fasta_iter(references):
    output.write("\nP\t")
    i = 0    
    while (i < len(ref)-k+1):
        kmer = ref[i:i+k]
        normalizedkmer = normalize(kmer)        
        unitigInfo = ref_kmers[normalizedkmer]
        if unitigInfo[0] != 0:
            exit("Error: this kmer is supposed to be a starting kmer, but it's not")
        output.write("%s%s," % (unitigInfo[1], "+-"[kmer==normalizedkmer]))
        i += (unitigInfo[2]-k+1)    


output.close()
print("done. result is in: %s.pufferized.fa" % unitigs)
print("to update unitig links in the fasta header (necessary to get a GFA file), run:")
print("mv %s.pufferized.fa %s" % (unitigs,unitigs))
prefix = "[prefix]" if ("unitigs.fa" not in unitigs) else (unitigs.split(".unitigs.fa")[0]+".h5")
script_dir = os.path.dirname(os.path.realpath(__file__))
bcalm_path = "bcalm" if not os.path.isfile("%s/../build/bcalm" % script_dir) else os.path.abspath("%s/../build/bcalm" %script_dir)
print("%s -in %s -skip-bcalm -skip-bglue -redo-links" % (bcalm_path,prefix))
print("%s/convertToGFA.py %s %s.gfa %d" % (script_dir,unitigs,unitigs,k))
