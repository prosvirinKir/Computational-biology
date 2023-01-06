# Computational-biology project

# Task

Implement algorithms that takes two **FASTA** files (reads) and outputs
a list of SNPs.

# Approach

## Description

1.  Count k-mer frequencies in both files.

2.  With high probability k-mers with low frequencies were gotten from
    reads with sequencing errors. We will not consider them.

3.  If k is high enough we can assume that each k-mer meets in genome
    only once. Therefore, k-mers with mutation should present only in
    one genome. So, we make two lists of k-mers which elements present
    in one file, but not present in another.

4.  Build the longest possible sequence from selected k-mers for each
    genome.

5.  Compare resulted parts of genomes to find SNPs.


### Execution: python3 contact_map_builder.py readsA.fasta readsB.fasta
