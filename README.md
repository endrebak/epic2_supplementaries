# epic paper

- speed
- memory usage (can use time to find it?)
- visualization
- diff. chip-seq

- last ned og ha SICER som en del av koden i pipelinen: http://home.gwu.edu/~wpeng/SICER_V1.1.tgz

# Reproduce:

Need to change PATHTO in SICER/SICER.sh

# New features

General
- accept bam files
- accept paired end data
- UCSC genomes included

Differential ChIP-Seq
- store matrix of read counts for all genomic bins (apply linear models to
  regions)

Ease of use
- require much fewer command line args
- has the effective genome size for any UCSC genome
    - can find the effective genome size for any fasta file
- automatically find effective genome size, based on read length and genome, no
  need for the user to input

Visualization
 - create individual bigwigs
 - create bigwigs, not wigs
 - create summary bed file containing all enriched regions for display in the
   UCSC genome browser.
