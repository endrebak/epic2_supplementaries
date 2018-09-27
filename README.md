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


chr1    34501400        34508399        1027    flank   3       1       1.0
chr1    34508400        34509599        1027    trunk   3       2       2.0
chr1    34509600        34510599        1027    valley  3       1       1.0
chr1    34510600        34516599        1027    trunk   3       2       2.0
chr1    34516600        34518199        1027    valley  3       1       1.0
chr1    34518200        34518799        1027    trunk   3       2       2.0
chr1    34518800        34520599        1027    valley  3       1       1.0
chr1    34520600        34529799        1027    trunk   3       2       2.8260869565217392
chr1    34529800        34531999        1027    valley  3       1       1.0
chr1    34532000        34541999        1027    trunk   3       2       2.76
chr1    34542000        34542999        1027    flank   3       1       1.0

34542999 - 34501400 = 41599
