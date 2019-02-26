# epic2 and sicer give the same results

By default, epic2 and SICER should give slightly different results. This is
because epic2 chooses an effective genome size that is slightly different than
SICER. Also, SICER contains a logical bug: when computing p-values, it also uses
those reads falling outside of the genome boundaries. epic2 has a flag that
allows you to run the original SICER algorithm, but it should only be used for
debugging, as this is unlikely to be the desired behavior. In the below tests
the `--original-algorithm` flag is used.

To illustrate that SICER and epic2 do indeed give the same results, I have done
three tests on publically available data. To avoid the problem of different
effective genome fractions, I have used the same effective genome fraction for
both SICER and epic2.

These are the three tests done:

1. The canonical test data which comes with SICER.
2. Aorta H3K27me3 data (both ChIP and input) from the epigenome roadmap
3. H3K27me3 and H3K4me3 data with input from a Hacat cell line

These results should be reproducible by running this Snakefile. You can also
test other files by updating sample_sheet.txt

The three tests I do are all on bed files, since that is the only format SICER
accepts. All tests are on ChIP + Input.

## Ensuring that the regions are the same

To ensure that the regions I get are exactly the same I used bedtools subtract on each pair.

```
bedtools subtract -a epic2_results.bed -b SICER_results.bed > locs_only_in_epic2.bed
bedtools subtract -b epic2_results.bed -a SICER_results.bed > locs_only_in_SICER.bed
```

This will tell me if there is any difference in the two files with respect to location.

## Ensuring the ordering is the same

To ensure that the regions have the exact same rank, I take the top 1000 lowest
FDR-values from each SICER/epic2 run and do an intersection. The end result
should contain 1000 regions if the rank is the same.

## Parameters used in the below tests

* Egf: 0.85
* Genome: hg38
* Remove all duplicate reads
* Bin size: 200
* Max gaps: 3 (or 600 nucleotides)
* FDR-cutoff: 0.05

## Test 1: SICER example data, without input

ChIP: https://raw.githubusercontent.com/biocore-ntnu/epic2/master/examples/test.bed

Input: https://raw.githubusercontent.com/biocore-ntnu/epic2/master/examples/control.bed

Here SICER finds 166 islands and so does SICER.

They are exactly the same.

## Test 2: epigenome roadmap aorta data

Chip: ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/Histone_H3K27me3/Aorta/UCSD.Aorta.H3K27me3.STL003.bed.gz

Input: ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/ChIP-Seq_Input/Aorta/UCSD.Aorta.Input.STL002.bed.gz

SICER and epic2 find 66011 regions, and they look exactly the same, down to the final decimal:

```
head -3 /mnt/scratch/projects/epic_same_results/data/sicer_results/aorta/aorta_chip-W200-G600-islands-summary-FDR1.0
chr1    10000   10599   23      13      3.9069431024123045e-05  2.47942226836   4.3105309354482855e-05
chr1    567400  568599  38      58      1       0.918166866996  1
chr1    569800  569999  25      52      1       0.673756051186  1
biocore-home ~/c/e/w/show_same_results (master DU=) head -3 /mnt/scratch/projects/epic_same_results/data/epic2_results/aorta/fdr_list.csv
#Chromosome     Start   End     PValue  Score   Strand  ChIPCount       InputCount      FDR     log2FoldChange
chr1    10000   10599   3.9069431024123045e-05  131.00039965542217      .       23      13      4.3105309354482855e-05  1.3100039965542218
chr1    567400  568599  1.0     -12.31717230456854      .       38      58      1.0     -0.12317172304568541
```

The order is also the same, as the intersection of the top 1000 results from each has a length of a 1000:

```
wc -l /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/*_aorta_significant.bed
 1000 45630 /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/epic2_aorta_significant.bed
 1000 45632 /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/sicer_aorta_significant.bed
```

# Test 3: Hacat dataset

In this test we use two ChIP-files and one Input file. We are mixing H3K27me3
and H3K4me3, which is something you should not do in practice, but this does not
matter when ensuring that the programs produce the same results.

Collection and citation info: https://zenodo.org/record/2548491/

H3K27me3.bed.gz: https://zenodo.org/record/2548491/files/Satrom-H3K27me3.bed.gz
H3K4me3.bed.gz: https://zenodo.org/record/2548491/files/Satrom-H3K4me3.bed.gz
Input.bed.gz: https://zenodo.org/record/2548491/files/Satrom-Input.bed.gz

(They were aligned with `hisat2 --threads 24 -x /mnt/cargo/genomes/hisat2/hg38/genome --no-spliced-alignment -k 1 --no-discordant --no-mixed  -U Satrom-chIP-05-Input_TGACCA_L002_R1_001.fastq`)

The islands are the same in number (342411) and exactly equal down to the last decimal:

```
head -3 /mnt/scratch/projects/epic_same_results/data/epic2_results/satrom/fdr_list.csv
#Chromosome     Start   End     PValue  Score   Strand  ChIPCount       InputCount      FDR     log2FoldChange
chr1    28400   30399   3.354391451012975e-120  268.6885321202972       .       265     29      4.9490571214294755e-119 2.6868853212029724
chr1    38000   51399   3.76376192797827e-10    28.25109196539736       .       1046    606     5.874535492271546e-10   0.2825109196539736
biocore-home ~/c/e/w/show_same_results (master DU=) head -3 /mnt/scratch/projects/epic_same_results/data/sicer_results/satrom/satrom_chip-W200-G600-islands-summary-FDR1.0
chr1    28400   30399   265     29      3.354391451012975e-120  6.43921723185   4.9490571214294755e-119
chr1    38000   51399   1046    606     3.76376192797827e-10    1.21630995399   5.874535492271546e-10
chr1    52200   73799   1516    1000    0.005218518673214944    1.06827828823   0.006389013758256891
```

The intersection of the top 1000 results from each are also the same:

```
wc -l /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/*_satrom_significant.bed
 1000 27944 /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/epic2_satrom_significant.bed
 1000 27944 /mnt/scratch/projects/epic_same_results/data/sicer_results/same_order/sicer_satrom_significant.bed
```
