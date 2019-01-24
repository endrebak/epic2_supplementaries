# epic2 and sicer give the same results

By default, epic2 and SICER should give slightly different results. This is
because epic2 chooses an effective genome size that is slightly different than
SICER. However, the regions found by SICER and epic2 are still the same in
location, length and poisson-score. They are also ranked the same way. The only
difference is where the cutoff for significance is placed, i.e. SICER might
accept the best ranked 60000 regions, while epic2 might only accept the 59000
first.

To illustrate that they do indeed give the same results, I have done three tests
on publically available data. To avoid the problem of different effective genome
fractions, I have used the same effective genome fraction for both SICER and
epic2.

These are the three tests done:

1. The canonical test data which comes with SICER.
2. Aorta H3K27me3 data (both ChIP and input) from the epigenome roadmap
3. H3K27me3 and H3K4me3 data without input from a Hacat cell line

These results should be reproducible by running this Snakefile. You can also
test other files by updating sample_sheet.txt

The three tests I do are all on bed files, since that is the only format SICER
accepts. The first test is on ChIP only, the last two on ChIP + Input.

In the below, there are some teensy differences, which I explain by pointing out
small bugs in SICER that I believe are largely insignificant.

## Ensuring that the regions are the same

To ensure that the regions I get are exactly the same I used bedtools subtract on each pair.

```
bedtools subtract -a epic2_results.bed -b SICER_results.bed > locs_only_in_epic2.bed
bedtools subtract -b epic2_results.bed -a SICER_results.bed > locs_only_in_SICER.bed
```

This will tell me if there is any difference in the two files with respect to location.

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

Here SICER finds 65652 significant regions and epic2 finds 65658.

The regions unique to epic2 are the six regions below:

```
chr2    163593600       163597399       0.04802850633859634     36.085391998291016      .       42      46      0.04828893765807152     0.36085391044616705
chr3    12212400        12217199        0.04802850633859634     36.085391998291016      .       42      46      0.04828893765807152     0.36085391044616705
chr4    27298600        27302599        0.04802850633859634     36.085391998291016      .       42      46      0.04828893765807152     0.36085391044616705
chr5    113599200       113604599       0.049435064196586616    32.85995864868164       .       50      56      0.04970160499215126     0.32859960198402405
chr10   54095000        54101199        0.04947330430150032     30.56852149963379       .       58      66      0.0497392937541008      0.3056852221488953
chr17   21561400        21565199        0.0473208948969841      27.528705596923828      .       74      86      0.04757966101169586     0.2752870619297028
```

As you can see, they all have a p-value just below 0.05. The discrepancy is likely due to this bug in SICER:
https://github.com/dariober/SICERpy/issues/10

These are the two regions unique to SICER:

```
chr8    145138199       145138599       112     44      2.0115472034441752e-29  3.5672320382800002      1.6437762125099466e-28
chr14   107040800       107043399       39      20      1.7864682548281746e-08  2.73275454361   2.592447672916509e-08
```

Investigating these two regions more closely in Python shows that neither should be an island:

```python
# Code to show that this should not be an island:

import pandas as pd
import pyranges as pr

# file is ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/Histone_H3K27me3/Aorta/UCSD.Aorta.H3K27me3.STL003.bed.gz
df = pd.read_table("/mnt/scratch/projects/epic_same_results/data/download/aorta/chip/aorta_chip.bed.gz", header=None)
df.columns = "Chromosome Start End Name Score Strand".split()
gr = pr.PyRanges(df)

def preprocess(gr, chromosome, start, end):
  # preprocess reads like epic/sicer does
  subset = gr[chromosome, start:end]
  ss = subset.df
  ss = ss.sort_values("Start End".split())
  ss = ss.drop_duplicates("Start End".split())
  ss.loc[ss.Strand == "+", "Start"] += 75
  ss.loc[ss.Strand == "-", "End"] -= 75
  bins = pd.concat([ss[ss.Strand == "-"].End, ss[ss.Strand == "+"].Start])
  bins = bins - (bins % 200)
  bins = bins[(bins < end) & (bins >= start)]
  return bins.value_counts().sort_index()

c14 = preprocess(gr, "chr14", 107040800, 107043399)
c8 = preprocess(gr, "chr8", 145131600, 145138599)

print(c14)

107040800     3
107041200     2
107041400     5
107041600     2
107041800     2
107042000     1
107042200     4
107042400     2
107042600     3
107042800     7
107043200     7
107043400     1
```

In the above you can see that there is a gap of 2200 - 1400 = 800 between the
two significant bins with 5 and 4 counts. A similar analysis shows a gap of 800
for c8 too. Since this is a bug that only affects the last potential island on a
chromosome I consider it insignificant, but it explains the discrepancy.

# Test 3: Hacat dataset

In this test we use two ChIP-files and one Input file. We are mixing H3K27me3
and H3K4me3, which is something you should not do in practice, but this does not
matter when ensuring that the programs produce the same results. Think of it as
throwing a curveball to both.

Collected bed files:

(They were aligned with `hisat2 --threads 24 -x /mnt/cargo/genomes/hisat2/hg38/genome --no-spliced-alignment -k 1 --no-discordant --no-mixed  -U Satrom-chIP-05-Input_TGACCA_L002_R1_001.fastq`)

Here epic2 finds 291239 significant regions, while SICER finds 291241.

The three regions unique to SICER are:

```
chr2    103063600       103065799       147     90      0.042517301272690655    1.1509594574200002      0.04998746438944515
chr8    18218800        18221199        147     90      0.042517301272690655    1.1509594574200002      0.04998746438944515
chr9    112048600       112049999       147     90      0.042517301272690655    1.1509594574200002      0.04998746438944515
```

As these all contain the same number of ChIP and Input reads, they seem like some weird artifact, and I take it as a good thing that epic2 did not consider them enriched.
