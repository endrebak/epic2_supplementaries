from subprocess import check_output

prefix = "/mnt/scratch/projects/epic_bencmarks"

print(expand("{prefix}/data/sicer_results/hiya_single/test_chip-W200-G600-islands-summary-FDR1.0", prefix=prefix))
rule all:
    input:
        expand("{prefix}/data/sicer_results/hiya_single/ex/test-W200-G600-E{eval}.scoreisland", prefix=prefix, eval=[1, 10, 100, 1000, 10000]),
        expand("{prefix}/data/sicer_results/hiya/ex/test-W200-G600-islands-summary-FDR1.0", prefix=prefix)

rule run_sicer:
    input:
        chip = "../../SICER/ex/test.bed",
        input = "../../SICER/ex/control.bed",
    output:
        result = "{prefix}/data/sicer_results/hiya/ex/test_chip-W200-G600-islands-summary-FDR1.0",
        # memory = "{prefix}/data/sicer_results/{subset}/memory.csv",
    resources:
        instances = 1
    run:
        pwd = check_output("echo $(dirname $(dirname `pwd`))", shell=True).decode().strip()

        # how to solve python2 requirements?

        shell("sh ../../SICER/SICER.sh ../../SICER ex/test.bed ex/control.bed {prefix}/data/sicer_results/hiya hg19 1 200 150 0.74 600 1 {pwd}")


rule run_sicer_no_control:
    input:
        chip = "../../SICER/ex/test.bed",
    output:
        result = "{prefix}/data/sicer_results/hiya_single/ex/test-W200-G600-E{eval}.scoreisland",
        # memory = "{prefix}/data/sicer_results/{subset}/memory.csv",
    resources:
        instances = 1
    run:
        pwd = check_output("echo $(dirname $(dirname `pwd`))", shell=True).decode().strip()

        # how to solve python2 requirements?

        shell("sh ../../SICER/SICER-rb.sh ../../SICER ex/test.bed {prefix}/data/sicer_results/hiya_single hg19 1 200 150 0.74 600 {wildcards.eval} {pwd}")

# Find candidate islands exhibiting clustering ...
# /mnt/work/endrebak/software/anaconda/envs/py27/bin/python /home/endrebak/code/epic_paper/SICER/src/find_islands_in_pr.py -s hg19 -b /mnt/scratch/projects/epic_bencmarks/data/sicer_results/hiya/ex/test-W200.graph -w 200 -g 600 -t 0.74 -e 1000 -f /mnt/scratch/projects/epic_bencmarks/data/sicer_results/hiya/ex/test-W200-G600.scoreisland
# Species:  hg19
# Window_size:  200
# Gap size:  600
# E value is: 1000.0
# Total read count: 9904.0
# Genome Length:  3095693983
# Effective genome Length:  2290813547
# Window average: 0.000864670982322
# Window pvalue: 0.2
# Minimum num of tags in a qualified window:  1
# Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows
# Determine the score threshold from random background
# The score threshold is:  7.055
# Make and write islands
# 	chrY does not have any islands meeting the required significance
# 	chr21 does not have any islands meeting the required significance
# Total number of islands:  166

# Find significant islands with E-value 1.0 for ex/test...
# python /home/endrebak/code/epic_paper/SICER/src/find_islands_in_pr.py -s hg19 -b /mnt/scratch/projects/epic_bencmarks/data/sicer_results/hiya_single/ex/test-W200.graph -w 200 -g 600 -t 0.74 -e 1.0 -f /mnt/scratch/projects/epic_bencmarks/data/sicer_results/hiya_single/ex/test-W200-G600-E1.0.scoreisland
# Species:  hg19
# Window_size:  200
# Gap size:  600
# E value is: 1.0
# Total read count: 9904.0
# Genome Length:  3095693983
# Effective genome Length:  2290813547
# Window average: 0.000864670982322
# Window pvalue: 0.2
# Minimum num of tags in a qualified window:  1
# Generate the enriched probscore summary graph and filter the summary graph to get rid of ineligible windows
# Determine the score threshold from random background
# The score threshold is:  14.801
# Make and write islands
# 	chrY does not have any islands meeting the required significance
# 	chrX does not have any islands meeting the required significance
# 	chr11 does not have any islands meeting the required significance
# 	chr10 does not have any islands meeting the required significance
# 	chr17 does not have any islands meeting the required significance
# 	chr15 does not have any islands meeting the required significance
# 	chr14 does not have any islands meeting the required significance
# 	chr19 does not have any islands meeting the required significance
# 	chr18 does not have any islands meeting the required significance
# 	chr22 does not have any islands meeting the required significance
# 	chr20 does not have any islands meeting the required significance
# 	chr21 does not have any islands meeting the required significance
# 	chr6 does not have any islands meeting the required significance
# 	chr5 does not have any islands meeting the required significance
# 	chr4 does not have any islands meeting the required significance
# 	chr3 does not have any islands meeting the required significance
# 	chr2 does not have any islands meeting the required significance
# 	chr9 does not have any islands meeting the required significance
# 	chr8 does not have any islands meeting the required significance
# Total number of islands:  6
