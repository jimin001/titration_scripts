# Titration Scripts
Scripts used to generate titration bams for DeepSomatic analysis. 

# Coverage
### coverage_titration.sh
This script takes in two bam files (`normal.bam`, `tumor.bam`) and downsamples the `tumor.bam` to 30x, 60x, 90x coverage and downsamples the `normal.bam` to 20x, 25x, 30x coverage.
If input bam file has a coverage less 30x, 60x, 90x (tumor) or 20x, 25x, 30x (tumor), the script will automatically skip downsampling to that coverage.
This script uses samtools with 30 threads.

run locally:
```
Usage: ./coverage_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> 
```

### downsample_bam.sh
This script is a generalized version of `coverage_titration.sh`, allowing user to specify which coverages to downsample to, and does not differentiate between tumor or normal bams.
This script uses samtools with 30 threads.

run locally:
```
Usage: ./downsample_bam.sh \
   -b <bam> \
   -c <starting_coverage> \
   -p <output_prefix> \
   -l <goal_coverage_list> \
   -o <output_directory>
```
# Purity

### purity_titration.sh
This script takes in two bam files  (`normal.bam`, `tumor.bam`) and creates purity titration bams.
Creates purity titration bams at : 60%, 70%, 80%, 90% tumor purity and 90%, 95% normal purity.
This script uses samtools with 30 threads.

run locally:
```
Usage: ./purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -g <tumor_goal_total_coverage> \
   -x <normal_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```

`tumor_purity_titration.sh` and `normal_purity_titration.sh` are almost identical, only difference is the naming convention for output files.

### tumor_purity_titration.sh
This script is an adaptation of `purity_titration.sh`, and only creates tumor purity titrations and allows user to specify which tumor purity percentages to titrate input bams to. 
This script uses samtools with 60 threads.

run locally:
```
Usage: ./tumor_purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -g <tumor_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> \
   -l <tumor_percent_list>
```
### normal_purity_titration.sh
This script is an adaptation of `purity_titration.sh`, and only creates normal purity titrations and allows user to specify which normal purity percentages to titrate input bams to. 
This script uses samtools with 30 threads.

run locally:
```
Usage: ./normal_purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -x <normal_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> \
   -l <normal_percent_list>
```
# Split Bams
These scripts were used as a pre-processing step before performing purity titrations, in order to ensure that there would not be duplicate reads in our titration set and evaluation set. 

These scripts with split a bam file into two sub-bams, so that there are no overlapping reads in each of the two split sub-bams. 
Two scripts are identical except for naming conventions specified for "normal" or "tumor" bams. These scripts use samtools with 60 threads and require large disk space and memory due to handling files in sam format during intermediate steps.

### split_bam_normal.sh
run locally:
```
Usage: ./split_bam_normal.sh \
   -n <normal_bam> \
   -q <normal_coverage> \
   -g <normal_goal_coverage> \ # Desired coverage for one of the two sub-bams. The other sub-bam will be the remainder of coverage from <normal_coverage>.
   -e <normal_evaluation_bam> \ # One of the two sub-bams, if already present, if not already present use "None". Coverage of <normal_evaluation_bam> must match <normal_goal_coverage>
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```

### split_bam_tumor.sh
run locally:
```
Usage: ./split_bam_tumor.sh \
   -n <tumor_bam> \
   -q <tumor_coverage> \
   -g <tumor_goal_coverage> \ # Desired coverage for one of the two sub-bams. The other sub-bam will be the remainder of coverage from <tumor_coverage>.
   -e <tumor_evaluation_bam> \ # One of the two sub-bams, if already present, if not already present use "None". Coverage of <tumor_evaluation_bam> must match <tumor_goal_coverage>
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```











