while getopts b:c:p:l:o: flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        c) starting_coverage=${OPTARG};;
        p) output_prefix=${OPTARG};;
        l) goal_coverage_list=${OPTARG};;
        o) output_directory=${OPTARG};;
    esac
done

echo "This simple script takes an input bam and downsamples it to desired coverage using samtools view -s, and checks coverages of downsampled bams with samtools depth."
echo "---------------------------------------------------------------------"
echo "INPUT VARIABLES:"
echo "Path to input bam: $bam";
echo "Coverage of input bam: $starting_coverage";
echo "Prefix to name downsampled bams: $output_prefix";
echo "List of coverages to downsample to, separated by spaces: $goal_coverage_list";
echo "Output directory: $output_directory";

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

echo "---------------------------------------------------------------------"
echo "Creating downsampled bams at coverages: $goal_coverage_list"
echo "-------------------------------------"

echo
echo "Downsampling bams..."
for goal_coverage in $goal_coverage_list
do
    echo "Downsample to: ${goal_coverage}x"

    # if original bam is within 3x of the goal coverage, then skip downsampling
    if [ $(($starting_coverage - $goal_coverage)) -lt 3 ]
    then
        echo "Downsampling unnecessary"
    else
        echo "Downsampling with samtools..."
        echo "Using 32 as seed for samtools view -s"
        fraction=$( bc <<< "scale=3;$goal_coverage/$starting_coverage + 32")
        echo "Fraction to downsample by: $fraction"

        samtools view -@30 -b -s ${fraction} ${bam} > ${output_directory}/${output_prefix}.${goal_coverage}x.bam
        samtools index -@30 ${output_directory}/${output_prefix}.${goal_coverage}x.bam

        echo "Downsampled bam name: ${output_prefix}.${goal_coverage}x.bam"

        echo "Verifying downsampling with samtools depth..."
        samtools depth ${output_directory}/${output_prefix}.${goal_coverage}x.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
        echo "-------------------------------------"
        echo
    fi
done
