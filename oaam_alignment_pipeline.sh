#!/bin/bash
#SBATCH--job-name=oaam_alignment
#SBATCH--output=/projects/foster_lab/john/birds/oaam_pop_gen/samples_2021/alignment%a.log
#SBATCH--time=10:00:00
#SBATCH--chdir=/projects/foster_lab/john/birds/oaam_pop_gen/samples_2021/
#SBATCH--cpus-per-task=10
#SBATCH--mem=70G
#SBATCH--array=1-9

module load bwa/0.7.17
module load samtools/1.11
module load picard/2.24.1
#.txt is a list of directories/files for SLURM array. Task ID will keep an index of each line
line_N=$( awk "NR==$SLURM_ARRAY_TASK_ID" third_run.txt)
REF=/projects/foster_lab/john/birds/oaam/oaam.v1.simple.fa
IDX=/projects/foster_lab/john/birds/oaam/oaam.bwa.index
READS_DIR=/projects/foster_lab/john/birds/oaam_pop_gen/samples_2021
prefix=$(echo $line_N | awk -F "," '{print $1}')
SM=$(echo $line_N | awk -F "," '{print $4}')
ID=$(echo $line_N | awk -F "," '{print $2}')
LB=$(echo $line_N | awk -F "," '{print $6}')
PU=$(echo $line_N | awk -F "," '{print $3}')
PL=$(echo $line_N | awk -F "," '{print $5}')

srun bwa mem -M -t 18 -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$PL" $IDX \
        $READS_DIR/${prefix}_trimmed_1P.fq.gz \
        $READS_DIR/${prefix}_trimmed_2P.fq.gz |
        samtools view -@ 10 -f 3 -bS - | samtools sort -@ 10 -O bam -o ${prefix}.bam

srun samtools index -@10 ${prefix}.bam
#mark and remove duplicates
srun java -Xmx90G -jar $PICARD MarkDuplicates \
  INPUT=${prefix}.bam \
  OUTPUT=${prefix}.second.bam \
  METRICS_FILE=${prefix}.metrics.txt \
  ASSUME_SORTED=true \
  REMOVE_DUPLICATES=true

srun samtools index -@10 ${prefix}.second.bam
rm ${prefix}.bam
rm ${prefix}.bam.bai
#load gatk3 for indel re-aligning
module purge
module load anaconda3
conda activate gatk3_env

srun gatk -T RealignerTargetCreator \
  -R $REF \
  -I ${prefix}.second.bam -o ${prefix}.target.intervals

srun gatk -T IndelRealigner \
-R $REF \
  -targetIntervals ${prefix}.target.intervals \
  -I ${prefix}.second.bam -o ${prefix}.indel.realigner.bam
conda deactivate
rm ${prefix}.second.bam
rm ${prefix}.second.bam.bai
