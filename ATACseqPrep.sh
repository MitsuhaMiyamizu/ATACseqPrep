#!/bin/bash
wkpath1=/home/path/to/clean 
wkpath2=/home/path/to/clean/bowtie2
bowtie2_index=/home/path/refgenome/grch38/grch38
ref=/home/path/blacklist/ENCFF356LFX.bed
temp=/home/path/to/clean/temp
for i in $(cat ${wkpath1}/sra.id)
do
        sample_name=`basename $i`
        bowtie2 -p 6 --very-sensitive -X 2000 -x $bowtie2_index \
                -1 ${wkpath1}/${sample_name}_1_val_1.fq.gz -2 ${wkpath1}/${sample_name}_2_val_2.fq.gz \
                -S ${wkpath2}/raw/${sample_name}.raw.sam
        samtools sort -@ 6 -o ${wkpath2}/raw/${sample_name}.raw.bam ${wkpath2}/raw/${sample_name}.raw.sam
        samtools index -@ 6 ${wkpath2}/raw/${sample_name}.raw.bam
        samtools flagstat -@ 6 ${wkpath2}/raw/${sample_name}.raw.bam > ${wkpath2}/raw/${sample_name}.raw.stat
        ## bedtools bamtobed -i ${wkpath2}/${sample_name}.raw.bam > ${wkpath2}/${sample_name}.raw.bed

        sambamba markdup --overflow-list-size 600000 --tmpdir=${temp} \
        -r ${wkpath2}/raw/${sample_name}.raw.bam ${wkpath1}/rmdup/${sample_name}.rmdup.bam
        samtools index -@ 6 ${wkpath1}/rmdup/${sample_name}.rmdup.bam 

        mtReads=$(samtools idxstats ${wkpath1}/rmdup/${sample_name}.rmdup.bam | grep -E 'chrM|chrMT' | cut -f 3)
        totalReads=$(samtools idxstats ${wkpath1}/rmdup/${sample_name}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
        echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%' > ${wkpath1}/mtDNA/${sample_name}.txt

        samtools flagstat -@ 6 ${wkpath1}/rmdup/${sample_name}.rmdup.bam > ${wkpath1}/rmdup/${sample_name}.rmdup.stat
        samtools view -h -f 2 -q 30 ${wkpath1}/rmdup/${sample_name}.rmdup.bam | grep -E -v "chrM|chrMT" | samtools sort -O bam -@ 6 -o - > ${wkpath1}/last/${sample_name}.last.bam
        samtools index ${wkpath1}/last/${sample_name}.last.bam
        samtools flagstat ${wkpath1}/last/${sample_name}.last.bam > ${wkpath1}/last/${sample_name}.last.stat
        bedtools bamtobed -i ${wkpath1}/last/${sample_name}.last.bam > ${wkpath1}/last/${sample_name}.last.bed
        bedtools intersect -v -a ${wkpath1}/last/${sample_name}.last.bed -b $ref > ${wkpath1}/blacklistfilter/${sample_name}.filtered.bed
done
