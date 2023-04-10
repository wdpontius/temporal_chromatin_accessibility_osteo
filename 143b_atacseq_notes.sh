#perform alignment to human (hg19)
for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
bwa mem -t 4 /mnt/rds/genetics01/ScacheriLab/ars51/genomes/hg19/hg19.fa "$file"_R1.fastq.gz "$file"_R2.fastq.gz | samtools sort -@ 4 -O bam -T "$file".tmp -o "$file".bam -
samtools index "$file".bam
bamCoverage -b "$file".bam -o "$file"_nodup_normed.bw --ignoreDuplicates --normalizeTo1x 2451960000
done

#perform alignment to mouse (mm9) for XenofilteR
species=mm9
for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
bwa mem -t 4 /mnt/rds/genetics01/ScacheriLab/ars51/genomes/mm9/mm9.fa "$file"_R1.fastq.gz "$file"_R2.fastq.gz | samtools sort -@ 4 -O bam -T "$file"_"$species".tmp -o "$file"_"$species".bam -
samtools index "$file"_"$species".bam
done

for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_3; do
bwa mem -t 4 /mnt/rds/genetics01/ScacheriLab/ars51/genomes/hg19/hg19.fa "$file"_R1.fastq.gz "$file"_R2.fastq.gz | samtools sort -@ 4 -O bam -T "$file".tmp -o "$file".bam -
samtools index "$file".bam
bamCoverage -b "$file".bam -o "$file"_nodup_normed.bw --ignoreDuplicates --normalizeTo1x 2451960000
done

# change to python 3
bamCoverage -b 143b_in_vitro_ATAC_1.bam -o 143b_in_vitro_ATAC_1_nodup_normed.bw --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2451960000
bamCoverage -b 143b_in_vitro_ATAC_3.bam -o 143b_in_vitro_ATAC_3_nodup_normed.bw --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2451960000

# make biwigs for all of them
for file in 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
bamCoverage -b "$file"_sorted_Filtered.bam -o "$file"_sorted_Filtered_nodup_normed.bw --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2451960000
done


#did not initially work, giving error about wrong seq_length or something
#ran this command 'bwa index mm9.fa && echo "OK. Done."' to remake the mm9 index, and will see if this helps.
#People say this is a problem with using different version of BWA to do alignment and make genome index
#yay it worked!!

#use this loop to make the sample.list file for xenofilter
for file in 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
echo "$file"_sorted.bam,"$file"_mm9.bam >> sample_list_xenofilter.csv
done

for file in 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1; do
echo "$file"_sorted.bam,"$file"_sorted_mm9.bam >> sample_list_xenofilter_2.csv
done

#needs to be sorted by chromosome, ugh why does this take so long

for file in 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
samtools sort -m 5G "$file".bam -o "$file"_sorted.bam
done

for file in 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
samtools sort -m 5G "$file"_mm9.bam -o "$file"_sorted_mm9.bam
done

#ran XenofilteR
library(XenofilteR)

sample.list <- read.csv("sample_list_xenofilter_2.csv", header = FALSE)
bp.param <- SnowParam(workers = 4, type = "SOCK")
XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param)

#make bigwigs from filtered bams
for file in 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1; do
bamCoverage -b "$file"_sorted_Filtered.bam -o "$file"_sorted_Filtered_nodup_normed.bw --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2451960000
done

for file in in-vivo_d1_1_S1 in-vivo_d1_2_S2 in-vivo_d1_3_S3 in-vivo_d22_1_S4 in-vivo_d22_2_S5 in-vivo_d22_3_S6 in-vivo_d22_4_S7 in-vivo_d22_5_S8; do
bamCoverage -b "$file"_sorted_Filtered.bam -o "$file"_sorted_Filtered_nodup_normed.bw --ignoreDuplicates --normalizeTo1x 2451960000
done

#call peaks - genrich requires sorting bams by query name, so have to resort...ugh
#sorting didn't take as long because bams are somewhat small after filtering
for file in 143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
samtools sort -n -m 20G "$file"_sorted_Filtered.bam -o "$file"_sorted-for-peaks_Filtered.bam
~/programs/Genrich/Genrich -t "$file"_sorted-for-peaks_Filtered.bam -o "$file"_sorted_Filtered.narrowPeak -j -r -e chrM
done


for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3; do
samtools sort -m 5G "$file".bam -o "$file"_sorted.bam
done

for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3; do
#samtools sort -n -m 20G "$file"_sorted.bam -o "$file"_sorted-for-peaks_Filtered.bam
~/programs/Genrich/Genrich -t "$file"_sorted-for-peaks_Filtered.bam -o "$file"_sorted_Filtered.narrowPeak -j -r -e chrM
done

# Remove duplicates
for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3; do
samtools rmdup "$file"_sorted.bam "$file"_sorted_nodup.bam
done

for file in  143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
samtools rmdup "$file"_sorted_Filtered.bam "$file"_sorted_Filtered_nodup.bam
done

for file in  143b_in_vivo_Day_1_ATAC_1 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3; do
samtools rmdup "$file"_sorted_Filtered.bam "$file"_sorted_Filtered_nodup.bam
done

# Downstream RPKM based analysis --------------------------------

#combine all 3 replicate in vitro peaks, and in vivo peaks (not the 95% mouse one)


#calc rpkm
for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3; do
perl ~/scripts/calc_rpkm.pl Filtered_bams_round1/143b_invitro-invivo_universe_autosome_noblacklist.bed "$file"_sorted_nodup.bam "$file" &
done

for file in 143b_in_vivo_Day_1_ATAC_combined 143b_in_vivo_Day_1_ATAC_2 143b_in_vivo_Day_1_ATAC_3 143b_in_vivo_Day_22_ATAC_1 143b_in_vivo_Day_22_ATAC_2 143b_in_vivo_Day_22_ATAC_3 143b_in_vivo_Day_22_ATAC_4 143b_in_vivo_Day_22_ATAC_5; do
perl ~/scripts/calc_rpkm.pl 143b_invitro-invivo_universe_autosome_noblacklist.bed "$file"_sorted_Filtered_nodup.bam "$file" &
done

# D1 looked like shit. Try combining all reps and making bigwig out of that.

# Combine D1 replicates ------------------------------------------------

samtools merge 143b_in_vivo_Day_1_ATAC_combined.bam 143b_in_vivo_Day_1_ATAC_1_sorted_Filtered_nodup.bam 143b_in_vivo_Day_1_ATAC_2_sorted_Filtered_nodup.bam 143b_in_vivo_Day_1_ATAC_3_sorted_Filtered_nodup.bam

bamCoverage -b 143b_in_vivo_Day_1_ATAC_combined.bam -o 143b_in_vivo_Day_1_ATAC_combined.bw --ignoreDuplicates --normalizeUsing RPGC --effectiveGenomeSize 2451960000

samtools sort -n -m 20G 143b_in_vivo_Day_1_ATAC_combined.bam -o 143b_in_vivo_Day_1_ATAC_combined_for-peaks.bam
~/programs/Genrich/Genrich -t 143b_in_vivo_Day_1_ATAC_combined_for-peaks.bam -o 143b_in_vivo_Day1_ATAC_combined_sorted_Filtered.narrowPeak -j -r -e chrM
done

cat 143b_in_vivo_Day1_ATAC_combined_sorted_Filtered.narrowPeak 143b_in_vivo_Day_22_ATAC_1_sorted_Filtered.narrowPeak 143b_in_vivo_Day_22_ATAC_2_sorted_Filtered.narrowPeak 143b_in_vivo_Day_22_ATAC_3_sorted_Filtered.narrowPeak 143b_in_vivo_Day_22_ATAC_4_sorted_Filtered.narrowPeak 143b_in_vivo_Day_22_ATAC_5_sorted_Filtered.narrowPeak > invivo_peaks.bed
cat invivo_peaks.bed ../143b_in_vitro_ATAC_1_sorted_Filtered.narrowPeak ../143b_in_vitro_ATAC_2_sorted_Filtered.narrowPeak ../143b_in_vitro_ATAC_3_sorted_Filtered.narrowPeak

# get rid of non chr 1-23 in R

for file in 143b_in_vitro_ATAC_1 143b_in_vitro_ATAC_2 143b_in_vitro_ATAC_3; do
perl ~/scripts/calc_rpkm.pl Filtered_bams_round1/143b_invitro_invivo_day1-combined_universe_noblacklist_autosome.bed "$file"_sorted_nodup.bam "$file"_day1_combined &
done

for file in 143b_in_vivo_Day_1_ATAC_combined 143b_in_vivo_Day_22_ATAC_1_sorted_Filtered_nodup 143b_in_vivo_Day_22_ATAC_2_sorted_Filtered_nodup 143b_in_vivo_Day_22_ATAC_3_sorted_Filtered_nodup 143b_in_vivo_Day_22_ATAC_4_sorted_Filtered_nodup 143b_in_vivo_Day_22_ATAC_5_sorted_Filtered_nodup; do
perl ~/scripts/calc_rpkm.pl 143b_invitro_invivo_day1-combined_universe_noblacklist_autosome.bed "$file".bam "$file"_day1_combined &
done
