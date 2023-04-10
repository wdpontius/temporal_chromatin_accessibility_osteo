# Step 1 - align to both genomes -------------------
for file in in_vitro_143b_3 in_vitro_143b_4 in_vitro_143b_5 d1_1_143b d1_2_143b d1_3_143b d22_1_143b d22_2_143b d22_3_143b d22_4_143b d22_5_143b; do
hisat2 -x ~/hg19/hg19 -1 "$file"_R1.fastq.gz -2 "$file"_R2.fastq.gz -p 4 | samtools view -Sbh > "$file"_hg19.bam
done

for file in in_vitro_143b_3 in_vitro_143b_4 in_vitro_143b_5 d1_1_143b d1_2_143b d1_3_143b d22_1_143b d22_2_143b d22_3_143b d22_4_143b d22_5_143b; do
hisat2 -x ~/mm9/mm9 -1 "$file"_R1.fastq.gz -2 "$file"_R2.fastq.gz -p 4 | samtools view -Sbh > "$file"_mm9.bam
done

# Step 2 - run xenofilter --------------------------
# First need to sort by coordinates
for file in in_vitro_143b_3 in_vitro_143b_4 in_vitro_143b_5 d1_1_143b d1_2_143b d1_3_143b d22_1_143b d22_2_143b d22_3_143b d22_4_143b d22_5_143b; do
samtools sort -m 5G "$file"_hg19.bam -o "$file"_hg19_sorted.bam
done

for file in in_vitro_143b_3 in_vitro_143b_4 in_vitro_143b_5 d1_1_143b d1_2_143b d1_3_143b d22_1_143b d22_2_143b d22_3_143b d22_4_143b d22_5_143b; do
samtools sort -m 5G "$file"_mm9.bam -o "$file"_mm9_sorted.bam
done

# This is in R
d1_1_143b_hg19_sorted.bam,d1_1_143b_mm9_sorted.bam
d1_2_143b_hg19_sorted.bam,d1_2_143b_mm9_sorted.bam
d1_3_143b_hg19_sorted.bam,d1_3_143b_mm9_sorted.bam
d22_1_143b_hg19_sorted.bam,d22_1_143b_mm9_sorted.bam
d22_2_143b_hg19_sorted.bam,d22_2_143b_mm9_sorted.bam
d22_3_143b_hg19_sorted.bam,d22_3_143b_mm9_sorted.bam
d22_4_143b_hg19_sorted.bam,d22_4_143b_mm9_sorted.bam
d22_5_143b_hg19.bam,d22_5_143b_mm9.bam


library(XenofilteR)

sample.list <- read.csv("xenofilter_list.csv", header = FALSE)
bp.param <- SnowParam(workers = 4, type = "SOCK")
XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param)


# Step 3 - run stringtie 
#run stringtie first
for file in `ls | grep .Filtered.bam$`; do stringtie "$file" -e -G ~/ensembl/hg19.ncbiRefSeq.gtf -p 16 -o "$file".gtf; done
for file in `ls | grep in_vitro | grep hg19_sorted.bam$`; do stringtie "$file" -e -G ~/ensembl/hg19.ncbiRefSeq.gtf -p 16 -o "$file".gtf; done


#file for prepDE.py
invitro_1	in_vitro_143b_3_hg19_sorted.bam.gtf
invitro_2	in_vitro_143b_4_hg19_sorted.bam.gtf
invitro_3	in_vitro_143b_5_hg19_sorted.bam.gtf
invivo_d1_1	d1_1_143b_hg19_sorted_Filtered.bam.gtf
invivo_d1_2	d1_2_143b_hg19_sorted_Filtered.bam.gtf
invivo_d1_3	d1_3_143b_hg19_sorted_Filtered.bam.gtf
invivo_d22_1	d22_1_143b_hg19_sorted_Filtered.bam.gtf
invivo_d22_2	d22_2_143b_hg19_sorted_Filtered.bam.gtf
invivo_d22_3	d22_3_143b_hg19_sorted_Filtered.bam.gtf
invivo_d22_4	d22_4_143b_hg19_sorted_Filtered.bam.gtf
invivo_d22_5	d22_5_143b_hg19_sorted_Filtered.bam.gtf

python2 ~/scripts/prepDE.py -i prep_de_file.txt





stringtie invivo_d1_3_h19.sorted -p 8 -o invivo_d1_3_h19.sorted.gtf

#merge gtf with stringtie merge
stringtie --merge prep_de_file.txt -o merged.gtf

#redo stringtie with merged gtf reference
for file in `ls | grep .sorted$`; do stringtie  $file -e -G merged.gtf -p 16 -o "$file".gtf.2 ; done
for file in `ls | grep .sorted.bam`; do stringtie -e -G ~/ensembl/hg19.ncbiRefSeq.gtf -p 16 -A "$file".tab ; done

#merge gtf with stringtie merge gtf_list.txt is just a list of all the gtfs that need to be merged
stringtie --merge gtf_list.txt -o merged.gtf


#redo stringtie with merged gtf reference
for file in `ls | grep .sorted$`; do stringtie -e -G merged.gtf -p 8 $file -o "$file".gtf.2 ; done
for file in `ls | grep .sorted$`; do stringtie -A -e -G merged.gtf -p 16 $file -o "$file".tab ; done





# ----------- MAKE RNASEQ BIGWIGS
# in vitro first
for file in in_vitro_143b_3_hg19 in_vitro_143b_4_hg19 in_vitro_143b_5_hg19; do
samtools index "$file"_sorted.bam
#bamCoverage -b "$file".bam -o "$file".bw
done

# have to activate python 3 (conda activate py37) first
for file in in_vitro_143b_3_hg19 in_vitro_143b_4_hg19 in_vitro_143b_5_hg19; do
bamCoverage -b "$file"_sorted.bam -o "$file".bw --ignoreDuplicates -p max/2 --normalizeUsing RPKM
done

# now in vivo - have to change directories to filtered directory
for file in `ls | grep bam$`; do
file2=${file%.bam}
bamCoverage -b "$file2".bam -o "$file2".bw --ignoreDuplicates -p max/2 --normalizeUsing RPKM
done

# the day 1 files kind of look like trash. going to try combining day 1 as I did with
# the ATAC seq

samtools merge d1_combined_143b_hg19_sorted_Filtered.bam d1_1_143b_hg19_sorted_Filtered.bam d1_2_143b_hg19_sorted_Filtered.bam d1_3_143b_hg19_sorted_Filtered.bam
samtools index d1_combined_143b_hg19_sorted_Filtered.bam
bamCoverage -b d1_combined_143b_hg19_sorted_Filtered.bam -o d1_combined_143b_hg19_sorted_Filtered.bw --ignoreDuplicates -p max/2 --normalizeUsing RPKM

