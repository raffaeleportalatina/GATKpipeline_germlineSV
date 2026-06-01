#!/bin/bash
#_________________________________________________________________________________________________________________________________________________________________________________________________________
#0)setting directories
#_____________________________________________________
samples=$(cat /home/worker/CARTELLA_DI_LAVORO/scripts\&nano/listasigle_positive_control)
#
wd="/home/worker/CARTELLA_DI_LAVORO/macao/o_GATK"
#wd="/media/worker/9d1935d4-054a-4514-ac05-5faba96d0e6a/macao/macao/o_GATK" # HardDisk
#gk="/media/worker/9d1935d4-054a-4514-ac05-5faba96d0e6a/macao/GATK/GATK_variant_calling"

ref="/home/worker/CARTELLA_DI_LAVORO/macao/output_rmasker_on_genomic/output_Rmasker_july"

#bam="/media/worker/OS/BU_linux/BAMs_uflag/bams_native"
bam="${wd}/bams_native"
mkdir -p ${wd}/dictionary ${wd}/true_vcfs ${wd}/results ${wd}/bootstrap_dir ${wd}/filter_files ${bam}/recal_tables ${bam}/recal_bams
#dict="${wd}/dictionary"
#reads="/home/worker/CARTELLA_DI_LAVORO/GATK/GATK_variant_calling/reads"
true="${wd}/true_vcfs"
results="${wd}/results"
boots="${wd}/bootstrap_dir"
filt="${wd}/filter_files"
mkdir -p ${filt}/last_filt ${filt}/gvcf_filt
zfilt="${filt}/last_filt"

#/media/worker/9d1935d4-054a-4514-ac05-5faba96d0e6a/macao/macao/output_rmasker_on_genomic/output_Rmasker_july
#_________________________________________________________________________________________________________________________________________________________________________________________________________
#1)create ref dictionary
#_____________________________________________________

for i in A F T ; do
samtools faidx ${ref}/${i}1G1_merged.masked.fasta \
 --fai-idx ${ref}/${i}1G1_merged.masked.fasta.fai
gatk CreateSequenceDictionary \
 -R ${ref}/${i}1G1_merged.masked.fasta \
 -O ${ref}/${i}1G1_merged.masked.dict
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#2)HaplotypeCaller(gVCF mode) and SelectVariants
#_______________________________________________________

for i in $samples ; do
#(GVCF MODE)
gatk --java-options "-Xmx26g" HaplotypeCaller \
 -I ${bam}/${i}.sort.mark.nat.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${filt}/gvcf_filt/${i}.raw.g.vcf.gz \
 --ERC GVCF \
 --sample-ploidy 1
done

#preparing to filter

for i in $samples ; do
gatk --java-options "-Xmx26g" GenotypeGVCFs \
 -V ${filt}/gvcf_filt/${i}.raw.g.vcf.gz \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${filt}/gvcf_filt/${i}.raw.vcf.gz \
 --sample-ploidy 1
done

for i in $samples ; do
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${filt}/gvcf_filt/${i}.raw.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.raw.snp.vcf.gz \
 --select-type-to-include SNP
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${filt}/gvcf_filt/${i}.raw.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.raw.indel.vcf.gz \
 --select-type-to-include INDEL
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${filt}/gvcf_filt/${i}.raw.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.raw.other.vcf.gz \
 --select-type-to-exclude SNP \
 --select-type-to-exclude INDEL
done

#line check (all=snp+indel+other)

for i in $samples ; do
echo ${i} >> ${filt}/gvcf_filt/vcf_line_count.txt
zcat ${filt}/gvcf_filt/${i}.raw.vcf.gz | grep -v "^#" | wc -l >> ${filt}/gvcf_filt/vcf_line_count.txt
zcat ${filt}/gvcf_filt/${i}.raw.snp.vcf.gz | grep -v "^#" | wc -l >> ${filt}/gvcf_filt/vcf_line_count.txt
zcat ${filt}/gvcf_filt/${i}.raw.indel.vcf.gz | grep -v "^#" | wc -l >> ${filt}/gvcf_filt/vcf_line_count.txt
zcat ${filt}/gvcf_filt/${i}.raw.other.vcf.gz | grep -v "^#" | wc -l >> ${filt}/gvcf_filt/vcf_line_count.txt
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#6)Hard-Filtering raw.vcfs to create "fil.vcf"
#______________________________________________

#SEE ${filt}/gvcf_filt/script_filter_on_gatk_standards.sh

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#building true sets
#________________________________________________


for i in $samples ; do

#select snps that pass the filters
gatk SelectVariants \
 -V ${filt}/gvcf_filt/${i}.fil.snp.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.fil.exc.snp.vcf.gz \
 --exclude-filtered &&

#select indels that pass the filters
gatk SelectVariants \
 -V ${filt}/gvcf_filt/${i}.fil.indel.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.fil.exc.indel.vcf.gz \
 --exclude-filtered &&

#merging ,sorting,indexing and normalizing filtered vcf files
gatk MergeVcfs \
 -I ${filt}/gvcf_filt/${i}.fil.exc.snp.vcf.gz \
 -I ${filt}/gvcf_filt/${i}.fil.exc.indel.vcf.gz \
 -O ${filt}/gvcf_filt/${i}.fil.vcf.gz

done

#sorting vcfs and merging multiallelics sites (plus tabix.idxs)

for i in $samples; do

bcftools sort ${filt}/gvcf_filt/${i}.fil.vcf.gz \
 -O z \
 -o ${filt}/gvcf_filt/${i}.fil.sort.vcf.gz &&
tabix -p vcf ${filt}/gvcf_filt/${i}.fil.sort.vcf.gz &&
bcftools norm ${filt}/gvcf_filt/${i}.fil.sort.vcf.gz \
 -m+both \
 -O z \
 -f ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -o ${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz &&
tabix -p vcf ${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz

done


#_________________________________________________________________________________________________________________________________________________________________________________________________________
#BOOTSTRAPPING (CYCLES 1to4...)
#_________________________________________________________

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE1) modeling and appling the model + recall
#_________________________________________________
#mkdir -p ${bam}/recal_tables ${bam}/recal_bams

for i in $samples ; do
gatk --java-options "-Xmx26g" BaseRecalibrator \
 -I ${bam}/${i}.sort.mark.nat.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --known-sites ${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz \
 -O ${bam}/recal_tables/${i}.recal_1.table &&

gatk --java-options "-Xmx26g" ApplyBQSR \
 -I ${bam}/${i}.sort.mark.nat.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --bqsr-recal-file ${bam}/recal_tables/${i}.recal_1.table \
 -O ${bam}/recal_bams/${i}.recal_1.bam &&

gatk --java-options "-Xmx26g" HaplotypeCaller \
 -I ${bam}/recal_bams/${i}.recal_1.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${filt}/gvcf_filt/${i}.cycle_1.g.vcf.gz \
 --ERC GVCF &&
gatk --java-options "-Xmx26g" GenotypeGVCFs \
 -V ${filt}/gvcf_filt/${i}.cycle_1.g.vcf.gz \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${boots}/${i}.cycle_1.vcf.gz
done

for i in $samples; do
#sort, norm, idxs
bcftools sort ${boots}/${i}.cycle_1.vcf.gz \
 -O z \
 -o ${boots}/${i}.cycle_1.sort.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_1.sort.vcf.gz &&
bcftools norm ${boots}/${i}.cycle_1.sort.vcf.gz \
 -m+both \
 -O z \
 -f ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -o ${boots}/${i}.cycle_1.sort.norm.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_1.sort.norm.vcf.gz
done


#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE1) compare the vcfs with VariantsEval: --eval: ${boots}/${i}.cycle_1.sort.norm.vcf.gz |VS| --comp: ${boots}/${i}.fil_1.sort.norm.vcf.gz
#_____________________________________________________________________________________________________________________

for i in $samples ; do
cp ${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.g* ${boots} &&

#VariantsEval(evaluating convergence 1)
gatk VariantEval \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${boots}/${i}.cycle1_eval.grp \
 --eval ${boots}/${i}.cycle_1.sort.norm.vcf.gz \
 --comp ${boots}/${i}.fil.sort.norm.vcf.gz \
 --ploidy 1 \
 --EV CompOverlap
done
#eval convergence script
> ${boots}/cycle_1.summary.txt
for i in $samples ; do
#1)concordant_rate
echo ${i}_concordant_rate >> ${boots}/cycle_1.summary.txt
awk '$1 == "CompOverlap"{print$NF}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#2)comp_rate
echo ${i}_comp_rate >> ${boots}/cycle_1.summary.txt
awk '$1 == "CompOverlap"{print$9}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#3)number_of_eval_var
echo ${i}_n_eval_var >> ${boots}/cycle_1.summary.txt
awk '$1 == "CompOverlap"{print$6}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#3')number_of_var_at_comp
echo ${i}_n_var_at_comp >> ${boots}/cycle_1.summary.txt
awk '$1 == "CompOverlap"{print$8}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#4)ti/tv
echo ${i}_ti/tv_ratio >> ${boots}/cycle_1.summary.txt
awk '$1 == "TiTvVariantEvaluator"{print$8}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#5)sensitivity
echo ${i}_sensitivity >> ${boots}/cycle_1.summary.txt
awk '$1 == "ValidationReport"{print$11}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#5')FDR
echo ${i}_FDR >> ${boots}/cycle_1.summary.txt
awk '$1 == "ValidationReport"{print$14}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#6) True Positives
echo ${i}_TP >> ${boots}/cycle_1.summary.txt
awk '$1 == "ValidationReport"{print$7}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#7) False Positives
echo ${i}_FP >> ${boots}/cycle_1.summary.txt
awk '$1 == "ValidationReport"{print$8}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
#8) False Negatives
echo ${i}_FN >> ${boots}/cycle_1.summary.txt
awk '$1 == "ValidationReport"{print$9}' ${boots}/${i}.cycle1_eval.grp >> ${boots}/cycle_1.summary.txt
done


#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CICLE2) modeling and appling the model + recall

#__________________________________________

#BQSR
for i in $samples ; do
gatk --java-options "-Xmx26g" BaseRecalibrator \
 -I ${bam}/recal_bams/${i}.recal_1.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --known-sites ${boots}/${i}.cycle_1.sort.norm.vcf.gz \
 -O ${bam}/recal_tables/${i}.recal_2.table &&

gatk --java-options "-Xmx26g" ApplyBQSR \
 -I ${bam}/recal_bams/${i}.recal_1.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --bqsr-recal-file ${bam}/recal_tables/${i}.recal_2.table \
 -O ${bam}/recal_bams/${i}.recal_2.bam &&

#Re-Calling
gatk --java-options "-Xmx26g" HaplotypeCaller \
 -I ${bam}/recal_bams/${i}.recal_2.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${filt}/gvcf_filt/${i}.cycle_2.g.vcf.gz \
 --ERC GVCF &&

gatk --java-options "-Xmx26g" GenotypeGVCFs \
 -V ${filt}/gvcf_filt/${i}.cycle_2.g.vcf.gz \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${boots}/${i}.cycle_2.vcf.gz
done

for i in $samples; do
#sort, norm, idxs
bcftools sort ${boots}/${i}.cycle_2.vcf.gz \
 -O z \
 -o ${boots}/${i}.cycle_2.sort.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_2.sort.vcf.gz &&
for i in $samples; do
bcftools norm ${boots}/${i}.cycle_2.sort.vcf.gz \
 -m+both \
 -O z \
 -f ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -o ${boots}/${i}.cycle_2.sort.norm.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_2.sort.norm.vcf.gz
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE2) compare the vcfs with VariantsEval: --eval: ${boots}/${i}.cycle_2.sort.norm.vcf.gz |VS| --comp: ${boots}/${i}.cycle_1.sort.norm.vcf.gz
#________________________________________________

for i in $samples ; do
#VariantsEval
gatk VariantEval \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${boots}/${i}.cycle2_eval.grp \
 --eval ${boots}/${i}.cycle_2.sort.norm.vcf.gz \
 --comp ${boots}/${i}.cycle_1.sort.norm.vcf.gz \
 --ploidy 1 \
 --EV CompOverlap
done

#eval convergence 2
> ${boots}/cycle_2.summary.txt &&
for i in $samples ; do
#1)concordant_rate
echo ${i}_concordant_rate >> ${boots}/cycle_2.summary.txt
awk '$1 == "CompOverlap"{print$NF}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#2)comp_rate
echo ${i}_comp_rate >> ${boots}/cycle_2.summary.txt
awk '$1 == "CompOverlap"{print$9}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#3)number_of_eval_var
echo ${i}_n_eval_var >> ${boots}/cycle_2.summary.txt
awk '$1 == "CompOverlap"{print$6}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#3')number_of_var_at_comp
echo ${i}_n_var_at_comp >> ${boots}/cycle_2.summary.txt
awk '$1 == "CompOverlap"{print$8}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#4)ti/tv
echo ${i}_ti/tv_ratio >> ${boots}/cycle_2.summary.txt
awk '$1 == "TiTvVariantEvaluator"{print$8}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#5)sensitivity
echo ${i}_sensitivity >> ${boots}/cycle_2.summary.txt
awk '$1 == "ValidationReport"{print$11}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#5')FDR
echo ${i}_FDR >> ${boots}/cycle_2.summary.txt
awk '$1 == "ValidationReport"{print$14}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#6) True Positives
echo ${i}_TP >> ${boots}/cycle_2.summary.txt
awk '$1 == "ValidationReport"{print$7}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#7) False Positives
echo ${i}_FP >> ${boots}/cycle_2.summary.txt
awk '$1 == "ValidationReport"{print$8}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
#8) False Negatives
echo ${i}_FN >> ${boots}/cycle_2.summary.txt
awk '$1 == "ValidationReport"{print$9}' ${boots}/${i}.cycle2_eval.grp >> ${boots}/cycle_2.summary.txt
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE3) modeling and appling model + recall
#__________________________________________


#BQSR
for i in $samples ; do
gatk --java-options "-Xmx26g" BaseRecalibrator \
 -I ${bam}/recal_bams/${i}.recal_2.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --known-sites ${boots}/${i}.cycle_2.sort.norm.vcf.gz \
 -O ${bam}/recal_tables/${i}.recal_3.table &&

gatk --java-options "-Xmx26g" ApplyBQSR \
 -I ${bam}/recal_bams/${i}.recal_2.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --bqsr-recal-file ${bam}/recal_tables/${i}.recal_3.table \
 -O ${bam}/recal_bams/${i}.recal_3.bam &&

#Re-Calling
gatk --java-options "-Xmx26g" HaplotypeCaller \
 -I ${bam}/recal_bams/${i}.recal_3.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${filt}/gvcf_filt/${i}.cycle_3.g.vcf.gz \
 --ERC GVCF &&

gatk --java-options "-Xmx26g" GenotypeGVCFs \
 -V ${filt}/gvcf_filt/${i}.cycle_3.g.vcf.gz \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${boots}/${i}.cycle_3.vcf.gz &&

#sort, norm, idxs
bcftools sort ${boots}/${i}.cycle_3.vcf.gz \
 -O z \
 -o ${boots}/${i}.cycle_3.sort.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_3.sort.vcf.gz &&
bcftools norm ${boots}/${i}.cycle_3.sort.vcf.gz \
 -m+both \
 -O z \
 -f ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -o ${boots}/${i}.cycle_3.sort.norm.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_3.sort.norm.vcf.gz
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE3: compare the vcfs with VariantsEval: --eval: ${boots}/${i}.cycle_3.sort.norm.vcf.gz |VS| --comp: ${boots}/${i}.cycle_2.sort.norm.vcf.gz
#__________________________________________

for i in $samples ; do
#VariantsEval
gatk VariantEval \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${boots}/${i}.cycle3_eval.grp \
 --eval ${boots}/${i}.cycle_3.sort.norm.vcf.gz \
 --comp ${boots}/${i}.cycle_2.sort.norm.vcf.gz \
 --ploidy 1 \
 --EV CompOverlap
done


#eval convergence 3 (example command for evaluating convergence between a cycle and the subsequent \
#: for i in 2 3 ; do echo ${i}_other_cycle ; more bootstrap_dir/cycle_${i}.summary.txt | grep -i -A 4 "nEval" ; done)
> ${boots}/cycle_3.summary.txt &&
for i in $samples ; do
#1)concordant_rate
echo ${i}_concordant_rate >> ${boots}/cycle_3.summary.txt
awk '$1 == "CompOverlap"{print$NF}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#2)comp_rate
echo ${i}_comp_rate >> ${boots}/cycle_3.summary.txt
awk '$1 == "CompOverlap"{print$9}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#3)number_of_eval_var
echo ${i}_n_eval_var >> ${boots}/cycle_3.summary.txt
awk '$1 == "CompOverlap"{print$6}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#3')number_of_var_at_comp
echo ${i}_n_var_at_comp >> ${boots}/cycle_3.summary.txt
awk '$1 == "CompOverlap"{print$8}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#4)ti/tv
echo ${i}_ti/tv_ratio >> ${boots}/cycle_3.summary.txt
awk '$1 == "TiTvVariantEvaluator"{print$8}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#5)sensitivity
echo ${i}_sensitivity >> ${boots}/cycle_3.summary.txt
awk '$1 == "ValidationReport"{print$11}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#5')FDR
echo ${i}_FDR >> ${boots}/cycle_3.summary.txt
awk '$1 == "ValidationReport"{print$14}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#6) True Positives
echo ${i}_TP >> ${boots}/cycle_3.summary.txt
awk '$1 == "ValidationReport"{print$7}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#7) False Positives
echo ${i}_FP >> ${boots}/cycle_3.summary.txt
awk '$1 == "ValidationReport"{print$8}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
#8) False Negatives
echo ${i}_FN >> ${boots}/cycle_3.summary.txt
awk '$1 == "ValidationReport"{print$9}' ${boots}/${i}.cycle3_eval.grp >> ${boots}/cycle_3.summary.txt
done
#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE4) modeling and appling model + recall.....to complete
#__________________________________________


#BQSR
for i in $samples ; do
gatk --java-options "-Xmx26g" BaseRecalibrator \
 -I ${bam}/recal_bams/${i}.recal_3.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --known-sites ${boots}/${i}.cycle_3.sort.norm.vcf.gz \
 -O ${bam}/recal_tables/${i}.recal_4.table &&

gatk --java-options "-Xmx26g" ApplyBQSR \
 -I ${bam}/recal_bams/${i}.recal_3.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --bqsr-recal-file ${bam}/recal_tables/${i}.recal_4.table \
 -O ${bam}/recal_bams/${i}.recal_4.bam &&

#Re-Calling
gatk --java-options "-Xmx26g" HaplotypeCaller \
 -I ${bam}/recal_bams/${i}.recal_4.bam \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${filt}/gvcf_filt/${i}.cycle_4.g.vcf.gz \
 --ERC GVCF &&

gatk --java-options "-Xmx26g" GenotypeGVCFs \
 -V ${filt}/gvcf_filt/${i}.cycle_4.g.vcf.gz \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 --ploidy 1 \
 -O ${boots}/${i}.cycle_4.vcf.gz
done

for i in $samples; do
#sort, norm, idxs
bcftools sort ${boots}/${i}.cycle_4.vcf.gz \
 -O z \
 -o ${boots}/${i}.cycle_4.sort.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_4.sort.vcf.gz &&
bcftools norm ${boots}/${i}.cycle_4.sort.vcf.gz \
 -m+both \
 -O z \
 -f ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -o ${boots}/${i}.cycle_4.sort.norm.vcf.gz &&
tabix -p vcf ${boots}/${i}.cycle_4.sort.norm.vcf.gz
done

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#CYCLE4: compare the vcfs with VariantsEval: --eval: ${boots}/${i}.cycle_4.sort.norm.vcf.gz |VS| --comp: ${boots}/${i}.cycle_3.sort.norm.vcf.gz
#__________________________________________

for i in $samples ; do
#VariantsEval
gatk VariantEval \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -O ${boots}/${i}.cycle4_eval.grp \
 --eval ${boots}/${i}.cycle_4.sort.norm.vcf.gz \
 --comp ${boots}/${i}.cycle_3.sort.norm.vcf.gz \
 --ploidy 1 \
 --EV CompOverlap
done


#eval convergence 4 (example command for evaluating convergence between a cycle and the subsequent \
#: for i in 2 3 ; do echo ${i}_other_cycle ; more bootstrap_dir/cycle_${i}.summary.txt | grep -i -A 4 "nEval" ; done)

> ${boots}/cycle_4.summary.txt &&
for i in $samples ; do
#1)concordant_rate
echo ${i}_concordant_rate >> ${boots}/cycle_4.summary.txt
awk '$1 == "CompOverlap"{print$NF}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#2)comp_rate
echo ${i}_comp_rate >> ${boots}/cycle_4.summary.txt
awk '$1 == "CompOverlap"{print$9}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#3)number_of_eval_var
echo ${i}_n_eval_var >> ${boots}/cycle_4.summary.txt
awk '$1 == "CompOverlap"{print$6}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#3')number_of_var_at_comp
echo ${i}_n_var_at_comp >> ${boots}/cycle_4.summary.txt
awk '$1 == "CompOverlap"{print$8}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#4)ti/tv
echo ${i}_ti/tv_ratio >> ${boots}/cycle_4.summary.txt
awk '$1 == "TiTvVariantEvaluator"{print$8}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#5)sensitivity
echo ${i}_sensitivity >> ${boots}/cycle_4.summary.txt
awk '$1 == "ValidationReport"{print$11}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#5')FDR
echo ${i}_FDR >> ${boots}/cycle_4.summary.txt
awk '$1 == "ValidationReport"{print$14}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#6) True Positives
echo ${i}_TP >> ${boots}/cycle_4.summary.txt
awk '$1 == "ValidationReport"{print$7}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#7) False Positives
echo ${i}_FP >> ${boots}/cycle_4.summary.txt
awk '$1 == "ValidationReport"{print$8}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
#8) False Negatives
echo ${i}_FN >> ${boots}/cycle_4.summary.txt
awk '$1 == "ValidationReport"{print$9}' ${boots}/${i}.cycle4_eval.grp >> ${boots}/cycle_4.summary.txt
done

#reaching the convergence!!
#_________________________________________________________________________________________________________________________________________________________________________________________________________
#VARIANTS SELECTION
#__________________________________________


for i in $samples ; do
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${boots}/${i}.cycle_4.sort.norm.vcf.gz \
 -O ${zfilt}/${i}.cy4.so.no.snp.vcf.gz \
 --select-type-to-include SNP
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${boots}/${i}.cycle_4.sort.norm.vcf.gz \
 -O ${zfilt}/${i}.cy4.so.no.indel.vcf.gz \
 --select-type-to-include INDEL
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${boots}/${i}.cycle_4.sort.norm.vcf.gz \
 -O ${zfilt}/${i}${i}.cy4.so.no.other.vcf.gz \
 --select-type-to-exclude SNP \
 --select-type-to-exclude INDEL
done

#extractring dp thr w ~/CARTELLA_DI_LAVORO/macao/o_GATK/bams_native/recal_bams/dp_stats/script_dp_cycle4.R to filter sv 

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#HARD-FILTERING :
#__________________________________________

#excluding filtered: select the snp that have pass the filters

for i in $samples ; do
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${zfilt}/${i}.fil.cy4.so.no.snp.vcf.gz \
 -O ${zfilt}/${i}.exc.fil.cy4.so.no.snp.vcf.gz \
 --exclude-filtered &&
#select the indel that have pass the filters
gatk SelectVariants \
 -R ${ref}/${i:0:1}1G1_merged.masked.fasta \
 -V ${zfilt}/${i}.fil.cy4.so.no.indel.vcf.gz \
 -O ${zfilt}/${i}.exc.fil.cy4.so.no.indel.vcf.gz \
 --exclude-filtered &&

#_________________________________________________________________________________________________________________________________________________________________________________________________________
#Merging cycle_3.fil.exc.snp e cycle_3.fil.exc.indel
#________________________________________________

gatk MergeVcfs \
 -I ${zfilt}/${i}.exc.fil.cy4.so.no.snp.vcf.gz \
 -I ${zfilt}/${i}.exc.fil.cy4.so.no.indel.vcf.gz \
 -O ${true}/${i}.ready.vcf.gz
done


#_________________________________________________________________________________________________________________________________________________________________________________________________________
#merging for multisample vcfs
#__________________________________________
#Merge
for i in A F T ; do
bcftools merge \
 ${true}/${i:0:1}1G2.ready.vcf.gz \
 ${true}/${i:0:1}15G.ready.vcf.gz \
 ${true}/${i:0:1}10G.ready.vcf.gz \
 -Ov \
 -o ${true}/${i}.ready.vcf
done
