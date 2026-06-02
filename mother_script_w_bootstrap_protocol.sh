#!/bin/env bash

#_________________________________________________________________________
# 0) SETTING DIRECTORIES
#_________________________________________________________________________
samples=$(cat /home/worker/CARTELLA_DI_LAVORO/scripts\&nano/listasigle_positive_control)
wd="/home/worker/CARTELLA_DI_LAVORO/macao/o_GATK"
ref="/home/worker/CARTELLA_DI_LAVORO/macao/output_rmasker_on_genomic/output_Rmasker_july"
bam="${wd}/bams_native"

mkdir -p \
    "${wd}/dictionary" \
    "${wd}/true_vcfs" \
    "${wd}/results" \
    "${wd}/bootstrap_dir" \
    "${wd}/filter_files" \
    "${bam}/recal_tables" \
    "${bam}/recal_bams"

dict="${wd}/dictionary"
true="${wd}/true_vcfs"
results="${wd}/results"
boots="${wd}/bootstrap_dir"
filt="${wd}/filter_files"

mkdir -p "${filt}/last_filt" "${filt}/gvcf_filt"
zfilt="${filt}/last_filt"

#_________________________________________________________________________
# 1) CREATE REF DICTIONARY
#_________________________________________________________________________
for i in A F T ; do
    samtools faidx "${ref}/${i}1G1_merged.masked.fasta" \
        --fai-idx "${ref}/${i}1G1_merged.masked.fasta.fai"
    gatk CreateSequenceDictionary \
        -R "${ref}/${i}1G1_merged.masked.fasta" \
        -O "${ref}/${i}1G1_merged.masked.dict"
done

#_________________________________________________________________________
# 2) HaplotypeCaller (gVCF mode) + GenotypeGVCFs + SelectVariants
#_________________________________________________________________________
for i in $samples ; do
    gatk --java-options "-Xmx26g" HaplotypeCaller \
        -I "${bam}/${i}.sort.mark.nat.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${filt}/gvcf_filt/${i}.raw.g.vcf.gz" \
        --ERC GVCF \
        --sample-ploidy 1
done

for i in $samples ; do
    gatk --java-options "-Xmx26g" GenotypeGVCFs \
        -V "${filt}/gvcf_filt/${i}.raw.g.vcf.gz" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${filt}/gvcf_filt/${i}.raw.vcf.gz" \
        --sample-ploidy 1
done

for i in $samples ; do
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${filt}/gvcf_filt/${i}.raw.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.raw.snp.vcf.gz" \
        --select-type-to-include SNP
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${filt}/gvcf_filt/${i}.raw.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.raw.indel.vcf.gz" \
        --select-type-to-include INDEL
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${filt}/gvcf_filt/${i}.raw.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.raw.other.vcf.gz" \
        --select-type-to-exclude SNP \
        --select-type-to-exclude INDEL
done

# line check (all = snp + indel + other)
for i in $samples ; do
    echo "${i}" >> "${filt}/gvcf_filt/vcf_line_count.txt"
    zcat "${filt}/gvcf_filt/${i}.raw.vcf.gz"   | grep -v "^#" | wc -l >> "${filt}/gvcf_filt/vcf_line_count.txt"
    zcat "${filt}/gvcf_filt/${i}.raw.snp.vcf.gz"   | grep -v "^#" | wc -l >> "${filt}/gvcf_filt/vcf_line_count.txt"
    zcat "${filt}/gvcf_filt/${i}.raw.indel.vcf.gz" | grep -v "^#" | wc -l >> "${filt}/gvcf_filt/vcf_line_count.txt"
    zcat "${filt}/gvcf_filt/${i}.raw.other.vcf.gz" | grep -v "^#" | wc -l >> "${filt}/gvcf_filt/vcf_line_count.txt"
done

#_________________________________________________________________________
# 6) HARD-FILTERING → "fil.vcf" (see script_filter_on_gatk_standards.sh)
#_________________________________________________________________________

# building true sets: exclude filtered variants, merge SNP+INDEL
for i in $samples ; do
    gatk SelectVariants \
        -V "${filt}/gvcf_filt/${i}.fil.snp.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.fil.exc.snp.vcf.gz" \
        --exclude-filtered &&
    gatk SelectVariants \
        -V "${filt}/gvcf_filt/${i}.fil.indel.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.fil.exc.indel.vcf.gz" \
        --exclude-filtered &&
    gatk MergeVcfs \
        -I "${filt}/gvcf_filt/${i}.fil.exc.snp.vcf.gz" \
        -I "${filt}/gvcf_filt/${i}.fil.exc.indel.vcf.gz" \
        -O "${filt}/gvcf_filt/${i}.fil.vcf.gz"
done

# sort, normalize, index
for i in $samples ; do
    bcftools sort "${filt}/gvcf_filt/${i}.fil.vcf.gz" \
        -O z \
        -o "${filt}/gvcf_filt/${i}.fil.sort.vcf.gz" &&
    tabix -p vcf "${filt}/gvcf_filt/${i}.fil.sort.vcf.gz" &&
    bcftools norm "${filt}/gvcf_filt/${i}.fil.sort.vcf.gz" \
        -m+both \
        -O z \
        -f "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -o "${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz" &&
    tabix -p vcf "${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz"
done

#_________________________________________________________________________
# BOOTSTRAPPING — CYCLES 0 → N
# convergence check: for i in 2 3 ; do
#   echo ${i}_other_cycle
#   more bootstrap_dir/cycle_${i}.summary.txt | grep -i -A 4 "nEval"
# done
#_________________________________________________________________________

#_________________________________________________________________________
# CYCLE 0 → BaseRecalibrator + ApplyBQSR + HaplotypeCaller + GenotypeGVCFs
#_________________________________________________________________________
for i in $samples ; do
    gatk --java-options "-Xmx26g" BaseRecalibrator \
        -I "${bam}/${i}.sort.mark.nat.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --known-sites "${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz" \
        -O "${bam}/recal_tables/${i}.recal_1.table" &&
    gatk --java-options "-Xmx26g" ApplyBQSR \
        -I "${bam}/${i}.sort.mark.nat.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --bqsr-recal-file "${bam}/recal_tables/${i}.recal_1.table" \
        -O "${bam}/recal_bams/${i}.recal_1.bam" &&
    gatk --java-options "-Xmx26g" HaplotypeCaller \
        -I "${bam}/recal_bams/${i}.recal_1.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${filt}/gvcf_filt/${i}.cycle_1.g.vcf.gz" \
        --ERC GVCF &&
    gatk --java-options "-Xmx26g" GenotypeGVCFs \
        -V "${filt}/gvcf_filt/${i}.cycle_1.g.vcf.gz" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${boots}/${i}.cycle_1.vcf.gz"
done

# sort, normalize, index
for i in $samples ; do
    bcftools sort "${boots}/${i}.cycle_1.vcf.gz" \
        -O z \
        -o "${boots}/${i}.cycle_1.sort.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_1.sort.vcf.gz" &&
    bcftools norm "${boots}/${i}.cycle_1.sort.vcf.gz" \
        -m+both \
        -O z \
        -f "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -o "${boots}/${i}.cycle_1.sort.norm.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_1.sort.norm.vcf.gz"
done

#_________________________________________________________________________
# CYCLE 1 — VariantEval: cycle_1 vs fil (true set)
#_________________________________________________________________________
for i in $samples ; do
    cp "${filt}/gvcf_filt/${i}.fil.sort.norm.vcf.gz"* "${boots}/" &&
    gatk VariantEval \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${boots}/${i}.cycle1_eval.grp" \
        --eval "${boots}/${i}.cycle_1.sort.norm.vcf.gz" \
        --comp "${boots}/${i}.fil.sort.norm.vcf.gz" \
        --ploidy 1 \
        --EV CompOverlap
done

> "${boots}/cycle_1.summary.txt"
for i in $samples ; do
    echo "${i}_concordant_rate" >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "CompOverlap"{print$NF}'          "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_comp_rate"       >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "CompOverlap"{print$9}'           "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_n_eval_var"      >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "CompOverlap"{print$6}'           "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_n_var_at_comp"   >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "CompOverlap"{print$8}'           "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_ti/tv_ratio"     >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "TiTvVariantEvaluator"{print$8}'  "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_sensitivity"     >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "ValidationReport"{print$11}'     "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_FDR"             >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "ValidationReport"{print$14}'     "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_TP"              >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "ValidationReport"{print$7}'      "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_FP"              >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "ValidationReport"{print$8}'      "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
    echo "${i}_FN"              >> "${boots}/cycle_1.summary.txt"
    awk '$1 == "ValidationReport"{print$9}'      "${boots}/${i}.cycle1_eval.grp" >> "${boots}/cycle_1.summary.txt"
done

#_________________________________________________________________________
# CYCLE 2 — BQSR + recall
#_________________________________________________________________________
for i in $samples ; do
    gatk --java-options "-Xmx26g" BaseRecalibrator \
        -I "${bam}/recal_bams/${i}.recal_1.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --known-sites "${boots}/${i}.cycle_1.sort.norm.vcf.gz" \
        -O "${bam}/recal_tables/${i}.recal_2.table" &&
    gatk --java-options "-Xmx26g" ApplyBQSR \
        -I "${bam}/recal_bams/${i}.recal_1.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --bqsr-recal-file "${bam}/recal_tables/${i}.recal_2.table" \
        -O "${bam}/recal_bams/${i}.recal_2.bam" &&
    gatk --java-options "-Xmx26g" HaplotypeCaller \
        -I "${bam}/recal_bams/${i}.recal_2.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${filt}/gvcf_filt/${i}.cycle_2.g.vcf.gz" \
        --ERC GVCF &&
    gatk --java-options "-Xmx26g" GenotypeGVCFs \
        -V "${filt}/gvcf_filt/${i}.cycle_2.g.vcf.gz" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${boots}/${i}.cycle_2.vcf.gz"
done

# sort, normalize, index
for i in $samples ; do
    bcftools sort "${boots}/${i}.cycle_2.vcf.gz" \
        -O z \
        -o "${boots}/${i}.cycle_2.sort.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_2.sort.vcf.gz" &&
    bcftools norm "${boots}/${i}.cycle_2.sort.vcf.gz" \
        -m+both \
        -O z \
        -f "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -o "${boots}/${i}.cycle_2.sort.norm.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_2.sort.norm.vcf.gz"
done

#_________________________________________________________________________
# CYCLE 2 — VariantEval: cycle_2 vs cycle_1
#_________________________________________________________________________
for i in $samples ; do
    gatk VariantEval \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${boots}/${i}.cycle2_eval.grp" \
        --eval "${boots}/${i}.cycle_2.sort.norm.vcf.gz" \
        --comp "${boots}/${i}.cycle_1.sort.norm.vcf.gz" \
        --ploidy 1 \
        --EV CompOverlap
done

> "${boots}/cycle_2.summary.txt"
for i in $samples ; do
    echo "${i}_concordant_rate" >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "CompOverlap"{print$NF}'          "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_comp_rate"       >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "CompOverlap"{print$9}'           "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_n_eval_var"      >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "CompOverlap"{print$6}'           "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_n_var_at_comp"   >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "CompOverlap"{print$8}'           "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_ti/tv_ratio"     >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "TiTvVariantEvaluator"{print$8}'  "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_sensitivity"     >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "ValidationReport"{print$11}'     "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_FDR"             >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "ValidationReport"{print$14}'     "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_TP"              >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "ValidationReport"{print$7}'      "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_FP"              >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "ValidationReport"{print$8}'      "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
    echo "${i}_FN"              >> "${boots}/cycle_2.summary.txt"
    awk '$1 == "ValidationReport"{print$9}'      "${boots}/${i}.cycle2_eval.grp" >> "${boots}/cycle_2.summary.txt"
done

#_________________________________________________________________________
# CYCLE 3 — BQSR + recall + sort/norm/index
#_________________________________________________________________________
for i in $samples ; do
    gatk --java-options "-Xmx26g" BaseRecalibrator \
        -I "${bam}/recal_bams/${i}.recal_2.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --known-sites "${boots}/${i}.cycle_2.sort.norm.vcf.gz" \
        -O "${bam}/recal_tables/${i}.recal_3.table" &&
    gatk --java-options "-Xmx26g" ApplyBQSR \
        -I "${bam}/recal_bams/${i}.recal_2.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --bqsr-recal-file "${bam}/recal_tables/${i}.recal_3.table" \
        -O "${bam}/recal_bams/${i}.recal_3.bam" &&
    gatk --java-options "-Xmx26g" HaplotypeCaller \
        -I "${bam}/recal_bams/${i}.recal_3.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${filt}/gvcf_filt/${i}.cycle_3.g.vcf.gz" \
        --ERC GVCF &&
    gatk --java-options "-Xmx26g" GenotypeGVCFs \
        -V "${filt}/gvcf_filt/${i}.cycle_3.g.vcf.gz" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${boots}/${i}.cycle_3.vcf.gz" &&
    bcftools sort "${boots}/${i}.cycle_3.vcf.gz" \
        -O z \
        -o "${boots}/${i}.cycle_3.sort.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_3.sort.vcf.gz" &&
    bcftools norm "${boots}/${i}.cycle_3.sort.vcf.gz" \
        -m+both \
        -O z \
        -f "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -o "${boots}/${i}.cycle_3.sort.norm.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_3.sort.norm.vcf.gz"
done

#_________________________________________________________________________
# CYCLE 3 — VariantEval: cycle_3 vs cycle_2
#_________________________________________________________________________
for i in $samples ; do
    gatk VariantEval \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${boots}/${i}.cycle3_eval.grp" \
        --eval "${boots}/${i}.cycle_3.sort.norm.vcf.gz" \
        --comp "${boots}/${i}.cycle_2.sort.norm.vcf.gz" \
        --ploidy 1 \
        --EV CompOverlap
done

> "${boots}/cycle_3.summary.txt"
for i in $samples ; do
    echo "${i}_concordant_rate" >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "CompOverlap"{print$NF}'          "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_comp_rate"       >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "CompOverlap"{print$9}'           "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_n_eval_var"      >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "CompOverlap"{print$6}'           "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_n_var_at_comp"   >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "CompOverlap"{print$8}'           "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_ti/tv_ratio"     >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "TiTvVariantEvaluator"{print$8}'  "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_sensitivity"     >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "ValidationReport"{print$11}'     "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_FDR"             >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "ValidationReport"{print$14}'     "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_TP"              >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "ValidationReport"{print$7}'      "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_FP"              >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "ValidationReport"{print$8}'      "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
    echo "${i}_FN"              >> "${boots}/cycle_3.summary.txt"
    awk '$1 == "ValidationReport"{print$9}'      "${boots}/${i}.cycle3_eval.grp" >> "${boots}/cycle_3.summary.txt"
done

#_________________________________________________________________________
# CYCLE 4 — BQSR + recall
#_________________________________________________________________________
for i in $samples ; do
    gatk --java-options "-Xmx26g" BaseRecalibrator \
        -I "${bam}/recal_bams/${i}.recal_3.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --known-sites "${boots}/${i}.cycle_3.sort.norm.vcf.gz" \
        -O "${bam}/recal_tables/${i}.recal_4.table" &&
    gatk --java-options "-Xmx26g" ApplyBQSR \
        -I "${bam}/recal_bams/${i}.recal_3.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --bqsr-recal-file "${bam}/recal_tables/${i}.recal_4.table" \
        -O "${bam}/recal_bams/${i}.recal_4.bam" &&
    gatk --java-options "-Xmx26g" HaplotypeCaller \
        -I "${bam}/recal_bams/${i}.recal_4.bam" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${filt}/gvcf_filt/${i}.cycle_4.g.vcf.gz" \
        --ERC GVCF &&
    gatk --java-options "-Xmx26g" GenotypeGVCFs \
        -V "${filt}/gvcf_filt/${i}.cycle_4.g.vcf.gz" \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        --ploidy 1 \
        -O "${boots}/${i}.cycle_4.vcf.gz"
done

# sort, normalize, index
for i in $samples ; do
    bcftools sort "${boots}/${i}.cycle_4.vcf.gz" \
        -O z \
        -o "${boots}/${i}.cycle_4.sort.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_4.sort.vcf.gz" &&
    bcftools norm "${boots}/${i}.cycle_4.sort.vcf.gz" \
        -m+both \
        -O z \
        -f "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -o "${boots}/${i}.cycle_4.sort.norm.vcf.gz" &&
    tabix -p vcf "${boots}/${i}.cycle_4.sort.norm.vcf.gz"
done

#_________________________________________________________________________
# CYCLE 4 — VariantEval: cycle_4 vs cycle_3
#_________________________________________________________________________
for i in $samples ; do
    gatk VariantEval \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -O "${boots}/${i}.cycle4_eval.grp" \
        --eval "${boots}/${i}.cycle_4.sort.norm.vcf.gz" \
        --comp "${boots}/${i}.cycle_3.sort.norm.vcf.gz" \
        --ploidy 1 \
        --EV CompOverlap
done

# convergence reached
> "${boots}/cycle_4.summary.txt"
for i in $samples ; do
    echo "${i}_concordant_rate" >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "CompOverlap"{print$NF}'          "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_comp_rate"       >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "CompOverlap"{print$9}'           "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_n_eval_var"      >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "CompOverlap"{print$6}'           "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_n_var_at_comp"   >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "CompOverlap"{print$8}'           "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_ti/tv_ratio"     >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "TiTvVariantEvaluator"{print$8}'  "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_sensitivity"     >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "ValidationReport"{print$11}'     "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_FDR"             >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "ValidationReport"{print$14}'     "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_TP"              >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "ValidationReport"{print$7}'      "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_FP"              >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "ValidationReport"{print$8}'      "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
    echo "${i}_FN"              >> "${boots}/cycle_4.summary.txt"
    awk '$1 == "ValidationReport"{print$9}'      "${boots}/${i}.cycle4_eval.grp" >> "${boots}/cycle_4.summary.txt"
done

#_________________________________________________________________________
# VARIANTS SELECTION — SelectVariants on cycle_4
#_________________________________________________________________________
for i in $samples ; do
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${boots}/${i}.cycle_4.sort.norm.vcf.gz" \
        -O "${zfilt}/${i}.cy4.so.no.snp.vcf.gz" \
        --select-type-to-include SNP
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${boots}/${i}.cycle_4.sort.norm.vcf.gz" \
        -O "${zfilt}/${i}.cy4.so.no.indel.vcf.gz" \
        --select-type-to-include INDEL
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${boots}/${i}.cycle_4.sort.norm.vcf.gz" \
        -O "${zfilt}/${i}.cy4.so.no.other.vcf.gz" \
        --select-type-to-exclude SNP \
        --select-type-to-exclude INDEL
done
# dp thresholds extracted via script_dp_cycle4.R before hard-filtering

#_________________________________________________________________________
# HARD-FILTERING — exclude filtered variants, merge SNP+INDEL → ready.vcf
#_________________________________________________________________________
for i in $samples ; do
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${zfilt}/${i}.fil.cy4.so.no.snp.vcf.gz" \
        -O "${zfilt}/${i}.exc.fil.cy4.so.no.snp.vcf.gz" \
        --exclude-filtered &&
    gatk SelectVariants \
        -R "${ref}/${i:0:1}1G1_merged.masked.fasta" \
        -V "${zfilt}/${i}.fil.cy4.so.no.indel.vcf.gz" \
        -O "${zfilt}/${i}.exc.fil.cy4.so.no.indel.vcf.gz" \
        --exclude-filtered &&
    gatk MergeVcfs \
        -I "${zfilt}/${i}.exc.fil.cy4.so.no.snp.vcf.gz" \
        -I "${zfilt}/${i}.exc.fil.cy4.so.no.indel.vcf.gz" \
        -O "${true}/${i}.ready.vcf.gz"
done

#_________________________________________________________________________
# MULTI-SAMPLE MERGE — one merged VCF per species (A, F, T)
#_________________________________________________________________________
for i in A F T ; do
    bcftools merge \
        "${true}/${i}1G2.ready.vcf.gz" \
        "${true}/${i}15G.ready.vcf.gz" \
        "${true}/${i}10G.ready.vcf.gz" \
        -Ov \
        -o "${true}/${i}.ready.vcf"
done

# END
