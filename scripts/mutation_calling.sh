mypwd=`pwd`
ref_genome=./pre_load_data/hg19_all.fa
known_vcf1=./pre_load_data/dbsnp_137.hg19.vcf
known_vcf2=./pre_load_data/Mills_and_1000G_gold_standard.indels.hg19.vcf
known_vcf3=./pre_load_data/1000G_phase1.indels.hg19.vcf

:<<B1
function merge_bam {
    cd $1
    #bamlist=$(cat ./*_bam_list.txt)
    #samtools merge ${1}.bam $bamlist
    #samtools sort ${1}.bam -o ${1}.sorted.bam
    #samtools index ${1}.sorted.bam
    #gatk ApplyBQSR -R $ref_genome -I ${1}.sorted.bam --bqsr-recal-file recal_data.table -O ${1}.sorted.recal.bam
    cd ..
}
export -f merge_bam
parallel -j 60 merge_bam ::: normal normalHR progenitor tumor

for i in normalHR normal progenitor tumor
do 
cd $i
gatk --java-options "-XX:ParallelGCThreads=50" AddOrReplaceReadGroups -I $i\.sorted.bam -O $i\.sorted.RG.bam -LB library -PL illumina -PU unit -SM $i
gatk --java-options "-XX:ParallelGCThreads=50" BaseRecalibrator -I $i\.sorted.RG.bam -R $ref_genome --known-sites $known_vcf1 --known-sites $known_vcf2 --known-sites $known_vcf3 -O recal_data.table
gatk --java-options "-XX:ParallelGCThreads=50" ApplyBQSR -R $ref_genome -I $i\.sorted.RG.bam --bqsr-recal-file recal_data.table -O $i\.sorted.recal.bam
cd ..
done


mkdir -p mutect
gatk --java-options "-Xmx120g -XX:ParallelGCThreads=20" Mutect2 --native-pair-hmm-threads 20 -R $ref_genome -I ./BCIS74T_merged_bam/Ancestral_c13.bam -I ./BCIS74T_merged_bam/Normal_LumHR.bam -normal Normal_LumHR -O ./mutect/progenitor_vs_normalHR_somatic.vcf.gz
gatk --java-options "-Xmx120g -XX:ParallelGCThreads=30" Mutect2 --native-pair-hmm-threads 20 -R $ref_genome -I ./BCIS74T_merged_bam/Tumor_c11.bam -I ./BCIS74T_merged_bam/Ancestral_c13.bam -normal Ancestral_c13 -O ./mutect/tumor_vs_progenitor_somatic.vcf.gz
gatk --java-options "-Xmx120g -XX:ParallelGCThreads=30" Mutect2 --native-pair-hmm-threads 20 -R $ref_genome -I ./BCIS74T_merged_bam/Tumor_c11.bam -I ./BCIS74T_merged_bam/Normal_LumHR.bam -normal Normal_LumHR -O ./mutect/tumor_vs_normal_somatic.vcf.gz
#---Two tumors one normal
gatk --java-options "-Xmx200g -XX:ParallelGCThreads=40" Mutect2 --native-pair-hmm-threads 30 -R $ref_genome -I ./BCIS74T_merged_bam/Tumor_c11.bam -I ./BCIS74T_merged_bam/Ancestral_c13.bam -I ./BCIS74T_merged_bam/Normal_LumHR.bam -normal Normal_LumHR -O ./mutect/tumor_progenitor_vs_normal_somatic.vcf.gz

for i in progenitor_vs_normal progenitor_vs_normalHR tumor_vs_progenitor tumor_vs_normal tumor_progenitor_vs_normal
do
#---filter vcf files
gatk --java-options "-Xmx120g -XX:ParallelGCThreads=30" FilterMutectCalls -V ./mutect/$i\_somatic.vcf.gz -R $ref_genome -O ./mutect/$i\_somatic_filtered.vcf.gz
#---seperate snp and indel
gatk --java-options "-Xmx120g -XX:ParallelGCThreads=30" SelectVariants -R $ref_genome  -V ./mutect/$i\_somatic_filtered.vcf.gz -O ./mutect/$i\_somatic_filtered_SNP.vcf.gz --exclude-filtered TRUE --select-type-to-include SNP
done
#---further filter SNPs to make sure there are at least three reads support REF and ALT alleles.
#---filter vcf files (for 1 tumor and 1 normal)
for i in progenitor_vs_normal progenitor_vs_normalHR tumor_vs_progenitor tumor_vs_normal
do
zcat ./mutect/$i\_somatic_filtered_SNP.vcf.gz | grep '^chr' | awk 'BEGIN{FS=OFS="\t"}{split($9,a,":"); split($10,b,":"); split($11,c,":"); for(i=1;i<=length(a);i++){if(a[i]=="FAD"){j=i;break}} print $0"\t"b[j], c[j]}' | awk 'BEGIN{FS=OFS="\t"} {split($NF,a,","); split($(NF-1),b,","); if (b[1]>2 && a[2]>2) print $0}' > ./mutect/$i\_somatic_filtered_SNP_ALT3reads.vcf
done
#---filter vcf files (for 2 tumor and 1 normal)
zcat ./mutect/tumor_progenitor_vs_normal_somatic_filtered_SNP.vcf.gz | grep '^chr' | awk 'BEGIN{FS=OFS="\t"}{split($9,a,":"); split($10,b,":"); split($11,c,":"); split($12,d,":"); for(i=1;i<=length(a);i++){if(a[i]=="FAD"){j=i;break}} print $0"\t"b[j], c[j], d[j]}' | awk 'BEGIN{FS=OFS="\t"} {split($NF,a,","); split($(NF-1),b,","); split($(NF-2),c,","); if (c[1] > 2 && b[2]>2 && a[2]>2) print $0}' > ./mutect/tumor_progenitor_vs_normal_somatic_filtered_SNP_ALT3reads.vcf
#cut -f1,2,4,5,13,14,15 ./mutect/tumor_progenitor_vs_normal_somatic_filtered_SNP_ALT3reads.vcf | awk '$1!~/_/' | sed 's/,/\t/g' > tumor_progenitor_vs_normal_somatic_filtered_SNP_ALT3reads_clean.txt
