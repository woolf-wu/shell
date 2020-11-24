########################## STEP 1 MarkDuplicatesSpark ##########################
cat all_sample | while read id; do
sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Mapping
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
temp=${dir}/${sample}/gatk/log/tmp

mkdir -p ${temp}
cd ${dir}/${sample}/gatk

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" MarkDuplicatesSpark -I ${dir}/${sample}/${sample}.nodup.bam -M ${sample}.metrics -O ${sample}_marked.bam -R ${reference} --spark-master local[20] --remove-sequencing-duplicates true" > ${dir}/${sample}/gatk/step1_${sample}_MarkDuplicatesSpark.sh
echo "rm -rf ${temp}" >> ${dir}/${sample}/gatk/step1_${sample}_MarkDuplicatesSpark.sh

cd ${dir}/${sample}/gatk/log
qsub -cwd -l vf=40G,p=10 ${dir}/${sample}/gatk/step1_${sample}_MarkDuplicatesSpark.sh

done


########################## STEP 2 BaseRecalibrator ##########################
## BaseRecalibrator运行很快，不建议用BaseRecalibratorSpark，报错是打开文件过多，WES偶尔报错，WGS一定报错，核心设置在20和5
cat all_sample | while read id; do
sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Mapping
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
## WES
snp=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/dbsnp/dbSNP_151.b37.exon.vcf
## WGS
#snp=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/dbsnp/00-common_all.vcf
indel_file1=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/1000G_phase1.indels.b37.vcf
indel_file2=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf
temp=${dir}/${sample}/gatk/log/tmp

mkdir -p ${temp}
cd ${dir}/${sample}/gatk

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" BaseRecalibrator -R ${reference} -I ${sample}_marked.bam --known-sites ${snp} --known-sites ${indel_file1} --known-sites ${indel_file2} -O ${sample}_recal.table" > ${dir}/${sample}/gatk/step2_${sample}_BaseRecalibratorSpark.sh
echo "sleep 5m" >> ${dir}/${sample}/gatk/step2_${sample}_BaseRecalibratorSpark.sh
echo "rm -rf ${temp}" >> ${dir}/${sample}/gatk/step2_${sample}_BaseRecalibratorSpark.sh

cd ${dir}/${sample}/gatk/log
qsub -cwd -l vf=20G,p=10 ${dir}/${sample}/gatk/step2_${sample}_BaseRecalibratorSpark.sh

done



########################## STEP 3 ApplyBQSRSpark ##########################
cat all_sample | while read id; do
sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Mapping
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
temp=${dir}/${sample}/gatk/log/tmp

mkdir -p ${temp}
cd ${dir}/${sample}/gatk

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" ApplyBQSRSpark -R ${reference} -I ${sample}_marked.bam -bqsr ${sample}_recal.table -O ${dir}/${sample}/gatk/${sample}_bqsr.bam --spark-master local[10]" > ${dir}/${sample}/gatk/step3_${sample}_ApplyBQSRSpark.sh
echo "sleep 5m" >> ${dir}/${sample}/gatk/step3_${sample}_ApplyBQSRSpark.sh
echo "rm -rf ${temp}" >> ${dir}/${sample}/gatk/step3_${sample}_ApplyBQSRSpark.sh

cd ${dir}/${sample}/gatk/log
qsub -cwd -l vf=200G,p=20  ${dir}/${sample}/gatk/step3_${sample}_ApplyBQSRSpark.sh

done




########################## STEP 4 Mutect2 ##########################

########################## STEP 4.1 PON   ##########################

cat sample_N | while read id; do
sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
dir_pon=${dir}/PON

mkdir -p ${dir_pon}/Step_1/tmp
echo "${gatk} --java-options \"-Djava.io.tmpdir=${dir_pon}/Step_1/tmp -Xmx50G\" Mutect2 -R ${reference} -I ${dir}/Mapping/${sample}/gatk/${sample}_bqsr.bam --max-mnp-distance 0 -O ${dir_pon}/Step_1/${sample}.vcf.gz" > ${dir_pon}/Step_1/${sample}.tumor-only.sh
# running time 18h
qsub -cwd -l vf=25G,p=2 ${dir_pon}/Step_1/${sample}.tumor-only.sh

done


for chr in {1..22} X Y MT; do
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
dir_pon=${dir}/PON

mkdir -p ${dir_pon}/Step_2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${dir_pon}/Step_2/tmp -Xmx50G\" GenomicsDBImport -R ${reference} -L ${chr} --genomicsdb-workspace-path ${dir_pon}/Step_2/pon_db_${chr} $(ls ${dir_pon}/Step_1/*vcf.gz | xargs -i echo "-V "{} | tr "\n" "\t")" > ${dir_pon}/Step_2/Step_2_GenomicsDBImport_${chr}.sh

cd ${dir_pon}/Step_2/
qsub -cwd -l vf=20G,p=2 ${dir_pon}/Step_2/Step_2_GenomicsDBImport_${chr}.sh

done

for chr in {1..22} X Y MT; do
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
af_gnomad=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/af-only-gnomad.raw.sites.b37.vcf.gz
dir_pon=${dir}/PON

mkdir -p ${dir_pon}/Step_3/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${dir_pon}/Step_3/tmp -Xmx50G\" CreateSomaticPanelOfNormals -R ${reference} --germline-resource ${af_gnomad} -V gendb://pon_db_${chr} -O ${dir_pon}/Step_3/pon_${chr}.vcf.gz" > ${dir_pon}/Step_3/Step_3_CreateSomaticPanelOfNormals_${chr}.sh

cd ${dir_pon}/Step_2/
qsub -cwd -l vf=20G,p=2 ${dir_pon}/Step_3/Step_3_CreateSomaticPanelOfNormals_${chr}.sh

done

gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
dir_pon=${dir}/PON
echo "${gatk} --java-options \"-Djava.io.tmpdir=${dir_pon}/tmp -Xmx50G\" GatherVcfs $(for i in {1..22} X Y ;do echo "-I ${dir_pon}/Step_3/pon_${i}.vcf.gz"; done | tr "\n" "\t") -O ${dir_pon}/pon.vcf.gz" > ${dir_pon}/Step_4_GatherVcfs.sh
echo "${gatk} IndexFeatureFile -I ${dir_pon}/pon.vcf.gz" >> ${dir_pon}/Step_4_GatherVcfs.sh





########################## STEP 4.2 Mutect2 paired  ##########################
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed
af_gnomad=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/af-only-gnomad.raw.sites.b37.vcf.gz
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

mkdir -p ${temp}

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" Mutect2 -R ${reference} -I ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam -I ${dir}/Mapping/${N_sample}/gatk/${N_sample}_bqsr.bam -normal ${N_sample}  --genotype-germline-sites true --genotype-pon-sites true --germline-resource ${af_gnomad} -L ${bed} -O ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_unfiltered.vcf.gz --f1r2-tar-gz ${dir}/Somatic/${T_sample}/Mutect2/f1r2.tar.gz" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_Mutect2.sh
# --panel-of-normals ${dir}/PON/pon.vcf.gz
cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=20G,p=10 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_Mutect2.sh
done


########################## STEP 4.2 Mutect2 single  ##########################
cat PAAD_single_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed
af_gnomad=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/af-only-gnomad.raw.sites.b37.vcf.gz
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

mkdir -p ${temp}

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" Mutect2 -R ${reference} -I ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --genotype-germline-sites true --genotype-pon-sites true --germline-resource ${af_gnomad} -L ${bed} --panel-of-normals ${dir}/PON/pon.vcf.gz -O ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_unfiltered.vcf.gz --f1r2-tar-gz ${dir}/Somatic/${T_sample}/Mutect2/f1r2.tar.gz" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_Mutect2.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=60G,p=10 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_Mutect2.sh
done


########################## STEP 4.3 Mutect2 Filter  ##########################

cat all_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" LearnReadOrientationModel -I ${dir}/Somatic/${T_sample}/Mutect2/f1r2.tar.gz -O ${dir}/Somatic/${T_sample}/Mutect2/read-orientation-model.tar.gz" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S1_LearnReadOrientationModel.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=4G,p=2 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S1_LearnReadOrientationModel.sh
done


cat all_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
exac_common=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/small_exac_common_3_b37.vcf.gz
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" GetPileupSummaries -I ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam -V ${exac_common} -L ${exac_common} -O ${dir}/Somatic/${T_sample}/Mutect2/getpileupsummaries.table" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S2_GetPileupSummaries.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=4G,p=2 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S2_GetPileupSummaries.sh
done


cat all_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" CalculateContamination -I ${dir}/Somatic/${T_sample}/Mutect2/getpileupsummaries.table -tumor-segmentation ${dir}/Somatic/${T_sample}/Mutect2/segments.table -O ${dir}/Somatic/${T_sample}/Mutect2/calculatecontamination.table" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S3_CalculateContamination.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=1G,p=1 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S3_CalculateContamination.sh
done


cat all_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" FilterMutectCalls -V ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_unfiltered.vcf.gz --reference ${reference} --tumor-segmentation ${dir}/Somatic/${T_sample}/Mutect2/segments.table --contamination-table ${dir}/Somatic/${T_sample}/Mutect2/calculatecontamination.table --ob-priors ${dir}/Somatic/${T_sample}/Mutect2/read-orientation-model.tar.gz -O ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_filtered.vcf.gz" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S4_FilterMutectCalls.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=2G,p=1 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S4_FilterMutectCalls.sh
done


########################## STEP 4.4 ffpe sample               ##########################
########################## STEP 4.4.1 FilterByOrientationBias ##########################
cat all_sample | while read id; do
T_sample=${id}
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
temp=${dir}/Somatic/${T_sample}/Mutect2/tmp

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" CollectSequencingArtifactMetrics -R ${reference} -I ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --FILE_EXTENSION \".txt\" -O tumor_artifact" > ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S5_CollectSequencingArtifactMetrics.sh

echo "${gatk} --java-options \"-Djava.io.tmpdir=${temp} -Xmx50G\" FilterByOrientationBias -V ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_filtered.vcf.gz --artifact-modes 'C/T' -P ${dir}/Somatic/${T_sample}/Mutect2/tumor_artifact.pre_adapter_detail_metrics.txt -O ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}.mutect2_filtered_ffpe.vcf.gz" >> ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S5_CollectSequencingArtifactMetrics.sh

cd ${dir}/Somatic/${T_sample}/Mutect2/
qsub -cwd -l vf=8G,p=2 ${dir}/Somatic/${T_sample}/Mutect2/${T_sample}_filter_S5_CollectSequencingArtifactMetrics.sh
done

########################## STEP 4.4.2 Pisces            ##########################
安装失败



########################## STEP 4.5 ANNOVAR Mutect2 ##########################

cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Mutect2
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

#zcat ${dir}/${T_sample}.mutect2_filtered.vcf.gz | grep "#" > ${dir}/${T_sample}.mutect2.vcf
#zcat ${dir}/${T_sample}.mutect2_filtered.vcf.gz | grep PASS | awk -F "\t" '{if ($1!="hs37d5" && $1!~/NC/ && $1!~/GL/) print $0}' >> ${dir}/${T_sample}.mutect2.vcf

#python /PUBLIC/software/CANCER/Module/CancerGenome/FFPE/CT.Filter.py -i /ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Mutect2/${T_sample}.filter.list -t 'C>T' -o /ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Mutect2/

#${gatk}  SelectVariants -R ${reference} -V ${dir}/${T_sample}.mutect2.vcf --select-type-to-include SNP -O ${dir}/${T_sample}.mutect2.somatic.snv.vcf
#${gatk}  SelectVariants -R ${reference} -V ${dir}/${T_sample}.mutect2.vcf --select-type-to-include INDEL -O ${dir}/${T_sample}.mutect2.somatic.indel.vcf

echo "cd ${dir}" > ${dir}/${T_sample}.mutect2.anno.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b SNP -z ${database} -d ${ANN} ${dir}/${T_sample}.mutect2.somatic.snv.vcf ${T_sample}" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "bcftools view -s ${T_sample},${N_sample} ${T_sample}.mutect2.somatic.snv.reformated.vcf.gz -Ov | bgzip > tmp.vcf.gz" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "mv tmp.vcf.gz ${T_sample}.mutect2.somatic.snv.reformated.vcf.gz" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "bcftools index -t ${T_sample}.mutect2.somatic.snv.reformated.vcf.gz"  >> ${dir}/${T_sample}.mutect2.anno.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.mutect2.somatic.snv.reformated.vcf.gz -o ${dir}/${T_sample}.muTect.somatic.snv.maf -t ${T_sample} -s WES_agi -p Hiseq -m muTect -f avsnp150 -x CADD,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b INDEL -z ${database} -d ${ANN} ${dir}/${T_sample}.mutect2.somatic.indel.vcf ${T_sample}" > ${dir}/${T_sample}.mutect2.anno.sh
echo "cd ${dir}" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "bcftools view -s ${T_sample},${N_sample} ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz -Ov | bgzip > tmp.vcf.gz" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "mv tmp.vcf.gz ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz" >> ${dir}/${T_sample}.mutect2.anno.sh
echo "bcftools index -t ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz"  >> ${dir}/${T_sample}.mutect2.anno.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.mutect2.somatic.indel.reformated.vcf.gz -o ${dir}/${T_sample}.mutect2.somatic.indel.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m muTect -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/${T_sample}.mutect2.anno.sh

cd ${dir}
qsub -cwd -l vf=2G ${dir}/${T_sample}.mutect2.anno.sh

done


########################## STEP 5 INDEL Strelka2 ##########################

########################## STEP 5.1 manta       ##########################
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
manta=/ifs/TJPROJ3/CANCER/share/software/manta-1.6.0/bin/configManta.py
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
#bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed.gz

mkdir -p ${dir}/Somatic/${T_sample}/manta/
${manta} --normalBam ${dir}/Mapping/${N_sample}/gatk/${N_sample}_bqsr.bam --tumorBam ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam  --referenceFasta ${reference} --runDir ${dir}/Somatic/${T_sample}/manta/
# --exome --callRegions=${bed}
echo "${dir}/Somatic/${T_sample}/manta/runWorkflow.py -m local -j 10" > ${dir}/Somatic/${T_sample}/manta/manta_${T_sample}.sh
cd ${dir}/Somatic/${T_sample}/manta/
qsub -cwd -l vf=30G -q mepi1.q@tjcompute038.hpc,cancer1.q@tjcompute135.hpc,cancer1.q@tjcompute136.hpc,cancer1.q@tjcompute137.hpc ${dir}/Somatic/${T_sample}/manta/manta_${T_sample}.sh
done


########################## STEP 5.2 Strelka2       ##########################
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
StrelkaSomatic=/ifs/TJPROJ3/CANCER/share/software/strelka-2.9.2/bin/configureStrelkaSomaticWorkflow.py
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
#bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed.gz

mkdir -p ${dir}/Somatic/${T_sample}/Strelka2_manta/
${StrelkaSomatic} --normalBam ${dir}/Mapping/${N_sample}/gatk/${N_sample}_bqsr.bam --tumorBam ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --referenceFasta ${reference} --runDir ${dir}/Somatic/${T_sample}/Strelka2_manta/ --indelCandidates ${dir}/Somatic/${T_sample}/manta/results/variants/candidateSmallIndels.vcf.gz
# --exome --callRegions=${bed}
echo "${dir}/Somatic/${T_sample}/Strelka2_manta/runWorkflow.py -m local -j 10" > ${dir}/Somatic/${T_sample}/Strelka2_manta/Strelka2_${T_sample}.sh
cd ${dir}/Somatic/${T_sample}/Strelka2_manta/
qsub -cwd -l vf=20G,p=2 -q mepi1.q@tjcompute038.hpc,cancer1.q@tjcompute135.hpc,cancer1.q@tjcompute136.hpc,cancer1.q@tjcompute137.hpc ${dir}/Somatic/${T_sample}/Strelka2_manta/Strelka2_${T_sample}.sh
done


########################## STEP 5.2 Strelka2 training  ##########################

########################## STEP 5.2.1 Strelka2-germline ##########################
#cat all_sample | while read id; do
#T_sample=${id}
#dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
#StrelkaGermline=/ifs/TJPROJ3/CANCER/share/software/strelka-2.9.2/bin/configureStrelkaGermlineWorkflow.py
#reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
#bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed.gz
#
#mkdir -p ${dir}/Germline/${T_sample}/Strelka2/
#${StrelkaGermline} --bam ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --referenceFasta ${reference} --runDir ${dir}/Germline/${T_sample}/Strelka2/ --exome --callRegions=${bed}
#echo "${dir}/Germline/${T_sample}/Strelka2/runWorkflow.py -m local -j 20" > ${dir}/Germline/${T_sample}/Strelka2/Strelka2_${T_sample}.sh
#cd ${dir}/Germline/${T_sample}/Strelka2/
#qsub -cwd -l vf=6G,p=4 ${dir}/Germline/${T_sample}/Strelka2/Strelka2_${T_sample}.sh
#done


########################## STEP 5.2.2 Strelka2 anno    ##########################
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Strelka2_manta
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

#zcat ${dir}/somatic.indels.vcf.gz | awk -F "\t" '{if ($1!="hs37d5" && $1!~/NC/ && $1!~/GL/) print $0}' | awk -F "\t" '{if ($7!~/Low/) print $0}' > ${dir}/${T_sample}.Strelka2.somatic.indel.vcf
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b INDEL -z ${database} -d ${ANN} -v vcf4old ${dir}/${T_sample}.Strelka2.somatic.indel.vcf ${N_sample},${T_sample}" > ${dir}/anno_Strelka2_${T_sample}.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.Strelka2.somatic.indel.reformated.vcf.gz -o ${dir}/${T_sample}.Strelka2.somatic.indel.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m Strelka -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/anno_Strelka2_${T_sample}.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b SNP -z ${database} -d ${ANN} -v vcf4old ${dir}/${T_sample}.Strelka2.somatic.snv.vcf ${N_sample},${T_sample}" >> ${dir}/anno_Strelka2_${T_sample}.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.Strelka2.somatic.snv.reformated.vcf.gz -o ${dir}/${T_sample}.Strelka2.somatic.snv.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m Strelka -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/anno_Strelka2_${T_sample}.sh
cd ${dir}
qsub -cwd -l vf=4G ${dir}/anno_Strelka2_${T_sample}.sh
done


cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Strelka2_manta
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

#zcat ${dir}/somatic.snvs.vcf.gz | awk -F "\t" '{if ($1!="hs37d5" && $1!~/NC/ && $1!~/GL/) print $0}' | awk -F "\t" '{if ($7!~/Low/) print $0}' > ${dir}/${T_sample}.Strelka2.somatic.snv.vcf
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b SNP -z ${database} -d ${ANN} -v vcf4old ${dir}/${T_sample}.Strelka2.somatic.snv.vcf ${N_sample},${T_sample}" > ${dir}/anno_Strelka2_${T_sample}.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.Strelka2.somatic.snv.reformated.vcf.gz -o ${dir}/${T_sample}.Strelka2.somatic.snv.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m Strelka -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/anno_Strelka2_${T_sample}.sh
cd ${dir}
qsub -cwd -l vf=4G ${dir}/anno_Strelka2_${T_sample}.sh
done

########################## STEP 5.2.1 Platypus       ##########################

cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
python=/PROJ/HUMAN/share/Cancer/python/Python-2.7.11/bin/python
Platypus=/ifs/TJPROJ3/CANCER/share/software/Platypus-master
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
#bed=/PAN1/CANCER/ONT_test/wululu/Database/Reference/GRCh37/GATK/agilent_region.B37.bed

mkdir -p ${dir}/Somatic/${T_sample}/Platypus/
echo "${python} ${Platypus}/bin/Platypus.py callVariants --filterDuplicates=1 --nCPU=10 --bamFiles=${dir}/Mapping/${N_sample}/gatk/${N_sample}_bqsr.bam,${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --refFile=${reference} --output=${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.vcf" > ${dir}/Somatic/${T_sample}/Platypus/Platypus_${T_sample}.sh
# --regions=${bed}
echo "${python} ${Platypus}/extensions/Cancer/somaticMutationDetector.py --minPosterior 5  --tumourSample=${T_sample} --normalSample=${N_sample} --inputVCF ${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.vcf --outputVCF ${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.somatic.vcf" >> ${dir}/Somatic/${T_sample}/Platypus/Platypus_${T_sample}.sh
echo "perl /ifs/TJPROJ3/CANCER/share/software/CGR/bin/PlatypusSomaticIndelFilter.pl ${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.somatic.vcf INDEL ${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.somatic.indel.vcf" >> ${dir}/Somatic/${T_sample}/Platypus/Platypus_${T_sample}.sh
cd ${dir}/Somatic/${T_sample}/Platypus/
qsub -cwd -l vf=25G ${dir}/Somatic/${T_sample}/Platypus/Platypus_${T_sample}.sh
done



########################## STEP 5.2.12 Platypus-ann    ##########################
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Platypus
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b INDEL -z ${database} -d ${ANN} ${dir}/${T_sample}.platypus.somatic.indel.vcf ${N_sample},${T_sample}" > ${dir}/anno_Platypus_${T_sample}.sh
echo "python /ifs/TJPROJ3/CANCER/share/software/CGR/bin/vcf2maf_platypus.py -v ${dir}/${T_sample}.platypus.somatic.indel.reformated.vcf.gz -o ${dir}/${T_sample}.platypus.somatic.indel.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m platypus -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score -rt ${T_sample} -rn ${N_sample}" >> ${dir}/anno_Platypus_${T_sample}.sh

cd ${dir}
qsub -cwd -l vf=2G -q mepi1.q@tjcompute038.hpc,cancer1.q@tjcompute135.hpc,cancer1.q@tjcompute136.hpc,cancer1.q@tjcompute137.hpc ${dir}/anno_Platypus_${T_sample}.sh
done


cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/Platypus
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb
python=/PROJ/HUMAN/share/Cancer/python/Python-2.7.11/bin/python
Platypus=/ifs/TJPROJ3/CANCER/share/software/Platypus-master

echo "perl /ifs/TJPROJ3/CANCER/share/software/CGR/bin/PlatypusSomaticIndelFilter.pl ${dir}/Somatic/${T_sample}/Platypus/${T_sample}.platypus.somatic.vcf SNP ${dir}/${T_sample}.platypus.somatic.snv.vcf" > ${dir}/anno_Platypus_${T_sample}.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b SNP -z ${database} -d ${ANN} ${dir}/${T_sample}.platypus.somatic.snv.vcf ${N_sample},${T_sample}" >> ${dir}/anno_Platypus_${T_sample}.sh
echo "python /ifs/TJPROJ3/CANCER/share/software/CGR/bin/vcf2maf_platypus.py -v ${dir}/${T_sample}.platypus.somatic.snv.reformated.vcf.gz -o ${dir}/${T_sample}.platypus.somatic.snv.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m platypus -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score -rt ${T_sample} -rn ${N_sample}" >> ${dir}/anno_Platypus_${T_sample}.sh
cd ${dir}
qsub -cwd -l vf=2G -q mepi1.q@tjcompute038.hpc,cancer1.q@tjcompute135.hpc,cancer1.q@tjcompute136.hpc,cancer1.q@tjcompute137.hpc ${dir}/anno_Platypus_${T_sample}.sh
done

########################## STEP 5.3 intersect Mutect2 Strelka2  ##########################

cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
m_dir=${dir}/Somatic/${T_sample}/Mutect2
s_dir=${dir}/Somatic/${T_sample}/Strelka2_manta/results/variants
indel_dir=${dir}/Somatic/${T_sample}/indel
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

mkdir -p ${indel_dir}

names=`cat ${m_dir}/${T_sample}.mutect2.somatic.indel.vcf | grep "#CHROM" | awk '{print $(NF-1)","$NF}'`

echo "bedtools intersect -a ${m_dir}/${T_sample}.mutect2.somatic.indel.vcf -b ${s_dir}/${T_sample}.Strelka2.somatic.indel.vcf -header > ${indel_dir}/${T_sample}.mutect2.somatic.indel.vcf" > ${indel_dir}/${T_sample}.anno.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b INDEL -z ${database} -d ${ANN} ${indel_dir}/${T_sample}.mutect2.somatic.indel.vcf ${names}" >> ${indel_dir}/${T_sample}.anno.sh
echo "cd ${indel_dir}" >> ${indel_dir}/${T_sample}.anno.sh
echo "bcftools view -s ${T_sample},${N_sample} ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz -Ov | bgzip > tmp.vcf.gz" >> ${indel_dir}/${T_sample}.anno.sh
echo "mv tmp.vcf.gz ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz" >> ${indel_dir}/${T_sample}.anno.sh
echo "bcftools index -t ${T_sample}.mutect2.somatic.indel.reformated.vcf.gz"  >> ${indel_dir}/${T_sample}.anno.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${indel_dir}/${T_sample}.mutect2.somatic.indel.reformated.vcf.gz -o ${indel_dir}/${T_sample}.mutect2.somatic.indel.maf -t ${T_sample} -n ${N_sample} -s WES_agi -p Hiseq -m muTect -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${indel_dir}/${T_sample}.anno.sh

cd ${indel_dir}
qsub -cwd -l vf=2G ${indel_dir}/${T_sample}.anno.sh

done





## Platypus-germline
ls -d W* | while read id; do
sample=${id}
dir=/NJPROJ3/CANCER/Proj/WES.H101SC19080677.164_Human.BeiZhong_X101SC19080677-Z01-J001-B1-40.20191222
Platypus=/NJPROJ2/CANCER/TEST/wululu/tools/Platypus-master
reference=/NJPROJ2/CANCER/share/database/Genome/human/b37/human_g1k_v37_decoy.fasta
python=/NJPROJ2/CANCER/share/software/Python-2.7.11/bin/python

#mkdir -p ${dir}/Mutation/Platypus/${sample}
echo "${python} ${Platypus}/bin/Platypus.py callVariants --filterDuplicates=1 --nCPU=8 --bamFiles=${dir}/Mapping/${sample}/${sample}_bqsr.bam --refFile=${reference} --output=${dir}/Mutation/Platypus/${sample}/${sample}.platypus.vcf" > ${dir}/Mutation/Platypus/${sample}/Platypus_${sample}.sh
cd ${dir}/Mutation/Platypus/${sample}/
qsub -cwd -V -l vf=25G Platypus_${sample}.sh
done








## Facet
cat sample_list | awk '{print $3}' | awk -F'_' '{print $1}' | sort -u | while read id; do
sample=${id}_1
facets=/NJPROJ2/CANCER/TEST/wululu/miniconda3/envs/facets/bin/cnv_facets.R
dir=/NJPROJ3/CANCER/Proj/WES.H101SC19080677.164_Human.BeiZhong_X101SC19080677-Z01-J001-B1-40.20191222
vcf_file=/NJPROJ2/CANCER/TEST/wululu/database/00-All.vcf.gz
bed=/NJPROJ2/CANCER/share/database/Genome/human/b37/SureSelectXT.Human.All.Exon.V6/agilent_region.B37.bed
mkdir -p ${dir}/Somatic/${sample}/Facets/
echo "${facets} -t ${dir}/Mapping/${sample}/${sample}_bqsr.bam -n ${dir}/Mapping/${id}_b/${id}_b_bqsr.bam -vcf ${vcf_file} -o ${sample}.facets --gbuild hg19 --snp-mapq 15 --snp-baq 20 --depth 100 10000 --snp-nprocs 10 --rnd-seed 1234 --cval 25 150 --targets ${bed} --nbhd-snp auto" > ${dir}/Somatic/${sample}/Facets/Step1_snp-pileup.sh
cd ${dir}/Somatic/${sample}/Facets/
qsub -cwd -V -l vf=20G ${dir}/Somatic/${sample}/Facets/Step1_snp-pileup.sh
done





## DeepVariant -- sudo权限 预想是否可以尝试rust和conda构建非sudo环境



########################## STEP 5.4 Lancet WGS ##########################
for chr in {1..22} X Y MT; do
cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
lancet=/ifs/TJPROJ3/CANCER/share/software/lancet/lancet
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta

mkdir -p ${dir}/Somatic/${T_sample}/lancet/chr_${chr}
echo "${lancet} --normal ${dir}/Mapping/${N_sample}/gatk/${N_sample}_bqsr.bam --tumor ${dir}/Mapping/${T_sample}/gatk/${T_sample}_bqsr.bam --ref ${reference} --reg ${chr} --num-threads 8 > ${dir}/Somatic/${T_sample}/lancet/${T_sample}.lancet_${chr}.vcf" > ${dir}/Somatic/${T_sample}/lancet/chr_${chr}/lancet_${T_sample}_${chr}.sh
cd ${dir}/Somatic/${T_sample}/lancet/chr_${chr}
qsub -cwd -l vf=8G ${dir}/Somatic/${T_sample}/lancet/chr_${chr}/lancet_${T_sample}_${chr}.sh

done
done


cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
lancet=/ifs/TJPROJ3/CANCER/share/software/lancet/lancet
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/lancet/
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta

vcf-concat ${dir}/*vcf > ${dir}/${id}.raw.vcf

done


cat paired_sample | while read id; do
T_sample=`echo $id | awk '{print $1}'`
N_sample=`echo $id | awk '{print $2}'`
dir=/ifs/TJPROJ3/CANCER/Proj/WES.H101SC20061371.11_Human.HuaXiYiYuan_X101SC20061371-Z01-J001-B1-40.20201110/Somatic/${T_sample}/lancet
gatk=/ifs/TJPROJ3/CANCER/share/software/gatk-4.1.9.0/gatk
reference=/PUBLIC/software/CANCER/Database/Genome/human/b37/b37_gatk/human_g1k_v37_decoy.fasta
database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568
ANN=/PUBLIC/software/CANCER/Database/ANNOVAR/humandb

#cat ${dir}/${T_sample}.raw.vcf | grep "#" > ${dir}/${T_sample}.lancet_snv.vcf
#cat ${dir}/${T_sample}.raw.vcf | grep PASS | awk -F';' '{if($3=="TYPE=snv") print $0}' >> ${dir}/${T_sample}.lancet_snv.vcf
#cat ${dir}/${T_sample}.raw.vcf | grep "#" > ${dir}/${T_sample}.lancet_indel.vcf
#cat ${dir}/${T_sample}.raw.vcf | grep PASS | awk -F';' '{if($3!="TYPE=snv") print $0}' >> ${dir}/${T_sample}.lancet_indel.vcf

echo "cd ${dir}" > ${dir}/${T_sample}.lancet_snv.anno.sh
echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b SNP -z ${database} -d ${ANN} ${dir}/${T_sample}.lancet_snv.vcf ${T_sample}" >> ${dir}/${T_sample}.lancet_snv.anno.sh
echo "bcftools view -s ${T_sample},${N_sample} ${T_sample}.lancet_snv.reformated.vcf.gz -Ov | bgzip > tmp.vcf.gz" >> ${dir}/${T_sample}.lancet_snv.anno.sh
echo "mv tmp.vcf.gz ${T_sample}.lancet_snv.reformated.vcf.gz" >> ${dir}/${T_sample}.lancet_snv.anno.sh
echo "bcftools index -t ${T_sample}.lancet_snv.reformated.vcf.gz"  >> ${dir}/${T_sample}.lancet_snv.anno.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.lancet_snv.reformated.vcf.gz -o ${dir}/${T_sample}.lancet.somatic.snv.maf -t ${T_sample} -s WGS -p Hiseq -m muTect -f avsnp150 -x CADD,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/${T_sample}.lancet_snv.anno.sh


database=GeneName,refGene,Genecode27,wgRna,cytoBand,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp150,cosmic,clinvar_20170905,gwasCatalognew,1000g2015aug_Chinese,1000g2015aug_eas,1000g2015aug_all,esp6500siv2_all,exac03_ALL_EAS,ljb30_sift,ljb30_pp2hvar,ljb30_pp2hdiv,ljb30_mt,gerp++gt2,caddgt10,NovoDb_WES_2573,NovoDb_WGS_568

echo "/PUBLIC/software/CANCER/Module/CancerGenome/Script/Annotation/Var_annotation_cancer_ANNOVAR2015Mar22.20180307.sh -a somatic -r ${reference} -u hg19 -b INDEL -z ${database} -d ${ANN} ${dir}/${T_sample}.lancet_indel.vcf ${T_sample}" > ${dir}/${T_sample}.lancet_indel.anno.sh
echo "bcftools view -s ${T_sample},${N_sample} ${T_sample}.lancet_indel.reformated.vcf.gz -Ov | bgzip > tmp.vcf.gz" >> ${dir}/${T_sample}.lancet_indel.anno.sh
echo "mv tmp.vcf.gz ${T_sample}.lancet_indel.reformated.vcf.gz" >> ${dir}/${T_sample}.lancet_indel.anno.sh
echo "bcftools index -t ${T_sample}.lancet_indel.reformated.vcf.gz"  >> ${dir}/${T_sample}.lancet_indel.anno.sh
echo "python /PUBLIC/software/CANCER/Module/CancerGenome/Advance/vcf2maf.py -v ${dir}/${T_sample}.lancet_indel.reformated.vcf.gz -o ${dir}/${T_sample}.lancet.somatic.indel.maf -t ${T_sample} -n ${N_sample} -s WGS -p Hiseq -m muTect -f avsnp150 -x caddindel,1000g2015aug_all,ExAC_ALL,SIFT_score,Polyphen2_HVAR_score,Polyphen2_HDIV_score" >> ${dir}/${T_sample}.lancet_indel.anno.sh


cd ${dir}
qsub -cwd -l vf=2G ${dir}/${T_sample}.lancet_snv.anno.sh
qsub -cwd -l vf=2G ${dir}/${T_sample}.lancet_indel.anno.sh

done

