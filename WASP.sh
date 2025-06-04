#!/bin/bash
#SBATCH -J MH_pre_SYF_186
#SBATCH -p LiuLiHe
##SBATCH --nodelist=b[01]
##SBATCH --array=1
#SBATCH --cpus-per-task=16
#SBATCH -o ./log/%x.o
#SBATCH -e ./log/%x.e

source /work/home/hxxgroup01/miniconda3/bin/activate /work/home/hxxgroup01/miniconda3
#================================================================================
######Files######
WorkDir='/work/home/hxxgroup01/muhongmei/MH/QC/step2_star/08_WASP'  
# DataDir='/work/home/hxxgroup01/muhongmei/MH/QC/step2_star/08_WASP/00pre/data'
IndexDir='/work/home/hxxgroup01/muhongmei/database/chicken_7b'
ID='MH_pre_SYF_186'
Thread='16'
PhasedSNPs='/work/home/hxxgroup01/muhongmei/MH/jiaozheng/panel/genomics/all_chr/0401_MH_wgs_all_chr.phased_rename.vcf.gz'
#================================================================================
#conda create -n WASP python=3.7 numpy scipy pysam hdf5 pytables
#wget https://github.com/bmvdgeijn/WASP/archive/master.zip
export LD_LIBRARY_PATH=/work/home/hxxgroup01/miniconda3/envs/WASP/lib:$LD_LIBRARY_PATH
WASPdir=/work/home/hxxgroup01/muhongmei/MH/QC/step2_star/08_WASP/env/WASP-master
Snp2h5=$WASPdir/snp2h5/snp2h5
Find_intersecting_snps=$WASPdir/mapping/find_intersecting_snps.py
Filter_remapped_reads=$WASPdir/mapping/filter_remapped_reads.py
# Samtools=~/miniconda3/envs/WASP/bin/samtools
# mamba activate WASP
#================================================================================
######Running######
#1.Extracting the WGS snps from a sample by RNA-seq sites
Ind=$(echo $ID | awk -F'_' '{print $1"_"$2"_"$4}') #cut tissue
#cat $WorkDir/Run/06-Rediting/DnaRna_*/outTable_* | \
                #awk 'FS="\t" {if ($8!="-" && $10>=10 && $13!="-" && substr($1,1,2)!="NW" && $1!="Region") print $1"\t"$2}' \
                #> $WorkDir/Run/08-WASP.candidates.txt

less $WorkDir/Run/07-Rediting/${ID}_candidates.gff | cut -f 1,4 > $WorkDir/Run/08-WASP/${ID}.candidates.txt
vcftools --gzvcf $PhasedSNPs \
                --indv $Ind \
                --positions $WorkDir/Run/08-WASP/${ID}.candidates.txt \
                --recode \
                --recode-INFO-all \
                --stdout | \
                grep -v "0|0" | \
                bgzip > $WorkDir/Run/08-WASP/${ID}.vcf.gz
tabix -p vcf $WorkDir/08-WASP/${ID}.vcf.gz

# mkdir $WorkDir/Run/08-WASPvcf
Chrom=$(zcat $WorkDir/Run/08-WASP/${ID}.vcf.gz | grep -v "#" | cut -f 1 | uniq | tr '\n' ' ')
for i in $Chrom
do
vcftools --gzvcf $WorkDir/Run/08-WASP/${ID}.vcf.gz \
                --recode \
                --recode-INFO-all \
                --stdout \
                --chr $i | \
                bgzip > $WorkDir/Run/08-WASP/${ID}.$i.vcf.gz
done

mkdir -p $WorkDir/Run/10-WASP
cp $WorkDir/Run/08-WASP/${ID}.vcf.gz $WorkDir/Run/10-WASP/${ID}.vcf.gz
cp $WorkDir/Run/08-WASP/${ID}.vcf.gz.tbi $WorkDir/Run/10-WASP/${ID}.vcf.gz.tbi
#================================================================================
#2.Creating HDF5 SNP files
$Snp2h5 --chrom $IndexDir/Chicken.chromInfo.txt \
                --format vcf \
                --haplotype $WorkDir/Run/10-WASP/${ID}.haplotypes.h5 \
                --snp_index $WorkDir/Run/10-WASP/${ID}.snp_index.h5 \
                --snp_tab $WorkDir/Run/10-WASP/${ID}.snp_tab.h5 \
                $WorkDir/Run/08-WASP/${ID}.*.vcf.gz
#================================================================================
###¹¹½¨Ë÷Òý
samtools index /work/home/hxxgroup01/muhongmei/MH/QC/step2_star/results/${ID}/Aligned.sortedByCoord.out.bam
#================================================================================
#3. to identify reads that may have mapping biases
mkdir -p $WorkDir/Run/08-WASP/${ID}

source /work/home/hxxgroup01/miniconda3/bin/activate WASP
python $Find_intersecting_snps \
                --is_paired_end \
                --is_sorted \
                --output_dir $WorkDir/Run/08-WASP/${ID} \
                --snp_tab $WorkDir/Run/10-WASP/${ID}.snp_tab.h5 \
                --snp_index $WorkDir/Run/10-WASP/${ID}.snp_index.h5 \
                --haplotype $WorkDir/Run/10-WASP/${ID}.haplotypes.h5 \
                /work/home/hxxgroup01/muhongmei/MH/QC/step2_star/results/${ID}/Aligned.sortedByCoord.out.bam
conda deactivate

samtools sort --threads $Thread \
                $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.to.remap.bam \
                -o $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.to.remap.sort.bam
samtools index $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.to.remap.sort.bam
#================================================================================
#4.Remapping the sequencing reads to reference genomes
STAR --runThreadN $Thread \
                --genomeDir $IndexDir/chicken7b_STAR \
                --sjdbGTFfile $IndexDir/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf \
                --runMode alignReads \
                --readFilesIn $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.remap.fq1.gz $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.remap.fq2.gz \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix $WorkDir/Run/08-WASP/${ID}/ \
                --outBAMsortingThreadN $Thread \
                --twopassMode Basic \
                --outSAMunmapped Within \
                --outFilterType BySJout \
                --outSAMattrRGline ID:$ID SM:$ID PL:DNBSEQ \
                --chimSegmentMin 10 \
                --twopass1readsN -1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverLmax 0.03 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --chimOutType SeparateSAMold \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --sjdbOverhang 100
samtools index $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.bam
#================================================================================
#5.filter out reads where one or more of the allelic versions of the reads fail to map back to the same location as the original read.
source /work/home/hxxgroup01/miniconda3/bin/activate WASP
python $Filter_remapped_reads \
                $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.to.remap.sort.bam \
                $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.bam \
                $WorkDir/Run/08-WASP/${ID}/${ID}.keep.bam
# conda deactivate
#================================================================================
#6.for Splicing. Merge bam for a complete set of mappability-filtered aligned reads.
samtools merge --threads $Thread \
                -f \
                -o $WorkDir/Run/08-WASP/${ID}/${ID}.keep.merge.bam \
                $WorkDir/Run/08-WASP/${ID}/${ID}.keep.bam \
                $WorkDir/Run/08-WASP/${ID}/Aligned.sortedByCoord.out.keep.bam
samtools sort --threads $Thread \
                $WorkDir/Run/08-WASP/${ID}/${ID}.keep.merge.bam \
                -o $WorkDir/Run/09-WASP/${ID}.forSplicing.bam
samtools index -@ $Thread $WorkDir/Run/09-WASP/${ID}.forSplicing.bam
#================================================================================
#7.for ASE. filter and duplicate reads.
samtools view --threads $Thread \
                -bq 255 \
                -o $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.bam \
                $WorkDir/Run/08-WASP/${ID}/${ID}.keep.bam
samtools sort --threads $Thread \
                $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.bam \
                -o $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.bam
samtools index -@ $Thread $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.bam
# #==================================================
mkdir -p $WorkDir/Run/10-WASP/${ID}
/work/home/hxxgroup01/muhongmei/software/sentieon-genomics-202503/bin/sentieon driver -t $Thread -i $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.bam --algo LocusCollector --fun score_info $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.score.gz

/work/home/hxxgroup01/muhongmei/software/sentieon-genomics-202503/bin/sentieon driver -t $Thread -i $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.bam --algo Dedup --rmdup --score_info $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.score.gz --metrics $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.dedup_metrics.txt $WorkDir/Run/10-WASP/${ID}/${ID}.keep.filter.ASE.dedup.bam
# #==================================================

# source /share/apps/anaconda3/bin/activate PanCattle
# picard -Xmx20g AddOrReplaceReadGroups \
#                 --TMP_DIR $WorkDir/Run/temp \
#                 -I $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.bam \
#                 -O $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.Groups.bam \
#                 --RGID 4 \
#                 --RGLB lib1 \
#                 --RGPL illumina \
#                 --RGPU run \
#                 --RGSM 20 \
#                 --CREATE_INDEX true \
#                 --VALIDATION_STRINGENCY SILENT \
#                 --SORT_ORDER coordinate 
# picard -Xmx20g MarkDuplicates \
#                 --TMP_DIR $WorkDir/Run/temp \
#                 -I $WorkDir/Run/08-WASP/${ID}/${ID}.keep.filter.sort.Groups.bam \
#                 -O $WorkDir/Run/10-WASP/${ID}.forASE.bam \
#                 -M $WorkDir/Run/08-WASP/${ID}/${ID}.keep.marked_dup_metrics.txt \
#                 --CREATE_INDEX true \
#                 --VALIDATION_STRINGENCY SILENT \
#                 --REMOVE_DUPLICATES true
#================================================================================
#8.Delete big temporary files
rm -rf $WorkDir/Run/08-WASP/${ID}*  #注意不要删除其他的内容，比如1，就会把11，12都删掉
# rm -r $WorkDir/Run/08-WASP/${ID}/
#================================================================================
