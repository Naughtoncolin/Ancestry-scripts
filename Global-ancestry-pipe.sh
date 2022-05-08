# Author: Colin Naughton, inspired by Varsha Bhat
# May 2022
# Global Ancestry prediction pipeline.


##################### Phase autosomes with Beagle #########################################
# FUTURE: Check for beagle map files. If not present, download. Delete after phasing
# FUTURE: Check for reference files. If not present, download. Is this necessary if no imputation is done?
# FUTURE: Check if "phase_out" directory exists, if not, then make
#Note: Expects 128GB of memory total per node
beagle_out=phase_out
find . -maxdepth 2 -name '*chr[1-9]*[mapgz]' | sort -k2 -t 'h' | parallel -N2 --sshloginfile $PBS_NODEFILE -j1 \
java -Xmx128g -jar /storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
gt=$PBS_O_WORKDIR/{2} \
map=$PBS_O_WORKDIR/{1} \
out=$PBS_O_WORKDIR/$beagle_out/{=1 's/\.\/beagle_maps\/plink\.//;s/.GRCh38.map//' =} \
impute=false \

#################### Get sample names specific to each query population ##################################
# Run from base directory
#cut -f1 igsr-ceu.tsv.tsv | tail -n +2  > CEU_samples_nogender.txt
#cut -f1 igsr-yri.tsv.tsv | tail -n +2  > YRI_samples_nogender.txt
ls igsr* | xargs -I{} tail -n +2 {} | cut -f1 >> query_population_sampleNames_noGender

# Include gender with sample names. NOTE: This was causing an issue for downstream analysis, so the "nogender" version was used.
#sed 's/female/F/' igsr-ceu.tsv.tsv | sed  's/male/M/' | cut -f1,2 |tail -n +2 > CEU_samples.txt 
#sed 's/female/F/' igsr-yri.tsv.tsv | sed  's/male/M/' | cut -f1,2 |tail -n +2 > YRI_samples.txt
ls igsr* | xargs -I{} tail -n +2 {} | sed 's/female/F/;s/male/M/' | cut -f1,2 >> query_population_sampleNames_withGender

############# Subset 1000 Genomes VCF based on population-specific sample names #########################
# Version2: For PBS scripts. Spreads jobs over multiple nodes
# Change later to include variable for conda env path
# TODO: Add output dir path
# TODO: Add variable associated with reference population acronyms; add to output path
# TODO: Check for output directory; if not present, make
ref_sub_out=ref/subset
find ref/ALL_POP/ -name *ALL.chr[0-9]*.gz | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools view \
--threads 12 \
--force-samples \
-S $PBS_O_WORKDIR/query_population_sampleNames_noGender \
-Oz \
-o $PBS_O_WORKDIR/$ref_sub_out/{=1 's/ref\/ALL_POP\/ALL\.//;s/.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz//' =}.CEU-YRI.vcf.gz \
$PBS_O_WORKDIR/{1} #Input

#################### Index Phased Autosomes & Subsetted Ref VCFs#############################################
find ref/subset phase_out -name "*chr*gz" | parallel -N1 --sshloginfile $PBS_NODEFILE -j1 \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools index \
--threads 12 \
$PBS_O_WORKDIR/{1}

#################### Intersect VCFs ##########################################################
# Intersect all chromosome vcfs between St. Jude patient data and CEU-YRI 1000 Genomes data
# Note: VCFs must be indexed
# Batch rename
ls CEU-YRI.ALL.* | xargs -P8 -I{} mv {} "$(echo {} | sed 's/.ALL//')" #Didn't work
ls *.gz  | sed 'p;s/.gz//' | xargs -n2 -P8 mv
ls | sed 'p;s/\.shapeit2\_integrated\_snvindels\_v2a\_27022019\.GRCh38\.phased\.vcf\.gz\_out//' | xargs -n2 -P8 mv # This works

intersect_out=intersect_out
find ref/subset phase_out -name '*chr*.gz' | sort  -k2 -t 'c' | parallel -N2 -j1 --sshloginfile $PBS_NODEFILE \
-Oz \
-p $intersect_out/{=2 's/phase_out\///;s/\.vcf\.gz//' =} \
$PBS_O_WORKDIR/{2} $PBS_O_WORKDIR/{1} #This worked

################## Add variant names to vcf #################################################
# The following loop works; the subsequent oneliners don't 
for file in intersect_out/chr*/0002*gz; do
	out_name=$(echo $file | sed 's/intersect_out\///;s/\/0002//')
	/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools annotate \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	-Oz \
	--threads 8 \
	-o $annotation_out/$out_name \
	$file
done

for file in intersect_out/chr*/0003*gz; do
	out_name=$(echo $file | sed 's/intersect_out\///;s/\/0003//')
	/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools annotate \
	--set-id +'%CHROM:%POS:%REF:%FIRST_ALT' \
	-Oz \
	--threads 12 \
	-o $annotation_out/CEU-YRI_$out_name \
	$file
done
#Test
find intersect_out/chr22 -maxdepth 2 -name '0002*gz' | parallel -N1 -j1 echo {1}

find intersect_out/chr22 -maxdepth 2 -name '0002*gz' | parallel -N1  \
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {1} | tail

find phase_out -maxdepth 2 -name '*chr22*gz' | parallel -N1 -j1 \
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {1} | tail
#######

find intersect_out/chr22 -maxdepth 2 -name '0002*gz' | xargs -I{} \
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {} | tail

bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' intersect_out/chr22/0002.vcf.gz | tail

annotation_out=annotation_out
find intersect_out -maxdepth 2 -name '0002*gz' | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools annotate \
--force \
--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
-Oz \
--threads 12 \
-o $PBS_O_WORKDIR/$annotation_out/{=1 's/intersect_out\///;s/\/0002.vcf.gz//' =}.vcf.gz \
$PBS_O_WORKDIR/{1}

annotation_out=annotation_out
find intersect_out -maxdepth 2 -name '0002*gz' | xargs -I{} \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools annotate \
--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
-Oz \
--threads 12 \
-o $PBS_O_WORKDIR/$annotation_out/{=1 's/intersect_out\///;s/\/0002.vcf.gz//' =}.vcf.gz \
$PBS_O_WORKDIR/{1}

###################### Convert to plink file.###############################################
# NOTE: All samples did not have sex recorded in the output; they came up as "ambiguous".
# Originally run in the beagle directory.
# Used Plink 1.9
# FUTURE: Make edits for use in PBS script
plink_out=plink_out
find annotation_out/ -name '*.gz' | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
plink --vcf {1} --out $plink_out/{=1 's/annotation_out\///;s/\.vcf\.gz//' =} --make-bed

#For ref
plink_out=plink_out
find intersect_out -maxdepth 2 -name '0003*gz' | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/plink \
--vcf $PBS_O_WORKDIR/{1} \
--out $PBS_O_WORKDIR/$plink_out/CEU-YRI_{=1 's/intersect_out\///;s/\/0003\.vcf\.gz//' =} \
--make-bed \
--threads 12
# The following should make the plink file and set the variant IDs in the process, but the $1 and $2 were not functioning as described in the manual. May need to use Plink2 or rename vcf in prior step.
find non-impute/ -maxdepth 3 -name '*.gz' | parallel -j4 plink --vcf {1} --out plink-out/{=1 's/non-impute\///;s/\.vcf\.gz//' =} --make-bed --set-missing-var-ids @:#:\$1,\$2
 
# Merge chromosomal plink files for downstream ancestry analysis with 'admixture'
ls plink_out/*[bf]* | xargs -n3 echo > mergedList.txt
plink --out plink_merge-all_out --merge-list mergedList.txt # Got warnings for multiallelic sites.


################## Associate ancestry with sample names####################################
# Creates file with ancestry acronyms on each line in the case of ref samples, or '-' in the case of non-ref samples
# TODO: Add output name variable.
plink_fam_file=plink_merge-all_out.fam
# Read through sample names in the plink fam file
cat $plink_fam_file | while read line; do
	sample_name=$(echo $line | cut -f1 -d ' ')

	# Check if sample name has prefix of non-reference samples.
	sample_prefix=$(echo $sample_name | grep -o ^..)
	if [[ $sample_prefix == 'SJ' ]]; then
		echo '-' >> plink_merge-all_out.pop
		continue
	fi
	# Loop through ref metadata to associate ancestry with sample name
	cat igsr* | while read line; do
		pop=$(echo $line | cut -f4 -d ' ') # Population acrononym e.g. "CEU"
		pop_sample=$(echo $line | cut -f1 -d ' ') # Reference sample name
		#Check if sample name from plink fam file is same as sample name from ref metadata
		if [[ $sample_name == $pop_sample ]]; then
			echo $pop >> plink_merge-all_out.pop
			break
		fi
	done
done &

################### Run admixture #####################################################
admixture ../plink_merge-all_out.bed 2 --supervised -j1


# R Code For plotting admixture results
tbl=read.table("chr22.2.Q")
newtbl <- tbl[order(tbl$V1),]
barplot(t(as.matrix(newtbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
newtbl2 <- tbl[order(tbl$V2),]
barplot(t(as.matrix(newtbl2)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)