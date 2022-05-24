# Author: Colin Naughton, inspired by Varsha Bhat
# May 2022
# Global Ancestry prediction pipeline.

# TODO: Get variable name for number of threads

##################### Phase autosomes with Beagle #########################################
# TODO: Check for beagle map files. If not present, download. Delete after phasing
# TODO: Check for reference files. If not present, download. Is this necessary if no imputation is done?
# TODO: Check if "phase_out" directory exists, if not, then make
# TODO: Check for total memory available on each node; set as variable
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
ls igsr* | xargs -I{} tail -n +2 {} | cut -f1 >> query_population_sampleNames_noGender

# Include gender with sample names. NOTE: This was causing an issue for downstream analysis, so the "nogender" version was used.
ls igsr* | xargs -I{} tail -n +2 {} | sed 's/female/F/;s/male/M/' | cut -f1,2 >> query_population_sampleNames_withGender

############# Subset 1000 Genomes VCF based on population-specific sample names #########################
# TODO: Include variable for conda env path
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
# Batch rename (Extra if needed)
#ls *.gz  | sed 'p;s/.gz//' | xargs -n2 -P8 mv

intersect_out=intersect_out
find ref/subset phase_out -name '*chr*.gz' | sort  -k2 -t 'c' | parallel -N2 -j1 --sshloginfile $PBS_NODEFILE \
-Oz \
-p $intersect_out/{=2 's/phase_out\///;s/\.vcf\.gz//' =} \
$PBS_O_WORKDIR/{2} $PBS_O_WORKDIR/{1} #This worked

################## Add variant names to vcf #################################################
# The following loops work; the subsequent oneliners don't 
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
#find intersect_out/chr22 -maxdepth 2 -name '0002*gz' | parallel -N1 -j1 echo {1}

#find intersect_out/chr22 -maxdepth 2 -name '0002*gz' | parallel -N1  \
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {1} | tail

#find phase_out -maxdepth 2 -name '*chr22*gz' | parallel -N1 -j1 \
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {1} | tail
#######
# Other attempts

#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' intersect_out/chr22/0002.vcf.gz | tail

#annotation_out=annotation_out
#find intersect_out -maxdepth 2 -name '0002*gz' | xargs -I{} \
#/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools annotate \
#--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
#-Oz \
#--threads 12 \
#-o $PBS_O_WORKDIR/$annotation_out/{=1 's/intersect_out\///;s/\/0002.vcf.gz//' =}.vcf.gz \
#$PBS_O_WORKDIR/{1}

###################### Convert to plink file.###############################################
# NOTE: All samples did not have sex recorded in the output; they came up as "ambiguous"..
# Used Plink 1.9
# FUTURE: Make edits for use in PBS script
# Query data conversion
plink_out=plink_out
find annotation_out/ -name '*.gz' | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
plink --vcf {1} \
--out $plink_out/{=1 's/annotation_out\///;s/\.vcf\.gz//' =} \
--make-bed \
--threads 12

# Ref data conversion
plink_out=plink_out
find intersect_out -maxdepth 2 -name '0003*gz' | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/plink \
--vcf $PBS_O_WORKDIR/{1} \
--out $PBS_O_WORKDIR/$plink_out/CEU-YRI_{=1 's/intersect_out\///;s/\/0003\.vcf\.gz//' =} \
--make-bed \
--threads 12
 
# Merge chromosomal plink files for downstream ancestry analysis with 'admixture'
ls plink_out/*[bf]* | xargs -n3 echo > mergedList.txt
plink --out plink_merge-all_out --merge-list mergedList.txt # Got warnings for multiallelic sites. Still merged when variant names were unique...I think.


################## Associate ancestry with sample names####################################
# Creates file with ancestry acronyms on each line in the case of ref samples, or '-' in the case of non-ref samples
# Takes a few minutes.
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
paste plink_merge-all_out.fam plink_merge-all_out.pop admixture_out/plink_merge-all_out.2.Q > admixture_out/ancestry_estimates.txt # Associate global ancestry esitmates with sample names.

grep SJC* ancestry_estimates.txt | cut -f8,9 > query_ancestry_estimates-noNames.tsv
# R Code For plotting admixture results
# TODO: Run from bash script
# TODO: Add population legend
tbl=read.table("query_ancestry_estimates-noNames.tsv")
newtbl <- tbl[order(tbl$V1),]
barplot(t(as.matrix(newtbl)), col=c("darkblue","red"), fill=c("darkblue","red"),
main="Global Ancestry Proportions",
xlab="Individual #", ylab="Ancestry", border=NA) +
  legend("top", 
         legend = c("European (CEU)", "African (YRI)"), 
         fill = c("darkblue", "red"))

################### Local Ancestry Inference with RFMix2 ##########################
# TODO: Isolate local ancestry inference pipeline in separate script
# TODO: Check for out folder
rfmix_out=rfmix_out

# Make file mapping sample names to ancestry required by RFMix2
ls igsr* | xargs -I{} tail -n +2 {} | cut -f1,4 > query_population_sampleNames_ancestry

# Edit genomic map files to format required by RFMix2
for file in beagle_maps/*chr[1-9]*; do
	out_name=$rfmix_out/$(echo $file | sed 's/beagle_maps\/plink/rfmix/')
	awk '{print $1"\t"$4"\t"$3}' $file > $out_name
done

# Run RFMix2
# Note: Took ~1.5hours with 22 nodes at 24 cores/node.
# Note: The developers recommended VCFs be joined prior for some reason.
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/rfmix \
-f $PBS_O_WORKDIR/intersect_out/chr{1}/0002.vcf.gz \
-r $PBS_O_WORKDIR/intersect_out/chr{1}/0003.vcf.gz \
-m $PBS_O_WORKDIR/query_population_sampleNames_ancestry \
-g $PBS_O_WORKDIR/$rfmix_out/rfmix.chr{1}.GRCh38.map \
-o $PBS_O_WORKDIR/$rfmix_out/chr{1} \
--chromosome={1} \
--n-threads=24

################### Tractor Analysis ######################################
git clone https://github.com/Atkinson-Lab/Tractor

# Step 1: Recovering tracts
# Required numpy installation
# TODO: Submit pull request for fixes
# Identify and correct switch errors in local ancestry calls
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
python3 $PBS_O_WORKDIR/Tractor/UnkinkMSPfile.py \
--msp $PBS_O_WORKDIR/rfmix_out/chr{1}

# Unzip phased VCFs
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
gunzip -c $PBS_O_WORKDIR/phase_out/chr{1}.vcf.gz '>' $PBS_O_WORKDIR/phase_out/chr{1}.vcf

# Correcting switch errors in genotype data
# Should this step be done on the subset VCF? Or full VCF?
# TODO: Un-hardcode python.
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/python \
$PBS_O_WORKDIR/Tractor/UnkinkGenofile.py \
--switches $PBS_O_WORKDIR/rfmix_out/chr{1}.switches.txt \
--genofile $PBS_O_WORKDIR/annotation_out/chr{1}.vcf

# Step 2: Extracting tracts & ancestral dosages
# TODO: Wasn't working on multiple nodes, additionally, a single node only supports 15 parallel jobs at one time regardless of whether there is enough cores
seq 1 22 | parallel -N1 -j22 --sshloginfile $PBS_NODEFILE \
python3 $PBS_O_WORKDIR/Tractor/ExtractTracts.py \
--msp $PBS_O_WORKDIR/rfmix_out/chr{1} \
--vcf $PBS_O_WORKDIR/tractor_out/unkinked/chr{1}.unkinked \
--num-ancs 2