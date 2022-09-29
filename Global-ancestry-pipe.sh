# Author: Colin Naughton, inspired by Varsha Bhat
# May 2022
# Global Ancestry prediction pipeline.

# TODO: Get variable name for number of threads

##################### Phase autosomes with Beagle #########################################
# The presence of a "chr" prefix on chromosomes on either the map of vcf file but not the other causes an error.
# Changed the map files to include the "chr" prefix via:
seq 1 22 | xargs -I{} sed -i 's/^/chr/' beagle_maps/plink.chr{}.GRCh38.map

# TODO: Check for beagle map files. If not present, download. Delete after phasing
# TODO: Check for reference files. If not present, download. Is this necessary if no imputation is done?
# TODO: Check if "phase_out" directory exists, if not, then make
# TODO: Check for total memory available on each node; set as variable
#Note: Expects 96GB of memory total per node with 12 threads per node 
mkdir $PBS_O_WORKDIR/phase_out
beagle_maps=beagle_maps
raw_data=raw/genotype
beagle_out=phase_out
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
java -Xmx96g -jar /storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/software/beagle.22Jul22.46e.jar \
gt=$PBS_O_WORKDIR/$raw_data/chr{1}_Emory-alloimmunization.vcf.gz \
map=$PBS_O_WORKDIR/$beagle_maps/plink.chr{1}.GRCh38.map \
out=$PBS_O_WORKDIR/$beagle_out/chr{1}_phased_Emory-alloimmunization \
impute=false \
nthreads=12

# Remove "chr" prefixes from vcf entries.
# TODO: Do this before phasing!
seq 1 22 | parallel \
"zcat phase_out/chr{1}_phased_Emory-alloimmunization.vcf.gz | sed 's/^chr//' | bgzip > phase_out/chr{1}_phased_Emory-alloimmunization.fixedChrPrefix.vcf.gz"

# Remove multiallelic sites & add missing variant IDs
# sort removes duplicates based on location; 
# awk replaces variants IDs with CHR_POS_REF_ALT from the same entry
mkdir $PBS_O_WORKDIR/rename_dedup_out
echo "Remove multiallelic sites & add missing variant IDs" > rename_dedup_out/README
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
"zcat $PBS_O_WORKDIR/phase_out/chr{1}_phased_Emory-alloimmunization.fixedChrPrefix.vcf.gz | (head -n 9 && tail -n +10 | \
sort -n -u -t $'\t' -k2 | \
awk -v OFS='\t' '{\$3=\$1\"_\"\$2\"_\"\$4\"_\"\$5;print \$0}') | \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bgzip > \
$PBS_O_WORKDIR/rename_dedup_out/chr{1}_phased_Emory-alloimmunization.fixedChrPrefix.biallelic.renamed.vcf.gz"

#Renaming ref variant IDs requires "head -n 22 && tail -n +23"

#################### Get sample names specific to each query population ##################################
# Originally used data from /~/gridftp/1000g/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
# Now using updated data from 
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
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools index \
$PBS_O_WORKDIR/rename_dedup_out/chr{1}_phased_.vcf.gz

seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools index \
$PBS_O_WORKDIR/rename_dedup_out/chr{1}.CEU-YRI.vcf.gz

#################### Intersect VCFs ##########################################################
# Intersect all chromosome vcfs between St. Jude patient data and CEU-YRI 1000 Genomes data
# Note: VCFs must be indexed
# Batch rename (Extra if needed)
#ls *.gz  | sed 'p;s/.gz//' | xargs -n2 -P8 mv

mkdir $PBS_O_WORKDIR/intersect_out
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bEmory-alloimmunizationcftools isec \
-Oz \
-p $PBS_O_WORKDIR/intersect_out/chr{1} \
$PBS_O_WORKDIR/rename_dedup_out/chr{1}_phased_Emory-alloimmunization.fixedChrPrefix.biallelic.renamed.vcf.gz \
$PBS_O_WORKDIR/rename_dedup_out/subset/chr{1}.CEU-YRI.vcf.gz 

################## Add variant names to vcf #################################################
# Remove?
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

###################### Convert to & merge plink files.###############################################
# NOTE: All samples did not have sex recorded in the output; they came up as "ambiguous"..
# Used Plink 1.9
# FUTURE: Make edits for use in PBS script
# Query data conversion
mkdir $PBS_O_WORKDIR/plink_out
seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/plink \
--vcf $PBS_O_WORKDIR/intersect_out/chr{1}/0002.vcf.gz \
--out $PBS_O_WORKDIR/plink_out/Emory-alloimmunization_chr{1} \
--make-bed \
--threads 12

seq 1 22 | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/plink \
--vcf $PBS_O_WORKDIR/intersect_out/chr{1}/0003.vcf.gz \
--out $PBS_O_WORKDIR/plink_out/CEU-YRI_chr{1} \
--make-bed \
--threads 12

# Merge chromosomal plink files for downstream ancestry analysis with 'admixture'
ls plink_out/*[bf]* | xargs -n3 echo > mergedList.txt
plink --make-bed --out plink_merge-all_out --merge-list mergedList.txt

################## Associate ancestry with sample names####################################
# Creates file with ancestry acronyms on each line in the case of ref samples, or '-' in the case of non-ref samples
# Takes a few minutes.
# TODO: Add output name variable.
# TODO: Add sample prefix variable.
plink_fam_file=plink_merge-all_out.fam
# Read through sample names in the plink fam file
cat $plink_fam_file | while read line; do
	sample_name=$(echo $line | cut -f1 -d ' ')

	# Check if sample name has prefix of non-reference samples.
	sample_prefix=$(echo $sample_name | grep -o ^..)
	if [[ $sample_prefix == 'NW' ]]; then
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

# Make .tsv with admixture ancestry estimates including header.
# Note: SJID column is needed for merging files in PRS calculation (as of 6/18/2022)
grep SJC* ancestry_estimates.txt | cut -f1,2,8,9 | \
awk 'BEGIN {print "FID\tIID\t%CEU\t%YRI\tSJID"} {print $0, $1_$2}' | \
sed 's/\s\+/\t/g' > query_ancestry_estimates.tsv

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