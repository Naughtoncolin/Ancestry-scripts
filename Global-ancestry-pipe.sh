# Author: Colin Naughton, inspired by Varsha Bhat
# May 2022
# Global Ancestry prediction pipeline.


##################### Phase autosomes with Beagle #########################################
# FUTURE: Check for beagle map files. If not present, download. Delete after phasing
# FUTURE: Check for reference files. If not present, download. Is this necessary if no imputation is done?

java -jar ~/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
gt="../intersect/CEU-YRI.chr22/0003.vcf.gz" \
ref="../ref/CEU-YRI.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" \
map="plink.chr22.GRCh38.map" \
out="beagle_stuff/non-impute_no-ref/chr22_beagleTest" \
impute=false \
nthreads=4

java -jar ~/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
gt="patient_data/baylorOnly.chr22.vqsr.decompose_norm_uniq.PASS_only.vcf.gz" \
map="beagle_stuff/plink.chr22.GRCh38.map" \
out="test/chr22_beagle-1st-step" \
impute=false \
nthreads=4

#find . -maxdepth 2 -name '*chr[1-9]*[mapgz]' | sort -k2 -t 'h' | parallel -N3 -j13 \
#java -jar ~/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
#gt={3} \
#ref={2} \
#map={1} \
#out=beagle_stuff/{=1 's/..beagle_stuff.plink\.//;s/.GRCh38.map//' =} \
#impute=false \

#For use in PBS script
#--sshloginfile $PBS_NODEFILE #List of nodes available
# With ref
find . -maxdepth 2 -name '*chr[1-9]*[mapgz]' | sort -k2 -t 'h' | parallel -N3 --sshloginfile $PBS_NODEFILE -j1 \
java -jar /storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
gt=$PBS_O_WORKDIR/{3} \
ref=$PBS_O_WORKDIR/{2} \
map=$PBS_O_WORKDIR/{1} \
out=$PBS_O_WORKDIR/beagle_stuff/{=1 's/..beagle_stuff.plink\.//;s/.GRCh38.map//' =} \
impute=false \

# Without ref
# FUTURE: Check if "phase_out" directory exists, if not, then make
find . -maxdepth 2 -name '*chr[1-9]*[mapgz]' | sort -k2 -t 'h' | parallel -N3 --sshloginfile $PBS_NODEFILE -j1 \
java -jar /storage/home/hcoda1/0/cnaughton7/p-ggibson3-0/software/beagle.19Apr22.7c0.jar \
gt=$PBS_O_WORKDIR/{3} \
map=$PBS_O_WORKDIR/{1} \
out=$PBS_O_WORKDIR/phase_out/{=1 's/..beagle_stuff.plink\.//;s/.GRCh38.map//' =} \
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
# Includes sex chromosomes
# Note: Using file with extra gender column didn't work
# Version1: 
find ref/ALL_POP/ -name *ALL.chr[0-9X]*.gz | parallel -N1 -j1 \
bcftools view \
--threads 4 \
--force-samples \
-S ../query_population_sampleNames_noGender \
-o CEU-YRI{=1 's/ref\/ALL_POP\/ALL\.//;s/.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz//' =} \
{1} #Input

# Version2: For PBS scripts. Spreads jobs over multiple nodes, each assumed to have 8 cores.
# Change later to include variable for conda env path
# TODO: Add output dir path
find ref/ALL_POP/ -name *ALL.chr[0-9X]*.gz | parallel -N1 -j1 --sshloginfile $PBS_NODEFILE \
/storage/home/hcoda1/0/cnaughton7/.conda/envs/ancestry-env1/bin/bcftools view \
--threads 8 \
--force-samples \
-S $PBS_O_WORKDIR/query_population_sampleNames_noGender \
-o $PBS_O_WORKDIR/{=1 's/ref\/ALL_POP\/ALL\.//;s/.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz//' =}".CEU-YRI" \
$PBS_O_WORKDIR/{1} #Input


#################### Intersect VCFs ##########################################################
# Intersect all chromosome vcfs between St. Jude patient data and CEU-YRI 1000 Genomes data
# Batch rename
ls CEU-YRI.ALL.* | xargs -P8 -I{} mv {} "$(echo {} | sed 's/.ALL//')" #Didn't work
ls *.gz  | sed 'p;s/.gz//' | xargs -n2 -P8 mv
ls | sed 'p;s/\.shapeit2\_integrated\_snvindels\_v2a\_27022019\.GRCh38\.phased\.vcf\.gz\_out//' | xargs -n2 -P8 mv # This works

# The loop below didn't work, nor was it intended to run multithreaded.
#for i in $(seq 1 23); do
#	if [[ $i -eq 23 ]]; then
#		i="X"
#	fi
#	out="intersect/chr${i}"
#	mkdir $out
#	find . -maxdepth 2 -name '*chr*.gz' | grep "chr${i}[YM\.]" | xargs -n2 bcftools isec -Oz -p out 
#done

find . -maxdepth 2 -name '*chr*.gz' | sort  -k2 -t 'h' | parallel -N2 \
bcftools isec -Oz -p intersect/{1}_out {1} {2} #This worked

################## Add variant names to vcf #################################################
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' file.vcf


###################### Convert to plink file.###############################################
# NOTE: All samples did not have sex recorded in the output; they came up as "ambiguous".
# Originally run in the beagle directory.
# Used Plink 1.9
# FUTURE: Make edits for use in PBS script
find non-impute/ -maxdepth 3 -name '*.gz' | parallel -j4 plink --vcf {1} --out plink-out/{=1 's/non-impute\///;s/\.vcf\.gz//' =} --make-bed
# The following should make the plink file and set the variant IDs in the process, but the $1 and $2 were not functioning as described in the manual. May need to use Plink2 or rename vcf in prior step.
find non-impute/ -maxdepth 3 -name '*.gz' | parallel -j4 plink --vcf {1} --out plink-out/{=1 's/non-impute\///;s/\.vcf\.gz//' =} --make-bed --set-missing-var-ids @:#:\$1,\$2

# Merge chromosomal plink files for downstream ancestry analysis with 'admixture'
ls plink-out/*[bf]* | xargs -n3 ls > mergedList.txt
plink --merge-list mergedList.txt #Error occured due to variant IDs all being set to "."; went back to the plink file conversion step which can rename the variant IDs...




# R Code For plotting admixture results
tbl=read.table("chr22.2.Q")
newtbl <- tbl[order(tbl$V1),]
barplot(t(as.matrix(newtbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
newtbl2 <- tbl[order(tbl$V2),]
barplot(t(as.matrix(newtbl2)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)