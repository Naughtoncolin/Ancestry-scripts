# Paint karyograms from RFMix2 output and Tagore
# Currently only processes a single sample at a time.
# Currently assumes only 2 populations
# Author: Colin Naughton

#Get column numbers for sample
sample_name=NWD117836
sample_position=$(tail -n +2 rfmix_out/chr22.msp.tsv | head -n1 | sed 's/\t/\n/g' | \
grep -n $sample_name | head -n1 | sed 's/:.*//g')

# Produce Tagore input file from RFMix2 output
# Subpopulation order/codes: CEU=0; YRI=1q
mkdir tagore_out
echo "Tagore input & output files related to sample ID '"$sample_name"'" > tagore_out/README
# Columns with NWD157077 estimates: 147 & 148 variable to awk command. Need to add +1 for position 2.
#TODO: Add $sample_position
seq 1 22 | xargs -I{} \
sed '/#/d' rfmix_out/chr{}.msp.tsv | \
awk '{
	if ($147 == 1) COLOR1="#0000ff";
	else COLOR1="#F4A500";
	if ($148 == 1) COLOR2="#0000ff";
	else COLOR2="#F4A500";
	print "chr"$1"\t"$2"\t"$3"\t0\t1\t"COLOR1"\t1";
	print "chr"$1"\t"$2"\t"$3"\t0\t1\t"COLOR2"\t2"
	}' > tagore_out/tagore_input.bedc

# Paint karyogram with Tagore
# Issues installing rsvg (required for Tagore) from conda; transferred Tagore input bed from PACE to the Gibson server and ran (it was already installed)
tagore --input tagore_out/tagore_input.bed --prefix tagore_out/tagore_out --build hg38
