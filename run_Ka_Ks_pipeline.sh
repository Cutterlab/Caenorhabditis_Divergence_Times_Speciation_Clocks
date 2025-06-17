# Load PRANK
export PATH=~/.local/bin/prank/bin/:$PATH

# Define the 2 species between which Ka and Ks will be calculated
species1=$1
species1_full=$2
species2=$3
species2_full=$4


# Go into the directory for this species pair
cd $SCRATCH/All_Pairs_Ka_Ks_Estimates
mkdir ${species1}_${species2}
cd ${species1}_${species2}

# Make the input tree for FitMG94
echo "(${species1},${species2});" > tree.nwk

# Get all the 1:1 orthologs between these 2 species
grep -v -e "," -e Orthogroup $SCRATCH/OrthoFinder/Results_Sep21/Orthologues/Orthologues_${species1_full}/${species1_full}__v__${species2_full}.tsv > orthologs.tsv


# Loop over each of those 1:1 orthologs
while read line
do

first_gene=`echo $line | sed 's/ /\t/g' | cut -f2`
second_gene=`echo $line | sed 's/ /\t/g' | cut -f3`

# Get the CDS for both of these genes and write them to a file
echo $first_gene | ~/programs/seqtk/seqtk subseq ../CDSs_All_Isoforms/${species1}.cds.fa - > unaligned.fa
echo $second_gene | ~/programs/seqtk/seqtk subseq ../CDSs_All_Isoforms/${species2}.cds.fa - >> unaligned.fa

# Run PRANK to align this pair of 1:1 orthologs in codon-alignment mode
prank -d=unaligned.fa -o=aligned.fa -codon

# Only run the subsequent steps if PRANK actually ran for this ortholog pair
if test -f aligned.fa.best.fas
then


# Make a copy of that alignment and add it to a big file that has all the individual pairwise alignments for all the 1:1 orthologs for this species pair
cat aligned.fa.best.fas >> all_alignments.fa

# For use with FitMG94 - rename the sequence names to just be the species names
sed -i "s/$first_gene/$species1/g" aligned.fa.best.fas
sed -i "s/$second_gene/$species2/g" aligned.fa.best.fas

# Replace all in-frame stop codons (all of which should be at the ends of sequences) with gaps
cd $SCRATCH/HyPhy/hyphy-develop/
./HYPHYMP $SCRATCH/HyPhy/hyphy-develop/res/TemplateBatchFiles/CleanStopCodons.bf Universal $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}/aligned.fa.best.fas "No/No" $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}/trimmed_alignment.fa
cd $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}


# If this cleaning step actually ran, run FitMG94 on the alignment and tree
if test -f trimmed_alignment.fa
then

cd $SCRATCH/HyPhy/hyphy-develop/
./HYPHYMP $SCRATCH/HyPhy/hyphy-analyses/FitMG94/FitMG94.bf --rooted No --alignment $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}/trimmed_alignment.fa --tree $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}/tree.nwk --type local --kill-zero-lengths Yes
cd $SCRATCH/All_Pairs_Ka_Ks_Estimates/${species1}_${species2}

# Extract the Ka and Ks values from this run of FitMG94 and write them to a file
Ka=`grep -A8 "${species2}\":{" trimmed_alignment.fa.FITTER.json | tail -n1 | sed 's/"dN"://' | sed 's/,//'`
Ks=`grep -A9 "${species2}\":{" trimmed_alignment.fa.FITTER.json | tail -n1 | sed 's/"dS"://' | sed 's/,//'`
echo $first_gene $second_gene $Ka $Ks >> all_Ka_Ks_values.txt

# Remove the cleaned alignment file (to handle subsequent iterations where the cleaning might not run)
rm trimmed_alignment.fa
fi

# Remove the initial alignment file (to handle subsequent iterations where PRANK might not run)
rm aligned.fa.best.fas
fi

done < orthologs.tsv 
