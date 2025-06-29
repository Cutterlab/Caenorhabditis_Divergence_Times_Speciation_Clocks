# Phylogenomic timetree-calibrated speciation clocks for _Caenorhabditis_ nematodes reveal slow but disproportionate accumulation of post-zygotic reproductive isolation
Supplementary data files of _Caenorhabditis_ divergence time and speciation clock calculations from Fusca et al. (2025). "_C._ sp. 24" (or abbreviations thereof) refers to _C. agridulce_ throughout.
## Description of files
### 205 Single Copy Ortholog Alignments
Directory of DNA coding sequence alignments for the 205 groups of orthologous genes that are single-copy in all 51 species used for our primary divergence time estimates (i.e. not including _C. auriculariae_ or _C. niphades_). Each file contains one sequence for each of the 51 species. These 205 FASTA files are subdivided into four non-overlapping sets, corresponding to the four gene subsets shown in Figure 2D.  Species are referred to using 5-letter abbreviations (e.g. "CELEG" for _C. elegans_).
### 692 Single Copy Orthologs Elegans Supergroup
Directory of DNA coding sequence alignments for the 692 groups of orthologous genes that are single-copy in the 31 species belonging to the Elegans supergroup clade, out of the 51 species that were used for our primary divergence time estimates (i.e. not including _C. niphades_). Each file contains one sequence for each of these 31 species. These 692 FASTA files are subdivided into four non-overlapping sets, corresponding to the four gene subsets shown in Figure 2E.  Species are referred to using 5-letter abbreviations (e.g. "CELEG" for _C. elegans_).
### Dated Species Trees
Directory of dated phylogeny files, in NEXUS format. <b>All branch lengths are in units of generations x 10<sup>7</sup> (i.e. tens of millions of generations).</b> In addition to the tree topology and all branch lengths, these phylogeny files also contain branch labels and node labels produced by BEAST (e.g. 95% HPD intervals for node ages and branch relative substitution rates). Species are referred to using 4-letter abbreviations (e.g. "ELEG" for _C. elegans_).

<ul>
<li><b>205_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: The primary set of divergence time estimates (i.e. those shown in Figure 1), using the 205 single-copy genes and our primary ASTRAL tree topology. </li>
<li><b>partition_1_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 1 of the 205 single-copy genes (red dots in Figure 2D) and our primary ASTRAL tree topology. </li>
<li><b>partition_2_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 2 of the 205 single-copy genes (blue dots in Figure 2D) and our primary ASTRAL tree topology. </li>
<li><b>partition_3_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 3 of the 205 single-copy genes (green dots in Figure 2D) and our primary ASTRAL tree topology. </li>
<li><b>partition_4_52_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 52 genes in Set 4 of the 205 single-copy genes (purple dots in Figure 2D) and our primary ASTRAL tree topology. </li>
<li><b>Elegans_supergroup_old_205_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the Elegans supergroup sequences of the 205 genes that were used in our primary date estimates (red dots in Figure 2E) and the topology of the Elegans supergroup from our primary ASTRAL tree topology. </li>
<li><b>Elegans_supergroup_set_1_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the Elegans supergroup sequences of the 163 genes in Set 1 of the 692 Elegans supergroup single-copy genes (blue dots in Figure 2E) and the topology of the Elegans supergroup from our primary ASTRAL tree topology. </li>
<li><b>Elegans_supergroup_set_2_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the Elegans supergroup sequences of the 162 genes in Set 2 of the 692 Elegans supergroup single-copy genes (green dots in Figure 2E) and the topology of the Elegans supergroup from our primary ASTRAL tree topology. </li>
<li><b>Elegans_supergroup_set_3_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the Elegans supergroup sequences of the 162 genes in Set 3 of the 692 Elegans supergroup single-copy genes (purple dots in Figure 2E) and the topology of the Elegans supergroup from our primary ASTRAL tree topology. </li>
<li><b>alternate_topology_1_partition_3_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 3 of the 205 single-copy genes and alternate tree topology #1 (red dots in Figure 2F). </li>
<li><b>alternate_topology_2_partition_3_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 3 of the 205 single-copy genes and alternate tree topology #2 (blue dots in Figure 2F). </li>
<li><b>alternate_topology_3_partition_3_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 3 of the 205 single-copy genes and alternate tree topology #3 (green dots in Figure 2F). </li>
<li><b>alternate_topology_4_partition_3_51_orthologs_gamma_priors_relaxed_clock_species_tree.tree</b>: Divergence time estimates produced using the 51 genes in Set 3 of the 205 single-copy genes and alternate tree topology #4 (purple dots in Figure 2F). </li>
</ul>

### BEAST XML Input Files
Directory of XML files used as input for runs of BEAST. Each XML input file corresponds to the BEAST run for one of the sets of divergence times listed above.
### run_Ka_Ks_pipeline.sh
Bash script used to calculate K<sub>A</sub> and K<sub>S</sub> for all 1:1 orthologs between a given pair of species. This script makes use of a file of known orthologs between the species pair that must be made beforehand, as well as files of the DNA coding sequences for each species. This script depends on the programs PRANK, Seqtk, and HYPHY. The command-line input to run this script are the abbreviated names and full names for both species being compared; for example, to calculate K<sub>A</sub> and K<sub>S</sub> at 1:1 orthologs between _C. elegans_ and _C. briggsae_, the command is "sh run_Ka_Ks_pipeline.sh CELEG Caenorhabditis_elegans CBRIG Caenorhabditis_briggsae".
### ASTRAL_Pro3_3130_BUSCO_species_tree.nwk
Newick tree file of the ASTRAL-Pro3 species tree for the 51 species used for our primary divergence time estimates (i.e. the tree whose topology was used in the primary divergence time estimates, after rooting with _C. monodelphis_ as the outgroup). Node labels indicate posterior probabilities. Branch lengths are in units of amino acid substitutions per site. 
### 3130_Single_Copy_BUSCO_Gene_Trees.tree
File of 3130 gene trees from IQ-TREE that were used as input to ASTRAL-Pro3 to create ASTRAL_Pro3_3130_BUSCO_species_tree.nwk. Each line corresponds to a separate gene tree. Prior to running ASTRAL-Pro3, all gene names were replaced with the name of the corresponding species. Branch lengths are in units of amino acid substitutions per site.
### niphades_Japonica_group_ASTRAL_tree.nwk
Newick tree file of the ASTRAL-Pro3 species tree for the 13 species in the Japonica group (including _C. niphades_), as well as _C. elegans_ and _C. briggsae_ which were rooted as the outgroup. Node labels indicate posterior probabilities. Branch lengths are in units of amino acid substitutions per site.
### niphades_Japonica_Group_Single_Copy_BUSCO_Gene_Trees.tree
File of 3128 gene trees from IQ-TREE that were used as input to ASTRAL-Pro3 to create niphades_Japonica_group_ASTRAL_tree.nwk. Each line corresponds to a separate gene tree. Prior to running ASTRAL-Pro3, all gene names were replaced with the name of the corresponding species. Branch lengths are in units of amino acid substitutions per site.
### nls_and_terminalstage_data.csv
File containing information used to generate logistic function fits for reproductive isolation and regression analyses for latest viable stage of F1 hybrids.
### Literature_RI_data.csv
File containing citations for reproductive isolation estimates obtained from the literature.
### nls_script.R
R script used to generate logistic function fits for components of reproductive isolation. Includes algorithm used to identify the maximum number of phylogenetically-independent species pairs that could be included in each model. 
### terminalstage_stats.R
R script used to run regression analysis on the latest viable stage of F1 hybrids from interspecies crosses. Includes algorithm used to select the maximum number of phylogenetically independent species for this analysis.
