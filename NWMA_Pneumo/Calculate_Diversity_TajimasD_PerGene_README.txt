
# USAGE
perl Calculate_Diversity_TajimasD_PerGene.pl SubDirectory

where "SubDirectory" is the name of the directory that contains the alignment files, one
for each gene. Ideally, these genes should be aligned with an aligner that used amino 
acid sequences and back translated to nucleotides, or a codon-aware aligner, which help 
preserve the reading frame of each gene. Other aligners, such as MAFFT, that ignore 
reading frames may artificially create frameshift mutations that affect inference of 
whether a mutation is synonymous or nonsynonymous. These gene alignment files should be in
FASTA format.

# OUTPUT
This script outputs four files, each with two columns. Two of these files contain 
Watterson's theta for each gene, one file containing this summary statistic only for 
4-fold degenerate sites, the other file for 0-fold degenerate sites. The other two output
files contain values of Tajima's D for each gene alignment, again for 4-fold and 0-fold
sites separately.