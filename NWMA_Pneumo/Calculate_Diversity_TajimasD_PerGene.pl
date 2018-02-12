#!usr/bin/perl
use warnings ;
use strict ;


my $wd_fa = $ARGV[0] ;

unless($wd_fa){
	print "please supply a subdirectory containing one fasta file per gene alignment!\n" ;
	exit ;
}

my @fasta_core = `ls ${wd_fa}/*` ;
my %core_ref_loci_seqs ; # Two-dimensional hash, stores the nucleotide sequences for each individual, for each alignment file (or gene) 
my $num_ind ;
my %Codon_Hash = Construct_Codon_hash () ; # Codon table to infer whether a nucleotide changes results in an amino acid difference
my %Individuals_to_include ;

open OUTPUT, ">./output.txt" ;
print OUTPUT "Total number of genes: ", scalar @fasta_core, "\n" ;

#############################################
# ASSEMBLE SEQUENCES PER GENE PER INDIVIDUAL
#############################################
LoadGeneSeqs(\@fasta_core, \%core_ref_loci_seqs) ;

print OUTPUT "Number of genes with seqs: " ;
print OUTPUT scalar keys %core_ref_loci_seqs, "\n" ;

#############################################
# FIND GENES WITH COHERENT READING FRAMES, to id nonsyn/syn SNPs
#############################################
FindGenesWithCoherentReadingFrames(\%core_ref_loci_seqs) ;

#############################################
# DELETE GENES WITH READING FRAME ERRORS, CLEAN UP %core_ref_loci_seqs
#############################################
DeleteGappedRegions(\%core_ref_loci_seqs) ;

my %Sample_Size_PerGene ;
foreach my $gene (keys %core_ref_loci_seqs){
    $Sample_Size_PerGene{$gene} = scalar keys %{$core_ref_loci_seqs{$gene}} ;
}

#############################################
## GET CONSENSUS BASE PER POSITION FOR CALCULATING DEGENERACY
#############################################
my %Consensus_Base ; # $Consensus_Base{gene}{site} = base
my %Degeneracy ; # $Degeneracy{gene}{site} = n-fold degeneracy
GetConsensusBase_DegeneracyPerSite(\%core_ref_loci_seqs, \%Codon_Hash, \%Consensus_Base, \%Degeneracy) ;

#############################################
## COLLECT SEGREGATING SITE POSITIONS
#############################################
my %core_ref_segsites ;
my %Functional_Effect ;
my %multiallelic_sites ;
my %multiallelic_sites_bases ;
my %biallelic_sites_freqs  ;
CollectSegSitePositions(\%core_ref_loci_seqs, \%Consensus_Base, \%Codon_Hash, \%core_ref_segsites, \%Functional_Effect, \%multiallelic_sites, \%multiallelic_sites_bases, \%biallelic_sites_freqs) ;
#Look at proportion of polymorphic sites 
my %Proportion_Variable_Sites ;
foreach my $gene (sort{$a cmp $b} keys %Degeneracy){
    foreach my $site (sort{$a <=> $b} keys %{$Degeneracy{$gene}}){
        $Proportion_Variable_Sites{"total"}{ $Degeneracy{$gene}{$site} } ++ ;
        if( exists $Functional_Effect{$gene}{$site} ){
            $Proportion_Variable_Sites{"biallelic"}{ $Degeneracy{$gene}{$site} } ++ ;
        }elsif(exists $multiallelic_sites{$gene}{$site}){
            $Proportion_Variable_Sites{"multiallelic"}{ $Degeneracy{$gene}{$site} } ++ ;
        }
    }
}
print OUTPUT "Degeneracy\tTotalNumSites\tBiAllelicDiversity\tMultiAllelicDiversity\n" ;
foreach my $key (sort{$a cmp $b} keys %{$Proportion_Variable_Sites{"total"}} ){
    print OUTPUT $key, "\t", $Proportion_Variable_Sites{"total"}{$key}, "\t", $Proportion_Variable_Sites{"biallelic"}{$key}/$Proportion_Variable_Sites{"total"}{$key}, "\t", $Proportion_Variable_Sites{"multiallelic"}{$key}/$Proportion_Variable_Sites{"total"}{$key},"\n" ;
}
#############################################
## CALCULATE DIVERSITY/Tajima's D
#############################################
my %SFS_4D ;
my %SFS_0D ;
CalculateDiversityPerGene(\%Degeneracy, \%Sample_Size_PerGene, \%biallelic_sites_freqs, \%multiallelic_sites_bases, \%Functional_Effect ) ;

	
exit ;




###################
# SUBROUTINES
###################


sub LoadGeneSeqs{

	my $fasta_genes= $_[0] ;
	my $core_ref_loci_seqs = $_[1] ;

	foreach my $file (@$fasta_genes){
		chomp $file ;
		my $gene = $file ;
		open IN, "<${file}" ;
		my $id ;
		my $ind ;
		while(<IN>){
			chomp $_ ;
			if($_ =~ m/^>/){
				$_ =~ s/>//g ;
				$ind = $_ ;
				my $x = <IN> ; chomp $x ;
				$core_ref_loci_seqs{$gene}{ $ind } = uc($x) ; #uppercase all DNA for consistency
				next ;
			}else{
				$_ = uc($_) ;
				$core_ref_loci_seqs{$gene}{ $ind } = $core_ref_loci_seqs{$gene}{ $ind }.${_} ;
			}
		}
		close IN ;
	}
}

sub FindGenesWithCoherentReadingFrames{

	my $core_ref_loci_seqs= $_[0] ;

	my %Start_codons ; # to make sure ROARY output genes in 5'-3'
	my %Gene_ReadingFrame_Info ;
	my %Genes_with_coherent_reading_frames ;
	my %PerInd_ReadFrame_Viol ;
	foreach my $gene (sort {$a cmp $b} keys %$core_ref_loci_seqs){
		foreach my $ind (sort {$a cmp $b} keys %{$$core_ref_loci_seqs{$gene}}){
			for (my $site = 0; $site <=length($$core_ref_loci_seqs{$gene}{$ind})-3; $site +=3 ){
				my $triplet = substr($$core_ref_loci_seqs{$gene}{$ind}, $site, 3) ;
				if( ($triplet =~ tr/[AGCTN]//) == 3 || ($triplet =~ tr/-//) == 3){ # check if reading frame intact
					next ;
				}else{ # there is some sort of frame-shift mutations (or alignment error) for this individual
					# This builds a reading-frame (RF) error profile per individual
					if( exists $PerInd_ReadFrame_Viol{$gene}{$ind} ){
						# sites exhibiting RF error are concatenated to make profile
						$PerInd_ReadFrame_Viol{$gene}{$ind} = $PerInd_ReadFrame_Viol{$gene}{$ind}."\.site_".$site ;
					}else{
						$PerInd_ReadFrame_Viol{$gene}{$ind} = "site_".$site ;
					}
				}
			}
		}
		## Double check sequence of start codons
		(my $test) = sort keys %{$$core_ref_loci_seqs{$gene}} ;
		my $dna = $$core_ref_loci_seqs{$gene}{$test} ;
		$dna =~ s/-//g ;
		$Start_codons{substr($dna, 0, 3)} ++ ;
		# Flag genes that contained frameshifts; reading frame cannot be accurately reconstructed
		if( !exists $PerInd_ReadFrame_Viol{$gene} ){
			$Genes_with_coherent_reading_frames{"coherent"}{$gene} ++ ;
		}else{
			$Genes_with_coherent_reading_frames{"incoherent"}{$gene} ++ ;
		}

	}
	open OUT, ">>./output.txt" ;
	print OUT "Checking first 3 bp of each alignment, which should be a recognizable start codon! Otherwise your DNA may be reverse-complemented and should be changed\n" ;
	print OUT "START CODONS:\n" ;
	foreach my $key (sort {$a cmp $b} keys %Start_codons){
		print OUT "\t", $key, "\t", $Start_codons{$key}, "\n" ;
	}
	foreach my $gene (sort {$a cmp $b} keys %$core_ref_loci_seqs){
		if(exists $Genes_with_coherent_reading_frames{"incoherent"}{$gene} ){
			print OUT "deleting ${gene}...\n" ;
			delete $$core_ref_loci_seqs{$gene} ;
		}
	}
	my $sum1 = 0;
	foreach my $gene (keys %$core_ref_loci_seqs){
		(my $first) = sort keys %{$$core_ref_loci_seqs{$gene}} ;
		if(exists $Genes_with_coherent_reading_frames{"coherent"}{$gene} ){
			$sum1 +=  length($$core_ref_loci_seqs{$gene}{$first}) ;
		}
	}
	 print OUT "Genes with coherent reading frames: ", scalar keys %{$Genes_with_coherent_reading_frames{"coherent"}}, "\n"  ;
	 print OUT "Genes with incoherent reading frames: ", scalar keys %{$Genes_with_coherent_reading_frames{"incoherent"}}, "\n"  ;
	 print OUT "Total coherent dna: ", $sum1, "\n" ;

	close OUT ;
}

sub DeleteGappedRegions{

	my $core_ref_loci_seqs= $_[0] ;

	#Remove sites with -'s; according to the filtering criteria above, these indels don't disrupt the reading frame
	foreach my $gene (keys %$core_ref_loci_seqs){
		my @sites_to_delete ;
		(my $first) = sort keys %{$$core_ref_loci_seqs{$gene}} ;
		# go thru ea ind, find each site, then delete, has to be 2 step process o.w. you modify string length in process
		SITELOOP: foreach my $site ( 0 .. length($$core_ref_loci_seqs{$gene}{$first})-1 ){
			foreach my $ind (keys %{$$core_ref_loci_seqs{$gene}}){
				my $base = substr($$core_ref_loci_seqs{$gene}{$ind}, $site, 1) ;
				if($base eq "-"){ # This finds -'s that occur in at least one individual, go to next site if found
					push @sites_to_delete, $site ;
					next SITELOOP ; 
				}
			}
		}
		foreach my $ind (keys %{$$core_ref_loci_seqs{$gene}}){
			my $counter = 0 ;
			foreach my $site (@sites_to_delete){
				# ($site-$counter) to account for fact that you modify string length and don't go overbounds
				substr($$core_ref_loci_seqs{$gene}{$ind}, ($site-$counter), 1) = "" ;
				$counter++ ;
			}
		}
	}
}

sub GetConsensusBase_DegeneracyPerSite{

	my $core_ref_loci_seqs= $_[0] ;
	my $Codon_Hash = $_[1] ;
	my $Consensus_Base = $_[2] ;
	my $Degeneracy = $_[3] ;

	foreach my $gene (keys %$core_ref_loci_seqs){
		(my $first) = sort keys %{$$core_ref_loci_seqs{$gene}} ;
		# go through each site and find consensus/majority base
		foreach my $site (0 .. length($$core_ref_loci_seqs{$gene}{$first})-1){ ## SOME INDIVIDUALS MAY HAVE DIFF LENGTHS
			my %bases ;
			foreach my $ind (keys %{$$core_ref_loci_seqs{$gene}}){
				my $base = substr($$core_ref_loci_seqs{$gene}{$ind}, $site, 1) ;
				$bases{$base}++ ;
			}
			(my $consensus) = sort {$bases{$b} <=> $bases{$a}} keys %bases ;
			$$Consensus_Base{$gene}{$site} = uc($consensus) ; #capitalize: codon table all caps
		}
	
		#go back through each site and calculate $Degeneracy with consensus bases
		foreach my $site (0 .. length($$core_ref_loci_seqs{$gene}{$first})-1){
			my $codon_position = ($site+1)%3 ; #add 1 to $site b/c starts at 0
			my $flank_site_1 ;
			my $flank_site_2 ;
			# Get positions of flanking sites in current codon; used to determine degeneracy of site
			if($codon_position == 1){ #site in 1st position
				$flank_site_1 = $site+1 ;
				$flank_site_2 = $site+2 ;
			}elsif($codon_position == 2){ #site in 2nd position
				$flank_site_1 = $site-1 ;
				$flank_site_2 = $site+1 ;
			}elsif($codon_position == 0){ #site in 3rd position
				$flank_site_1 = $site-2 ;
				$flank_site_2 = $site-1 ;
			}
			my $fb1_consensus = $$Consensus_Base{$gene}{$flank_site_1} ;
			my $fb2_consensus = $$Consensus_Base{$gene}{$flank_site_2} ;
			
			if($fb1_consensus !~ m/[Nn]/ && $fb2_consensus !~ m/[Nn]/){
				my %codon_degeneracy ; 
				foreach my $base ( "A", "G", "C", "T" ){
					if($codon_position == 1){ #current site is in 1st position
						$codon_degeneracy{ $$Codon_Hash{$base.$fb1_consensus.$fb2_consensus} } ++ ;
					}elsif($codon_position == 2){ #current site is in 2nd position
						$codon_degeneracy{ $$Codon_Hash{$fb1_consensus.$base.$fb2_consensus} } ++ ;
					}elsif($codon_position == 0){ #current site is in 3rd position
						$codon_degeneracy{ $$Codon_Hash{$fb1_consensus.$fb2_consensus.$base} } ++ ;
					}
				}
				if((scalar keys %codon_degeneracy) == 1){
					$$Degeneracy{$gene}{$site} = "4D" ; #n-fold degeneracy -> n nucleotides specify same aa
				}elsif((scalar keys %codon_degeneracy) == 2){
					(my $tmp) = sort {$codon_degeneracy{$a} <=> $codon_degeneracy{$b}} keys %codon_degeneracy ; #take lower freq aa
					if( $codon_degeneracy{$tmp} == 2){ # if 2 aa's equally frequent (2/4 ea.), 2D b/c 1 of 3 synonymous
						$$Degeneracy{$gene}{$site} = "2D" ;
					}elsif($codon_degeneracy{$tmp} == 1){
						$$Degeneracy{$gene}{$site} = "3D" ;
					}
				}elsif((scalar keys %codon_degeneracy) == 3){ # always 2D
					$$Degeneracy{$gene}{$site} = "2D" ;
				}elsif((scalar keys %codon_degeneracy) == 4){
					$$Degeneracy{$gene}{$site} = "0D" ;
				}
			}
		}
	}
}

sub CollectSegSitePositions{

	my $core_ref_loci_seqs= $_[0] ;
	my $Consensus_Base = $_[1] ;
	my $Codon_Hash = $_[2] ;
	my $core_ref_segsites = $_[3] ;
	my $Functional_Effect = $_[4] ;
	my $multiallelic_sites = $_[5] ;
	my $multiallelic_sites_bases = $_[6] ;
	my $biallelic_sites_freqs = $_[7] ;
	
	my $segsite_counter = 0 ;
	foreach my $gene (keys %$core_ref_loci_seqs){
		my $seg_site_index = 0 ;
		(my $first) = sort keys %{$$core_ref_loci_seqs{$gene}} ;
		# go thru ea ind, for each site
		foreach my $site (0 .. length($$core_ref_loci_seqs{$gene}{$first})-1){
			my %poly_bases ;
			my $indel = 0 ;
			foreach my $ind (keys %{$$core_ref_loci_seqs{$gene}}){
				my $base = substr($$core_ref_loci_seqs{$gene}{$ind}, $site, 1) ;
				if($base !~ m/-/ && $base !~ m/[Nn]/){
					$base = uc($base) ; #capitalize: codon table all caps
					$poly_bases{$base} ++ ;
				}else{
					$indel ++ ;
				}
			}
			if( (scalar keys %poly_bases >= 3) ){
				$$multiallelic_sites{$gene}{$site}++ ;
				foreach my $base ( keys %poly_bases ){
					push @{$$multiallelic_sites_bases{$gene}{$site}}, $base ;
				}
			}
			if( (scalar keys %poly_bases == 2) && ($indel==0) ){
				$$core_ref_segsites{$gene}{$seg_site_index} = $site  ;
				$seg_site_index++ ;
				$segsite_counter++ ;
				# find codon position and positions of other sites in same codon
				my $codon_position = ($site+1)%3 ; #add 1 to $site b/c starts at 0
				my $flank_site_1 ;
				my $flank_site_2 ;
				if($codon_position == 1){ #site in 1st position
					$flank_site_1 = $site+1 ;
					$flank_site_2 = $site+2 ;
				}elsif($codon_position == 2){ #site in 2nd position
					$flank_site_1 = $site-1 ;
					$flank_site_2 = $site+1 ;
				}elsif($codon_position == 0){ #site in 3rd position
					$flank_site_1 = $site-2 ;
					$flank_site_2 = $site-1 ;
				}
			
				my $fb1_consensus = $$Consensus_Base{$gene}{$flank_site_1} ;
				my $fb2_consensus = $$Consensus_Base{$gene}{$flank_site_2} ;
				#construct biallelic codons, only 2 according to if(condition) above
				my @biallelic_codon ;
				foreach my $base ( keys %poly_bases ){
					$$biallelic_sites_freqs{$gene}{$site}{$base} = $poly_bases{$base} ;
					if($codon_position == 1){ #site in 1st position
						push @biallelic_codon, $base.$fb1_consensus.$fb2_consensus ;
					}elsif($codon_position == 2){ #site in 2nd position
						push @biallelic_codon, $fb1_consensus.$base.$fb2_consensus ;
					}elsif($codon_position == 0){ #site in 3rd position
						push @biallelic_codon, $fb1_consensus.$fb2_consensus.$base ;
					}
				}
			
				#################
				if(!exists $$Codon_Hash{$biallelic_codon[0]} || !exists $$Codon_Hash{$biallelic_codon[1]}){
					print "DOESNT EXIST: ", "@biallelic_codon", , "\t", "in gene: ", $gene, " ", $site+1, "\n" ;
				}else{
					if( $$Codon_Hash{$biallelic_codon[0]} eq $$Codon_Hash{$biallelic_codon[1]} ){
						$$Functional_Effect{$gene}{$site} = "SYNONYMOUS" ;
					}else{
						$$Functional_Effect{$gene}{$site} = "NONSYNONYMOUS" ;
					}
				}
			}
		}
	}	

	open OUT, ">>./output.txt" ;
	print OUT "Found ${segsite_counter} biallelic segregating sites\n" ;
	close OUT ;
}

sub CalculateDiversityPerGene{

	my $Degeneracy = $_[0] ;
	my $Sample_Size_PerGene = $_[1] ;
	my $biallelic_sites_freqs = $_[2] ;
	my $multiallelic_sites_bases = $_[3] ;
	my $Functional_Effect = $_[4] ;

	my %Wattersons_Theta_4D ;
	my %Pi_Theta_4D ;
	my %Lengths_4D ;
	my $Total_4D = 0 ;
	my %SFS_4D ;
	my %Wattersons_Theta_0D ;
	my %Pi_Theta_0D ;
	my %Lengths_0D ;
	my $Total_0D = 0 ;
	my %SFS_0D ;
	
	
	my %SFS_Synonymous ;
	
	foreach my $gene (sort{$a cmp $b} keys %$Degeneracy){
		$Wattersons_Theta_4D{"biallelic"}{$gene} = 0 ; # Initialize Wattersons theta ONLY for biallelic sites
		$Wattersons_Theta_4D{"multiallelic"}{$gene} = 0 ; # Initialize Wattersons theta ONLY for multiallelic sites
		$Wattersons_Theta_4D{"finitesites"}{$gene} = 0 ; # Initialize Wattersons theta that incorporates both biallelic and multiallelic sites
		$Pi_Theta_4D{"biallelic"}{$gene} = 0 ; # Initialize Pi theta only for biallelic sites

		$Wattersons_Theta_0D{"biallelic"}{$gene} = 0 ; # Initialize Wattersons theta ONLY for biallelic sites
		$Wattersons_Theta_0D{"multiallelic"}{$gene} = 0 ; # Initialize Wattersons theta ONLY for multiallelic sites
		$Wattersons_Theta_0D{"finitesites"}{$gene} = 0 ; # Initialize Wattersons theta that incorporates both biallelic and multiallelic sites
		$Pi_Theta_0D{"biallelic"}{$gene} = 0 ; # Initialize Pi theta only for biallelic sites
		
		my $n = $$Sample_Size_PerGene{$gene} ;
		## CALCULATE DIVERSITY FOR 4D sites
		foreach my $site (sort{$a <=> $b} keys %{$$Degeneracy{$gene}}){
			if($$Degeneracy{$gene}{$site} eq "4D"){
				$Lengths_4D{$gene}++ ;
				$Total_4D ++ ;
				if( exists $$biallelic_sites_freqs{$gene}{$site} ){
					$Wattersons_Theta_4D{"biallelic"}{$gene}++ ;
					$Wattersons_Theta_4D{"finitesites"}{$gene}++ ;
					(my $first_base) = sort {$a cmp $b} keys %{$$biallelic_sites_freqs{$gene}{$site}} ;
					my $af = $$biallelic_sites_freqs{$gene}{$site}{$first_base} ;
					$Pi_Theta_4D{"biallelic"}{$gene} += ($af*($n-$af))/bc($n,2) ;
					$SFS_4D{$gene}{$site} = $af ;
				}else{
					$SFS_4D{$gene}{$site} = 0 ;
				}
			
				if( exists $multiallelic_sites{$gene}{$site} ){
					$Wattersons_Theta_4D{"multiallelic"}{$gene}++ ;
					$Wattersons_Theta_4D{"finitesites"}{$gene}+= ((scalar @{$$multiallelic_sites_bases{$gene}{$site}}) - 1) ;
				}
			}
			if(exists $$Functional_Effect{$gene}{$site}){
				if( $$Functional_Effect{$gene}{$site} eq "SYNONYMOUS" ){
					if( exists $$biallelic_sites_freqs{$gene}{$site} ){
						(my $first_base) = sort {$a cmp $b} keys %{$$biallelic_sites_freqs{$gene}{$site}} ;
						my $af = $$biallelic_sites_freqs{$gene}{$site}{$first_base} ;
						$SFS_Synonymous{$gene}{$site} = $af ;
					}else{
						$SFS_Synonymous{$gene}{$site} = 0 ;
					}        
				}
			}
			if($$Degeneracy{$gene}{$site} eq "0D"){
				$Lengths_0D{$gene}++ ;
				$Total_0D ++ ;
				if( exists $$biallelic_sites_freqs{$gene}{$site} ){
					$Wattersons_Theta_0D{"biallelic"}{$gene}++ ;
					$Wattersons_Theta_0D{"finitesites"}{$gene}++ ;
					(my $first_base) = sort {$a cmp $b} keys %{$$biallelic_sites_freqs{$gene}{$site}} ;
					my $af = $$biallelic_sites_freqs{$gene}{$site}{$first_base} ;
					$Pi_Theta_0D{"biallelic"}{$gene} += ($af*($n-$af))/bc($n,2) ;
					$SFS_0D{$gene}{$site} = $af ;
				}else{
					$SFS_0D{$gene}{$site} = 0 ;
				}
			
				if( exists $multiallelic_sites{$gene}{$site} ){
					$Wattersons_Theta_0D{"multiallelic"}{$gene}++ ;
					$Wattersons_Theta_0D{"finitesites"}{$gene}+= ((scalar @{$$multiallelic_sites_bases{$gene}{$site}}) - 1) ;
				}
			}
		}
		# finite sites estimator uses proportion of seg sites
		$Wattersons_Theta_4D{"finitesites"}{$gene} = $Wattersons_Theta_4D{"finitesites"}{$gene}/$Lengths_4D{$gene} ;
		$Wattersons_Theta_0D{"finitesites"}{$gene} = $Wattersons_Theta_0D{"finitesites"}{$gene}/$Lengths_0D{$gene} ;
		# constants for wattersons estimators, see Tajima 1996, genetics, pp 1457-1465
		my $a1 = 0 ;
		foreach ( 1 .. ($n-1) ){ # sample size may vary depending on $max_miss_data
			$a1 += (1/$_) ;
		}
		my $a1_square = 0 ;
		foreach ( 1 .. ($n-1) ){
			$a1_square += (1/$_**2) ;
		}
		my $a2 = (($a1**2) - $a1_square)/2 ;
		my $c2 = (4*$a1/3) - (7*$a2/(3*$a1)) ;
		## Wattersons theta per gene
		$Wattersons_Theta_4D{"biallelic"}{$gene} = $Wattersons_Theta_4D{"biallelic"}{$gene}/$a1 ;
		$Wattersons_Theta_4D{"multiallelic"}{$gene} = $Wattersons_Theta_4D{"multiallelic"}{$gene}/$a1 ;
		$Wattersons_Theta_4D{"finitesites"}{$gene} = $Wattersons_Theta_4D{"finitesites"}{$gene}/($a1 - ($c2*$Wattersons_Theta_4D{"finitesites"}{$gene})) ;
		
		$Wattersons_Theta_0D{"biallelic"}{$gene} = $Wattersons_Theta_0D{"biallelic"}{$gene}/$a1 ;
		$Wattersons_Theta_0D{"multiallelic"}{$gene} = $Wattersons_Theta_0D{"multiallelic"}{$gene}/$a1 ;
		$Wattersons_Theta_0D{"finitesites"}{$gene} = $Wattersons_Theta_0D{"finitesites"}{$gene}/($a1 - ($c2*$Wattersons_Theta_0D{"finitesites"}{$gene})) ;
		
	}
	## Normalize wattersons theta by gene length -> per nucleotide
	my $watt_avg_bi = 0;
	my $watt_avg_multi = 0 ;
	my $watt_avg_finite = 0 ;
	my $pi_avg_bi = 0 ;
	foreach my $gene (sort{$a cmp $b} keys %$Degeneracy){
		# note finite sites estimator already divided by gene length above
		$Wattersons_Theta_4D{"biallelic"}{$gene} = $Wattersons_Theta_4D{"biallelic"}{$gene}/($Lengths_4D{$gene}) ;
		$Wattersons_Theta_4D{"multiallelic"}{$gene} = $Wattersons_Theta_4D{"multiallelic"}{$gene}/($Lengths_4D{$gene}) ;
		$Pi_Theta_4D{"biallelic"}{$gene} = $Pi_Theta_4D{"biallelic"}{$gene}/($Lengths_4D{$gene}) ;

		$Wattersons_Theta_0D{"biallelic"}{$gene} = $Wattersons_Theta_0D{"biallelic"}{$gene}/($Lengths_0D{$gene}) ;
		$Wattersons_Theta_0D{"multiallelic"}{$gene} = $Wattersons_Theta_0D{"multiallelic"}{$gene}/($Lengths_0D{$gene}) ;
		$Pi_Theta_0D{"biallelic"}{$gene} = $Pi_Theta_0D{"biallelic"}{$gene}/($Lengths_0D{$gene}) ;
		###
		$watt_avg_bi += $Wattersons_Theta_4D{"biallelic"}{$gene}*($Lengths_4D{$gene}/$Total_4D) ;
		$watt_avg_multi += $Wattersons_Theta_4D{"multiallelic"}{$gene}*($Lengths_4D{$gene}/$Total_4D) ;
		$watt_avg_finite += $Wattersons_Theta_4D{"finitesites"}{$gene}*($Lengths_4D{$gene}/$Total_4D) ;
		$pi_avg_bi += $Pi_Theta_4D{"biallelic"}{$gene}*($Lengths_4D{$gene}/$Total_4D) ;
	}
	open OUT, ">>./output.txt" ;
	print OUT "Diversity estimates across all genes:\n " ;
	print OUT "\tWattersons_Theta_BiallelicSites: ",$watt_avg_bi, "\n" ;
	print OUT "\tWattersons_Theta_FiniteSitesModel: ", $watt_avg_finite, "\n" ;
	print OUT "\tPi_Theta_BiallelicSites: ",$pi_avg_bi, "\n" ;
	close OUT ;

	open OUT, ">./Wattersons_Theta_PerGene_4Dsites.txt" ;
	print OUT "GeneName\tTheta\n" ;
	foreach my $gene (sort {$a cmp $b} keys %{$Wattersons_Theta_4D{"finitesites"}} ){
		print OUT $gene, "\t", $Wattersons_Theta_4D{"finitesites"}{$gene}, "\n" ;
	}	
	close OUT ;
	open OUT, ">./Wattersons_Theta_PerGene_0Dsites.txt" ;
	print OUT "GeneName\tTheta\n" ;
	foreach my $gene (sort {$a cmp $b} keys %{$Wattersons_Theta_0D{"finitesites"}} ){
		print OUT $gene, "\t", $Wattersons_Theta_0D{"finitesites"}{$gene}, "\n" ;
	}	
	close OUT ;
	
	
	
	my %SumStat_Results_4D ;
	my %SumStat_Results_0D ;

	foreach my $gene (keys %SFS_4D){
		$SumStat_Results_4D{$gene}{"PI"} = theta_pi( \%{$SFS_4D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		$SumStat_Results_4D{$gene}{"W"} = theta_W( \%{$SFS_4D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		$SumStat_Results_4D{$gene}{"S"} = seg_sites( \%{$SFS_4D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		if( $SumStat_Results_4D{$gene}{"S"} > 2 ){
			$SumStat_Results_4D{$gene}{"TajD"} = Taj_D( $SumStat_Results_4D{$gene}{"PI"}, $SumStat_Results_4D{$gene}{"W"}, $SumStat_Results_4D{$gene}{"S"}, $Sample_Size_PerGene{$gene} ) ;
		}else{
			$SumStat_Results_4D{$gene}{"TajD"} = "NA" ;
		}
	}
	foreach my $gene (keys %SFS_0D){
		$SumStat_Results_0D{$gene}{"PI"} = theta_pi( \%{$SFS_0D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		$SumStat_Results_0D{$gene}{"W"} = theta_W( \%{$SFS_0D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		$SumStat_Results_0D{$gene}{"S"} = seg_sites( \%{$SFS_0D{$gene}}, $Sample_Size_PerGene{$gene} ) ;
		if( $SumStat_Results_0D{$gene}{"S"} > 2 ){
			$SumStat_Results_0D{$gene}{"TajD"} = Taj_D( $SumStat_Results_0D{$gene}{"PI"}, $SumStat_Results_0D{$gene}{"W"}, $SumStat_Results_0D{$gene}{"S"}, $Sample_Size_PerGene{$gene} ) ;
		}else{
			$SumStat_Results_0D{$gene}{"TajD"} = "NA" ;
		}
	}
	
	open OUT, ">./TajD_PerGene_4Dsites.txt" ;
	print OUT "Gene\tTajD\n" ;
	foreach my $gene(keys %SumStat_Results_4D ){
		print OUT $gene, "\t", $SumStat_Results_4D{$gene}{"TajD"}, "\n" ;
	}
	close OUT ;
	open OUT, ">./TajD_PerGene_0Dsites.txt" ;
	print OUT "Gene\tTajD\n" ;
	foreach my $gene(keys %SumStat_Results_0D ){
		print OUT $gene, "\t", $SumStat_Results_0D{$gene}{"TajD"}, "\n" ;
	}
	close OUT ;
		
}






sub bc {
    my ($n,$k) = @_;
    my $r=1;
    $r*=$n/($n-$k),$n--while$n>$k;
    $r;
}

sub Construct_Codon_hash{
    my %hash ;
    $hash{"ATT"} = "I" ;
    $hash{"ATC"} = "I" ;
    $hash{"ATA"} = "I" ;
    
    $hash{"CTT"} = "L" ;
    $hash{"CTC"} = "L" ;
    $hash{"CTA"} = "L" ;
    $hash{"CTG"} = "L" ;
    $hash{"TTA"} = "L" ;
    $hash{"TTG"} = "L" ;
    
    $hash{"GTT"} = "V" ;
    $hash{"GTC"} = "V" ;
    $hash{"GTA"} = "V" ;
    $hash{"GTG"} = "V" ;
    
    $hash{"TTT"} = "F" ;
    $hash{"TTC"} = "F" ;
    
    $hash{"ATG"} = "M" ;
    
    $hash{"TGT"} = "C" ;
    $hash{"TGC"} = "C" ;
    
    $hash{"GCT"} = "A" ;
    $hash{"GCC"} = "A" ;
    $hash{"GCA"} = "A" ;
    $hash{"GCG"} = "A" ;
    
    $hash{"GGT"} = "G" ;
    $hash{"GGC"} = "G" ;
    $hash{"GGA"} = "G" ;
    $hash{"GGG"} = "G" ;
    
    $hash{"CCT"} = "P" ;
    $hash{"CCC"} = "P" ;
    $hash{"CCA"} = "P" ;
    $hash{"CCG"} = "P" ;
    
    $hash{"ACT"} = "T" ;
    $hash{"ACC"} = "T" ;
    $hash{"ACA"} = "T" ;
    $hash{"ACG"} = "T" ;
    
    $hash{"TCT"} = "S" ;
    $hash{"TCC"} = "S" ;
    $hash{"TCA"} = "S" ;
    $hash{"TCG"} = "S" ;
    $hash{"AGT"} = "S" ;
    $hash{"AGC"} = "S" ;
    
    $hash{"TAT"} = "Y" ;
    $hash{"TAC"} = "Y" ;
    
    $hash{"TGG"} = "W" ;
    
    $hash{"CAA"} = "Q" ;
    $hash{"CAG"} = "Q" ;
    
    $hash{"AAT"} = "N" ;
    $hash{"AAC"} = "N" ;
    
    $hash{"CAT"} = "H" ;
    $hash{"CAC"} = "H" ;
    
    $hash{"GAA"} = "E" ;
    $hash{"GAG"} = "E" ;
    
    $hash{"GAT"} = "D" ;
    $hash{"GAC"} = "D" ;
    
    $hash{"AAA"} = "K" ;
    $hash{"AAG"} = "K" ;
    
    $hash{"CGT"} = "R" ;
    $hash{"CGC"} = "R" ;
    $hash{"CGA"} = "R" ;
    $hash{"CGG"} = "R" ;
    $hash{"AGA"} = "R" ;
    $hash{"AGG"} = "R" ;
    
    $hash{"TAA"} = "Stop" ;
    $hash{"TAG"} = "Stop" ;
    $hash{"TGA"} = "Stop" ;
    
    return (%hash) ;
}

sub seg_sites{
	#returns seg sites, #singletons
	my ($AFS, $n) = @_ ;
	my $S = 0 ;
	my $S_singleton = 0 ;
	foreach my $pos (keys ( %{$AFS} )){
		if(($$AFS{$pos} > 0) && ($$AFS{$pos} < $n)){
			$S += 1 ;
		}
		if($$AFS{$pos} == 1){
			$S_singleton ++ ;
		}
	}
	return($S) ;
}

sub theta_pi{
	#usage: theta_H(\%AFS, $n), where AFS is $AFS{scaff}{pos}
	my ($AFS, $n) = @_ ;
	my $theta_pi = 0 ;
	foreach my $pos (keys ( %{$AFS} )){
		if(($$AFS{$pos} > 0) && ($$AFS{$pos} < $n)){
			$theta_pi += (($$AFS{$pos})*($n-$$AFS{$pos}))/bc($n,2) ;
		}
	}
	return($theta_pi) ;
}

sub theta_W{
	my ($AFS, $n) = @_ ;
	my $theta_W = 0 ;
	my $a = 0 ;
	foreach ( 1 .. ($n - 1) ) { 
		$a = $a + 1/$_ ; 
	}
	foreach my $pos (keys ( %{$AFS} )){
		if(($$AFS{$pos} > 0) && ($$AFS{$pos} < $n)){
			$theta_W += 1/$a ;
		}
	}
	return($theta_W) ;
}

sub Taj_D{
	my($theta_pi, $theta_W, $S, $n) = @_ ;
	my $tajD = 0 ;
	my $a1 = 0 ; 
	my $a2 = 0 ;
	foreach ( 1 .. ($n - 1) ) { 
		$a1 = $a1 + 1/$_ ; 
		$a2 = $a2 + 1/($_)**2 ;
	}
	my $e1 = (1/$a1)*(($n+1)/(3*($n-1)) - (1/$a1)) ;
	my $e2 = (1/($a1**2 + $a2))*((2*(($n**2) + $n + 3))/(9*$n*($n-1)) - ($n+2)/($n*$a1) + ($a2)/($a1**2)) ;
	if ( (sqrt($S*$e1 + $S*($S-1)*$e2)) > 0 ) { 
		$tajD = ($theta_pi - $theta_W) / (sqrt($S*$e1 + $S*($S-1)*$e2)) ; 
		return ($tajD) ;
	}else{
		return ("NA") ;
	}
}

