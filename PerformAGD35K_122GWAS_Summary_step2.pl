#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);
use Scalar::Util ();

require "./SharedFunctionList.pl";
my $Population="EUR";#AFR
#my $TraitType="Phecode";#LabValue
my $TraitType="LabValue";
my $SignificantSNVListDir= "AGD35K_".$Population."_$TraitType";
#AGD35kAFR_789.regenie.significant.tsv
#1 169903875 chr1:169903875:A:G A G 0.00973865 24182 ADD 0.619723 0.109247 32.1795 7.85211 NA
  my @SignificantSNVList=get_allfileswithSpecificPostfix($SignificantSNVListDir,
          ".regenie.significant.tsv");
  my $NumberofFiles=scalar(@SignificantSNVList);
  print "NumberofSignifcantSNVFiles:$NumberofFiles\n";
  ##Create the file with trait  amd chr based count of significant SNVs 
my $chr_trait_linecountFile="$SignificantSNVListDir/chr_trait_linecount.txt"; 
count_sig_snvs_by_chr($SignificantSNVListDir, 
     \@SignificantSNVList, 
     $chr_trait_linecountFile);
#Create the trait name amd chr list that we need perform ld_clumping.
my ($ld_traitsRef, 
   $ld_trait_chrRef)=build_ld_clumping_TraitChr_hash($chr_trait_linecountFile);
my %ld_traits=%{$ld_traitsRef};
my %ld_trait_chr=%{$ld_trait_chrRef};
my $Size_ldClumpingTraits=keys %ld_traits;
print "Size_ldClumpingTraits:$Size_ldClumpingTraits\n";

my $GWASsstDir=$SignificantSNVListDir."_GWAS";
my $iTraitbeAnalyzed=0;
my $iTraitChrbeAnalyzed=0;
foreach my $trait (sort keys %ld_trait_chr) {
   $iTraitbeAnalyzed++;	
      #Reads a compressed .gz file
      #Extracts the ID column and LOG10P column
      #Converts LOG10P → P-value
      #Outputs a simple tab-delimited file with ID an d P
       my $InsstFile = "$GWASsstDir/${trait}.regenie.gz";
       #print "InsstFile:$InsstFile\n";
       foreach my $chr (sort keys %{ $ld_trait_chr{$trait} }) {
	next if($ld_trait_chr{$trait}{$chr} == 1);
	my $ConvertedsstFile = "$GWASsstDir/${trait}_chr$chr.txt";       
	convert_log10p_to_p_sstFile($InsstFile, $ConvertedsstFile,$chr);
	 my $SpecificChr=$chr;
         if($chr ==23)
         {
           $SpecificChr="X";
         }
	 #  print "Processing $trait chr$chr\n";
        # define files
        #Independent significant SNPs (r² < 0.6)
	my $outpref_1  = $GWASsstDir."_Clump/${trait}.chr${chr}IndependentSNVs";
 	#AFR_WGS
	$iTraitChrbeAnalyzed++;
	my $genotype_data = "/home/jupyter/AGD35K_GWAS122_2026/AGD35k_$Population";
	$genotype_data .= "/AGD35K_MAC100MR002QCed_$Population"."_chr${SpecificChr}";
	# run plink 
	 system(
                "bash", "-c",
                "/home/jupyter/AGD35K_GWAS122_2026/plink2 ".
                "--pfile $genotype_data ".
                "--chr ${chr} ".
                "--clump $ConvertedsstFile ".
                 "--clump-p1 5e-8 ".
                "--clump-p2 1e-5 ".
                "--clump-r2 0.1 ".
                "--clump-kb 1000 ".
                "--out $outpref_1"
                ) == 0 or warn "PLINK2 failed for $trait chr$chr.ld06\n";  
     }
  }     
print "iTraitbeAnalyzed:$iTraitbeAnalyzed\tiTraitChrbeAnalyzed:$iTraitChrbeAnalyzed\n"; 

#scans a trait–chromosome count hash and 
#  identifies chromosomes that contain exactly one significant SNV 
#   for each trait.
my %Trait_chrListCount1_RawGWAS=collectTraitChrWithSingleSNV(\%ld_trait_chr);

my $ClumpedSNVListDir= $SignificantSNVListDir."_GWAS_Clump";
#CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
#11      5015263 chr11:5015263:G:A       5.04197e-18     402     24      31      80      41      226
my @ClumpedSNVFileList=get_allfileswithSpecificPostfix($ClumpedSNVListDir,
          ".clumps");
 my %Trait_snvid=(); 
 #AGD35kAFR_300.1.chr20IndependentSNVs.clumps
for(my $ii=0;$ii<@ClumpedSNVFileList;$ii++)
{
  my $ClumpedSNPFile="$ClumpedSNVListDir/$ClumpedSNVFileList[$ii]";
  my $ClumpedSNVFileStr=$ClumpedSNVFileList[$ii];
  $ClumpedSNVFileStr =~ s/AGD35kAFR\_//ig;
  $ClumpedSNVFileStr =~ s/AGD35kEUR\_//ig; 
  $ClumpedSNVFileStr =~ s/IndependentSNVs\.clumps//ig;
  my @Trait_Chr=split(".chr",$ClumpedSNVFileStr);
  my $trait=$Trait_Chr[0];
  #print "$ClumpedSNVFileList[$ii]\t$trait\t$chr\n";
  my $Trait_snvidRef=clumpsBasedSNVsCollected(
	  $ClumpedSNPFile,
	  $trait,
	  \%Trait_snvid);
  %Trait_snvid=%{$Trait_snvidRef};
}


my $MergedAFRchr_trait_clumpSNVsFile=$SignificantSNVListDir."_GWAS_ClumpSNVs/";

$MergedAFRchr_trait_clumpSNVsFile.=$Population."_allmergeged_clumps.txt";
#Merge all Clumping SNVs and single SNVs/per Chr into a total file
merge_clumps_snvs_Trait(
     $SignificantSNVListDir,
     \@SignificantSNVList,
     $MergedAFRchr_trait_clumpSNVsFile,
     \%Trait_snvid,
     \%Trait_chrListCount1_RawGWAS);
=cut;
my $EUR_TypeDir="AGD35K_EUR_Phecode";
#my $EUR_TypeDir="AGD35K_EUR_LabValue";
my $EURPhecode_GWASDir= $EUR_TypeDir."_GWAS";
#AGD35kAFR_789.regenie.significant.tsv
#1 169903875 chr1:169903875:A:G A G 0.00973865 24182 ADD 0.619723 0.109247 32.17795 7.85211 NA
 my @EURGWASFileList=get_allfileswithSpecificPostfix($EURPhecode_GWASDir,
         ".regenie.gz");
  my $NumberofGWASFiles=scalar(@EURGWASFileList);
  print "NumberofEUR GWASFiles:$NumberofGWASFiles\n";
   my $MergedEURchr_trait_SNVsFile=$SignificantSNVListDir."_GWAS_ClumpSNVs/EUR_allmergeged_clumps.txt";
   print "write: $MergedEURchr_trait_SNVsFile\n";
   collectEUR_AFRSameSNVEffectFromGWASsstList(
  	 $MergedAFRchr_trait_clumpSNVsFile,
          $EURPhecode_GWASDir,
          \@EURGWASFileList,
          $MergedEURchr_trait_SNVsFile);
=cut;	  
