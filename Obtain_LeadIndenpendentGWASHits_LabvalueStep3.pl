#Merge AGD35K AFR all related data into an output file
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);

#Trait CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTR
#BASOAB 1 157834858 chr1:157834858:G:A G A 0.0746193 7880 ADD 0.00349444 0.000598392 34.1024 8.28161 NA
#A1FREQ	N	BETA	SE	CHISQ	LOG10P
sub ReadTraitBasedSignificantVariantInfo
{
  my ($SignificantSNVsFile)=@_;	
  my %Trait_chrPos_SNVID;
  my %Trait_chrPos_RefAlt;
  my %Trait_chrPos_BetaInfo;  
  my $iLineNum=0;
  print "$SignificantSNVsFile\n";
  open(my $fh, "<", $SignificantSNVsFile) or die "Cannot 
     open $SignificantSNVsFile: $!";
  while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    my @cols= split(/\s++/,$_);
    #A1FREQ N TEST BETA SE CHISQ LOG10P
   my $BetaInfo="$cols[6]\t$cols[7]\t$cols[9]\t$cols[10]\t$cols[11]\t$cols[12]";
   my $Trait=$cols[0];
   my $ChrPos="$cols[1]\t$cols[2]";
   my $A0A1="$cols[4]\t$cols[5]";
   $Trait_chrPos_RefAlt{$Trait}{$ChrPos}=$A0A1;
   $Trait_chrPos_SNVID{$Trait}{$ChrPos}=$cols[3];
   $Trait_chrPos_BetaInfo{$Trait}{$ChrPos}=$BetaInfo;
  }
  close $fh;
  return (\%Trait_chrPos_SNVID,\%Trait_chrPos_RefAlt,\%Trait_chrPos_BetaInfo);
}


# Subroutine to get all .regenie.gz files from a directory
# Params:
#   $dir => directory path
# Returns:
#   array of file names
sub get_allfileswithSpecificPostfix {
    my ($dir,$FilePostfix) = @_;
    my @files;

    opendir(my $dh, $dir) or die "Cannot open directory $dir: $!";
    while (my $file = readdir($dh)) {
        # Skip non-files
        next unless -f File::Spec->catfile($dir, $file);
        # Only include files ending with .regenie.gz
        push @files, $file if $file =~ /$FilePostfix$/;
    }
    closedir($dh);
    return @files;
}

sub convert_log10p_to_p_sstFile {
   my ($infile, $outfilePrefix) = @_;
   my $outfile="$outfilePrefix.txt";
   print "$infile\n";
   print "$outfile\n";
   open(my $out, ">", $outfile) or die "Cannot open $outfile: $!";
  # Print header
  print $out "SNP\tP\tCHR\tBP\n";
  my $iLineNumber=0;
  #CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
  #1 777135 chr1:777135:T:TC T TC 0.0997168 4237 ADD 0.252579 0.247327 1.04292 0.512659 NA

  # open compressed file using zcat
   open(my $fh, "-|", "zcat $infile") or die "Cannot open $infile: $!";
  # find column indices for ID and LOG10P
  my ($id_col, $log10p_col);
   while (<$fh>) {
    chomp;
    my @cols = split(/\s+/,$_);
    if($iLineNumber ==0)
    {
      $iLineNumber++;
      next;
    }
    last if($cols[0]>22);
    # convert LOG10P -> P
    my $pval = 10 ** (-$cols[11]);
    $iLineNumber++;
    print $out "$cols[2]\t$pval\t$cols[0]\t$cols[1]\n";
   }
   close $fh;
   close $out;
   print "Total used lines: $iLineNumber\n";
}

sub PerformLDClump_onGWASFile
{
    my ($sstCompressedFile,$trait,$OutFilePrefix,$GenotypePrefix)=@_;	
      #Reads a compressed .gz file
      #Extracts the ID column and LOG10P column
      #Converts LOG10P → P-value
      #Outputs a simple tab-delimited file with ID an d P
     convert_log10p_to_p_sstFile($sstCompressedFile,$OutFilePrefix);
     for(my $chr=1;$chr<=22;$chr++)
     {
        my $outpref_1  = $OutFilePrefix."chr$chr";
        #AFR_WGS
        my $genotype_data = $GenotypePrefix."_chr$chr";
         system(
                "bash", "-c",
                "./plink2 ".
                "--pfile $genotype_data ".
		"--maf 0.01 ".
                "--clump $OutFilePrefix.txt ".
		"--exclude range longLD_region_hg38.txt ".
                 "--clump-p1 5e-5 ".
                "--clump-p2 1e-2 ".
                "--clump-r2 0.1 ".
                "--clump-kb 1000 ".
                "--out $outpref_1"
                ) == 0 or warn "PLINK2 failed for $trait chr$chr.ld06\n";
     }
}

sub clumpsBasedSNVsCollected
{
  my ($ClumpedSNPFile,$trait,$chr,$type,$Trait_Type_Chr_snvid)=@_;
  my %Trait_Type_Chr_snvid=%{$Trait_Type_Chr_snvid};
  my $iLineNum=0;
  open(my $fh, "<", $ClumpedSNPFile) or die "Cannot open $ClumpedSNPFile: $!";
  while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    my @cols = split(/\t/,$_,5);
    if(!defined($Trait_Type_Chr_snvid{$trait}{$type}{$chr}))
    {
      $Trait_Type_Chr_snvid{$trait}{$type}{$chr}=$cols[2];
    }
    else
    {
     $Trait_Type_Chr_snvid{$trait}{$type}{$chr}.="\t$cols[2]";
    }
  }
  close $fh;
  return (%Trait_Type_Chr_snvid);
}


my $GWASsstDir_AFR="AGD35K_AFR_Phecode_GWAS";
my @AllAFRFWASFileList=get_allfileswithSpecificPostfix($GWASsstDir_AFR,
	".regenie.gz");
my $Size=scalar(@AllAFRFWASFileList);
print "NumberofFiles: $Size\n";
my $GenotypePrefix_AFR= "AGD35K_AFRWGS/AGD35K_MAC100MR002QCed_AFR";
my $GWASFilePrefix_AFR="Phecode_AFRAFR_";
my $GWASFilePostfix=".regenie.gz";
my $ClumpFilePrefix_AFR="Phecode_AFR_";
my $OutFilesDir_AFR="AGD35K_AFR_Phecode_GWAS_Corr_1";
my @PhenotypeList=();
for(my $ii=0;$ii<@AllAFRFWASFileList;$ii++)
{
   my $sstFile="$GWASsstDir_AFR/$AllAFRFWASFileList[$ii]";	
   my $Phenotype=$AllAFRFWASFileList[$ii];
   $Phenotype =~ s/\Q$GWASFilePrefix_AFR\E//ig;
   $Phenotype =~ s/\Q$GWASFilePostfix\E//ig;
   #print "$ii\t$AllAFRFWASFileList[$ii]\t$GenotypePrefix_AFR\t$Phenotype\n";
   my $OutFilePrefix= "$OutFilesDir_AFR/$ClumpFilePrefix_AFR$Phenotype";
   #PerformLDClump_onGWASFile($sstFile,$Phenotype,
   #	   $OutFilePrefix_AFR,$GenotypePrefix_AFR);
   push @PhenotypeList,$Phenotype;
}
my $GWASsstDir_EUR="AGD35K_EUR_Phecode_GWAS_Full";
my $GenotypePrefix_EUR= "AGD35K_EURWGS/AGD35K_MAC100MR002QCed_EUR";
my $GWASFilePrefix_EUR="Phecode_EUREUR_";
my $ClumpFilePrefix_EUR="Phecode_EUR_";
my $OutFilesDir_EUR="AGD35K_EUR_Phecode_GWAS_Corr_1";

#for(my $ii=0;$ii<@PhenotypeList;$ii++)
for(my $ii=0;$ii<12;$ii++)
{
   my $Phenotype=$PhenotypeList[$ii];
   my $sstFile="$GWASsstDir_EUR/$GWASFilePrefix_EUR$PhenotypeList[$ii].regenie.gz";
   my $OutFilePrefix_EUR= "$OutFilesDir_EUR/$ClumpFilePrefix_EUR$Phenotype";
   print "$ii\t$AllAFRFWASFileList[$ii]\t$GenotypePrefix_EUR\t$Phenotype\t$OutFilePrefix_EUR\n";
   #  PerformLDClump_onGWASFile($sstFile,$Phenotype,$OutFilePrefix_EUR,$GenotypePrefix_EUR);
}


my @postfixes =(".clump.clumps",".clumps");
my @clumpTypes = ("AFR","EUR");
my %Trait_Type_Chr_snvid=();
my %Trait_Chr_Type_snvCount=();
for(my $iT=0;$iT<@postfixes;$iT++)
{
   my @ClumpFiles=();
   my $
   if($iT ==0)
   {
    @ClumpFiles=get_allfileswithSpecificPostfix($OutFilesDir_AFR,$postfixes[$iT]);
   }
   else
   {
   @ClumpFiles=get_allfileswithSpecificPostfix($OutFilesDir_EUR,$postfixes[$iT]);
   }
    # Open directory
   for(my $jj=0;$jj<@ClumpFiles;$jj++)
   {
     my $file=$ClumpFiles[$jj];
     my $fullPath="$dir/$file";
     $file =~ s/\Q$postfixes[$iT]\E$//ig;
     $file =~ s/^\Q$FilePrefix\E//ig;
     my ($trait,$chr)=split(/\.chr/,$file);
     $chr =~ s/^chr//ig;
     #  print "$file\t$trait\t$chr\n";
     %Trait_Type_Chr_snvid=clumpsBasedSNVsCollected($fullPath,
             $trait,$chr,$clumpTypes[$iT],\%Trait_Type_Chr_snvid);
   }
}
=cut;
 my $SharedDir="AGD35K_NoveltyGWASHits";
 my $missingBetaInfo="NA\tNA\tNA\tNA\tNA\tNA";
 my $LabValueAFRGWASFile="$SharedDir/AGD35KAFR_LabValue_DiseaseIDsLDC0.1.txt";
 my ($Trait_chrPos_AFRSNVIDRef,$Trait_chrPos_AFRRefAltRef,
	 $Trait_chrPos_AFRBetaInfoRef)=ReadTraitBasedSignificantVariantInfo($LabValueAFRGWASFile);
 my $PheCodeAFRGWASFile="$SharedDir/AGD35KAFR_PheCode_DiseaseIDsLDC0.1.txt";
  my ($Phecode_chrPos_AFRSNVIDRef,$Phecode_chrPos_AFRRefAltRef,
	  $Phecode_chrPos_AFRBetaInfoRef)=ReadTraitBasedSignificantVariantInfo($PheCodeAFRGWASFile);

  my $LabValueEURGWASFile="$SharedDir/AGD35KAFR_LabValue_DiseaseIDsLDC0.1.txt.EURPairList.txt";
  my ($Trait_chrPos_EURSNVIDRef,$Trait_chrPos_EURRefAltRef, 
	   $Trait_chrPos_EURBetaInfoRef)=ReadTraitBasedSignificantVariantInfo($LabValueEURGWASFile);
  my $PheCodeEURGWASFile="$SharedDir/AGD35KAFR_PheCode_DiseaseIDsLDC0.1.txt.EURPairList.txt";
  my ($Phecode_chrPos_EURSNVIDRef,$Phecode_chrPos_EURRefAltRef, 
	  $Phecode_chrPos_EURBetaInfoRef)=ReadTraitBasedSignificantVariantInfo($PheCodeEURGWASFile);
my  %Trait_chrPos_AFRSNVID=%{$Trait_chrPos_AFRSNVIDRef};
my  %Trait_chrPos_AFRBetaInfo=%{$Trait_chrPos_AFRBetaInfoRef};
my  %Trait_chrPos_AFRRefAlt=%{$Trait_chrPos_AFRRefAltRef};
my  %Phecode_chrPos_AFRSNVID=%{$Phecode_chrPos_AFRSNVIDRef};
my  %Phecode_chrPos_AFRBetaInfo=%{$Phecode_chrPos_AFRBetaInfoRef};
my  %Phecode_chrPos_AFRRefAlt=%{$Phecode_chrPos_AFRRefAltRef};
my  %Trait_chrPos_EURSNVID=%{$Trait_chrPos_EURSNVIDRef};
my  %Trait_chrPos_EURBetaInfo=%{$Trait_chrPos_EURBetaInfoRef};
my  %Trait_chrPos_EURRefAlt=%{$Trait_chrPos_EURRefAltRef};
my  %Phecode_chrPos_EURSNVID=%{$Phecode_chrPos_EURSNVIDRef};
my  %Phecode_chrPos_EURBetaInfo=%{$Phecode_chrPos_EURBetaInfoRef};
my  %Phecode_chrPos_EURRefAlt=%{$Phecode_chrPos_EURRefAltRef};
#AGD35K_OTG_GWASCata_Disease_Labvalue.tsv
#AGD35K_OTG_GWASCata_Disease_Phecode.tsv
my $iLineNum=0;
#LabTrait_AGD35K diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO name_AGD35K
my %LabTrait_EFO=();
my %LabTrait_ReportedTrait=();
my %LabTrait_name_EFO=();
my %LabTrait_nameAGD35K=();

my $LabValueEFOFile="$SharedDir/AGD35K_OTG_GWASCata_Disease_Labvalue.tsv";
open(my $fh, "<", $LabValueEFOFile) or die "Cannot
     open $LabValueEFOFile: $!";
  while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      print "$_\n";	    
      $iLineNum++;
      next;
    }
    my @cols= split(/\t/,$_);
    if(!defined($LabTrait_EFO{$cols[0]}))
    {
       $LabTrait_EFO{$cols[0]}=$cols[1];
       $LabTrait_ReportedTrait{$cols[0]}=$cols[2];
       $LabTrait_name_EFO{$cols[0]}=$cols[4];
       $LabTrait_nameAGD35K{$cols[0]}=$cols[6];
    }
    else
    {
       $LabTrait_EFO{$cols[0]}.=":$cols[1]";
       $LabTrait_ReportedTrait{$cols[0]}.=":$cols[2]";
       $LabTrait_name_EFO{$cols[0]}.=":$cols[4]";
       $LabTrait_nameAGD35K{$cols[0]}.=":$cols[6]";    
    }
    #print "$_\n";
 }
 close $fh;
  foreach my $Trait (keys %LabTrait_EFO) {
   if(defined($Trait_chrPos_AFRSNVID{$Trait}))
   {	  

     print "$LabTrait_EFO{$Trait}\t"; #$LabTrait_ReportedTrait{$Trait}\t";
     my @Name_EFOList=split(/:/,  $LabTrait_name_EFO{$Trait});
     my %seen;
     my @unique_Str = grep { !$seen{$_}++ } @Name_EFOList;
     my $Name_EFOListStr=join(":",@unique_Str);
     $LabTrait_name_EFO{$Trait}=$Name_EFOListStr;
     
     my @nameAGD35KList=split(/:/,$LabTrait_nameAGD35K{$Trait});
     my %seen1;
     my @unique_nameAGD35KStr = grep { !$seen1{$_}++ } @nameAGD35KList;
     my $Name_AGD35KListStr=join(":",@unique_nameAGD35KStr);
     $LabTrait_nameAGD35K{$Trait}=$Name_AGD35KListStr;
     # print "$LabTrait_name_EFO{$Trait}\t$LabTrait_nameAGD35K{$Trait}\n";
    #print "$LabTrait_name_EFO{$Trait}\t$LabTrait_nameAGD35K{$Trait}\n";
   }
}

$iLineNum=0;
#phecode diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO phecode1.2_simpleLabel  phecode1.2_category
#LabTrait_AGD35K diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO name_AGD35K
my %Phecode_EFO=();
my %Phecode_ReportedTrait=();
my %Phecode_name_EFO=();
my %Phecode_nameAGD35K=();
my %Phecode_category=();
my $PhecodeEFOFile="$SharedDir/AGD35K_OTG_GWASCata_Disease_Phecode.tsv";
open(my $fh, "<", $PhecodeEFOFile) or die "Cannot
     open $PhecodeEFOFile: $!";
  while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      print "$_\n";
      $iLineNum++;
      next;
    }
    my @cols= split(/\t/,$_);
    if(!defined($Phecode_EFO{$cols[0]}))
    {
       $Phecode_EFO{$cols[0]}=$cols[1];
       $Phecode_ReportedTrait{$cols[0]}=$cols[2];
       $Phecode_name_EFO{$cols[0]}=$cols[4];
       $Phecode_nameAGD35K{$cols[0]}=$cols[6];
       $Phecode_category{$cols[0]}=$cols[7]
    }
    else
    {
       $Phecode_EFO{$cols[0]}.=":$cols[1]";
       $Phecode_ReportedTrait{$cols[0]}.=":$cols[2]";
       $Phecode_name_EFO{$cols[0]}.=":$cols[4]";
       $Phecode_nameAGD35K{$cols[0]}.=":$cols[6]";
       $Phecode_category{$cols[0]}=$cols[7];
    }
    #print "$_\n";
 }
 close $fh;
 
 foreach my $phecode (keys %Phecode_EFO) {
   next if(!defined($Phecode_chrPos_AFRSNVID{$phecode}));
   if(defined($Phecode_chrPos_AFRSNVID{$phecode}))
   {
     print "$Phecode_EFO{$phecode}\t"; #$LabTrait_ReportedTrait{$Trait}\t";
     my @Name_EFOList=split(/:/,  $Phecode_name_EFO{$phecode});
     my %seen;
     my @unique_Str = grep { !$seen{$_}++ } @Name_EFOList;
     my $Name_EFOListStr=join(":",@unique_Str);
     $Phecode_name_EFO{$phecode}=$Name_EFOListStr;

     my @nameAGD35KList=split(/:/,$Phecode_nameAGD35K{$phecode});
     my %seen1;
     my @unique_nameAGD35KStr = grep { !$seen1{$_}++ } @nameAGD35KList;
     my $Name_AGD35KListStr=join(":",@unique_nameAGD35KStr);
     $Phecode_nameAGD35K{$phecode}=$Name_AGD35KListStr;

     if(!defined($Phecode_category{$phecode}))
     {
	     $Phecode_category{$phecode}="NA";
     }
     else
     {
        my @PhecodeCategoryList=split(/:/,$Phecode_category{$phecode});
        my %seen2;
        my @unique_CategoryStr = grep { !$seen2{$_}++ } @PhecodeCategoryList;
        my $CategoryListStr=join(":",@unique_CategoryStr);
        $Phecode_category{$phecode}=$Name_AGD35KListStr;
     }
     #print "$phecode\t$Phecode_category{$phecode}\n";
     # print "$Phecode_name_EFO{$phecode}\t$Phecode_nameAGD35K{$phecode}\t$Phecode_category{$phecode}\n";
    #print "$LabTrait_name_EFO{$Trait}\t$LabTrait_nameAGD35K{$Trait}\n";
   }
}

my $OutFile="$SharedDir/AGD35K_NovelSNVResult_0208_2026.tsv";
open(my $OUT, ">", "$OutFile") or die "Cannot open $OutFile: $!";
print $OUT "Number\tID\tA0\tA1\tTrait/Phecode\tName\tCategory\tEFO\tName_EFO\t";
print $OUT "A1FREQ_AFR\tN_AFR\tBETA_AFR\tSE_AFR\tCHISQ_AFR\tLOG10P_AFR\t";
print $OUT "A1FREQ_EUR\tN_EUR\tBETA_EUR\tSE_EUR\tCHISQ_EUR\tLOG10P_EUR\n";
my $NumberofSNVs=0;
foreach my $Trait (keys %Trait_chrPos_AFRSNVID) {
    foreach my $chrpos (keys %{$Trait_chrPos_AFRSNVID{$Trait} })
     {
       if(defined($Trait_chrPos_EURSNVID{$Trait}{$chrpos}))
        {
	    if($Trait_chrPos_EURSNVID{$Trait}{$chrpos} ne $Trait_chrPos_AFRSNVID{$Trait}{$chrpos})
	    {
	     print "$Trait_chrPos_EURSNVID{$Trait}{$chrpos}\t$Trait_chrPos_AFRSNVID{$Trait}{$chrpos}\n";	    
	    }
	    $NumberofSNVs++;
	    print $OUT "$NumberofSNVs\t$Trait_chrPos_AFRSNVID{$Trait}{$chrpos}\t$Trait_chrPos_AFRRefAlt{$Trait}{$chrpos}\t$Trait\t";
	    print $OUT "$LabTrait_nameAGD35K{$Trait}\tNA\t$LabTrait_EFO{$Trait}\t$LabTrait_name_EFO{$Trait}\t";
	    print $OUT "$Trait_chrPos_AFRBetaInfo{$Trait}{$chrpos}\t$Trait_chrPos_EURBetaInfo{$Trait}{$chrpos}\n";
       }
       else
       {
	    $NumberofSNVs++;   
	      print $OUT "$NumberofSNVs\t$Trait_chrPos_AFRSNVID{$Trait}{$chrpos}\t$Trait_chrPos_AFRRefAlt{$Trait}{$chrpos}\t$Trait\t";
	      print $OUT "$LabTrait_nameAGD35K{$Trait}\tNA\t$LabTrait_EFO{$Trait}\t$LabTrait_name_EFO{$Trait}\t";
	      print $OUT "$Trait_chrPos_AFRBetaInfo{$Trait}{$chrpos}\t$missingBetaInfo\n";
       }       
    }
}
foreach my $phecode (keys %Phecode_chrPos_AFRSNVID) {
    foreach my $chrpos (keys %{$Phecode_chrPos_AFRSNVID{$phecode} })
     {
       if(defined($Phecode_chrPos_EURSNVID{$phecode}{$chrpos}))
         {
	    $NumberofSNVs++;	 
            if($Phecode_chrPos_EURSNVID{$phecode}{$chrpos} ne $Phecode_chrPos_AFRSNVID{$phecode}{$chrpos})
            {
	     #print "$Phecode_chrPos_EURSNVID{$chr}{$pos}\t$Phecode_chrPos_AFRSNVID{$chr}{$pos}\n";
            }
	     print $OUT "$NumberofSNVs\t$Phecode_chrPos_AFRSNVID{$phecode}{$chrpos}\t$Phecode_chrPos_AFRRefAlt{$phecode}{$chrpos}\t$phecode\t";
             my $phecode_num = 0 + $phecode; 
	     print $OUT "$Phecode_nameAGD35K{$phecode_num}\t$Phecode_category{$phecode_num}\t";
	     print $OUT "$Phecode_EFO{$phecode_num}\t$Phecode_name_EFO{$phecode_num}\t";
	      print $OUT "$Phecode_chrPos_AFRBetaInfo{$phecode}{$chrpos}\t$Phecode_chrPos_EURBetaInfo{$phecode}{$chrpos}\n";
         }
	 else
	 {
	   $NumberofSNVs++;	 
	     print $OUT "$NumberofSNVs\t$Phecode_chrPos_AFRSNVID{$phecode}{$chrpos}\t$Phecode_chrPos_AFRRefAlt{$phecode}{$chrpos}\t$phecode\t";
	     my $phecode_num = 0 + $phecode;
	     print $OUT "$Phecode_nameAGD35K{$phecode_num}\t$Phecode_category{$phecode_num}\t";
	     print $OUT "$Phecode_EFO{$phecode_num}\t$Phecode_name_EFO{$phecode_num}\t";
             print $OUT "$Phecode_chrPos_AFRBetaInfo{$phecode}{$chrpos}\t$missingBetaInfo\n";  
	 } 
    }
}
close $OUT;
# AGD35K_NoveltyGWASHits
#AGD35K_OTG_GWASCata_Disease_Labvalue.tsv
#AGD35K_OTG_GWASCata_Disease_Phecode.tsv
#AGD35KAFR_LabValue_DiseaseIDsLDC0.1.txt
#AGD35KAFR_LabValue_DiseaseIDsLDC0.1.txt.EURPairList.txt
#AGD35KAFR_PheCode_DiseaseIDsLDC0.1.txt
#AGD35KAFR_PheCode_DiseaseIDsLDC0.1.txt.EURPairList.txt
=cut;
