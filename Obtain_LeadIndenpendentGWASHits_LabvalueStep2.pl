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

