#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);
use Scalar::Util ();

require "./SharedFunctionList.pl";
my $Population="AFR";#AFR,EUR
#LabvalueAFR_MONOABClumpedSNVs.txt
#Phecode_AFRAFR_851ClumpedSNVs.txt
my @PhenotypeList=("Phecode","LabValue");
my @FilePrefixList=("Phecode_AFRAFR_","LabvalueAFR_");
my $FilePostfix="ClumpedSNVs.txt";
my %TraitChr_AlleleListFull_AFR=();
my %TraitChr_SNVIDListFull_AFR=();
my %TraitChrPos_BetaFreqFull_AFR=();
my %TraitChrPos_SNVIDListFull_AFR=();
my  %TraitChr_pos_id_log10List_AFR=();
for(my $iT=0;$iT<@PhenotypeList;$iT++)
{
  my $ClumpSNVDir="AGD35K_$Population"."_$PhenotypeList[$iT]_GWAS_ClumpSNVs";
  my @ClumpedSNVList=get_allfileswithSpecificPostfix($ClumpSNVDir,
	  "ClumpedSNVs.txt"); 
  my $NumberofFiles=scalar(@ClumpedSNVList);
  print "NumberofFiles:$NumberofFiles\n";
 
   for(my $jj=0;$jj<@ClumpedSNVList;$jj++)
   {
     my $Trait=$ClumpedSNVList[$jj];
     my $fullPath="$ClumpSNVDir/$ClumpedSNVList[$jj]";
     #  print "fullPath:$fullPath\n";
      $Trait =~ s/\Q$FilePostfix\E$//ig;
      $Trait =~ s/^\Q$FilePrefixList[$iT]\E//ig;
      # print "$Trait\n";
      my ($a0a1_AFRRef,
	  $SNVIDList_AFRRef,
          $betaFreq_AFRRef,
	  $posid_Log10PRef,
	  $SNVIDListbypos_AFRRef)=
         ReadGWASclumpedSNVListFile($fullPath,
		 $Trait,
         \%TraitChr_AlleleListFull_AFR,
	 \%TraitChr_SNVIDListFull_AFR,
         \%TraitChrPos_BetaFreqFull_AFR,
	 \%TraitChr_pos_id_log10List_AFR,
         \%TraitChrPos_SNVIDListFull_AFR);
      %TraitChr_AlleleListFull_AFR=%{$a0a1_AFRRef};
      %TraitChr_SNVIDListFull_AFR=%{$SNVIDList_AFRRef};
      %TraitChrPos_BetaFreqFull_AFR=%{$betaFreq_AFRRef};
      %TraitChr_pos_id_log10List_AFR=%{$posid_Log10PRef};
      %TraitChrPos_SNVIDListFull_AFR=%{$SNVIDListbypos_AFRRef};
    }
}

my $iTotalIndependentLoci=0;
my %MergedTraitChrPosList=();
my %RemovedTraitChrPosList=();
my $iRemoveSNVs=0;
my $iRelatedTrait=0;
my $iTrait=0;
my %Chr_RemoveSegments=();
$Chr_RemoveSegments{1}=0;
$Chr_RemoveSegments{11}=0;
$Chr_RemoveSegments{16}=0;
$Chr_RemoveSegments{22}=0;
foreach my $trait (keys %TraitChr_pos_id_log10List_AFR) {
    foreach my $chr (keys %{$TraitChr_pos_id_log10List_AFR{$trait}}) {
	    # print "$trait\t$chr\t$TraitChr_pos_id_log10List_AFR{$trait}{$chr}\n";	
    if(($trait =~ /^(NTAuto|WBC|288.2|MONOAB)$/ && $chr ==1) ||($chr ==16) ||($chr ==11) ||($chr ==22))
    {    
	    my $iStartPos=0;
	    my $iEndPos=0;
	    if($trait =~ /^(NTAuto|WBC|288.2|MONOAB)$/ && $chr ==1)
             {
                $iStartPos=90385392;
		$iEndPos=177814914;
             }
             elsif($chr ==11)#SCD 
             {
                 $iStartPos=2000000;
		 $iEndPos=8000000;
	     }
	     elsif($chr ==16)
             {
                 $iStartPos=1;
                 $iEndPos=3000000;
             }
             elsif($chr ==22)#apoi1
             {
                $iStartPos=35650526;
	        $iEndPos=37650526;
	     }
	    my ($SegmentKeptTraitChrPosListRef,
	    $SegmentRemovedTraitChrPosListRef,
	    $iRemoveSNVs_segment)=OverlappIndependentLocitoLargeSegments($trait,
	    $chr,
	    $TraitChr_pos_id_log10List_AFR{$trait}{$chr},
	    $iStartPos,
	    $iEndPos,
	    \%MergedTraitChrPosList,
	    \%RemovedTraitChrPosList);
	     %MergedTraitChrPosList=%{$SegmentKeptTraitChrPosListRef};
	    %RemovedTraitChrPosList=%{$SegmentRemovedTraitChrPosListRef};
           
	    my $iRawCount=scalar(split(/\t/,$TraitChr_pos_id_log10List_AFR{$trait}{$chr}));
	    my $iKeepVariants=$iRawCount-$iRemoveSNVs_segment;
	    #if($iStartPos ==0 && $iEndPos ==0)
	    {
		 print "Group: $iRawCount\t$iRemoveSNVs_segment\t$iKeepVariants\n";
	    }
            $iTotalIndependentLoci+=$iRawCount;
	    $iRemoveSNVs+=$iRemoveSNVs_segment;
	    $Chr_RemoveSegments{$chr}+=$iRemoveSNVs_segment;
	    #   print "iTotalIndependentLoci\t$trait\t$chr\t$iRemoveSNVs_segment\n";
       }
   }
}

for my $chr (keys %Chr_RemoveSegments) {
    print "$chr => $Chr_RemoveSegments{$chr}\n";
}
#print "iTotalIndependentLoci\t$iRemoveSNV\t$iTotalIndependentLoci\n";

#Construct collapsed Chr POS LIST
my %CollapedTraitChr_PosList_AFR=();
my %CollapedTraitChr_SNVIDList_AFR=();
my $iTotalTraitChr_PosList_AFR=0;
foreach my $trait (sort keys %TraitChr_AlleleListFull_AFR) {
   foreach my $chrpos (sort keys %{$TraitChr_AlleleListFull_AFR{$trait}}) {
       next if(defined($RemovedTraitChrPosList{$trait}{$chrpos}) &&
               !defined($MergedTraitChrPosList{$trait}{$chrpos}));
        my ($chr,$pos)=split(/\t/,$chrpos);
	if(!defined($CollapedTraitChr_PosList_AFR{$trait}{$chr}))
	{
          $CollapedTraitChr_PosList_AFR{$trait}{$chr}="$pos";
	  $CollapedTraitChr_SNVIDList_AFR{$trait}{$chr}
	    =$TraitChrPos_SNVIDListFull_AFR{$trait}{$chrpos};
	}
	else
	{
         $CollapedTraitChr_PosList_AFR{$trait}{$chr}.="\t$pos";
	 $CollapedTraitChr_SNVIDList_AFR{$trait}{$chr}
	           .= "\t$TraitChrPos_SNVIDListFull_AFR{$trait}{$chrpos}";
	}
	$iTotalTraitChr_PosList_AFR++;
     }
}
print "iTotalTraitChr_PosList_AFR:$iTotalTraitChr_PosList_AFR\n";


my $GroupHits_Type="DiseaseIDs";
my $OTGGWASHitsFile="AGD35K_NoveltyGWASHits/Known_LabValue".
        "_GWASHits_OTG_GWASCataby$GroupHits_Type"."_0204_2026.tsv";
my %TraitChr_Pos_KnownHits=();
%TraitChr_Pos_KnownHits=ReadKnownGWASHits_ByTraitFromSourceFile($OTGGWASHitsFile,\%TraitChr_Pos_KnownHits);
$OTGGWASHitsFile="AGD35K_NoveltyGWASHits/Known_PheCode".
        "_GWASHits_OTG_GWASCataby$GroupHits_Type"."_0204_2026.tsv";
%TraitChr_Pos_KnownHits=ReadKnownGWASHits_ByTraitFromSourceFile($OTGGWASHitsFile,\%TraitChr_Pos_KnownHits);

  my ($iOverlappedLocus,
     $iNovelLocus,
     $iTotalSignificantSNVs,
     $NovelSNVsListInfoRef)=
     DetectNovelityAssoctionBasedOnKnownGWASsHits(
          \%CollapedTraitChr_PosList_AFR,
          \%CollapedTraitChr_SNVIDList_AFR,
          \%TraitChr_Pos_KnownHits);
print "OTGKnown Hits Overlap: $iOverlappedLocus\tNovel: $iNovelLocus\tTotal:$iTotalSignificantSNVs\n";

my $GWASCatalogHitsFile="AGD35K_NoveltyGWASHits/Known_LabValue".
        "_GWASHits_GWASCataby$GroupHits_Type"."_0206_2026.tsv";
my %TraitChr_Pos_KnownHits_GWASCata=();
   %TraitChr_Pos_KnownHits_GWASCata=ReadKnownGWASHits_ByTraitFromSourceFile(
           $GWASCatalogHitsFile,
	   \%TraitChr_Pos_KnownHits_GWASCata);
$GWASCatalogHitsFile="AGD35K_NoveltyGWASHits/Known_PheCode".
        "_GWASHits_GWASCataby$GroupHits_Type"."_0206_2026.tsv";
%TraitChr_Pos_KnownHits_GWASCata=ReadKnownGWASHits_ByTraitFromSourceFile(
           $GWASCatalogHitsFile,
	   \%TraitChr_Pos_KnownHits_GWASCata);
    my ($TraitChr_PosRef,
        $TraitChr_SNVIDlistRef)= DetectNovelityAssoctionBasedOnKnownGWASsHits_2ndTime(
	\%TraitChr_Pos_KnownHits_GWASCata,
	$NovelSNVsListInfoRef);
   ($iOverlappedLocus,
    $iNovelLocus, 
    $iTotalSignificantSNVs,
    $NovelSNVsListInfoRef)=
    DetectNovelityAssoctionBasedOnKnownGWASsHits(
     $TraitChr_PosRef,
     $TraitChr_SNVIDlistRef,
     \%TraitChr_Pos_KnownHits_GWASCata);
     my @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
   print "GWAS catalogue Overlap: $iOverlappedLocus\tNovel: $iNovelLocus\tTotal:$iTotalSignificantSNVs\n";

   my $iTotalNoveltySNVs=0;
    my %TraitSNVID_Novelty=();
    my %Trait_ChrPosNovelty=();
    my %Trait_ChrPos_SNVIDList_Novel=();
     my %TraitChr_SNVIDList_Novelty=();
   for(my $ii=0;$ii<@NovelSNVsListInfo;$ii++)
    {
      my @cols=split(/\t/,$NovelSNVsListInfo[$ii]);
      $TraitSNVID_Novelty{$cols[0]}{$cols[3]}=1;
     if(!defined($Trait_ChrPosNovelty{$cols[0]}))
      {
         $Trait_ChrPosNovelty{$cols[0]}="$cols[1],$cols[2]";
      }
      else
      {
        $Trait_ChrPosNovelty{$cols[0]}.="\t$cols[1],$cols[2]";
      }
      $Trait_ChrPos_SNVIDList_Novel{$cols[0]}{$cols[1]}{$cols[2]}=$cols[3];
      #print "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]\n";
      $iTotalNoveltySNVs++;
   }

   #KnownGWAS_AGD35K_Phecode_EFO_OTGV25.12.tsv
   #KnownGWAS_AGD35K_Phecode_EFO_OTGV25.12.tsv
my $LabvalueEFOFile= "AGD35K_NoveltyGWASHits/AGD35K_OTG_GWASCata_Disease_Labvalue.tsv";
#LabTrait_AGD35K diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO name_AGD35K
#LYMPAB  EFO_0004587     Lymphocyte count (UKB data field 30120) OTG_Exact       lymphocyte count       A quantification of lymphocytes in blood.        Lymphocyte count

my ($LabTrait_EFORef,$LabTrait_EFONameRef)=ObtainLabValueBaseEFO_Name($LabvalueEFOFile);
my $PhecodeEFOFile= "AGD35K_NoveltyGWASHits/AGD35K_OTG_GWASCata_Disease_Phecode.tsv";
my ($Phecode_EFORef,$Phecode_EFONameRef)= ObtainPhecodeBaseEFO_Name($PhecodeEFOFile);
#phecode diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO phecode1.2_simpleLabel  phecode1.2_category
#8       MONDO_0000916   Intestinal infection    OTG_Exact       intestinal infectious disease   An infectious disease involving a pathogenic inflammatory response in the intestinal mucosa.    intestinal infection    infectious diseases
my %LabTrait_EFO=%{$LabTrait_EFORef};
my %LabTrait_EFOName=%{$LabTrait_EFONameRef};
my %Phecode_EFO=%{$Phecode_EFORef};
my %Phecode_EFOName=%{$Phecode_EFONameRef};
#Obtain All Trait_Chr INDependent List
 print "TraitChr_novel SNVs count:$iTotalNoveltySNVs\n";
 my $iTotalTrait_ChrList=0;
 my $iCollectedTrait_ChrList=0;
 my %AllTraitSNVIDInfo=();
  my $EURBasedIndependentLociInfoFile="AGD35K_NoveltyGWASHits/AllAFRrelatedEUR_IndLoci_GWASSubsetCollapse_0217_2026.txt";
  #open(my $fh, ">", $EURBasedIndependentLociInfoFile) or die "Cannot open $EURBasedIndependentLociInfoFile: $!";
  #close $fh;
 foreach my $trait (keys %CollapedTraitChr_SNVIDList_AFR) {
	 my %TraitBasedSNVIDList=();	 
         my $iIndependentLoci_AFR=0;
	 foreach my $chr (keys %{$CollapedTraitChr_SNVIDList_AFR{$trait}}) {
		 #print "$CollapedTraitChr_SNVIDList_AFR{$trait}{$chr}\n";
               my @SNVIDList=split(/\t/,$CollapedTraitChr_SNVIDList_AFR{$trait}{$chr});
               for(my $ii=0;$ii<@SNVIDList;$ii++)
	       {
                  $TraitBasedSNVIDList{$SNVIDList[$ii]}=1;
	       }
               my $iCountLoci=scalar(@SNVIDList);
	       $iIndependentLoci_AFR+=$iCountLoci;
               $iTotalTrait_ChrList+=$iCountLoci;
           }
	 my $GWASFilePrefix=""; 
	 if(Scalar::Util::looks_like_number($trait)) { 
	      $GWASFilePrefix="AGD35K_EUR_Phecode_GWAS_Full/Phecode_EUREUR_$trait";
          }
	  else
	  {
	     $GWASFilePrefix="AGD35K_EUR_LabValue_GWAS/LabvalueEUR_$trait";	  
	  } 
	  # print "$trait $GWASFilePrefix: AFR Count: $iIndependentLoci_AFR\n";
	  #   my $iObtainedList=SelectGWASSubsetbyID($GWASFilePrefix,
	  # \%TraitBasedSNVIDList,
	  # $trait, 
	  # $EURBasedIndependentLociInfoFile);	
	  # $iCollectedTrait_ChrList+=$iObtainedList; 
}
print "Number of total: $iTotalTrait_ChrList\t$iCollectedTrait_ChrList\n";


#LabTrait        GWAStrait
#Johns_AGD35K_LabTraits.txt
my $JohnLabTraitNameFile="./Johns_AGD35K_LabTraits.txt";
my $iLineNum=0;
my %Lab_FullName=();
open(my $fh, "<", $JohnLabTraitNameFile) or
  die "Cannot open $JohnLabTraitNameFile: $!";
while (my $line = <$fh>) {
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    chomp $line;
    my @cols=split(/\t/,$line);
    $Lab_FullName{$cols[0]}=$cols[1];
}
close($fh);

#phecode1.2.tsv
#phecode1.2_code phecode1.2_label        phecode1.2_simpleLabel  phecode1.2_category
my $Phecode1_2File="./phecode1.2.tsv";
$iLineNum=0;
my %Code_Label=();
open(my $fh, "<", $Phecode1_2File) or
  die "Cannot open $Phecode1_2File: $!";
while (my $line = <$fh>) {
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    chomp $line;
    my @cols=split(/\t/,$line);
    $Code_Label{$cols[0]}=$cols[1];
}
close($fh);

# #CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
my ($trait_snvid_eurRef, $trait_refalt_eurRef,$trait_betafreq_eurRef) =
       ReadTraitBasedSignificantVariantInfo($EURBasedIndependentLociInfoFile);
my %trait_snvid_eur=%{$trait_snvid_eurRef};
my %trait_refalt_eur=%{$trait_refalt_eurRef};
my %trait_betafreq_eur=%{$trait_betafreq_eurRef};

my $iSame=0;
my $log10P_threshold=(-1)*log(0.05)/log(10);
#print "log10P_threshold:$log10P_threshold\n";
#A1FREQ N BETA SE LOG10P
my $outsst_AGD35KAFRfile = "sst_AllIndependentloici_AGD35kAFR_Merged_0215_2026.txt";
open(my $out, ">", $outsst_AGD35KAFRfile) 
    or die "Cannot open $outsst_AGD35KAFRfile: $!";
 my $is_signal_Total=0;
 my $is_signal_005_Total=0;
 print $out "trait\tPhecode/abbreviation\tchr\tpos\tNovelStatus\tEFOList\tEFONameList\tSNVID\tA0\tA1";
 print $out "\tA1Freq_AFR\tN_AFR\tBeta_AFR\tSE_AFR\tLog10P_AFR\t";
  print $out "A1Freq_EUR\tN_EUR\tBeta_EUR\tSE_EUR\tLog10P_EUR\t";#A1FREQ N BETA SE LOG10P for EUR
 print $out "is_same_direction\tis_same_direction_005\n";#$is_signal\t$is_signal_005

foreach my $trait (sort keys %TraitChr_AlleleListFull_AFR) {
    foreach my $chrpos (sort keys %{$TraitChr_AlleleListFull_AFR{$trait}}) {
       next if(defined($RemovedTraitChrPosList{$trait}{$chrpos}) &&  
	       !defined($MergedTraitChrPosList{$trait}{$chrpos}));
       my ($chr,$pos)=split(/\t/,$chrpos);
       my $NovelStatus=0; 
       if(defined($Trait_ChrPos_SNVIDList_Novel{$trait}{$chr}{$pos}))
       {
       $NovelStatus=1;;
       }   
        my $TraitName="";
	my $ConvertedTrait=$trait;
	$ConvertedTrait =~ s/^0+//;
	if(defined($Lab_FullName{$trait}))
	{
           $TraitName=$Lab_FullName{$trait};
        }
        else
	{
	   $TraitName=$Code_Label{$trait};
         }
	 my $EFOList="NA";
	 my $EFONameList="NA";
	if(defined($LabTrait_EFO{$ConvertedTrait}))
	{
    	  $EFOList=$LabTrait_EFO{$ConvertedTrait};
          $EFONameList=$LabTrait_EFOName{$ConvertedTrait};	  
	}
	elsif(defined($Phecode_EFO{$ConvertedTrait}))
	{
         $EFOList=$Phecode_EFO{$ConvertedTrait};
	 $EFONameList=$Phecode_EFOName{$ConvertedTrait};
	}
       if(!defined($trait_betafreq_eur{$trait}{$chrpos}))
        {
          print $out "$TraitName\t$trait\t$chrpos\t$NovelStatus\t";
	  print $out "$EFOList\t$EFONameList\t$TraitChr_AlleleListFull_AFR{$trait}{$chrpos}";
	  print $out "\t$TraitChrPos_BetaFreqFull_AFR{$trait}{$chrpos}\t";
	  print $out "NA\tNA\tNA\tNA\tNA\t";#A1FREQ N BETA SE LOG10P for EUR
	  print $out "NA\tNA\n";#$is_signal\t$is_signal_005
	 }  
	else #if(defined($trait_betafreq_eur{$trait}{$chrpos}))
	{
         my @BetaFreq_AFR=split(/\t/,
		 $TraitChrPos_BetaFreqFull_AFR{$trait}{$chrpos});
	 my @BetaFreq_EUR=split(/\t/,
                 $trait_betafreq_eur{$trait}{$chrpos});
	 my $Beta_sign_AFR = $BetaFreq_AFR[2] < 0 ? -1 : 
	            ($BetaFreq_AFR[2] > 0 ? 1 : 0);
	 my $Beta_sign_EUR = $BetaFreq_EUR[2]< 0 ? -1 : 
	            ($BetaFreq_EUR[2] > 0 ? 1 : 0);
	 my $point05Status=$BetaFreq_EUR[4] >= $log10P_threshold ? 1 : 0;	
	 my $is_signal=$Beta_sign_AFR == $Beta_sign_EUR ? 1 : 0;
	 # Flag if the EUR p-value passes threshold AND 
	 # effect directions are concordant  
	 my $is_signal_005 = ($point05Status && $is_signal) ? 1 : 0;
          print $out "$TraitName\t$trait\t$chrpos\t$NovelStatus\t";
	  print $out "$EFOList\t$EFONameList\t$TraitChr_AlleleListFull_AFR{$trait}{$chrpos}\t";
	  print $out "$TraitChrPos_BetaFreqFull_AFR{$trait}{$chrpos}\t";
	  print $out "$trait_betafreq_eur{$trait}{$chrpos}\t";
	  print $out "$is_signal\t$is_signal_005\n";
	 $is_signal_Total+=$is_signal;
         $is_signal_005_Total+=$is_signal_005;
        }
     }
}
close $out;
#print "SAme alleles:$iSame\t$is_signal_Total\t$is_signal_005_Total\n";


        #collect all novelty SNVs into a file
	#for(my $jj=0;$jj<@LeadSNVsFiles;$jj++)
	# {
	#my $Trait=$LeadSNVsFiles[$jj];
	# my $fullPath="$IndependentLeadSNVsDir/$LeadSNVsFiles[$jj]";
	#$Trait =~ s/\Q$FilePostFix\E$//ig;
	#$Trait =~ s/^\Q$FilePrefix\E//ig;
	# WriteOutSignificantFile($fullPath,$Trait,\%TraitSNVID_Novelty,$OutNoveltyGWASResultFile);
	# }
