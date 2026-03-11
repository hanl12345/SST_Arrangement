#Merge AGD35K AFR Sample size,Case AND Control count
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);
require "./SharedFunctionList.pl";

 my $psam_AFRCompressedFile="AGD35K_Phecode122_0308_2026/AGD35K_MAC100MR002QCed_AFR_chr22.psam.gz";
 my ($SampleIDList_AFRRef,$TotalSampleSize_AFR)=readPsamSampleIDList(
	$psam_AFRCompressedFile);
 print "TotalSampleSize_AFR:$TotalSampleSize_AFR\n";
 my %SampleIDList_AFR=%{$SampleIDList_AFRRef};
#AGD35K_Phecode_summary_lambda.tsv
#popultaion      phecode TotalVariants   SignificantVariants     LambdaGC        SampleSize
#AFR     530.1   14507284        0       1.0461  24752
 my $PhecodeLambdaFile_AFR="/home/jupyter/AGD35K_GWAS122_2026/GWAS_AnalysisV0305_2026/AGD35K_AFR_Phecode";
 $PhecodeLambdaFile_AFR.="/AGD35K_Phecode_summary_lambda.tsv";
 my %Phecode_LambdaSampleSize= ObtainLambda_SampleSize_FromSummaryLambdaFile(
	$PhecodeLambdaFile_AFR);
 my $SampleSize_Phecode=keys %Phecode_LambdaSampleSize;
print "SampleSize_Phecode:$SampleSize_Phecode\n";
#AGD35K_AFR_Labvaluesummary_lambda.tsv
#File    TotalVariants   SignificantVariants     LambdaGC
#LabvalueAFR_RBC 12330123        390     1.0370  13173
 my $LabvalueLambdaFile_AFR="/home/jupyter/AGD35K_GWAS122_2026/GWAS_AnalysisV0305_2026/AGD35K_AFR_LabValue";
 $LabvalueLambdaFile_AFR.="/AGD35K_AFR_Labvaluesummary_lambda.tsv";
 my %Labvalue_LambdaSampleSize= ObtainLambda_SampleSize_FromSummaryLambdaFile(
        $LabvalueLambdaFile_AFR);
 my $SampleSize_Labvalue=keys %Labvalue_LambdaSampleSize;
print "SampleSize_Labvalue:$SampleSize_Labvalue\n";
my $PhecodeFile_AGD35K="AGD35K_Phecode122_0308_2026/AGD_35k_code_data_floor_binarized_wide_format_phecode_12_gwas.txt.gz";
my $OutPhecode122_AGD35KFile="AGD35K_Phecode122_0308_2026/AGD_35k_Phecode122_AGD35K.txt.gz";
my ($caseCount_AGD35K_122Ref,
     $controlCount_AGD35K_122Ref) =constructSubsetAndCountCaseControl(
        $PhecodeFile_AGD35K,
        $SampleIDList_AFRRef,
        \%Phecode_LambdaSampleSize,
        $OutPhecode122_AGD35KFile);
my $OutLambda_SampleSize_AGD35KFile="AGD35K_Phecode122_0308_2026/AGD_35k_LambdaSampleSize_AGD35K.txt.gz";
my $Phecode1_2AnnotationFile="AGD35K_Phecode122_0308_2026/phecode1.2.tsv.gz";
my $Phecode1_2AnnotationRef=LoadTrait_AnnotationHash($Phecode1_2AnnotationFile);
my $LabTrait_AnnotationFile="AGD35K_Phecode122_0308_2026/Johns_AGD35K_LabTraits.txt.gz";
my $LabTrait_AnnotationRef=LoadTrait_AnnotationHash($LabTrait_AnnotationFile);

OutSampleSize_Lambda_AGD35K
($caseCount_AGD35K_122Ref,
          $controlCount_AGD35K_122Ref,
          \%Phecode_LambdaSampleSize,,
          \%Labvalue_LambdaSampleSize,
          $Phecode1_2AnnotationRef,
          $LabTrait_AnnotationRef,
          $OutLambda_SampleSize_AGD35KFile);
#zcat Johns_AGD35K_LabTraits.txt.gz |head
#LabTrait        GWAStrait
#BASOAB  Basophil count
#AGD35K_Phecode122_0308_2026$ zcat phecode1.2.tsv.gz |head
#phecode1.2_code phecode1.2_label        phecode1.2_simpleLabel  phecode1.2_category
