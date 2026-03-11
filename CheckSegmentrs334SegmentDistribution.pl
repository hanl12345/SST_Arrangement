#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

# Subroutine to get all .regenie.gz files from a directory
# Params:
#   $dir => directory path
# Returns:
#   array of file names
sub get_regenie_SignificantSNVs_files {
    my ($dir) = @_;
    my @files;

    opendir(my $dh, $dir) or die "Cannot open directory $dir: $!";
    while (my $file = readdir($dh)) {
        # Skip non-files
        next unless -f File::Spec->catfile($dir, $file);
        # Only include files ending with .regenie.gz
        push @files, $file if $file =~ /\.regenie.significant.tsv$/;
    }
    closedir($dh);

    return @files;
}


sub get_ClumpedSNVs_files {
    my ($dir) = @_;
    my @files;

    opendir(my $dh, $dir) or die "Cannot open directory $dir: $!";
    while (my $file = readdir($dh)) {
        # Skip non-files
        next unless -f File::Spec->catfile($dir, $file);
        # Only include files ending with .regenie.gz
        push @files, $file if $file =~ /ClumpedSNVs.txt$/;
    }
    closedir($dh);

    return @files;
}
# Subroutine to construct EUR satisfy hash from gzipped file
# Params:
#   $input_file => EUR GWAS gzipped file
# Returns:
#   Hashref keyed by RSID with variant info
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
sub construct_EUR_satisfy_hash_gz {
    my ($input_file) = @_;
    my %SNVID_Info=();
    my $iNumberofLines=0;
    my $fh = IO::Uncompress::Gunzip->new($input_file)
        or die "Cannot open $input_file: $GunzipError\n";

    while (<$fh>) {
	next if($_ =~/^23/);
        chomp;
	if($iNumberofLines ==0)
	{
            $iNumberofLines++;
	    next;
	}
        my @cols = split( /\s++/,$_,7);
        # Skip INDELs
        next if (length($cols[3]) > 1 || length($cols[4]) > 1);
       next if($cols[0] ==6 && $cols[1]>25000000 && $cols[1]<34000000);
       	# Keep only MAF > 0.05
        next if ($cols[5] < 0.01 || $cols[5]>0.99);
	@cols = split( /\s++/,$_);
	my $p_value = 10 ** ((-1)*$cols[11]);
	my $Info="$cols[0]\t$cols[1]\t$cols[3]\t$cols[4]\t$cols[5]\t";
	$Info.="$cols[6]\t$cols[8]\t$cols[9]\t$p_value";
	$SNVID_Info{$cols[2]}=$Info;
    }
    $fh->close();
    my $Size=keys %SNVID_Info;
    print "Total EUR size: $Size\n";
    return (\%SNVID_Info);
}

sub construct_AFR_satisfy_hash_gz {
    my ($input_file,$EURSNVID_InfoRef) = @_;
    my %EURSNVID_Info=%{$EURSNVID_InfoRef};
    #my $SampleSize=keys %EURSNVID_Info;
    #print "EUR SampleSize:$SampleSize\n";
    my $outliers=0;
    my $iNumberofLines=0;
    my $fh = IO::Uncompress::Gunzip->new($input_file)
        or die "Cannot open $input_file: $GunzipError\n";
    while (<$fh>) {
	next if($_ =~ /^23/);    
        chomp;
	 if($iNumberofLines == 0)
	 {
             $iNumberofLines++;
	     next;
	 }
        my @cols = split( /\s++/,$_,7);
	next if(!defined($EURSNVID_Info{$cols[2]}));
	next if($cols[0] ==6 && $cols[1]>25000000 && $cols[1]<34000000);
	@cols = split(/\s++/,$_);
        # Skip INDELs
        next if (length($cols[3]) > 1 || length($cols[4]) > 1);
        
     	# Keep only MAF > 0.05
        next if ($cols[5] < 0.01 || $cols[5]>0.99);
	my $p_value = 10 ** ((-1)*$cols[11]);
	 my $Info="$cols[0]\t$cols[1]\t$cols[3]\t$cols[4]\t$cols[6]\t$cols[8]";
	$Info.="\t$cols[9]\t$p_value";
	 my @EURInfolist=split(/\t/,$EURSNVID_Info{$cols[2]});
	 if(abs($cols[6] -$EURInfolist[4])>0.2)
	 { 
      		 
	   $outliers++;	 
	 }

	$iNumberofLines++;
    }
    $fh->close();
    print "Total AFR size: $iNumberofLines\n";
    print "outliers:$outliers\n";
    # return (\%SNVID_A1_Freq,\%SNVID_Info);
}	    


#rs73885319,chr22:36,265,860
# rs334, chr11:5,225,464, 
#AGD35K_AFR_Phecode$ more chr_trait_linecount.txt
#CHR     TRAIT   LINE_COUNT
#11      Phecode_AFRAFR_038.3    29
#Find the file name that the chr11 contains many significANT SNVs
	#HR     TRAIT   LINE_COUNT
	#3       LabvalueAFR_BASOAB      30
my %SelectedPhecodes=();
my $SpecifiedChr=22;##22
my @DirList=("AGD35K_AFR_Phecode","AGD35K_AFR_LabValue");
my @TypePrefix=("Phecode_AFRAFR_","LabvalueAFR_");
for(my $iD=0;$iD<@DirList;$iD++)
{
my $iLineNumber = 0;
my $CountFiles=0;
my $Chr_TraitFile="$DirList[$iD]/chr_trait_linecount.txt";
open(my $fh, "<", $Chr_TraitFile) or die "Cannot open $Chr_TraitFile: $!";
while (my $line = <$fh>) {
    chomp $line;        # remove newline
     my @cols=split(/\s++/,$line);
    if($line =~/^$SpecifiedChr/)
    {
      if($cols[2] >3)
      {      
        my @cols=split(/\s++/,$line);	    
	$cols[1] =~ s/\Q$TypePrefix[$iD]\E//ig;
	$cols[1] =~ s/Phecode\_AFRAFR\_//ig;
	$cols[1] =~ s/LabvalueAFR\_//ig;
	$SelectedPhecodes{$cols[1]}=1;
	#print "$cols[0]\t$cols[1]\n";
	$CountFiles++;
       }
     }
 }

 close($fh);
print "CountFiles:$CountFiles\n";
}
#Collect the common 
#For estimating genetic effect correlations across populations, restrict to SNPs that are common in both populations
#and have MAF differences < 0.2. This reduces noise and gives more stable correlation estimates.
#Trait CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTR
#512.9 11 4611959 chr11:4611959:T:C T C 0.0292411 24640 ADD 0.443581 0.0676789 42.9575 10.2523 NA
#AGD35KAFR_LabValue_DiseaseIDsLDC0.1.txt
#AGD35KAFR_PheCode_DiseaseIDsLDC0.1.txt

#AGD35K_AFR_Phecode_GWAS_ClumpSNVs,AGD35K_AFR_LabValue_GWAS_ClumpSNVs
my @Types=("Phecode","LabValue");
my %hPhenotype_ClumpedPos=();
my @ClumpingDirList=("AGD35K_AFR_Phecode_GWAS_ClumpSNVs","AGD35K_AFR_LabValue_GWAS_ClumpSNVs");
for(my $iD=0;$iD<@ClumpingDirList;$iD++)
{
  my $ClumpedSNVsGWASDir="$ClumpingDirList[$iD]";
  my @LDClumpFiles=get_ClumpedSNVs_files($ClumpedSNVsGWASDir);
  my $SizeofLDClumpedFiles=scalar(@LDClumpFiles);
  print "SizeofLDClumpedFiles:$SizeofLDClumpedFiles\n";
  for(my $ii=0;$ii<@LDClumpFiles;$ii++)
  {
   my $ClumpingFile="$ClumpingDirList[$iD]/$LDClumpFiles[$ii]";
    my $Trait=$LDClumpFiles[$ii];          
    $Trait =~ s/Phecode\_AFRAFR\_//ig;
    $Trait =~ s/LabvalueAFR\_//ig;
    $Trait =~ s/ClumpedSNVs\.txt//ig;
    open(my $fh, "<", $ClumpingFile) or die "Cannot open $ClumpingFile: $!";
     while (my $line = <$fh>) {
     chomp $line;        # remove newline
     my @cols=split(/\s++/,$line);
     if($cols[0] == $SpecifiedChr)
     {
	  print "pos: $cols[1]\n";
	  if(!defined($hPhenotype_ClumpedPos{$Trait}))
         {
           $hPhenotype_ClumpedPos{$Trait}=$cols[1];
	 }
         else
         {
          $hPhenotype_ClumpedPos{$Trait}.=":$cols[1]";
	 }	 
     }
    }
   close($fh);
  } 
}
print "Phenotype_ClumpedPos\n";
foreach my $key (keys %hPhenotype_ClumpedPos) {
	#print "$key\n";
}

#my $Rs334LDPruningFile="AGD35K_NoveltyGWASHits/AGD35K_RS334_LDPruningPosCollection.txt";
my $Rs334LDPruningFile="AGD35K_NoveltyGWASHits/AGD35K_apol1_LDPruningPosCollection.txt";
open(my $outF, ">", $Rs334LDPruningFile) or die "Cannot open $Rs334LDPruningFile: $!";
print $outF "Phenotype\tStart_Pos\tEnd_Pos\tClumpedPosList\tOverlappedStatus_1Mb\tOverlappedStatus_550kb\n";
my $iCount_500Kb=0;
my $iCount_1M=0;
for(my $iD=0;$iD<@DirList;$iD++)
{
  my $AFRBasedSignificantSNVsGWASDir="$DirList[$iD]";
  my @AFRBasedGWASFiles=get_regenie_SignificantSNVs_files($AFRBasedSignificantSNVsGWASDir);
  my $SizeofAFRSignificantFiles=scalar(@AFRBasedGWASFiles);
  print "SizeofAFRFiles:$SizeofAFRSignificantFiles\n";
#Phecode_AFRAFR_425.1.regenie.significant.tsv
#10 33923116 chr10:33923116:A:G A G 0.0228291 26326 ADD 0.756711 0.135812 31.0446 7.5983 NA
  my $SignficantSNVFilePrefix=$TypePrefix[$iD];
  my $SignficantSNVFilePostfix=".regenie.significant.tsv";
  for(my $iA=0;$iA<@AFRBasedGWASFiles;$iA++)
  {
    my $Phenotype= $AFRBasedGWASFiles[$iA];
    $Phenotype =~ s/\Q$SignficantSNVFilePrefix\E//ig;
    $Phenotype =~ s/\Q$SignficantSNVFilePostfix\E//ig;
    next if(!defined($SelectedPhecodes{$Phenotype}));
    my $SignificantFile="$AFRBasedSignificantSNVsGWASDir/$AFRBasedGWASFiles[$iA]";
    my @SelectedPosList=();
    open(my $fh, "<", $SignificantFile) or die "Cannot open $SignificantFile: $!";
     while (my $line = <$fh>) {
     chomp $line;        # remove newline
     my @cols=split(/\s++/,$line,4);
     if($line =~/^$SpecifiedChr/)
     {
	push @SelectedPosList,$cols[1];     
	# print "$cols[0]\t$cols[1]\t$cols[2]\n";
     }
    }
   close($fh);
   my @sortedPos = sort { $a <=> $b } @SelectedPosList;
   my $first = $sortedPos[0];          # smallest
   my $last  = $sortedPos[-1];         # largest
    my $Range1_1M=35265860;
   my $Range2_1M=37265860; 
    my $Range1_500K=35765860;
   my $Range2_500K=36765860;
   #5,225,464,rs334
   #  my $Range1_1M=4225464;
   # my $Range2_1M=6225464;
   # my $Range1_500K=4725464;
   # my $Range2_500K=5725464;  

   if(($first < $Range1_1M && $last> $Range1_1M)||($first > $Range1_1M && $first< $Range2_1M)) 
    {
      print $outF "$Phenotype\t$first\t$last\t";
     # print "$Phenotype\t$first\t$last\t";
     # print "$hPhenotype_ClumpedPos{$Phenotype}\n";
      print $outF "$hPhenotype_ClumpedPos{$Phenotype}\t";
     if(!defined($hPhenotype_ClumpedPos{$Phenotype}))
     {
       print "short: $Phenotype\n";
     }
     #if($first>4227002 && $last< 6227002)
    if(($first < $Range1_1M && $last> $Range1_1M)||($first > $Range1_1M && $first< $Range2_1M))
    {
      $iCount_1M++;	    
      print $outF "1\t";
    }
    else
    {
	 print $outF "0\t";
    }
    if(($first < $Range1_500K && $last> $Range1_500K)||($first > $Range1_500K && $first< $Range2_500K))
    {
      $iCount_500Kb++;	    
        print $outF "1\n";
    }
    else
    {
	  print $outF "0\n";
    }
   }
 }
}
close $outF;
print "$iCount_500Kb\t$iCount_1M\n";


