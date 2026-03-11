#Remove the gwas files with number of gwas hits ==0
#: <<'COMMENT'
#for f in *.regenie.significant.tsv; do
#    # Count lines with awk
#    if [ "$(awk 'END{print NR}' "$f")" -eq 0 ]; then
#        echo "Deleting empty file: $f"
#        rm "$f"
#    fi
#done
#COMMENT
#Only read files with at least 2 lines
#For each file, count the number of lines per chromosome
#Record the result as chr, trait, line_count into an output file
# Output file
#: <<'COMMENT'
#out="chr_trait_linecount.txt"
#echo -e "CHR\tTRAIT\tLINE_COUNT" > "$out"

#for f in *.regenie.significant.tsv; do
#    # Skip files with less than 2 lines
#    if [ "$(wc -l < "$f")" -lt 2 ]; then
#        continue
#    fi
#
#    # Extract trait name from filename (remove .txt)
#    trait="${f%.regenie.significant.tsv}"
#
    # Count lines per chromosome, only keep counts > 1
#    awk -v trait="$trait" '
    #        {chr_count[$1]++}            # assuming first column is chromosome
    #        END {
    #        for (chr in chr_count)
    #            if (chr_count[chr] > 1)
    #                print chr "\t" trait "\t" chr_count[chr]
    #    }
    #' "$f" >> "$out"
    #done
    #COMMENT
#Reads chr_trait_linecount.txt
#Extracts the 2nd column (TRAIT = regenie file prefix)
#Skips the header
#Downloads the corresponding regenie GWAS files
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);

=cut;
#The goal of the part is used to
#obtain the chromsome based independent lead significant SNVs 
#The used method is LD clumping 
#Reads a compressed .gz file
#Extracts the ID column and LOG10P column
#Converts LOG10P → P-value
#Outputs a simple tab-delimited file with ID and P

sub convert_log10p_to_p_sstFile {	
   my ($infile, $outfile,$SelectedChrsRef) = @_;
   my %SelectedChrsHash=%{$SelectedChrsRef};

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
    my @cols = split(/\s+/,$_,2);
    next if(!defined($SelectedChrsHash{$cols[0]}));
    chomp;
    @cols = split(/\s+/,$_);
    if($iLineNumber ==0)
    {
      $iLineNumber++;
      next;
    }
    # convert LOG10P -> P
    my $pval = 10 ** (-$cols[11]);
    $iLineNumber++;
    print $out "$cols[2]\t$pval\t$cols[0]\t$cols[1]\n";
   }
   close $fh;
   close $out;   
   print "Total used lines: $iLineNumber\n";
}

#my $InputDir="AGD35K_EUR_Phecode";
#my $outdir = "AGD35K_EUR_Phecode_GWAS";

#my $InputDir="AGD35K_EUR_LabValue";
#my $outdir = "AGD35K_EUR_LabValue_GWAS";
#AGD35K_AFR_LabValue
#my $InputDir="AGD35K_AFR_LabValue";
#my $outdir = "AGD35K_AFR_LabValue_GWAS";
#AGD35K_AFR_Phecode

my $InputDir="AGD35K_AFR_Phecode";
my $outdir = "AGD35K_AFR_Phecode_GWAS";


my $input  = "$InputDir/chr_trait_linecount.txt";

mkdir $outdir unless -d $outdir;

# --------------------------------------------------
# Step 1: read file and build prefix hash
# --------------------------------------------------
my %prefix_hash;
my %prefix_Chr_hash;
my $iLineNum=0;
open(my $fh, "<", $input) or die "Cannot open $input: $!";
while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    my ($chr, $trait, $count) = split /\t/;
    next unless $count > 1;
    $prefix_hash{$trait} = 1;
    $prefix_Chr_hash{$trait}{$chr}=1;
}
close $fh;
# --------------------------------------------------
# Step 2: loop over prefixes and download files
# --------------------------------------------------

#EUR_Phecode
#my $BASE_URL = "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/RegenieWorkflowV1/Phecode_EURProb05";
#EUR_Labvalue
#my $BASE_URL = "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/RegenieWorkflowV1/LabValue_EURProb05";
#AFR_Labvalue
#my $BASE_URL = "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/RegenieWorkflowV1/LabValue_AFRProb05";
#AFR_Phecode
my $BASE_URL = "gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/RegenieWorkflowV1/Phecode_AFRProb05";

foreach my $prefix (sort keys %prefix_hash) {

    my $file = "$prefix.regenie.gz";
    my $outf = "$outdir/$file";

    my $url = "$BASE_URL/$file";
    print "Downloading $file\n";

    system("gsutil -m cp $url $outdir/") == 0
        or warn "Failed to download $file\n";
}


foreach my $trait (sort keys %prefix_Chr_hash) {
    # collect all chromosomes for this trait
    my @chrs = sort keys %{ $prefix_Chr_hash{$trait} };
    my %chr_hash=(); 
    # store array reference in hash
      for(my $iChr=0;$iChr<@chrs;$iChr++)
      {
       $chr_hash{$chrs[$iChr]}=1;
      }
      #Reads a compressed .gz file
      #Extracts the ID column and LOG10P column
      #Converts LOG10P → P-value
      #Outputs a simple tab-delimited file with ID an d P
       my $InsstFile = "$outdir/${trait}.regenie.gz";
       my $ConvertedsstFile = "$outdir/${trait}.txt";
     print "InsstFile:$InsstFile\n";
     convert_log10p_to_p_sstFile($InsstFile, $ConvertedsstFile,\%chr_hash);

      foreach my $chr (sort keys %{ $prefix_Chr_hash{$trait} }) {
         my $SpecificChr=$chr;
	 if($chr ==23)
	 {
           $SpecificChr="X";
	 }
        print "Processing $trait chr$chr\n"; 
        # define files
	#Independent significant SNPs (r² < 0.6) 
        my $outpref_1  = $outdir."_Clump06/${trait}.chr${chr}IndependentSNVs";
	#AFR_WGS
	my $genotype_data = "AGD35K_AFRWGS/AGD35K_MAC100MR002QCed_AFR_chr${SpecificChr}";
	#EUR WGS
	#my $genotype_data = "AGD35K_EURWGS/AGD35K_MAC100MR002QCed_EUR_chr${chr}";
        # run plink
	system(
	       	"bash", "-c",
  		"./plink2 ".
  		"--pfile $genotype_data ".
		"--chr ${chr} ".
  		"--clump $ConvertedsstFile ".
 		 "--clump-p1 5e-8 ".
  		"--clump-p2 1e-5 ".
  		"--clump-r2 0.1 ".
  		"--clump-kb 1000 ".
  		"--out $outpref_1"
		) == 0 or warn "PLINK2 failed for $trait chr$chr.ld06\n";
		#Lead SNPs (r² < 0.1)
		#my $outpref_2  = $outdir."_Clump06/${trait}.chr${chr}LeadSNVs";
		#system(
		# "bash", "-c",
		# "./plink2 ".
		#"--pfile $genotype_data ".
		#"--chr ${chr} ".
		#"--clump $outpref_1.clumps ".
		# "--clump-p1 5e-8 ".
		#"--clump-p2 1e-5 ".
		#"--clump-r2 0.1 ".
		#"--clump-kb 1000 ".
		#"--out $outpref_2"
		#) == 0 or warn "PLINK2 failed for $trait chr$chr.leadSNVs\n";	
    }
}
=cut;
#The following part is used to merge all independent leaded clumped association 
#into one file by trait and population
## LabvalueAFR_MONOAB.chr1IndependentSNVs.clumps
##CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
#LabvalueAFR_MONOAB.chr4LeadSNVs.clumps
#CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
#4       75726621        chr4:75726621:A:G       2.66373e-08     0       0       0       0       0       0       .
#read all files with the two different postfix into a vector
##CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
#4       75726621        chr4:75726621:A:G       2.66373e-08     7       0       0      

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

#Collect the trait with one significant SNVs information
#popultaion      phecode TotalVariants   SignificantVariants     LambdaGC        SampleSize
#AFR     250.22  14733100        11      1.0498  26315
sub collectTraitWithonly1SignificantSNVs
{
#AGD35K_AFR_Phecode/MEGA_AFR_summary_lambda.tsv	
   my ($summaryLambdaFile)=@_;
  my %Trait_1snvid=();
  # my $summaryLambdaFile="AGD35K_AFR_Phecode/MEGA_AFR_summary_lambda.tsv";
  my $iLineNum=0;
  open(my $fh, "<", $summaryLambdaFile) or die "Cannot open $summaryLambdaFile: $!";
  while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    my @cols = split(/\s++/,$_);
    if($cols[3] ==1)
    { 
      $Trait_1snvid{$cols[1]}=1;	    
    }
  }
  close $fh;
  return (%Trait_1snvid);
}


# Directory containing files
#AFR_Labvalue
my $SpecificDir="AGD35K_AFR_LabValue";#Labvalue
my $FilePrefix="LabvalueAFR_";#Labvalue
my $summary_lambdaFile= "$SpecificDir/AGD35K_AFR_Labvaluesummary_lambda.tsv";

#EUR_LabValue
#my $SpecificDir="AGD35K_EUR_LabValue";#Labvalue
#my $FilePrefix="LabvalueEUR_";#Labvalue
#my $summary_lambdaFile= "$SpecificDir/AGD35K_EUR_Labvaluesummary_lambda.tsv";

#my $SpecificDir="AGD35K_AFR_Phecode";#Phecode
#my $FilePrefix="Phecode_AFRAFR_";#Phecode
#my $summary_lambdaFile= "$SpecificDir/AGD35K_AFR_summary_lambda.tsv"; 

#use EUR_Phecode
#my $SpecificDir="AGD35K_EUR_Phecode";#Phecode
#my $FilePrefix="Phecode_EUREUR_";#Phecode
#my $summary_lambdaFile= "$SpecificDir/AGD35K_EUR_summary_lambda.tsv";

my $dir = $SpecificDir."_GWAS_Clump";
my $ClumpedSNVDir=$dir."SNVs";
#
#Number of variant in raw GWAS file
#AGD35K_AFR_LabValue
#chr_trait_linecount.txt
#CHR     TRAIT   LINE_COUNT
#3       LabvalueAFR_BASOAB      30
my %Trait_Chr_Count_RawGWAS=();#Trait-chr based sum of signidficant SNV
my %Trait_SumCount_RawGWAS=();#Trait based sum of signidficant SNVs
my %Trait_SumCount_1chr_RawGWAS=();#For the part only 1 significant variat/per chr
my %Trait_Count1_chrList_RawGWAS=();#collect the chr list with only one significant SNVs for each trait
my $RawSignificantSNVChrSumFile="$SpecificDir/chr_trait_linecount.txt"; 
my $iLineNum=0;
 open(my $fh, "<", $RawSignificantSNVChrSumFile) or die "Cannot open $RawSignificantSNVChrSumFile: $!";
 while (<$fh>) {
    chomp;
    if($iLineNum ==0)
    {
      $iLineNum++;
      next;
    }
    my @cols = split(/\t/,$_);
    my $Trait=$cols[1];
    $Trait =~ s/^\Q$FilePrefix\E//ig;
    $Trait_Chr_Count_RawGWAS{$Trait}{$cols[0]}=$cols[2];
 }
 close $fh;

 for my $trait (sort keys %Trait_Chr_Count_RawGWAS)
 {  
    $Trait_SumCount_RawGWAS{$trait}=0;
    $Trait_SumCount_1chr_RawGWAS{$trait}=0;
    for my $chr (sort keys %{$Trait_Chr_Count_RawGWAS{$trait}}) 
     {
       $Trait_SumCount_RawGWAS{$trait}+=$Trait_Chr_Count_RawGWAS{$trait}{$chr};
       if($Trait_Chr_Count_RawGWAS{$trait}{$chr} ==1)
       {
        $Trait_SumCount_1chr_RawGWAS{$trait}+=1;
         if(!defined($Trait_Count1_chrList_RawGWAS{$trait}))
        {
	  $Trait_Count1_chrList_RawGWAS{$trait}="$chr";     
        }
        else
       {
	  $Trait_Count1_chrList_RawGWAS{$trait}.="\t$chr";    
       }
      }
    }
 }
 my $SizeTrait_Count1_chrList_RawGWAS=keys %Trait_Count1_chrList_RawGWAS;
 print "SizeTrait_Count1_chrList_RawGWAS:$SizeTrait_Count1_chrList_RawGWAS\n";
 #for my $trait (sort keys %Trait_Count1_chrList_RawGWAS)
 # {
 #  print "$trait\t$Trait_Count1_chrList_RawGWAS{$trait}\n";
 # }
# Postfixes
 #my @postfixes = ("IndependentSNVs.clumps","LeadSNVs.clumps");
 #my @clumpTypes = ("Independent","Lead");
 #Only consider IndependentSNVs R2>0.1
my @postfixes = ("IndependentSNVs.clumps");
my @clumpTypes = ("Independent");

my %Trait_Type_Chr_snvid=();
my %Trait_Chr_Type_snvCount=();
for(my $iT=0;$iT<@postfixes;$iT++)
{
   my @ClumpFiles=();	
   # Open directory
   opendir(my $dh, $dir) or die "Cannot open directory $dir: $!";
   # Loop through files
   while (my $file = readdir($dh)) {
       next unless -f File::Spec->catfile($dir, $file);  # skip non-files
       push @ClumpFiles, $file if $file =~ /\Q$postfixes[$iT]\E$/;
    }
   closedir($dh);
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

my %Trait_Type_SNVCount=();
for(my $iT=0;$iT<@clumpTypes;$iT++)
{  	
   for my $trait (sort keys %Trait_Type_Chr_snvid) {
     $Trait_Type_SNVCount{$trait}{$clumpTypes[$iT]}=$Trait_SumCount_1chr_RawGWAS{$trait};	   
        for my $chr (sort keys %{ $Trait_Type_Chr_snvid{$trait}{$clumpTypes[$iT]} }) {
            my $snv_ref = $Trait_Type_Chr_snvid{$trait}{$clumpTypes[$iT]}{$chr};
            my @SNVIDs=split(/\t/,$snv_ref);
            my $count =scalar(@SNVIDs) ;
	    $Trait_Type_SNVCount{$trait}{$clumpTypes[$iT]}+=$count;
        }
	#  print "$trait\t$clumpTypes[$iT]\t$Trait_Type_SNVCount{$trait}{$clumpTypes[$iT]}\n";
  }
}
#$summary_lambdaFile
#collect the Trait with only 1 significant SNVs across WG data
 my %Trait_1snvid=collectTraitWithonly1SignificantSNVs($summary_lambdaFile);

my $ClumpedSNVs_Trait="$ClumpedSNVDir/Clumped_Raw_GWAS_SignificantSNVs.txt";
open(my $out, ">", "$ClumpedSNVs_Trait") or die "Cannot open $ClumpedSNVs_Trait: $!";
print $out "Trait\tRaw";
for(my $iT=0;$iT<@clumpTypes;$iT++)
{
   print $out "\t$clumpTypes[$iT]";	
}
print $out "\n";
my @SumSNVsCount=(0,0,0);
for my $trait (keys %Trait_Type_SNVCount) {
    print $out "$trait\t$Trait_SumCount_RawGWAS{$trait}";
    $SumSNVsCount[0]+=$Trait_SumCount_RawGWAS{$trait};	
    for(my $iT=0;$iT<@clumpTypes;$iT++) {
	 print $out "\t$Trait_Type_SNVCount{$trait}{$clumpTypes[$iT]}";
	 $SumSNVsCount[$iT+1]+=$Trait_Type_SNVCount{$trait}{$clumpTypes[$iT]};
    }
     print $out "\n";
     #   print "$Trait_Count1_chrList_RawGWAS{$trait}\n";  
}

my $Size_Withonly1snvid=keys %Trait_1snvid;
print "Size_Withonly1snvid:$Size_Withonly1snvid\n";
for my $trait (keys %Trait_1snvid) {   	
   print $out "$trait\t1\t1\n";
   $SumSNVsCount[0]+=1;
   $SumSNVsCount[1]+=1;
   #   $SumSNVsCount[2]+=1;
}

#Count the trait for all their SNV list only one
for my $trait (keys %Trait_Count1_chrList_RawGWAS) {
   next if(defined($Trait_Type_SNVCount{$trait}));
   my @ChrList=split(/\t/,$Trait_Count1_chrList_RawGWAS{$trait});
   my $iSumSNVList=scalar(@ChrList);
   $SumSNVsCount[0]+=$iSumSNVList;
   $SumSNVsCount[1]+=$iSumSNVList;
   #  $SumSNVsCount[2]+=$iSumSNVList;
   print "$trait\t$iSumSNVList\t$iSumSNVList\n";
   print $out "$trait\t$iSumSNVList\t$iSumSNVList\n";
}



print $out "Sum\t$SumSNVsCount[0]\t$SumSNVsCount[1]\n";
close $out;
 
 my @SignifcantSNVsFiles=();
 my $Postfix=".regenie.significant.tsv";

   # Open directory
   opendir(my $dh, $SpecificDir) or die "Cannot open directory $SpecificDir: $!";
   # Loop through files
   while (my $file = readdir($dh)) {
       next unless -f File::Spec->catfile($SpecificDir, $file);  # skip non-files
       push @SignifcantSNVsFiles, $file if $file =~ /\Q$Postfix\E$/;
    }
    my $iSignificantSNVsFiles=scalar(@SignifcantSNVsFiles);
    print "iSignificantSNVsFiles:$iSignificantSNVsFiles\n";
  for(my $jj=0;$jj<@SignifcantSNVsFiles;$jj++)
  {
     my $Trait=$SignifcantSNVsFiles[$jj];
     $Trait =~ s/\Q$Postfix\E$//ig;
     $Trait =~ s/^\Q$FilePrefix\E//ig;
     if(defined($Trait_1snvid{$Trait}))
     {
       my $fullPath="$SpecificDir/$SignifcantSNVsFiles[$jj]";
       my $outPath="$ClumpedSNVDir/$FilePrefix".$Trait."ClumpedSNVs.txt";
       #1 158863706 chr1:158863706:GCTT:G GCTT G 0.0139629 7878 ADD 0.36276 0.0647703 31.368 7.67066 NA
       open(my $out, ">", "$outPath") or die "Cannot open $outPath: $!";
       #print $out "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA\n";
       open(my $fh, "<", $fullPath) or die "Cannot open $fullPath: $!";
       while (<$fh>) {
         print $out "$_";
        }
        close $fh;
        close $out;
      }
     else
      {
        my %hChrList_SingleSignificantSNVs=();
        if(defined($Trait_Count1_chrList_RawGWAS{$Trait}))
        {
          my @ChrList_SingleSignificantSNVs=split(/\t/,$Trait_Count1_chrList_RawGWAS{$Trait});
          for(my $kk=0;$kk<@ChrList_SingleSignificantSNVs;$kk++)
          {
            $hChrList_SingleSignificantSNVs{$ChrList_SingleSignificantSNVs[$kk]}=1;	      
          }
        }
        my %ClumpedSNVIDList_Trait=();
        for my $chr (sort keys %{$Trait_Type_Chr_snvid{$Trait}{$clumpTypes[0]}})
        {
	  my @ChrBasedSNVIDList=split(/\t/,$Trait_Type_Chr_snvid{$Trait}{$clumpTypes[0]}{$chr});     
          for(my $kk=0;$kk<@ChrBasedSNVIDList;$kk++)
          {
            $ClumpedSNVIDList_Trait{$ChrBasedSNVIDList[$kk]}=1;
          }
         }
       my $fullPath="$SpecificDir/$SignifcantSNVsFiles[$jj]";
       my $outPath="$ClumpedSNVDir/$FilePrefix".$Trait."ClumpedSNVs.txt";
       #1 158863706 chr1:158863706:GCTT:G GCTT G 0.0139629 7878 ADD 0.36276 0.0647703 31.368 7.67066 NA
       open(my $out, ">", "$outPath") or die "Cannot open $outPath: $!";
       # print $out "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA\n";
       open(my $fh, "<", $fullPath) or die "Cannot open $fullPath: $!";
       while (<$fh>) {
          chomp;
          my @cols=split(/\s++/,$_);
          if(defined($hChrList_SingleSignificantSNVs{$cols[0]}))
          {
            print $out "$_\n";
	    next;
          }
          if(defined($ClumpedSNVIDList_Trait{$cols[2]}))
	  {
           print $out "$_\n";		 
	  }
         }
         close $fh;
         close $out;
      }
  }

 
#Merge all GWAS lead SNVs into one set

 my %IndependentPhecodeList=();
 my $AGD35KIndependentPhecodeListFile="IndependentPhecodeList_AGD35K.txt";
 open( $fh, "<", $AGD35KIndependentPhecodeListFile) or die "Cannot open $AGD35KIndependentPhecodeListFile: $!";
 while (<$fh>) {
   chomp;
   my @cols=split(/\s++/,$_);
   $cols[0] =~ s/^0+//;
   $cols[0] = $cols[0]+0;
   $IndependentPhecodeList{$cols[0]}=1;
 }
 close $fh;
my @Pops=("AFR","EUR"); 
my %IndPhecode_Count=();
my %Labvalues_Count=();
for(my $iPop=0;$iPop<@Pops;$iPop++)
{
  my $iFile="AGD35K_$Pops[$iPop]_Phecode_GWAS_ClumpSNVs/Clumped_Raw_GWAS_SignificantSNVs.txt";
  #Trait   Raw     Independent
  #288.2   905     208 
  my $iLineNumber=0;
  my $iNonZeroCount=0;
   open(my $fh, "<", $iFile) or die "Cannot open $iFile: $!";
   while (<$fh>) {
   chomp;
   if($iLineNumber == 0)
   {
     $iLineNumber++;
     next;
   }
   my @cols=split(/\s++/,$_,2);
   next if($cols[0] =~/^Sum/);
   $cols[0] =~ s/^0+//;
   $cols[0] = $cols[0]+0;
    $iNonZeroCount++;
   $IndPhecode_Count{$cols[0]}{$Pops[$iPop]}=$cols[1];
  }
  close $fh;
  print "$Pops[$iPop]\t$iNonZeroCount\n";
}
#AGD35K_AFR_LabValue_GWAS_ClumpSNVs
#Trait   Raw     Independent     Lead
#MCHC    1088    275     75
for(my $iPop=0;$iPop<@Pops;$iPop++)
{
  my $iFile="AGD35K_$Pops[$iPop]_LabValue_GWAS_ClumpSNVs/Clumped_Raw_GWAS_SignificantSNVs.txt";
  #Trait   Raw     Independent     Lead
  #288.2   905     208     45
  my $iLineNumber=0;
   open(my $fh, "<", $iFile) or die "Cannot open $iFile: $!";
   while (<$fh>) {
   chomp;
   if($iLineNumber == 0)
   {
     $iLineNumber++;
     next;
   }   
   my @cols=split(/\s++/,$_,2);
   next if($cols[0] =~/^Sum/);
   $Labvalues_Count{$cols[0]}{$Pops[$iPop]}=$cols[1];
  }
  close $fh;
}

my $outPath="AGD35K_NoveltyGWASHits/AGD35K_SignificantSNVsCountSummary.txt";
open(my $out, ">", "$outPath") or die "Cannot open $outPath: $!";
my @Populations=("AFR","EUR");
print $out "Phenotype\tiTotalSNVs_AFR\tiIndepenfentSNVs_AFR\t";
print $out "iTotalSNVs_EUR\tiIndepenfentSNVs_EUR\n";
for my $phecode (sort { $a <=> $b }  keys %IndependentPhecodeList) {
    print $out "$phecode";
    for(my $iPop=0;$iPop<@Populations;$iPop++)
    {    
	if(defined($IndPhecode_Count{$phecode}{$Pops[$iPop]}))
	{
          print $out "\t$IndPhecode_Count{$phecode}{$Pops[$iPop]}"; 		
        }
	else
	{
         print $out "\t0\t0";
	}
    }	 
    print $out "\n";
}
for my $lab (sort  keys %Labvalues_Count) {
    print $out "$lab";
    for(my $iPop=0;$iPop<@Populations;$iPop++)
    {
        if(defined($Labvalues_Count{$lab}{$Pops[$iPop]}))
        {
          print $out "\t$Labvalues_Count{$lab}{$Pops[$iPop]}";
        }
        else
        {
         print $out "\t0\t0";
        }
    }
    print $out "\n";
}
close $out;


#start from 0308_2026
#Use Known EFO/traitfromSource to obtain novel SNVs information 
#Detect the novel signficant SNVs for all LabValue based GWAS hits
#AGD35K_AFR_LabValue_GWAS_ClumpSNVs/LabvalueAFR_EOSIAClumpedSNVs.txt
#3 103831674 chr3:103831674:C:A C A 0.0224437 6171 ADD 0.0510826 0.00933651 29.9349 7.34988 NA
#Now using all independent SNVs

my $Population="AFR";#AFR,EUR
my $TraitType="LabValue";
my $TraitTypeAlias="LabValue";
#my $TraitType="PheCode";
#my $TraitTypeAlias="Phecode";
my $GroupHits_Type="DiseaseIDs";#"TraitFromSource";#"DiseaseIDs";#"TraitFromSource";#"DiseaseIDs";#"TraitFromSource"; 
my $EURBasedGWASDir="AGD35K_EUR_LabValue_GWAS";
my $FilePrefix="Labvalue$Population"."_";
my $GWASfilePrefix="$EURBasedGWASDir/LabvalueEUR_";
if($Population eq "AFR" && $TraitType eq "PheCode")
{
   $FilePrefix="Phecode_AFRAFR_";
   $EURBasedGWASDir="AGD35K_EUR_Phecode_GWAS_Full"; 
   $GWASfilePrefix="$EURBasedGWASDir/Phecode_EUREUR_";
}
elsif($Population eq "EUR" && $TraitType eq "PheCode")
{
   $FilePrefix="Phecode_EUREUR_";
   $EURBasedGWASDir= "AGD35K_EUR_Phecode_GWAS_Full";
   $GWASfilePrefix="$EURBasedGWASDir/Phecode_EUREUR_";
}
my $LDClumpingThreshold="0.1";
#Known_LabValue_GWASHits_OTG_GWASCatabyDiseaseIDs_0204_2026.tsv
#Known_LabValue_GWASHits_OTG_GWASCatabyTraitFromSource_0204_2026.tsv

my $IndependentLeadSNVsDir="AGD35K_$Population"."_$TraitTypeAlias"."_GWAS_ClumpSNVs";

my $OTGGWASHitsFile="AGD35K_NoveltyGWASHits/Known_$TraitType".
        "_GWASHits_OTG_GWASCataby$GroupHits_Type"."_0204_2026.tsv";
my $GWASCatalogHitsFile="AGD35K_NoveltyGWASHits/Known_$TraitType".
        "_GWASHits_GWASCataby$GroupHits_Type"."_0206_2026.tsv";

my $OutNoveltyGWASResultFile="AGD35K_NoveltyGWASHits/AGD35K"."$Population"."_$TraitType";
$OutNoveltyGWASResultFile.="_$GroupHits_Type"."LDC$LDClumpingThreshold.txt";
open(my $OUTN, ">", $OutNoveltyGWASResultFile) or die "Cannot open $OutNoveltyGWASResultFile: $!";
print $OUTN "Trait CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTR\n";
close $OUTN;

#1 157834858 chr1:157834858:G:A G A 0.0746193 7880 ADD 0.00349444 0.000598392 34.1024 8.28161 NA
sub ReadGWASclumpedSNVListFile
{
  my ($SNVFile,$Trait,$TraitChr_PosListRef,
         $TraitChr_SNVIDListRef)=@_;
   my %TraitChr_PosList=%{$TraitChr_PosListRef};
   my %TraitChr_SNVIDList=%{$TraitChr_SNVIDListRef};
   open(my $fh, "<", $SNVFile) or die "Cannot open $SNVFile: $!";
   while (<$fh>) {
      chomp;
      my @cols=split(/\s++/,$_,4);
      if(!defined($TraitChr_PosList{$Trait}{$cols[0]}))
      {
        $TraitChr_PosList{$Trait}{$cols[0]}=$cols[1];
        $TraitChr_SNVIDList{$Trait}{$cols[0]}=$cols[2];
     }
     else
      {
       $TraitChr_PosList{$Trait}{$cols[0]}.="\t$cols[1]";
       $TraitChr_SNVIDList{$Trait}{$cols[0]}.="\t$cols[2]";
      } 
   }
  close $fh;
  return (\%TraitChr_PosList,\%TraitChr_SNVIDList);
}

#$LabValueBasedGWASHitsFile
#traitFromSource diseaseIds.x    variantId       chromosome      position        LabTrait_AGD35K name_AGD35K       diseaseIds.y    MatchMethod
sub ReadKnownGWASHits_ByTraitFromSourceFile
{
   my ($OTGGWASHitsFile)=@_;
   my %TraitChr_Pos_KnownHits=();
   my $iLineNumber=0;
   open(my $fh, "<", $OTGGWASHitsFile) or die "Cannot open 
      $OTGGWASHitsFile: $!";
   while (<$fh>) {
      chomp;
      if($iLineNumber ==0)
      {
	$iLineNumber++;      
        next;
      }
      my @cols=split(/\t/,$_);
      #next if(scalar(@cols)<3); 
      if(!defined($TraitChr_Pos_KnownHits{$cols[0]}{$cols[1]}))
      {
        $TraitChr_Pos_KnownHits{$cols[0]}{$cols[1]}=$cols[2];
      }
      else
      {
       $TraitChr_Pos_KnownHits{$cols[0]}{$cols[1]}.="\t$cols[2]";
       }
      $iLineNumber++; 
    }
  close $fh;
  return (%TraitChr_Pos_KnownHits);
}


sub JudgeGWASHitsNoveltyBasedon_1McM
{
  my ($trait,$Chromosome,$PosList_AGD35KStr,$VariantIDList_AGD35KStr,
	  $KnownPosList_OTGStr,$iOverlappedLocus,
	  $iNovelLocus,$NovelSNVListRef)=@_;
  my @PosList_AGD35K=split(/\t/,$PosList_AGD35KStr);
  my @VariantIDList_AGD35K=split(/\t/,$VariantIDList_AGD35KStr);
  my @KnownPosList_OTG=split(/\t/,$KnownPosList_OTGStr);
  my @NovelSNVList=@{$NovelSNVListRef};
  for(my $ii=0;$ii<@PosList_AGD35K;$ii++)
  {
     my $iOverlapStatus=0;	  
     for(my $jj=0;$jj<@KnownPosList_OTG;$jj++)
     {
	 my $AbsDistance=abs($KnownPosList_OTG[$jj]-$PosList_AGD35K[$ii]);    
         if($AbsDistance<1000000)
         {
	   $iOverlapStatus=1;
	    last;
         }
     }
    if($iOverlapStatus == 0)
    {
      my @Result=($trait,$Chromosome,$PosList_AGD35K[$ii],$VariantIDList_AGD35K[$ii]);	
      my $ResultStr=join("\t",@Result); 
      push @NovelSNVList,$ResultStr; 
      $iNovelLocus++;  
    }
    else
    {
      $iOverlappedLocus++;   
    }
  }
  return ($iOverlappedLocus,$iNovelLocus,\@NovelSNVList);
}


#Collect novel GWAS association hits
sub CompleteNovelGWASHitsCollection
{
  my ($trait,$Chromosome,$PosList_AGD35KStr,$VariantIDList_AGD35KStr,
	  $iNovelLocus,$NovelSNVListRef)=@_;
  my @PosList_AGD35K=split(/\t/,$PosList_AGD35KStr);
  my @VariantIDList_AGD35K=split(/\t/,$VariantIDList_AGD35KStr);
  my @NovelSNVList=@{$NovelSNVListRef};
  for(my $ii=0;$ii<@PosList_AGD35K;$ii++)
  {
      my @Result=($trait,$Chromosome,$PosList_AGD35K[$ii],$VariantIDList_AGD35K[$ii]);
      my $ResultStr=join("\t",@Result);
      push @NovelSNVList,$ResultStr;
      $iNovelLocus++;
  }
  return ($iNovelLocus,\@NovelSNVList);
}
#Detect Novelity GWASHitsbasedon Known GWAS hits
sub DetectNovelityAssoctionBasedOnKnownGWASsHits
{
   my ($TraitChr_PosListFullRef,$TraitChr_SNVIDListFullRef,
	   $TraitChr_Pos_KnownHitsRef)=@_;
   my %TraitChr_PosListFull=%{$TraitChr_PosListFullRef};
   my %TraitChr_Pos_KnownHits=%{$TraitChr_Pos_KnownHitsRef}; 
   my %TraitChr_SNVIDListFull=%{$TraitChr_SNVIDListFullRef};
   my $iOverlappedLocus=0;
   my $iNovelLocus=0;
   my $iTotalSignificantSNVs=0;
   my @NovelSNVsListInfo=();
   my $NovelSNVsListInfoRef;
   foreach my $trait (sort keys %TraitChr_PosListFull) {
       foreach my $chr (sort {$a<=>$b} keys %{ $TraitChr_PosListFull{$trait} }) {
             my @PosList_AGD35K=split(/\t/,$TraitChr_PosListFull{$trait}{$chr});
               $iTotalSignificantSNVs+=scalar(@PosList_AGD35K);
                my $SpecificChr=$chr;
                if($chr ==23)
                {
                 $SpecificChr='X';
                }
		if(defined($TraitChr_Pos_KnownHits{$trait}{$SpecificChr}))
		{                  
		     ($iOverlappedLocus,$iNovelLocus,$NovelSNVsListInfoRef)= 
		   JudgeGWASHitsNoveltyBasedon_1McM($trait,$chr,
			$TraitChr_PosListFull{$trait}{$chr},$TraitChr_SNVIDListFull{$trait}{$chr},
                        $TraitChr_Pos_KnownHits{$trait}{$SpecificChr},$iOverlappedLocus,
                        $iNovelLocus,\@NovelSNVsListInfo);
                }
                else
                {
                   ($iNovelLocus,$NovelSNVsListInfoRef)=CompleteNovelGWASHitsCollection($trait,$chr,
                           $TraitChr_PosListFull{$trait}{$chr},
                           $TraitChr_SNVIDListFull{$trait}{$chr},
                           $iNovelLocus,\@NovelSNVsListInfo);
                }
                @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
           }
    }
   return ($iOverlappedLocus,$iNovelLocus, $iTotalSignificantSNVs,\@NovelSNVsListInfo);
}


sub DetectNovelityAssoctionBasedOnKnownGWASsHits_2ndTime
{
   my ($TraitChr_Pos_KnownHits_GWASCataRef,$NovelSNVsListInfoRef)=@_;
   my %TraitChr_Pos_KnownHits_GWASCata=%{$TraitChr_Pos_KnownHits_GWASCataRef};
   my @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
   my %TraitChr_Pos=();
   my %TraitChr_SNVIDlist=();
    for(my $ii=0;$ii<@NovelSNVsListInfo;$ii++)
    {
      my @cols=split(/\t/,$NovelSNVsListInfo[$ii]);
      if(!defined($TraitChr_Pos{$cols[0]}{$cols[1]}))
      {
         $TraitChr_Pos{$cols[0]}{$cols[1]}=$cols[2];
	 $TraitChr_SNVIDlist{$cols[0]}{$cols[1]}=$cols[3];
      }
      else
      {
          $TraitChr_Pos{$cols[0]}{$cols[1]}.="\t$cols[2]";
          $TraitChr_SNVIDlist{$cols[0]}{$cols[1]}.="\t$cols[3]";
      }
    }
    return (\%TraitChr_Pos,\%TraitChr_SNVIDlist)
 }


#Based on Novelty SNVs list,write out all significnat SNV GWAS signala
#whose format is same as GWAS regenie format,the additional column list
#is add trait column
sub WriteOutSignificantFile
{
   my ($Significant_SST_File,$Trait,
	  $TraitSNV_noveltyRef,$OutFile)=@_;
   my %TraitSNVID_Novelty=%{$TraitSNV_noveltyRef};
   open(my $OUTN, ">>", $OutFile) or die "Cannot open $OutFile: $!";
   open(my $fh, "<", $Significant_SST_File) or die "Cannot open $Significant_SST_File: $!";
   while (<$fh>) {
      chomp;
      my @cols=split(/\s++/,$_,4);
      if(defined($TraitSNVID_Novelty{$Trait}{$cols[2]}))
      {
	print $OUTN "$Trait $_\n";      
      }
   }
  close $fh;
  close $OUTN;
}
#Detect the novel signficant SNVs for all LabValue based GWAS hits
#AGD35K_AFR_LabValue_GWAS_ClumpSNVs/LabvalueAFR_EOSIAClumpedSNVs.txt
#3 103831674 chr3:103831674:C:A C A 0.0224437 6171 ADD 0.0510826 0.00933651 29.9349 7.34988 NA

my $FilePostFix="ClumpedSNVs.txt";
   my @LeadSNVsFiles=();
   # Open directory
   opendir(my $dh, $IndependentLeadSNVsDir) or die "Cannot open directory $IndependentLeadSNVsDir: $!";
   # Loop through files
   while (my $file = readdir($dh)) {
       next unless -f File::Spec->catfile($IndependentLeadSNVsDir, $file);  # skip non-files
       push @LeadSNVsFiles, $file if $file =~ /\Q$FilePostFix\E$/;
    }
   closedir($dh);
   my $NumberofFiles=scalar(@LeadSNVsFiles);
   print "IndependentLeadSNVsDir:$IndependentLeadSNVsDir: File Count:$NumberofFiles\n";
   my %TraitChr_PosListFull=();
   my %TraitChr_SNVIDListFull=();
   for(my $jj=0;$jj<@LeadSNVsFiles;$jj++)
   {
     my $Trait=$LeadSNVsFiles[$jj];
     my $fullPath="$IndependentLeadSNVsDir/$LeadSNVsFiles[$jj]";
     #  print "fullPath:$fullPath\n";
      $Trait =~ s/\Q$FilePostFix\E$//ig;
      $Trait =~ s/^\Q$FilePrefix\E//ig;
      # print "$Trait\n";
      my ($TraitChr_PosListFullRef,$TraitChr_SNVIDListFullRef)= 
         ReadGWASclumpedSNVListFile($fullPath,$Trait,
         \%TraitChr_PosListFull,\%TraitChr_SNVIDListFull);
      %TraitChr_PosListFull=%{$TraitChr_PosListFullRef};	
      %TraitChr_SNVIDListFull=%{$TraitChr_SNVIDListFullRef}; 	
    }
 #Obtain All Trait_Chr INDependent List
 print "TraitChr_SNVIDListFull count\n";
 my $iTtoalTrait_ChrList=0;
 foreach my $trait (keys %TraitChr_SNVIDListFull) {
	 foreach my $chr (keys %{$TraitChr_SNVIDListFull{$trait}}) {
		 #print "$TraitChr_SNVIDListFull{$trait}{$chr}\n";
		 my @IndependentList=split(/\t/,$TraitChr_SNVIDListFull{$trait}{$chr});
		 my $iCountLoci=scalar(@IndependentList);
		 $iTtoalTrait_ChrList+=$iCountLoci;
	   }
}
print "Number of total TraitChr based independent Loci: $iTtoalTrait_ChrList\n";

    

  my %TraitChr_Pos_KnownHits=ReadKnownGWASHits_ByTraitFromSourceFile($OTGGWASHitsFile);

  my ($iOverlappedLocus,$iNovelLocus, $iTotalSignificantSNVs,$NovelSNVsListInfoRef)=
  DetectNovelityAssoctionBasedOnKnownGWASsHits(\%TraitChr_PosListFull,
	  \%TraitChr_SNVIDListFull,\%TraitChr_Pos_KnownHits);
  my %TraitChr_Pos_KnownHits_GWASCata=ReadKnownGWASHits_ByTraitFromSourceFile($GWASCatalogHitsFile); 
  my ($TraitChr_PosRef,$TraitChr_SNVIDlistRef)= DetectNovelityAssoctionBasedOnKnownGWASsHits_2ndTime(\%TraitChr_Pos_KnownHits_GWASCata,$NovelSNVsListInfoRef); 
  ($iOverlappedLocus,$iNovelLocus, $iTotalSignificantSNVs,$NovelSNVsListInfoRef)=
   DetectNovelityAssoctionBasedOnKnownGWASsHits($TraitChr_PosRef,$TraitChr_SNVIDlistRef,
	   \%TraitChr_Pos_KnownHits_GWASCata); 
  my @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
  print "Overlap: $iOverlappedLocus\tNovel: $iNovelLocus\tTotal:$iTotalSignificantSNVs\n";
#Obtain Trait_SNVID list obtain these list SNVs effect size
    my %TraitSNVID_Novelty=();
    my %Trait_ChrPosNovelty=();
    my %Trait_ChrPos_SNVIDList_Novel=();
    #  my %TraitChr_SNVIDList_Novelty=();    
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
    }
    #collect all novelty SNVs into a file
   for(my $jj=0;$jj<@LeadSNVsFiles;$jj++)
   {
     my $Trait=$LeadSNVsFiles[$jj];
     my $fullPath="$IndependentLeadSNVsDir/$LeadSNVsFiles[$jj]";
      $Trait =~ s/\Q$FilePostFix\E$//ig;
      $Trait =~ s/^\Q$FilePrefix\E//ig;
      WriteOutSignificantFile($fullPath,$Trait,\%TraitSNVID_Novelty,$OutNoveltyGWASResultFile);
    }  
    #gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/RegenieWorkflowV1/Phecode_EURProb05/*manhattan.png  
    # for my $trait (keys %Trait_Novelty) {
    # my $cmd="gsutil cp ";
    # print "$key\t$value\n";
    #}

 #Based on AFR filtered phecode/lab value based novelity SNVs to find the corresponding EUR subset
 #AGD35K_EUR_Phecode_GWAS_Full/Phecode_EUREUR_585.4.regenie.gz
 #AGD35K_EUR_LabValue_GWAS/LabvalueEUR_BASOAB.regenie.gz
 # $OutNoveltyGWASResultFile, trait,chr,position to collect the information
 #
 #my $EURBasedGWASDir="AGD35K_EUR_Phecode_GWAS_Full";
 my %TraitChrPos=();
 my $OutFile= "$OutNoveltyGWASResultFile.EURPairList.txt";
 my $OutFile2= "$OutNoveltyGWASResultFile.EURPairList_Unavailable.txt";

 open(my $OUTN, ">", $OutFile) or die "Cannot open $OutFile: $!";
 open(my $OUTN2, ">", $OutFile2) or die "Cannot open $OutFile2: $!";
 for my $Trait (keys %Trait_ChrPosNovelty) {
     my %ChrPos_Novel=();
     my @ChrPosList=split(/\t/,$Trait_ChrPosNovelty{$Trait});
    my $iPredictedSNVs=0;
      for(my $ii=0;$ii<@ChrPosList;$ii++)
     {
      my ($chr,$pos)=split(/,/,$ChrPosList[$ii]);
      $ChrPos_Novel{$chr}{$pos}=$Trait_ChrPos_SNVIDList_Novel{$Trait}{$chr}{$pos};
      $iPredictedSNVs++;
     }
      my %ObtainedChrPos=();
      my $iObtainedSNVs=0;
      my $iLineNumber=0;
      my $GWASfile="$GWASfilePrefix$Trait.regenie.gz";
      open(my $fh, "-|", "zcat $GWASfile") or die "Cannot open $GWASfile: $!";
     while (<$fh>) {
      chomp;
      my @cols = split(/\s+/,$_,4);
      if($iLineNumber ==0)
      {
       $iLineNumber++;
       next;
      }
      if(defined($ChrPos_Novel{$cols[0]}{$cols[1]}))
      {
	$ObtainedChrPos{$cols[0]}{$cols[1]}=1;      
	$iObtainedSNVs++;      
       print $OUTN "$Trait\t$_\t$ChrPos_Novel{$cols[0]}{$cols[1]}\n";
      }     
     }
    close $fh;
    for my $chr (keys %ChrPos_Novel) {
       for my $pos (keys %{ $ChrPos_Novel{$chr} }) {
	if(!defined($ObtainedChrPos{$chr}{$pos}))
         {
          print $OUTN2 "$Trait\t$chr\t$pos\n";   		 
	 }	 
    }
}
    print "$Trait\tPredictedSNVs:$iPredictedSNVs\t"; 
    print  "Obtain:$iObtainedSNVs\n";
  }
  close $OUTN;
  close $OUTN2;
# AGD35K_AFR_LabValue/AGD35K_AFR_Labvaluesummary_lambda.tsv
# AGD35K_AFR_Phecode/AGD35K_AFR_summary_lambda.tsv
#AGD35K_EUR_LabValue/AGD35K_EUR_Labvaluesummary_lambda.tsv
#AGD35K_EUR_Phecode/AGD35K_EUR_summary_lambda.tsv 
 #File    TotalVariants   SignificantVariants     LambdaGC
 #LabvalueEUR_BASOAB      5109160 0       1.0024  1915

my %Phenotype_SampleSize_Lambda=();
my @Populations=("AFR","EUR");
my @TraitTypes=("Phecode","LabValue");
for(my $iPop=0;$iPop<@Populations;$iPop++)
{	
   my $LabvaluePrefix="Labvalue"."$Populations[$iPop]_"; 
   for(my $iT=0;$iT<@TraitTypes;$iT++)
   {   
     my $FilePrefix="AGD35K_$Populations[$iPop]_$TraitTypes[$iT]/AGD35K_$Populations[$iPop]_";
     my $MiddleStr="";
     if($iT ==1)
      {
      $MiddleStr= "Labvalue";
      }
      my $SummaryLambdaFile="$FilePrefix$MiddleStr"."summary_lambda.tsv";
      print "$SummaryLambdaFile\n";
       # Open the file for reading
       my $iLineNumber=0;
       open(my $fh, '<', $SummaryLambdaFile) or die "Could not open file '$SummaryLambdaFile' $!";
      # Read the file line by line
       while (my $line = <$fh>) {	      
	   if($iLineNumber ==0)
	   {
               $iLineNumber++;
	       next;
	   }    
            chomp $line;  # remove newline character
           my @cols=split(/\s++/,$line);
	   if($iT ==0)
	   {
                 $Phenotype_SampleSize_Lambda{$cols[1]}{$Populations[$iPop]}="$cols[4]\t$cols[5]";
	    }
            else
            {
		 $cols[0] =~ s/^\Q$LabvaluePrefix\E//ig;
		 $Phenotype_SampleSize_Lambda{$cols[0]}{$Populations[$iPop]}="$cols[3]\t$cols[4]";
	    }	    
 	   $iLineNumber++; 
      }
    # Close the file
      close($fh);
  }	  
}

my $OutLambdaFile="AGD35K_NoveltyGWASHits/AGD35K_Lambda_SampleSize.txt";
open(my $OUTLambda, ">", $OutLambdaFile) or die "Cannot open $OutLambdaFile: $!";
print $OUTLambda "Phenotype\tLambda_AFR\tSampleSize_AFR\tLambda_EUR\tSampleSize_EUR\n";
foreach my $pheno (keys %Phenotype_SampleSize_Lambda)
{
   print $OUTLambda "$pheno";	
   foreach my $pop (keys %{$Phenotype_SampleSize_Lambda{$pheno}})
   {
     print $OUTLambda "\t$Phenotype_SampleSize_Lambda{$pheno}{$pop}";
   }	   
   print $OUTLambda "\n";
}
close $OUTLambda;
