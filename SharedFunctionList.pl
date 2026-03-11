#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec::Functions qw(catfile);

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

#Description:
#   This subroutine reads a GWAS clumped SNV list file and
#   organizes SNVs by trait and chromosome.
#
#   For each trait, SNVs are grouped by chromosome.
#   If multiple SNVs exist on the same chromosome, their
#   positions and SNV IDs are concatenated using tab
#   separators.
#only collect the chr position and SNVID information
#my ($pos_AFRRef,$SNVIDList_AFRRef,
#         $betaFreq_AFRRef)=
#         ReadGWASclumpedSNVListFile($fullPath,$Trait,
#         \%TraitChr_PosListFull_AFR,\%TraitChr_SNVIDListFull_AFR,
#         \%TraitChr_BetaFreqFull_AFR);
 #CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
sub ReadGWASclumpedSNVListFile
{
  my ($SNVFile,
      $Trait,
      $TraitChr_AlleleListRef,
      $TraitChr_SNVIDListRef,
      $TraitChrPos_BetaFreqListRef,
      $TraitChr_pos_id_log10ListRef,
      $TraitChrPos_SNVIDListRef)=@_;
   my %TraitChr_AlleleList=%{$TraitChr_AlleleListRef};
   my %TraitChr_SNVIDList=%{$TraitChr_SNVIDListRef};
   my %TraitChrPos_BetaFreqList=%{$TraitChrPos_BetaFreqListRef};
   my %TraitChr_pos_id_log10List=%{$TraitChr_pos_id_log10ListRef};
   my %TraitChrPos_SNVIDList=%{$TraitChrPos_SNVIDListRef};
   open(my $fh, "<", $SNVFile) or die "Cannot open $SNVFile: $!";
   while (<$fh>) {
      chomp;
      my @cols=split(/\s++/,$_);
      my $ChrPos="$cols[0]\t$cols[1]";
      my $betaFreq  = join "\t", @cols[5,6,8,9,11];#A1FREQ N TEST BETA SE LOG10P
      $TraitChr_AlleleList{$Trait}{$ChrPos}="$cols[2]\t$cols[3]\t$cols[4]";
      $TraitChrPos_SNVIDList{$Trait}{$ChrPos}=$cols[2];
      $TraitChrPos_BetaFreqList{$Trait}{$ChrPos}=$betaFreq;
      if(!defined($TraitChr_SNVIDList{$Trait}{$cols[0]}))
      {
         $TraitChr_SNVIDList{$Trait}{$cols[0]}=$cols[2];
	 $TraitChr_pos_id_log10List{$Trait}{$cols[0]}="$cols[1],$cols[2],$cols[11]";
      }
      else
      {
        $TraitChr_SNVIDList{$Trait}{$cols[0]}.="\t$cols[2]";
	$TraitChr_pos_id_log10List{$Trait}{$cols[0]}.="\t$cols[1],$cols[2],$cols[11]";
      }
   }
  close $fh;
  return (\%TraitChr_AlleleList,
	  \%TraitChr_SNVIDList,
	  \%TraitChrPos_BetaFreqList,
	  \%TraitChr_pos_id_log10List,
          \%TraitChrPos_SNVIDList);
}


sub ReadGWASclumpedSNVListFileV1
{
  my ($SNVFile,
      $TraitChr_AlleleListRef,
      $TraitChr_SNVIDListRef,
      $TraitChrPos_BetaFreqListRef,
      $TraitChr_pos_id_log10ListRef,
      $TraitChrPos_SNVIDListRef)=@_;
   my %TraitChr_AlleleList=%{$TraitChr_AlleleListRef};
   my %TraitChr_SNVIDList=%{$TraitChr_SNVIDListRef};
   my %TraitChrPos_BetaFreqList=%{$TraitChrPos_BetaFreqListRef};
   my %TraitChr_pos_id_log10List=%{$TraitChr_pos_id_log10ListRef};
   my %TraitChrPos_SNVIDList=%{$TraitChrPos_SNVIDListRef};
   my $iLine=0;
   open(my $fh, "<", $SNVFile) or die "Cannot open $SNVFile: $!";
   while (<$fh>) {
      if($iLine ==0)
      {
        $iLine++;
	next;
      }      
      chomp;
      my @cols=split(/\s++/,$_);
       my $Trait=$cols[0];
      my $ChrPos="$cols[1]\t$cols[2]";
      my $betaFreq  = join "\t", @cols[6,7,9,10,12];#A1FREQ N TEST BETA SE LOG10P
      $TraitChr_AlleleList{$Trait}{$ChrPos}="$cols[3]\t$cols[4]\t$cols[5]";
      $TraitChrPos_SNVIDList{$Trait}{$ChrPos}=$cols[3];
      $TraitChrPos_BetaFreqList{$Trait}{$ChrPos}=$betaFreq;
      if(!defined($TraitChr_SNVIDList{$Trait}{$cols[1]}))
      {
         $TraitChr_SNVIDList{$Trait}{$cols[1]}=$cols[3];
         $TraitChr_pos_id_log10List{$Trait}{$cols[1]}="$cols[2],$cols[3],$cols[12]";
      }
      else
      {
        $TraitChr_SNVIDList{$Trait}{$cols[1]}.="\t$cols[3]";
        $TraitChr_pos_id_log10List{$Trait}{$cols[1]}.="\t$cols[2],$cols[3],$cols[12]";
      }
   }
  close $fh;
  return (\%TraitChr_AlleleList,
          \%TraitChr_SNVIDList,
          \%TraitChrPos_BetaFreqList,
          \%TraitChr_pos_id_log10List,
          \%TraitChrPos_SNVIDList);
}

# ------------------------------------------------------------------
# Subroutine: build_ld_clumping_hash
# Purpose: Build hash of traits and chromosomes that need LD clumping
# Input : trait_chr_count.txt file (CHR\tTRAIT\tCOUNT)
# Output:
#   $ld_traits_ref     -> hash ref: trait => 1
#   $ld_trait_chr_ref  -> hash ref: trait => { chr => 1, ... }
# ------------------------------------------------------------------
sub build_ld_clumping_TraitChr_hash {
    my ($input) = @_;

    die "Input file not defined" unless $input;
    die "Cannot open input file $input: $!" unless -e $input;

    my %ld_traits;
    my %ld_trait_chr;
    my $line_num = 0;

    open(my $fh, "<", $input) or die "Cannot open $input: $!";
    while (<$fh>) {
        chomp;
        $line_num++;

        # skip header
        next if $line_num == 1;

        my ($chr, $trait, $count) = split /\t/;
        next unless defined $chr && defined $trait && defined $count;
        $ld_traits{$trait} = 1;
        $ld_trait_chr{$trait}{$chr} = $count;
    }
    close $fh;

    return (\%ld_traits, \%ld_trait_chr);
}


#---------------------------------------------------------------
# Function: clumpsBasedSNVsCollected
#
# Description:
#   Read a clumped SNP file and collect SNV IDs from the third
#   column. The SNVs are stored in a nested hash structure
#   organized by trait, type, and chromosome.
#
# Input:
#   $ClumpedSNPFile         : Path to the clumped SNP file
#   $trait                  : Trait name (e.g., T2D, BMI)
#   $chr                    : Chromosome number
#   $type                   : Analysis type (e.g., EUR, AFR)
#   $Trait_Type_Chr_snvid   : Reference to a nested hash storing SNV IDs
#
# Data Structure:
#   $Trait_Type_Chr_snvid->{trait}->{type}->{chr} = [SNV1, SNV2, ...]
#
# Output:
#   Returns a reference to the updated nested hash containing
#   collected SNV IDs.
#
# Notes:
#   - The first line of the clumped SNP file is assumed to be a header.
#   - SNV IDs are taken from the third column of the file.
#CHROM  POS     ID      P       TOTAL   NONSIG  S0.05   S0.01   S0.001  S0.0001 SP2
#11      4932543 chr11:4932543:C:T       1.38931e-12     434     65      91      107     75      96      ch
sub clumpsBasedSNVsCollected {

    my ($ClumpedSNPFile, $trait, $Trait_snvidRef) = @_;

    my %Trait_snvid = %{$Trait_snvidRef};

    open(my $fh, "<", $ClumpedSNPFile) or die "Cannot open $ClumpedSNPFile: $!";

    my $header = <$fh>;   # Skip header line

    while (my $line = <$fh>) {
        chomp $line;
        my @cols = split(/\t/, $line);
	if(!defined($Trait_snvid{$trait}))
	{
          $Trait_snvid{$trait}=$cols[2];
	}
	else
	{
         $Trait_snvid{$trait}.="\t$cols[2]";
	}
    }
    close $fh;

    return \%Trait_snvid;
}

#The goal of the part is used to
#obtain the chromsome based independent lead significant SNVs
#The used method is LD clumping
#Reads a compressed .gz file
#Extracts the ID column and LOG10P column
#Converts LOG10P → P-value
#Outputs a simple tab-delimited file with ID and P
sub convert_log10p_to_p_sstFile {
   my ($infile, $outfile,$chr) = @_;
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
    chomp;
    @cols = split(/\s+/,$_);
    if($iLineNumber ==0)
    {
      $iLineNumber++;
      next;
    }
    next if($cols[0] != $chr);
    last if($cols[0]>$chr);
    # convert LOG10P -> P
    my $pval = 10 ** (-$cols[11]);
    $iLineNumber++;
    print $out "$cols[2]\t$pval\t$cols[0]\t$cols[1]\n";
   }
   close $fh;
   close $out;
   print "Total used lines: $iLineNumber\n";
}


#This Perl routine scans a set of GWAS summary statistic files that already contain only significant SNVs, counts the number of SNVs per chromosome for each file, and outputs a summary table.
#Directory containing the GWAS SST files (each file contains only significant SNVs).
#List of files (vector/array of file names) or all files in the directory.
#Output file path to save the summary table.
#Created it on 0306_2026
sub count_sig_snvs_by_chr {
    my ($dir, $files_ref, $outfile) = @_;
    my $iEmptyFiles=0;
    my $TotalSignficantSNVs=0;
    open(my $out, ">", $outfile) or die "Cannot write $outfile\n";
    print $out "CHR\tTRAIT\tLINE_COUNT\n";
    foreach my $file (@$files_ref) {
        my $path = "$dir/$file";
        # skip empty files
        if (-z $path) {
		#warn "Empty file skipped: $file\n";
	    $iEmptyFiles++;
            next;
        }

        # extract file prefix (remove extension)
        my $prefix = $file;
        $prefix =~ s/\.regenie\.significant\.tsv$//;#Only remove the target files

        open(my $fh, "<", $path) or die "Cannot open $path\n";

        my %chr_count;

        while (<$fh>) {
            chomp;
            my @cols = split(/\s++/);
            my $chr  = $cols[0];   # chromosome is in column 1
            $chr_count{$chr}++;
        }

        close $fh;

        foreach my $chr (sort {$a <=> $b} keys %chr_count) {
            print $out "$chr\t$prefix\t$chr_count{$chr}\n";
	    $TotalSignficantSNVs+=$chr_count{$chr};
        }
    }
    close $out;
    print "Empty Files count: $iEmptyFiles\tTotalSignficantSNVs:$TotalSignficantSNVs\n";
}


#10 104929961 chr10:104929961:C:T C T 0.00915789 24296 ADD 0.660795 0.118094 31.3094 7.65756 NA
sub merge_clumps_snvs_Trait {
    my ($dir, 
	$files_ref, 
	$outfile,
        $Trait_snvidRef,
        $Trait_chrListCount1_RawGWASRef) = @_;
    my $iEmptyFiles=0;
    my %Trait_snvid=%{$Trait_snvidRef};
    my %Trait_chrListCount1_RawGWAS=%{$Trait_chrListCount1_RawGWASRef};
    my $TotalSignficantSNVs=0;
    open(my $out, ">", $outfile) or die "Cannot write $outfile\n";
    print $out "TRAIT CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA\n";
    foreach my $file (@$files_ref) {
        my $path = "$dir/$file";
        # skip empty files
        if (-z $path) {
                #warn "Empty file skipped: $file\n";
            $iEmptyFiles++;
            next;
        }
        # extract file prefix (remove extension)
        my $trait = $file;
        $trait =~ s/\.regenie\.significant\.tsv$//;#Only remove the target files
	$trait =~ s/AGD35kAFR\_//ig;
	my $trait1=$trait;
	$trait1 =~ s/LabvalueAFR\_//ig;
	my  %SNVIDList_Selected=();
	if(defined($Trait_snvid{$trait}))
	{
           my @SNVIDList_Selected=split(/\t/,$Trait_snvid{$trait});
  	   for(my $ii=0;$ii<@SNVIDList_Selected;$ii++)
           {
            $SNVIDList_Selected{$SNVIDList_Selected[$ii]}=1;
	   }
        }
	my %ChrList_Selected=();
	if(defined($Trait_chrListCount1_RawGWAS{$trait}))
	{
	   print "Chr: $trait\t$Trait_chrListCount1_RawGWAS{$trait}\n";
           my @ChrList_Selected=split(/\t/,$Trait_chrListCount1_RawGWAS{$trait});
           for(my $ii=0;$ii<@ChrList_Selected;$ii++)
           {
            $ChrList_Selected{$ChrList_Selected[$ii]}=1;
           }    
	   # print "Chr: $trait\t$Trait_chrListCount1_RawGWAS{$trait}\n"; 
	}
	# print "$trait\n";
       	open(my $fh, "<", $path) or die "Cannot open $path\n";
        while (<$fh>) {
            chomp;
            my @cols = split(/\s++/);
            my $chr  = $cols[0];   # chromosome is in column 1
	    if(defined($ChrList_Selected{$chr}))
	    {
             print $out "$trait1 $_\n";
	     $TotalSignficantSNVs++;
	    }
	    else
	    {
	      #  print "$Trait_snvid{$trait}\n";	    
              if(defined($SNVIDList_Selected{$cols[2]}))
	      {
                print $out "$trait1 $_\n";
		$TotalSignficantSNVs++;
	       }
	    }
        }
        close $fh;
    }
    close $out;
    print "Empty Files count: $iEmptyFiles\tTotalSignficantSNVs:$TotalSignficantSNVs\n";
}	


###############################################################
# Subroutine: collectTraitChrWithSingleSNV
#
# Description:
#   This subroutine scans a trait–chromosome count hash and
#   identifies chromosomes that contain exactly one significant SNV
#   for each trait. The prefix "AGD35kAFR_" will be removed from the
#   trait name before reporting.
#
# Input:
#   $ld_trait_chr_ref  - Reference to a hash structured as:
#                       $ld_trait_chr{trait}{chr} = SNV_count
#
# Output:
#   Prints traits with chromosomes that contain exactly one SNV.
#   Returns a hash where:
#       key   = trait name (without prefix)
#       value = tab-separated chromosome list
#
# Example structure:
#   %ld_trait_chr = (
#       "AGD35kAFR_trait1" => { chr1 => 2, chr2 => 1 },
#       "AGD35kAFR_trait2" => { chr3 => 1 }
#   );
#
###############################################################

sub collectTraitChrWithSingleSNV
{
    my ($ld_trait_chr_ref) = @_;

    my %Trait_chrListCount1_RawGWAS = ();
    my $iSignificantSNVs_ChrTrait = 0;

    foreach my $trait (sort keys %{$ld_trait_chr_ref}) {
        my $trait1 = $trait;
        $trait1 =~ s/AGD35kAFR\_//ig;
        foreach my $chr (sort keys %{ $ld_trait_chr_ref->{$trait} }) {
            if ($ld_trait_chr_ref->{$trait}{$chr} == 1) {
                $iSignificantSNVs_ChrTrait++;
                if (!defined($Trait_chrListCount1_RawGWAS{$trait1})) {
                    $Trait_chrListCount1_RawGWAS{$trait1} = $chr;
                }
                else {
                    $Trait_chrListCount1_RawGWAS{$trait1} .= "\t$chr";
                }
            }
        }

        if (defined($Trait_chrListCount1_RawGWAS{$trait1})) {
            print "Trait $trait\t$Trait_chrListCount1_RawGWAS{$trait1}\n";
        }
    }
    return %Trait_chrListCount1_RawGWAS;
}

###############################################################
# Subroutine: collectEUR_AFRSameSNVEffectFromGWASsstList
#
# Description:
# This subroutine extracts SNP effects from EUR GWAS summary
# statistics files that correspond to the same trait–SNP pairs
# identified in an AFR clumped SNP list.
#
# Workflow:
# 1. Read the AFR clumped SNP file and construct a two-dimensional
#    hash structure:
#        %trait_snp{trait}{SNPID} = 1
#
# 2. Iterate through a list of EUR GWAS summary statistic files
#    (compressed .gz files). The filenames contain trait
#    information.
#
# 3. For each EUR GWAS file:
#       - Determine the trait from the filename.
#       - Scan the GWAS file.
#       - Keep rows whose SNPID appears in the AFR
#         trait–SNP hash.
#
# 4. Write the matching SNP rows into a merged output file.
#
# Input:
#   $MergedAFRchr_trait_clumpSNVsFile
#       File containing AFR clumped SNPs with at least:
#           trait   SNPID
#
#   $EURPhecode_GWASDir
#       Directory containing EUR GWAS summary statistic files.
#
#   $EURGWASFileListRef
#       Reference to an array containing GWAS filenames.
#
#   $MergedEURchr_trait_SNVsFile
#       Output file storing the matched EUR SNP records.
#
# Output:
#   A merged file containing:
#       trait   SNPID   GWAS_row
#
###############################################################

sub collectEUR_AFRSameSNVEffectFromGWASsstList
{
    my ($MergedAFRchr_trait_clumpSNVsFile,
        $EURPhecode_GWASDir,
        $EURGWASFileListRef,
        $MergedEURchr_trait_SNVsFile)=@_;
    my @EURGWASFileList=@{$EURGWASFileListRef};
    my %trait_snp;

    ###########################################################
    # Step 1: Read AFR trait-SNP list and build lookup hash
    ###########################################################
    open(my $fh,"<",$MergedAFRchr_trait_clumpSNVsFile)
        or die "Cannot open $MergedAFRchr_trait_clumpSNVsFile\n";

    my $header=<$fh>;

    while(<$fh>)
    {
        chomp;
        next if $_ eq "";
        my @arr=split(/\s++/);
        my $trait=$arr[0];
        my $snpid=$arr[3];
        $trait_snp{$trait}{$snpid}=1;
    }
    close $fh;

    ###########################################################
    # Step 2: Open output file
    ###########################################################
    open(my $out,">",$MergedEURchr_trait_SNVsFile)
        or die "Cannot write $MergedEURchr_trait_SNVsFile\n";
    print $out "TRAIT CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA\n";
#AGD35kEUR_585.3.regenie.gz
    ###########################################################
    # Step 3: Loop through EUR GWAS files
    ###########################################################
    for(my $ii=0;$ii<@EURGWASFileList;$ii++)
    #for(my $ii=0;$ii<1;$ii++)
    {
	my $file=$EURGWASFileList[$ii];    
        my $fullfile="$EURPhecode_GWASDir/$file";
        my $trait=$file;
	$trait =~ s/\.regenie\.gz//ig;
	$trait =~ s/AGD35kEUR\_//ig;
        $trait =~ s/LabvalueEUR\_//ig;
        next unless exists $trait_snp{$trait};
        open(my $gz,"gzip -dc $fullfile |")
            or die "Cannot open $fullfile\n";
        my $header=<$gz>;
	#my $iLine=0;
        chomp $header;
        while(<$gz>)
        {
            chomp;
            my @arr=split(/\s++/,$_,5);
            my $snpid=$arr[2];
	    #print "$trait\t$snpid\n";
            if(exists $trait_snp{$trait}{$snpid})
            {
                print $out "$trait $_\n";
		print "$trait $_\n";
            }
	    # $iLine++;
	     #last if($iLine>5);
        }
        close $gz;
    }
    close $out;
}


sub DetectNovelityAssoctionBasedOnKnownGWASsHits_2ndTime
{
   my ($TraitChr_Pos_KnownHits_GWASCataRef,$NovelSNVsListInfoRef)=@_;
   my %TraitChr_Pos_KnownHits_GWASCata=%{$TraitChr_Pos_KnownHits_GWASCataRef};   my @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
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
    return (\%TraitChr_Pos,\%TraitChr_SNVIDlist);
}   


#The function scans a .regenie.gz GWAS result file and retrieves only the variants whose SNV IDs match a provided list. It returns the matched GWAS records as a hash reference for downstream analysis.
sub SelectGWASSubsetbyID {
    my ($GWASFilePrefix, $SNVIDListRef,
	    $trait,$OutFile) = @_;
     my $iCollectSharedSNVs=0;
    my %SNVIDList     = %{$SNVIDListRef};
    my $GWASfile = "$GWASFilePrefix.regenie.gz";
    my $iLineNumber=0;
    open(my $Outfh, ">>", $OutFile) or die "Cannot open $OutFile: $!"; 
    open(my $fh, "-|", "zcat", $GWASfile)
        or die "Cannot open $GWASfile: $!";
    while (<$fh>) {
        chomp;
        next if $. == 1;   # skip header        
        my @cols = split(/\s+/, $_, 4);
        $iLineNumber++;	
        if (exists $SNVIDList{$cols[2]}) {
	    print $Outfh "$trait $_\n";	
	    $iCollectSharedSNVs++;
        }
    }

    close $fh;
    close $Outfh;
    print "$trait $iLineNumber\t$iCollectSharedSNVs\n";
    return ($iCollectSharedSNVs);
}

############################################################
# Function: ReadTraitBasedSignificantVariantInfo
#
# Description:
#   Read a significant SNVs result file (e.g., GWAS output)
#   and extract variant information grouped by Trait and Chr/Pos.
#
# Input:
#   $file  - Path to significant SNV result file
#
# Assumptions:
#   - First line is header and will be skipped
#   - Columns are whitespace separated
#   - At least 13 columns exist
#
#
sub ReadTraitBasedSignificantVariantInfo {
    my ($file,
	$trait_snvidRef,
        $trait_refaltRef,
        $trait_betaRef) = @_;

    die "No input file provided\n" unless defined $file;
    die "Cannot open $file: $!" unless open(my $fh, "<", $file);

    my (%trait_snvid, %trait_refalt, %trait_beta);
    my $line_number = 0;

    print "Reading file: $file\n";

    while (my $line = <$fh>) {
        chomp $line;

        # Skip header line
        if ($line_number++ == 0) {
            next;
        }

        next if $line =~ /^\s*$/;   # skip empty lines

        my @cols = split /\s+/, $line;
        # Column assignment
        my $trait  = $cols[0];
        my $chr    = $cols[1];
        my $pos    = $cols[2];
        my $snvid  = $cols[3];
        my $ref    = $cols[4];
        my $alt    = $cols[5];

        my $chrpos     = "$chr\t$pos";#5,6,8,9,11
        my $refalt     = "$ref\t$alt";
        my $beta_info  = join "\t", @cols[6,7,9,10,12];

        # Store into nested hashes
        $trait_snvid{$trait}{$chrpos} = $snvid;
        $trait_refalt{$trait}{$chrpos} = $refalt;
        $trait_beta{$trait}{$chrpos}   = $beta_info;
    }

    close $fh;
    return (\%trait_snvid, 
	    \%trait_refalt, 
	    \%trait_beta);
}

# Purpose:
#   Read known GWAS hits from a source file (tab-delimited),
#   and store them in a hash organized by:
#       Trait -> "chr:pos" -> KnownHitInfo
sub ReadKnownGWASHits_ByTraitFromSourceFile
{
   my ($OTGGWASHitsFile,$TraitChr_Pos_KnownHitsRef)=@_;
   my %TraitChr_Pos_KnownHits=%{$TraitChr_Pos_KnownHitsRef};
   my $iLineNumber=0;
   open(my $fh, "-|", "zcat", $OTGGWASHitsFile)
          or die "Cannot open $OTGGWASHitsFile: $!";
   while (<$fh>) {
      chomp;
      if($iLineNumber ==0)
      {
        $iLineNumber++;
        next;
      }
      my @cols=split(/\t/,$_);
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



#my $LabvalueEFOFile= "AGD35K_NoveltyGWASHits/AGD35K_OTG_GWASCata_Disease_Labvalue.tsv";
#LabTrait_AGD35K diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO name_AGD35K
#LYMPAB  EFO_0004587     Lymphocyte count (UKB data field 30120) OTG_Exact       lymphocyte countA quantification of lymphocytes in blood.        Lymphocyte count

sub ObtainLabValueBaseEFO_Name
{
   my ($LabvalueEFOFile)=@_;
   my %LabTrait_EFO=();
   my %LabTrait_EFOName=();
   my $iLineNumber=0;
   open(my $fh, "-|", "zcat", $LabvalueEFOFile)
            or die "Cannot open $LabvalueEFOFile: $!";
   while (<$fh>) {
      chomp;
      if($iLineNumber ==0)
      {
        $iLineNumber++;
        next;
      }
      my @cols=split(/\t/,$_);
      if(!defined($LabTrait_EFO{$cols[0]}))
      {
	   $LabTrait_EFO{$cols[0]}=$cols[1];
	   $LabTrait_EFOName{$cols[0]}=$cols[4];   
      }
      else
      {  
	    my @EFOlist=split(/:/,$LabTrait_EFO{$cols[0]});
	    my %ExistedEFO=();
	    for(my $kk=0;$kk<@EFOlist;$kk++)
	    {
              $ExistedEFO{$EFOlist[$kk]}=1;
            }		    
	   if(!defined($ExistedEFO{$cols[1]}))
	   {   
            $LabTrait_EFO{$cols[0]}.=":$cols[1]";
            $LabTrait_EFOName{$cols[0]}.=":$cols[4]";
          }
      }
    }
    close $fh;  
    return (\%LabTrait_EFO,\%LabTrait_EFOName);
}      


#my $PhecodeEFOFile= "AGD35K_NoveltyGWASHits/AGD35K_OTG_GWASCata_Disease_Phecode.tsv";
#phecode diseaseIds      trait_OTG       MatchMethod     name_EFO        description_EFO phecode1.2_simpleLabel  phecode1.2_category
sub ObtainPhecodeBaseEFO_Name
{
   my ($PhecodeEFOFile)=@_;
   my %Phecode_EFO=();
   my %Phecode_EFOName=();
   my $iLineNumber=0;
   open(my $fh, "-|", "zcat", $PhecodeEFOFile)
     or die "Cannot open $PhecodeEFOFile: $!";
   while (<$fh>) {
      chomp;
      if($iLineNumber ==0)
      {
        $iLineNumber++;
        next;
      }
      my @cols=split(/\t/,$_);
     if(!defined($Phecode_EFO{$cols[0]}))
      {
           $Phecode_EFO{$cols[0]}=$cols[1];
           $Phecode_EFOName{$cols[0]}=$cols[4];
      }
      else
      {
	   my @EFOlist=split(/:/,$Phecode_EFO{$cols[0]});
            my %ExistedEFO=();
            for(my $kk=0;$kk<@EFOlist;$kk++)
            {
              $ExistedEFO{$EFOlist[$kk]}=1;
            }
            if(!defined($ExistedEFO{$cols[1]}))
           {
             $Phecode_EFO{$cols[0]}.=":$cols[1]";
             $Phecode_EFOName{$cols[0]}.=":$cols[4]";
          }
      }
    }
    close $fh;
    return (\%Phecode_EFO,\%Phecode_EFOName);
}

# This subroutine evaluates the novelty of GWAS significant
# variants identified in the AGD35K dataset by comparing their
# genomic positions against previously reported loci (e.g., OTG).
# A distance-based criterion is applied: if a variant lies within
# ±1 Mb (1,000,000 bp) of any known locus, it is considered
# overlapped; otherwise, it is classified as a novel locus.
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
      my @Result=($trait,$Chromosome,$PosList_AGD35K[$ii],
	      $VariantIDList_AGD35K[$ii]);
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

#This function appends newly identified GWAS variants (from AGD35K) to an existing collection of novel loci for a given trait.
#It:
#Takes a trait name and chromosome.
#Parses position and variant ID lists (tab-delimited strings).
#Combines them into structured records:
#Collect novel GWAS association hits
sub CompleteNovelGWASHitsCollection
{
   my (
        $trait,
        $Chromosome,
        $PosList_AGD35KStr,
        $VariantIDList_AGD35KStr,
        $iNovelLocus,
        $NovelSNVListRef
    ) = @_;

    # Convert tab-delimited strings to arrays
    my @PosList_AGD35K       = split(/\t/, $PosList_AGD35KStr);
    my @VariantIDList_AGD35K = split(/\t/, $VariantIDList_AGD35KStr);

    # Safety check: ensure equal length
    if (@PosList_AGD35K != @VariantIDList_AGD35K)
    {
        die "Error: Position list and VariantID list have different lengths!";
    }
        # Use original array reference directly (no unnecessary copy)
    my @NovelSNVList = @{$NovelSNVListRef};

    for my $idx (0 .. $#PosList_AGD35K)
    {
        my $pos = $PosList_AGD35K[$idx];
        my $vid = $VariantIDList_AGD35K[$idx];
        # Skip undefined or empty entries
        next unless defined $pos && defined $vid;
        next if $pos eq '' || $vid eq '';

        # Build output line
        my $ResultStr = join("\t", $trait, $Chromosome, $pos, $vid);
        push @NovelSNVList, $ResultStr;
        $iNovelLocus++;
    }
    return ($iNovelLocus, \@NovelSNVList);
}

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
       foreach my $chr (sort {$a<=>$b} keys %{ $TraitChr_PosListFull{$trait} })
{
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
                      $TraitChr_PosListFull{$trait}{$chr},
		      $TraitChr_SNVIDListFull{$trait}{$chr},
                        $TraitChr_Pos_KnownHits{$trait}{$SpecificChr},
			$iOverlappedLocus,
                        $iNovelLocus,\@NovelSNVsListInfo);
                }
		else
                {
                   ($iNovelLocus,$NovelSNVsListInfoRef)=
		   CompleteNovelGWASHitsCollection($trait,$chr,
                           $TraitChr_PosListFull{$trait}{$chr},
                           $TraitChr_SNVIDListFull{$trait}{$chr},
                           $iNovelLocus,\@NovelSNVsListInfo);
                }
                @NovelSNVsListInfo=@{$NovelSNVsListInfoRef};
           }
    }
   return ($iOverlappedLocus,$iNovelLocus, 
	   $iTotalSignificantSNVs,\@NovelSNVsListInfo);
}


############################################################
#  Function: MergeIndexSNVs_500kb
#
#  Description:
#    This function merges genome-wide significant index SNVs
#    within a 500 kb window for a single chromosome and
#    selects the most significant SNV (largest -log10P)
#    as the lead SNV for each locus.
#
#    The function is designed for African ancestry (AFR)
#    GWAS analyses, where linkage disequilibrium (LD)
#    blocks are generally smaller and a 500 kb window
#    better reflects local LD structure compared to 1 Mb.

sub MergeIndexSNVs_500kb_fromString
{
    my ($chr, $chrSNPs_string,$window) = @_;

    # Array to hold parsed SNP info
    my @chrSNPs;
    my @ChrBasedSNVs=split(/\t/,$chrSNPs_string);
    my $iNumberofSNVs=scalar(@ChrBasedSNVs);
    if($iNumberofSNVs == 1 )
    {
       my ($pos, $id, $log10P) = split(/,/, $chrSNPs_string);     
       return ($chrSNPs_string,$pos);
    }
    else
    {
      #print "Raw: $chrSNPs_string\n";    
      # Parse each string line into array [pos, id, log10P]
     foreach my $line (@ChrBasedSNVs)
     {
        my ($pos, $id, $log10P) = split(/,/, $line);
        push @chrSNPs, [$pos, $id, $log10P];
     }

    # Sort by position
    @chrSNPs = sort { $a->[0] <=> $b->[0] } @chrSNPs;

    my @MergedLeadSNVs;
    my @MergedLeadPOS;
    #my $window = 250000;
    #my $window = 50000;
    my $locus_end = $chrSNPs[0]->[0] + $window;
    my $lead_pos  = $chrSNPs[0]->[0];
    my $lead_id   = $chrSNPs[0]->[1];
    my $lead_logP = $chrSNPs[0]->[2];

    for(my $i = 1; $i < @chrSNPs; $i++)
    {
        my ($pos, $id, $logP) = @{$chrSNPs[$i]};
        if($pos <= $locus_end)
        {
            $locus_end = $pos + $window if $pos + $window > $locus_end;
            # larger -log10P = more significant
            if($logP > $lead_logP)
            {
                $lead_pos  = $pos;
                $lead_id   = $id;
                $lead_logP = $logP;
            }
        }
        else
        {
            my $LeadSNVsStr="$lead_pos,$lead_id,$lead_logP";	    	    
	    push @MergedLeadSNVs,$LeadSNVsStr;
	    push @MergedLeadPOS,$lead_pos;
            $locus_end = $pos + $window;
            $lead_pos  = $pos;
            $lead_id   = $id;
            $lead_logP = $logP;
        }
    }
    # push last locus
     my $LeadSNVsStr="$lead_pos,$lead_id,$lead_logP";
     push @MergedLeadSNVs,$LeadSNVsStr;
     push @MergedLeadPOS,$lead_pos;
     $chrSNPs_string=join("\t",@MergedLeadSNVs);
     my $MergedPosStr=join("\t", @MergedLeadPOS);
     #print "$MergedPosStr\n";
     return ($chrSNPs_string,$MergedPosStr);
   } 
}

sub OverlappIndependentLocitoLargeSegments
 {
     my ($trait,
     $chr,
     $Pos_id_log10ListStr,
     $iStartPos,
     $iEndPos,
     $KeptTraitChrPosListRef,
     $RemovedTraitChrPosListRef)=@_;
     my $iRemoveSNVs=0;
     my $iWindow=$iEndPos-$iStartPos;
     my %RemovedTraitChrPosList=%{$RemovedTraitChrPosListRef};
     my %KeptTraitChrPosList=%{$KeptTraitChrPosListRef};
     my @RawChrPosList=split(/\t/,$Pos_id_log10ListStr);
      my $iNumberofRawChrPosList=scalar(@RawChrPosList);
      my $FilteredStr="";  
      my $iCount=0;
      for(my $jj=0;$jj<@RawChrPosList;$jj++)
      {
           my ($pos,$snvid,$Log10P)=split(/,/,$RawChrPosList[$jj]);
           if($pos>=$iStartPos &&  $pos<=$iEndPos)
           {
             if($iCount == 0)
             {
               $FilteredStr=$RawChrPosList[$jj];
             }
             else
             {
              $FilteredStr.="\t$RawChrPosList[$jj]";
             }
             $iCount++;
             my $Chr_PosStr="$chr\t$pos";
             $RemovedTraitChrPosList{$trait}{$Chr_PosStr}=1;
           }
        }
	 $iRemoveSNVs+=$iCount;
   return (\%KeptTraitChrPosList,
           \%RemovedTraitChrPosList,
           $iRemoveSNVs);
}

sub MergeIndependentLocitoLargeSegments
{
     my ($trait,
     $chr,
     $Pos_id_log10ListStr,
     $iStartPos,
     $iEndPos,
     $KeptTraitChrPosListRef,
     $RemovedTraitChrPosListRef)=@_;
     my $iRemoveSNVs=0;
     my $iWindow=$iEndPos-$iStartPos;
     my %RemovedTraitChrPosList=%{$RemovedTraitChrPosListRef};
     my %KeptTraitChrPosList=%{$KeptTraitChrPosListRef};

      #print "$trait\t$chr\t$TraitChr_pos_id_log10List_AFR{$trait}{$chr}\n";
      my @RawChrPosList=split(/\t/,$Pos_id_log10ListStr);
      my $iNumberofRawChrPosList=scalar(@RawChrPosList);
      my $FilteredStr="";
      if($iNumberofRawChrPosList >1)
      {
          my $iCount=0;
          for(my $jj=0;$jj<@RawChrPosList;$jj++)
          {
           my ($pos,$snvid,$Log10P)=split(/,/,$RawChrPosList[$jj]);
           if($pos>=$iStartPos &&  $pos<=$iEndPos)
           {
             if($iCount == 0)
             {
               $FilteredStr=$RawChrPosList[$jj];
             }
             else
             {
              $FilteredStr.="\t$RawChrPosList[$jj]";
             }
             $iCount++;
              my $Chr_PosStr="$chr\t$pos";
              $RemovedTraitChrPosList{$trait}{$Chr_PosStr}=1;
            }
          }
          $iRemoveSNVs+=$iCount;
	  #print "$iTrait\t$iNumberofRawChrPosList\t$iKeptLociWithinSpecificRegion\n";
	  # print "FilteredStr:$FilteredStr\n";
         my  ($chrSNPs_string,
	      $MergedPosStr)=MergeIndexSNVs_500kb_fromString(
	      $chr,
              $FilteredStr,
	      $iWindow);
         my $Chr_PosStr="$chr\t$MergedPosStr";
	   $KeptTraitChrPosList{$trait}{$Chr_PosStr}=1;
	 $iRemoveSNVs--;
      }
   return (\%KeptTraitChrPosList,
	   \%RemovedTraitChrPosList,
	   $iRemoveSNVs);
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
      chomp;      my @cols=split(/\s++/,$_,4);
      if(defined($TraitSNVID_Novelty{$Trait}{$cols[2]}))
      {
        print $OUTN "$Trait $_\n";
      }
   }
  close $fh;
  close $OUTN;
}


############################################################
# Subroutine: readPsamSampleIDList
# Description:
#   Read a PLINK .psam file and collect Sample IDs into a hash.
#   The first line (header) will be skipped.
#   Sample IDs are taken from the second column.
#
# Input:
#   $psamFile_Pop  - Path to the .psam file
#
# Output:
#   \%SampleIDList_Pop  - Hash reference containing sample IDs
#                         (SampleID => 1)
#   $TotalSampleSize    - Total number of unique samples
#
# Usage Example:
#   my ($SampleHashRef, $TotalSampleSize) =
#       readPsamSampleIDList($psamFile);
#
##FID    IID     SEX
#R223471825      R223471825      NA
############################################################
sub readPsamSampleIDList
{
    my ($psamCompressedFile_Pop) = @_;
    my %SampleIDList_Pop = ();
    my $iLineNumber = 0;
     open(my $fh, "-|", "zcat", $psamCompressedFile_Pop)
            or die "Cannot open $psamCompressedFile_Pop: $!";
    while(<$fh>)
    {
        if($iLineNumber == 0)
        {
            $iLineNumber++;
            next;
        }
        chomp;
        my @cols = split(/\s+/, $_);
        $SampleIDList_Pop{$cols[1]} = 1;
    }
    close $fh;
    my $TotalSampleSize = keys %SampleIDList_Pop;
    return (\%SampleIDList_Pop, $TotalSampleSize);
}

#-----------------------------------------------------------
# Subroutine: ObtainLambda_SampleSize_FromSummaryLambdaFile
#
# Description:
#   Read the phecode GWAS summary lambda file and extract
#   LambdaGC and SampleSize for each phecode.
#
# Input:
#   $summaryLambdaFile
#       Path to AGD35K_Phecode_summary_lambda.tsv
#
# File format:
# population  phecode  TotalVariants  SignificantVariants  LambdaGC  SampleSize
# AFR         530.1    14507284       0                    1.0461    24752
#
# Output:
#   Returns a hash reference:
#       key   = phecode
#       value = "LambdaGC\tSampleSize"
#
# Example:
#   $hash{"530.1"} = "1.0461\t24752"
#-----------------------------------------------------------
sub ObtainLambda_SampleSize_FromSummaryLambdaFile
{
    my ($summaryLambdaFile) = @_;
    my %Trait_LambdaSampleSize;
    my $iLineNumber = 0;
    open(my $fh, "<", $summaryLambdaFile)
        or die "Cannot open file $summaryLambdaFile\n";

    while(<$fh>)
    {
        chomp;

        if($iLineNumber == 0){
            $iLineNumber++;
            next;
        }
        my @cols = split(/\t/, $_);
        my $trait     = $cols[1];
        my $lambda      = $cols[4];
        my $sample_size = $cols[5];
        $Trait_LambdaSampleSize{$trait} = "$lambda\t$sample_size";
    }

    close($fh);

    return %Trait_LambdaSampleSize;
}


############################################################
# Subroutine: constructSubsetAndCountCaseControl
#
# Description:
#   Read a compressed phenotype matrix (.gz), extract a
#   subset of samples and phenotype columns, and count
#   case/control numbers for each phenotype column.
#
# Input:
#   $pheFile          : phenotype gz file
#   $sampleHashRef    : hash ref of selected sample IDs
#   $phecodeHashRef   : hash ref of selected phecode columns
#   $subsetOutFile    : output subset file
#
# Output:
#   1. Subset phenotype file
#   2. Case/control counts per phecode
#
############################################################

sub constructSubsetAndCountCaseControl
{
    my ($pheFile, 
	$sampleHashRef, 
	$phecodeHashRef,
	$subsetOutFile) = @_;

    my %colIndex=();
    my %caseCount=();
    my %controlCount=();
   my %phecodeHash=%{$phecodeHashRef};
    open(my $IN, "gunzip -c $pheFile |") or die "Cannot open $pheFile\n";
    open(my $OUT, "| gzip -c > $subsetOutFile") 
	    or die "Cannot write $subsetOutFile\n";
    my $lineNum=0;

    while(<$IN>)
    {
        chomp;
        my @cols = split(/\s+/, $_);

        # header
        if($lineNum == 0)
        {
            for(my $i=0; $i<@cols; $i++)
            {
                if(exists $phecodeHashRef->{$cols[$i]})
                {
                    $colIndex{$cols[$i]} = $i;
                }
            }

            print $OUT "FID\tIID";
            foreach my $phe (sort keys %colIndex)
            {
                print $OUT "\t$phe";
            }
            print $OUT "\n";

            $lineNum++;
            next;
        }

        my $fid=$cols[0];
        my $iid=$cols[1];

        next unless exists $sampleHashRef->{$iid};

        print $OUT "$fid\t$iid";

        foreach my $phe (sort keys %colIndex)
        {
            my $idx = $colIndex{$phe};
            my $val = $cols[$idx];

            print $OUT "\t$val";

            next if(!defined $val || $val eq "NA");

            if($val == 1 || $val == 1.0)
            {
                $caseCount{$phe}++;
            }
            elsif($val == 0 || $val == 0.0)
            {
                $controlCount{$phe}++;
            }
        }
        print $OUT "\n";
    }

    close($IN);
    close($OUT);
    return (\%caseCount,\%controlCount);
}

###############################################################
# Subroutine: LoadLabTraitGWAStraitHash
#
# Description:
#   Read a gzipped trait mapping file and construct a hash that
#   maps LabTrait IDs to GWAS trait names.
#
# Input File Format (tab-delimited):
#   LabTrait        GWAStrait
#   BASOAB          Basophil count
#
# Input:
#   $LabTraitFile   - gzipped file (e.g. Johns_AGD35K_LabTraits.txt.gz)
#
# Output:
#   Reference to hash:
#       key   = LabTrait
#       value = GWAStrait
#
# Usage Example:
#   my $TraitHashRef = LoadLabTraitGWAStraitHash($file);
#   my %TraitHash = %{$TraitHashRef};
#
###############################################################

sub LoadTrait_AnnotationHash
{
    my ($TraitFile) = @_;

    my %Trait_Annotation = ();
    my $iLine = 0;
   open(my $fh, "gunzip -c $TraitFile |") or die "Cannot open $TraitFile\n";
    while(<$fh>)
    {
        chomp;

        # skip header
        if($iLine == 0)
        {
            $iLine++;
            next;
        }

        my @cols = split(/\t/,$_);

        $Trait_Annotation{$cols[0]} = $cols[1];
    }

    close($fh);

    return \%Trait_Annotation;
}


###############################################################
# Subroutine: OutSampleSize_Lambda_AGD35K
#
# Description:
#   Output case/control counts, effective sample size,
#   lambda GC, and sample size for AGD35K phecode traits
#   and lab traits in a single table.
#
# Input:
#   $caseCount_AGD35K_122Ref
#   $controlCount_AGD35K_122Ref
#   $phecodeHashRef
#   $labvalue_LambdaSampleSizeRef
#   $Phecode1_2AnnotationRef
#   $LabTrait_AnnotationRef
#   $OutLambda_SampleSizeFile
#my $Phecode1_2AnnotationFile="AGD35K_Phecode122_0308_2026/phecode1.2.tsv.gz";
#AGD35K_Phecode122_0308_2026$ zcat phecode1.2.tsv.gz |head
#phecode1.2_code phecode1.2_label        phecode1.2_simpleLabel  phecode1.2_category
sub OutSampleSize_Lambda_AGD35K
{
      my ($caseCount_AGD35K_122Ref,
          $controlCount_AGD35K_122Ref,	      
	  $phecodeHashRef,
          $labvalue_LambdaSampleSizeRef,
	  $Phecode1_2AnnotationRef,
	  $LabTrait_AnnotationRef,
          $OutLambda_SampleSizeFile) = @_;
    my %caseCount=%{$caseCount_AGD35K_122Ref};
    my %controlCount=%{$controlCount_AGD35K_122Ref};
    my %phecodeHash=%{$phecodeHashRef};
    my %labvalue_LambdaSampleSize=%{$labvalue_LambdaSampleSizeRef};
    my %Phecode1_2Annotation=%{$Phecode1_2AnnotationRef};
    my %LabTrait_Annotation=%{$LabTrait_AnnotationRef};
    open(my $OUTLAB, "| gzip -c > $OutLambda_SampleSizeFile")
            or die "Cannot write $OutLambda_SampleSizeFile\n";

    print $OUTLAB "Trait\tPhecode1.2\tLambda_GC\tTotal_Sample_Size\tNumber_of_Cases\t";
    print $OUTLAB "Number_of_Controls\tCase:Control_Ratio\tEffective_Sample_Size\n";
    foreach my $phe (sort keys %caseCount)
    {
        my $case = $caseCount{$phe} // 0;
        my $ctrl = $controlCount{$phe} // 0;
        my $fRatio=sprintf("%.4f",$case/$ctrl);
        my $iEffective_Sample_Size=(4*$case*$ctrl)/($case+$ctrl);
        $iEffective_Sample_Size=sprintf("%.0f",$iEffective_Sample_Size);
        print $OUTLAB "$Phecode1_2Annotation{$phe}\t$phe\t$phecodeHash{$phe}\t$case\t$ctrl\t$fRatio\t$iEffective_Sample_Size\n";
    }
    foreach my $lab (keys %labvalue_LambdaSampleSize) {
     my ($lambda,$SampleSize)=split(/\s++/,$labvalue_LambdaSampleSize{$lab});
      print $OUTLAB "$LabTrait_Annotation{$lab}\t$lab\t$labvalue_LambdaSampleSize{$lab}\tNA\tNA\tNA\t$SampleSize\n";
    }
   close ($OUTLAB);
 }  

1;   # IMPORTANT: module must return true
