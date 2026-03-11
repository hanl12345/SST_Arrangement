#!/usr/bin/env perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Basename;

# ==========================
# User settings
my $Population="AFR";
#my $outdir = "/home/jupyter/AGD35K_GWAS122_2026/GWAS_AnalysisV0305_2026/AGD35K_$Population"."_Phecode";
my $outdir  = "/home/jupyter/AGD35K_GWAS122_2026/GWAS_AnalysisV0305_2026/AGD35K_$Population"."_LabValue";
my $indir  =$outdir."_GWAS";

my $log10p_thresh = 7.30103;   # -log10(5e-8)
my $chisq_median_expected = 0.4549364;  # median of chi-square(1 df)

mkdir $outdir unless -d $outdir;

# Summary output
open my $SUM, ">", "$outdir/summary_lambda.tsv"
  or die "Cannot write summary file\n";
print $SUM join("\t",
    qw(popultaion phecode TotalVariants SignificantVariants LambdaGC SampleSize)
), "\n";

# ==========================
# Process each .gz file
# ==========================
opendir my $DH, $indir or die "Cannot open directory $indir\n";
my $iFiles=0;
while (my $file = readdir($DH)) {
    next unless $file =~ /\.gz$/;
    $iFiles++;
    #last if($iFiles>1);
    #next unless $file =~ /840/;      # only filenames containing 840
    my $path = "$indir/$file";
   print "$path\n";

    open my $fh, "-|", "zcat $path" or die $!;
    my $outfile = "$outdir/" . basename($file, ".gz") . ".significant.tsv";
    open my $OUT, ">", $outfile or die "Cannot write $outfile\n";

    my (%idx, @chisq);
    my ($total, $sig,$SampleSize) = (0, 0, 0);
    # ==========================
# Read header (whitespace-delimited)
# ==========================
    my $iLineNum=0;
    #CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
    #1 176482 chr1:176482:T:A T A 0.0143804 25034 ADD 0.043795 0.12248 0.127856 0.142266 NA
    # Process rows
    while (my $line = <$fh>) {
       if($iLineNum ==0)
       {
	     $iLineNum++;
	       next;
       }
    	chomp $line;
        my @f = split /\s++/, $line;
	next if ($f[10] eq "NA" || $f[11] eq "NA");
        $total++;
	 if($f[6] >= $SampleSize)
	 {
          $SampleSize=$f[6];
	 }
	 push @chisq, $f[10];
	 if ($f[11] >= $log10p_thresh) 
	 {
	    print $OUT $line, "\n";
	     $sig++;
	  }
	  #$iLineNum++;
    }

    close $OUT;
    close $fh;
    #print "Total read lines: $iLineNum\n";
    # ==========================
    # Calculate lambda GC
    # ==========================
    my $lambda = "NA";
    if (@chisq) {
        @chisq = sort { $a <=> $b } @chisq;
        my $n = scalar @chisq;
        my $median = ($n % 2)
            ? $chisq[int($n/2)]
            : ($chisq[$n/2 - 1] + $chisq[$n/2]) / 2;

        $lambda = sprintf("%.4f", $median / $chisq_median_expected);
    }
   my $Filestr=$file;
   $Filestr =~ s/\.regenie\.gz//ig;
   $Filestr =~ s/AGD35kAFR\_//ig;
   $Filestr =~ s/AGD35kEUR\_//ig;
   $Filestr =~ s/LabvalueAFR\_//ig;
   $Filestr =~ s/LabvalueEUR\_//ig;
   print $SUM join("\t","$Population", $Filestr, $total, $sig, $lambda,$SampleSize), "\n";
   #print $SUM join("\t","AFR", $Filestr, $total, $sig, $lambda,$SampleSize), "\n";
     $iFiles++;
     #  last if($iFiles>=5);
}

closedir $DH;
close $SUM;

print "Done. Results in: $outdir\n";
