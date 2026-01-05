#!/usr/bin/env perl
use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Basename;

# ==========================
# User settings
# ==========================
#my $indir  = "PhecodeAFR_GWASRegenie";
#my $outdir = "PhecodeAFR_GWASRegenie_SignificantResult";

#my $indir  = "PhecodeEUR_GWASRegenie";
#my $outdir = "PhecodeEUR_GWASRegenie_SignificantResult";

#my $indir  = "LabvalueAFR_GWASRegenie";
#my $outdir = "LabvalueAFR_GWASRegenie_SignificantResult";  

my $indir  = "LabvalueEUR_GWASRegenie";
my $outdir = "LabvalueEUR_GWASRegenie_SignificantResult";

my $log10p_thresh = 7.30103;   # -log10(5e-8)
my $chisq_median_expected = 0.4549364;  # median of chi-square(1 df)

mkdir $outdir unless -d $outdir;

# Summary output
open my $SUM, ">", "$outdir/summary_lambda.tsv"
  or die "Cannot write summary file\n";
print $SUM join("\t",
    qw(File TotalVariants SignificantVariants LambdaGC)
), "\n";

# ==========================
# Process each .gz file
# ==========================
opendir my $DH, $indir or die "Cannot open directory $indir\n";
my $iFiles=0;
while (my $file = readdir($DH)) {

    next unless $file =~ /\.gz$/;

    my $path = "$indir/$file";
   print "$path\n";
    my $z = IO::Uncompress::Gunzip->new($path)
        or die "Cannot open $path: $GunzipError\n";

    my $outfile = "$outdir/" . basename($file, ".gz") . ".significant.tsv";
    open my $OUT, ">", $outfile or die "Cannot write $outfile\n";

    my (%idx, @chisq);
    my ($total, $sig,$SampleSize) = (0, 0, 0);
    # ==========================
# Read header (whitespace-delimited)
# ==========================
  my $header;
  while (my $line = <$z>) {
    next if $line =~ /^\s*$/;   # skip empty lines
    $header = $line;
    last;
  }

   die "Empty file or missing header in $file\n" unless defined $header;
   chomp $header;
    
   my @cols = split /\s++/, $header;

    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }

    die "LOG10P column not found in $file\n" unless exists $idx{LOG10P};
    die "CHISQ column not found in $file\n"  unless exists $idx{CHISQ};

    print $OUT $header, "\n";

    # Process rows
    while (my $line = <$z>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split /\s++/, $line;
        $total++;
        my $variantSize=$f[ $idx{N} ];
        my $log10p = $f[ $idx{LOG10P} ];
        my $chisq  = $f[ $idx{CHISQ} ];

	 if($variantSize >= $SampleSize)
	 {
          $SampleSize=$variantSize;
	 }

        push @chisq, $chisq if defined $chisq && $chisq ne "NA";

        if (defined $log10p && $log10p ne "NA" && $log10p > $log10p_thresh) {
            print $OUT $line, "\n";
            $sig++;
        }
    }

    close $OUT;
    $z->close();

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
    print $SUM join("\t", $Filestr, $total, $sig, $lambda,$SampleSize), "\n";
     $iFiles++;
     #  last if($iFiles>=5);
}

closedir $DH;
close $SUM;

print "Done. Results in: $outdir\n";
