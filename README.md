This script ObtainSignificantSNV_Inflator.pl is used to detect number of significant SNVs,Sample size, and lambda genome inflator, draw out number of significant SNVs
These two codes makeRegeneie_GWASWorkflow_LabValue_SubmitTable_allXall.sh makeRegeneie_GWASWorkflow_Phecode_SubmitTable_allXall.sh are used to create submitted table files for the workflow analysis 
01_20_2026
How can make use Parallel::ForkManager; be used in Terra parallel compuation
This bypasses CPAN config entirely.
Step 1: Install cpanminus locally
curl -L https://cpanmin.us | perl - App::cpanminus
Step 2: Install Parallel::ForkManager into your home directory
~/perl5/bin/cpanm --local-lib=~/perl5 Parallel::ForkManager
Step 3: Set environment variables (VERY important)
mkdir -p ~/perl5
echo 'export PERL5LIB=$HOME/perl5/lib/perl5:$PERL5LIB' >> ~/.bashrc
echo 'export PATH=$HOME/perl5/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
Step 4: Verify
perl -MParallel::ForkManager -e 'print "OK\n"'
<img width="468" height="228" alt="image" src="https://github.com/user-attachments/assets/e9b2cd19-d809-425b-9425-79b5ee023f6f" />
