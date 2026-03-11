Prepared Data: Covariate list including age, gender and PC1-PC20
gs://working-set-redeposit/ica-agd/agd35k/AGD35KcovariateList.txt
Phenotype data for GWAS analysis
gs://working-set-redeposit/ica-agd/agd35k/AGD_35k_code_data_floor_binarized_wide_format_phecode_12_gwas.txt
#AGD35K EUR per-chromosome PGEN files 
#gs://working-set-redeposit/ica-agd/agd35k/eligible_v5_202410/AGD35K_RegenieSubset/EUR/ChrBased/AGD35K_MAC100MR002QCed_EUR_chr*.pgen
#AGD35K AFR per-chromosome PGEN files 
#gs://working-set-redeposit/ica-agd/agd35k/eligible_v5_202410/AGD35K_RegenieSubset/AFR/ChrBased/AGD35K_MAC100MR002QCed_AFR_chr*.pgen
#AGD35K 122 phecode based gwas SST files
gs://fc-efa2b59f-4130-48b7-b894-2d859978701a/Phecode-based_GWAS_for_African_and_European_ancestry_cohorts_regenie_gz_files.tar.gz
Total process
1.Estimate the genomic inflation factor (λGC) and sample size from the Regenie GWAS SST files
#Require: gwas SST files as input
ObtainSignificantSNV_Inflator_122_Step1.pl

2.Collect the significant SNVs with p-value < 5e-8 from the SST files.
Perform LD clumping on the significant SNVs using the study samples as the reference panel.
Extract SNVs from the EUR-based SST files that correspond to the leading independent SNVs identified in the AFR SST files.
PerformAGD35K_122GWAS_Summary_step2.pl

3. Construct a reference set of known GWAS loci 
based on the GWAS Catalog and OTG by phecode and labvalue independently
PhecodeMatch_0308_2026_Step3_1.r
LabValueMatchV2_0308_2026_Step3_2.r 
4.Obtain lead independent loci by removing known regions and 
identify novel loci based on a reference set of known GWAS loci using distance based method. 
ReOrganizedIndependentLoci_0308_2026_step4.pl
5.Obtain the case and control counts, 
construct the AFR GWAS phecode subset, 
and  so that we can estimate the correlations among these phenotypes.
 Estimate_agd35K_SampleSize_Step5.pl

