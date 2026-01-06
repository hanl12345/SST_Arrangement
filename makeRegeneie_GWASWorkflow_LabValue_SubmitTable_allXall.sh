out_tsv="Regenie_GWAS35K_LabValue_Table.tsv"

# header
echo -e "idx\tcovarColList\tCovFile\tSuperPop\tstep1_Options\tstep2_Options\toutBucket_prefix\tkeepFile\tphenoFile\tphenoColList_Array" > "${out_tsv}"

idx=1

# build PC string once
PC_STR=$(printf "PC%d," {1..20} | sed 's/,$//')

for pop in AFR EUR; do
  for t in BASOAB EOSIA Hgb LYMPAB MCH MCHC MCV MONOAB MPV NTAuto PCV Plt-Ct RBC RDWSD WBC; do

    out_bucket="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/LabValue_${pop}Prob05/Phecode"

    keep_file="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_AGD35KSampleClass/AGD35KQV2_SupervisedProb0.5_${pop}_Subset.txt"
   phenotype_File="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_PCsTest/LabValue_agd35k_outpatient_only_exclusions_pervisit.txt"
    covar_str="Sex,Age_${t},Age2_${t},${PC_STR}"
    covFile_str="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_PCsTest/Covariate_ProjectedPCs_agd35k_outpatient_only_exclusions_pervisit.txt"
    step1_OptionsStr="--bsize 1000 --cv 10 --qt"
    step2_OptionsStr="--bsize 1000 --qt"
   # echo -e "${idx}\t${covar_str}\t${covFile_str}\t${pop}\t${step1_OptionsStr}\t${step2_OptionsStr}\t${out_bucket}\t${keep_file}" >> "${out_tsv}"
    echo -e "${idx}\t\"${covar_str}\"\t\"${covFile_str}\"\t\"${pop}\"\t\"${step1_OptionsStr}\"\t\"${step2_OptionsStr}\"\t\"${out_bucket}\"\t\"${keep_file}\"\t\"${phenotype_File}\"\t[\"${t}\"]" >> "${out_tsv}"

    ((idx++))
  done
done

