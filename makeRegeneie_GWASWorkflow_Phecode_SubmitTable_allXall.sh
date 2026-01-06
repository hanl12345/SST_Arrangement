#keep files:
#AFR:N=27807,EUR: 5143
#gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_AGD35KSampleClass/AGD35KQV2_SupervisedProb0.5_AFR_Subset.txt
#gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_AGD35KSampleClass/AGD35KQV2_SupervisedProb0.5_EUR_Subset.txt
#Split total file into 219 files
#split -l 20 -d --additional-suffix=.txt IndependentPhecodeList_AGD35K.txt IndPhecodeList_AGD35K_Part_
#phenoColList_File
#gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/PhecodeList_AGD35K_Sub/IndPhecodeList_AGD35K_Part_18.txt
#Outfile _Prefix
#gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/Phecode_AFRProb05/Phecode
#gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/Phecode_EURProb05/Phecode

out_tsv="Regenie_GWAS35K_Phecode_Table.tsv"

echo -e "idx\tSuperPop\tphenoColList_File\toutBucket_prefix\tkeepFile" > "${out_tsv}"

idx=1

for pop in AFR EUR; do
  for i in $(seq -w 0 17); do
    pheno_file="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/PhecodeList_AGD35K_Sub/IndPhecodeList_AGD35K_Part_${i}.txt"

    out_bucket="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_GWAS35KResultV2/Phecode_${pop}Prob05/Phecode"

    keep_file="gs://fc-secure-540f27be-97ea-4ffd-adb7-c195458eb278/Regenie_allxall/AllxAll_AGD35KSampleClass/AGD35KQV2_SupervisedProb0.5_${pop}_Subset.txt"

    echo -e "${idx}\t${pop}\t${pheno_file}\t${out_bucket}\t${keep_file}" >> "${out_tsv}"

    idx=$((idx + 1))
  done
done


