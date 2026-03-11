#install.packages("arrow")
#install.packages("R.utils")
library(arrow)
library(dplyr)
library(data.table)
library(tidyr)
library(R.utils);
library(readr);
library(stringr)
# -----------------------------
# 1. Check input path
# -----------------------------

GWASCatalogueDF <- fread("AGD35K_Phecode122_0308_2026/gwas-catalog-download-associations-alt-full.tsv.gz", 
			 sep = "\t",quote = "");
GWASCatalogueDF  =as.data.frame(GWASCatalogueDF);
GWASCatalogueDF=GWASCatalogueDF[,c(8,36)];
colnames(GWASCatalogueDF)=c("traitFromSource","EFO_URL");
print(str(GWASCatalogueDF));
print(nrow(GWASCatalogueDF));
# filter rows containing "Phecode"
GWASCatalogueDF <- GWASCatalogueDF[
  grepl("Phecode", GWASCatalogueDF$traitFromSource, ignore.case = TRUE),
]
GWASCatalogueDF=unique(GWASCatalogueDF);

GWASCatalogueDF$trait <- sub("\\s*\\(.*\\)", "", GWASCatalogueDF$traitFromSource)
GWASCatalogueDF$phecode <- sub(".*\\(|\\)", "", GWASCatalogueDF$traitFromSource)
GWASCatalogueDF$phecode <- gsub("\\)", "", GWASCatalogueDF$phecode)
GWASCatalogueDF$phecode <- gsub("^PheCode\\s*", "", GWASCatalogueDF$phecode);
GWASCatalogueDF$phecode  =as.numeric(GWASCatalogueDF$phecode);
GWASCatalogueDF$diseaseIds <- sub(".*/", "", GWASCatalogueDF$EFO_URL);
GWASCatalogueDF=GWASCatalogueDF[,c("diseaseIds","trait","phecode")];
GWASCatalogueDF=unique(GWASCatalogueDF);
print(str(GWASCatalogueDF));
#diseaseIds","traitFromSource","phecode"

AGD35K_AFR_summaryDF <- fread("/home/jupyter/AGD35K_GWAS122_2026/GWAS_AnalysisV0305_2026/AGD35K_AFR_Phecode/summary_lambda.tsv");
AGD35K_AFR_summaryDF =as.data.frame(AGD35K_AFR_summaryDF);
#AGD35K_AFR_summaryDF=AGD35K_AFR_summaryDF[which(AGD35K_AFR_summaryDF$SignificantVariants != 0),];
#print("AGD35K_AFR_summaryDF");
#print(str(AGD35K_AFR_summaryDF));
#AGD35K_EUR_summaryDF <- fread("AGD35K_EUR_summary_lambda.tsv");
#AGD35K_EUR_summaryDF =as.data.frame(AGD35K_EUR_summaryDF);
#AGD35K_EUR_summaryDF=AGD35K_EUR_summaryDF[which(AGD35K_EUR_summaryDF$SignificantVariants != 0),]
#print("AGD35K_EUR_summaryDF");
#print(str(AGD35K_EUR_summaryDF));
#AGD35K_AFR_EUR_summaryDF=rbind(AGD35K_AFR_summaryDF,AGD35K_EUR_summaryDF);
AGD35K_AFR_EUR_summaryDF=AGD35K_AFR_summaryDF;
rownames(AGD35K_AFR_EUR_summaryDF)=NULL;
AGD35K_AFR_EUR_summaryDF=as.data.frame(AGD35K_AFR_EUR_summaryDF);
AGD35K_AFR_EUR_summaryDF$phecode=as.numeric(AGD35K_AFR_EUR_summaryDF$phecode);
UniquePhecode_AGD35K=unique(AGD35K_AFR_EUR_summaryDF$phecode);
print(c("number of UniquePhecode_AGD35K",length(UniquePhecode_AGD35K)));

OTGDVIDDF <- fread("AGD35K_Phecode122_0308_2026/disease_variant_study_Merged_OTGV25.12.tsv.gz",
				       	sep = "\t",quote = "");
#print(str(OTGDVIDDF));
OTGDVIDDF=OTGDVIDDF[,c("traitFromSource","diseaseIds")];
print(c("Number of OTGDVIDDF",nrow(OTGDVIDDF)));
OTGDVIDDF <- OTGDVIDDF[
  grepl("Phecode", OTGDVIDDF$traitFromSource, ignore.case = TRUE),
];
OTGDVIDDF=unique(OTGDVIDDF);
OTGDVIDDF$trait <- sub("\\s*\\(.*\\)", "", OTGDVIDDF$traitFromSource)
OTGDVIDDF$phecode <- sub(".*\\(|\\)", "", OTGDVIDDF$traitFromSource)
OTGDVIDDF$phecode <- gsub("\\)", "", OTGDVIDDF$phecode)
OTGDVIDDF$phecode <- gsub("^PheCode\\s*", "", OTGDVIDDF$phecode);
OTGDVIDDF$phecode <- gsub("^Phecode:\\s*", "", OTGDVIDDF$phecode); 
#print(unique(OTGDVIDDF$phecode));
OTGDVIDDF$phecode  =as.numeric(OTGDVIDDF$phecode);
OTGDVIDDF=OTGDVIDDF[,c("diseaseIds","trait","phecode")];
OTGDVIDDF=OTGDVIDDF[OTGDVIDDF$phecode %in% UniquePhecode_AGD35K,];
OTGDVIDDF=unique(OTGDVIDDF);
OTGDVIDDF$MatchMethod="OTG_Exact";
print("OTGDVIDDF_Part");
print(str(OTGDVIDDF));
#OTGDVIDDF1=OTGDVIDDF[,c("diseaseIds","phecode")];
#Remove the duplicate line with same diseaseIds and phecode  that have existed in OTGDVIDDF 
#At first keep the lines with its phecode from agd369 phecode

GWASCatalogueDF=GWASCatalogueDF[which(GWASCatalogueDF$phecode %in% UniquePhecode_AGD35K),];
print("GWASCataDF 122phecodes  subset");
print(str(GWASCatalogueDF));
GWASCatalogueDF <- GWASCatalogueDF %>%
  anti_join(OTGDVIDDF, by = c("diseaseIds", "phecode"));
print("GWASCataDF_after removing OTGset");
print(str(GWASCatalogueDF));
GWASCatalogueDF$MatchMethod="GWASCata_Exact";
GWASCatalogue_OTGDVIDDF= rbind(GWASCatalogueDF,OTGDVIDDF);
rownames(GWASCatalogue_OTGDVIDDF)=NULL;
GWASCatalogue_OTGDVIDDF =as.data.frame(GWASCatalogue_OTGDVIDDF);
print(c("MergedData", nrow(GWASCatalogue_OTGDVIDDF)));
GWASCatalogue_OTGDVIDDF=unique(GWASCatalogue_OTGDVIDDF);
print(c(" unique MergedData", nrow(GWASCatalogue_OTGDVIDDF)));
GWASCatalogue_OTGDVIDDF=GWASCatalogue_OTGDVIDDF[GWASCatalogue_OTGDVIDDF$phecode %in% UniquePhecode_AGD35K,];
GWASCatalogue_OTGDVIDDF=unique(GWASCatalogue_OTGDVIDDF[,
			       c("diseaseIds",
				 "trait",
				 "phecode",
				 "MatchMethod")]);
print(c("AGD35K MergedData", nrow(GWASCatalogue_OTGDVIDDF)));
print("Phecode Table");
tab <- table(GWASCatalogue_OTGDVIDDF$phecode);
tab[tab > 1];

print("GWAS Catalogue Overlaps");
overlap <- intersect(UniquePhecode_AGD35K, GWASCatalogueDF$phecode);
print(length(unique(overlap)));
# Elements of v1 NOT in v2
print(sort(UniquePhecode_AGD35K[!UniquePhecode_AGD35K %in% overlap]));

print("OTG data overlaps");
overlap <- intersect(UniquePhecode_AGD35K, OTGDVIDDF$phecode);
print(length(unique(overlap)));

# Elements of v1 NOT in v2
print(sort(UniquePhecode_AGD35K[!UniquePhecode_AGD35K %in% overlap]));
# Combine phecodes from GWAS catalogue and OTG list
print("OTG _JOIN GWAS catalogue");
UnitedPhecodes <- c(GWASCatalogueDF$phecode, OTGDVIDDF$phecode);

# Find overlap with AGD35K phecodes
overlap <- intersect(UniquePhecode_AGD35K, UnitedPhecodes)

# Number of overlapping phecodes
print(length(overlap))

# Elements in UniquePhecode_AGD35K NOT in overlap (i.e., not in GWAS/OTG)
missing_phecodes <- sort(UniquePhecode_AGD35K[!UniquePhecode_AGD35K %in% overlap])
print(missing_phecodes)

OTG_Phecode1.2_ManuallDF <- read.table("AGD35K_Phecode122_0308_2026/PairofEFO_Phecode1.2_ManuallyOTG.txt.gz",
                 sep = "\t",
		 header = TRUE,   # use TRUE if the first line has column names
                 fill = TRUE,     # fill missing values
                 strip.white = TRUE)  # remove leading/trailing spaces);
colnames(OTG_Phecode1.2_ManuallDF)[colnames(OTG_Phecode1.2_ManuallDF) == "UnfoundPhecode"] <- "phecode";
colnames(OTG_Phecode1.2_ManuallDF)[colnames(OTG_Phecode1.2_ManuallDF) == "EFO"] <- "diseaseIds";
colnames(OTG_Phecode1.2_ManuallDF)[colnames(OTG_Phecode1.2_ManuallDF) == "Trait"] <- "phenotype";
colnames(OTG_Phecode1.2_ManuallDF)[colnames(OTG_Phecode1.2_ManuallDF) == "ReportedTraits"] <- "trait";
OTG_Phecode1.2_ManuallDF$MatchMethod="OTG_Manually";
OTG_Phecode1.2_ManuallDF$MatchMethod[which(is.na(OTG_Phecode1.2_ManuallDF$phenotype))]="ChatGPT";
OTG_Phecode1.2_ManuallDF=OTG_Phecode1.2_ManuallDF[,c("diseaseIds","trait","phecode","MatchMethod")];

OTG_Phecode1.2_ManuallDF <-
OTG_Phecode1.2_ManuallDF[rowSums(is.na(OTG_Phecode1.2_ManuallDF)) <
                            ncol(OTG_Phecode1.2_ManuallDF), ];
OTG_GWASCataDF=rbind(GWASCatalogue_OTGDVIDDF,OTG_Phecode1.2_ManuallDF);
rownames(OTG_GWASCataDF)=NULL;
OTG_GWASCataDF=as.data.frame(OTG_GWASCataDF);
print("OTG_GWAS_set Total Lines");
print(nrow(OTG_GWASCataDF));
print("Number of phecodes");
print(length(unique(OTG_GWASCataDF$phecode)));
ObtainedPhedcode=unique(OTG_GWASCataDF$phecode);
print("set difference");
print(setdiff(ObtainedPhedcode,UniquePhecode_AGD35K));

OTG_DiseaseDF <- read.table("AGD35K_Phecode122_0308_2026/disease_OTGV25.12.tsv.gz",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
  );

OTG_DiseaseDF=OTG_DiseaseDF[,c("id","name","description")];
colnames(OTG_DiseaseDF)[colnames(OTG_DiseaseDF) == "id"] <- "diseaseIds";
OTG_GWASCata_DiseaseDF=merge(OTG_GWASCataDF,OTG_DiseaseDF,
			     by="diseaseIds", all.x = TRUE);
print("OTG_GWASCata_DiseaseDF");
print(str(OTG_GWASCata_DiseaseDF));
print("OTG_GWASCata_DiseaseDF Phecode");
OTG_GWASCata_DiseaseDF_Phecode=unique(OTG_GWASCata_DiseaseDF$phecode);

print("final lost phecode:")
missing_phecodes <- setdiff(
  UniquePhecode_AGD35K,
  OTG_GWASCata_DiseaseDF_Phecode
)
print(sort(missing_phecodes));
#print(length(unique(OTG_GWASCata_DiseaseDF$phecode)));
phecode1.2DF <- read.table("AGD35K_Phecode122_0308_2026/phecode1.2.tsv.gz",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
  );

phecode1.2DF$phecode1.2_label=NULL;

#print(str(phecode1.2DF));
OTG_GWASCata_Disease_PhecodeDF=merge(OTG_GWASCata_DiseaseDF, phecode1.2DF,
				     by.x="phecode",by.y="phecode1.2_code");
OTG_GWASCata_Disease_PhecodeDF=as.data.frame(OTG_GWASCata_Disease_PhecodeDF);
#OTG_GWASCata_Disease_PhecodeDF$phecode=NULL;
print(c("OTGDVIDDF",str(OTGDVIDDF)));
print(c("OTG_GWASCata_Disease_PhecodeDF", 
	str(OTG_GWASCata_Disease_PhecodeDF)));
MatchedDF_byEFO=merge(OTGDVIDDF,
		      OTG_GWASCata_Disease_PhecodeDF,
		      by="diseaseIds");
print("Final disease_variant pair size by EFO:");
print(nrow(MatchedDF_byEFO));
print("Final phecode size by EFO:");
print(length(unique(MatchedDF_byEFO$phecode.y)));

MatchedDF_byPhecode=merge(OTGDVIDDF,
			  OTG_GWASCata_Disease_PhecodeDF,
			  by="phecode");
print("Final disease_variant pair size by Phecode:");
print(nrow(MatchedDF_byPhecode));
print("Final phecode size by Phecode:");
print(length(unique(MatchedDF_byPhecode$phecode)));

colnames(OTG_GWASCata_Disease_PhecodeDF)[colnames(OTG_GWASCata_Disease_PhecodeDF) == "trait"] <- "trait_OTG";
colnames(OTG_GWASCata_Disease_PhecodeDF)[colnames(OTG_GWASCata_Disease_PhecodeDF) == "name"] <- "name_EFO";
colnames(OTG_GWASCata_Disease_PhecodeDF)[colnames(OTG_GWASCata_Disease_PhecodeDF) == "description"] <- "description_EFO";
#OTG_GWASCata_Disease_PhecodeDF$trait_OTG=NULL;

#print(str(OTG_GWASCata_Disease_PhecodeDF));

fwrite(
   OTG_GWASCata_Disease_PhecodeDF,
  "AGD35K_Phecode122_0308_2026/AGD35K_OTG_GWASCata_Disease_Phecode_0308_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );
print(str(OTG_GWASCata_Disease_PhecodeDF));
colnames(OTG_GWASCata_Disease_PhecodeDF)[colnames(OTG_GWASCata_Disease_PhecodeDF) == "phecode"] <- "Phecode1.2";

OTG_StudysetDF <- read.table("AGD35K_Phecode122_0308_2026/disease_variant_study_Merged_OTGV25.12.tsv.gz",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
  );
print("OTG_GWASCata_Disease_PhecodeDF");
print(colnames(OTG_GWASCata_Disease_PhecodeDF));
OTG_StudysetDF=OTG_StudysetDF[,c("diseaseIds","variantId",
                "chromosome","position","traitFromSource")];
print(colnames(OTG_StudysetDF));
FinalMergedData=merge(OTG_GWASCata_Disease_PhecodeDF,OTG_StudysetDF,by="diseaseIds");
print(colnames(FinalMergedData));
FinalMergedData=FinalMergedData[,c("Phecode1.2","chromosome","position")];

#Known_LabValue_GWASHits_OTG_GWASCatabyDiseaseIDs_0204_2026.tsv
#KnownGWASHits_OTG_ByEFO_0202_2026.tsv
fwrite(
   FinalMergedData,
  "AGD35K_Phecode122_0308_2026/Known_PheCode_GWASHits_OTG_GWASCatabyDiseaseIDs_0308_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );
#KnownGWASHits_OTG_ByTraitFromSource_0202_2026.tsv
OTG_StudysetDF$trait_OTG <- sub("\\s*\\(.*\\)", "", OTG_StudysetDF$traitFromSource);
FinalMergedData=merge(OTG_GWASCata_Disease_PhecodeDF,OTG_StudysetDF,by="trait_OTG");
print(colnames(FinalMergedData));
FinalMergedData=FinalMergedData[,c("Phecode1.2","chromosome","position")];
fwrite(
   FinalMergedData,
   "AGD35K_Phecode122_0308_2026/Known_PheCode_GWASHits_OTG_GWASCatabyTraitFromSource_0308_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );

#gwas-catalog-download-associations-alt-full.tsv
GWASCatalog_StudysetDF <- read.table("AGD35K_Phecode122_0308_2026/gwas-catalog-download-associations-alt-full.tsv.gz",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
  );
print("GWASCatalog_StudysetDF");
print(dim(GWASCatalog_StudysetDF));
#[1] "DATE.ADDED.TO.CATALOG"      "PUBMEDID"
GWASCatalog_StudysetDF  =as.data.frame(GWASCatalog_StudysetDF);
GWASCatalog_StudysetDF=GWASCatalog_StudysetDF[,c("DISEASE.TRAIT",
		"CHR_ID","CHR_POS","MAPPED_TRAIT_URI")];

#GWASCatalog_KnownHitsDF <- GWASCatalog_KnownHitsDF[
#  grepl("Phecode", GWASCatalog_KnownHitsDF$DISEASE.TRAIT, ignore.case = TRUE),
#];
#print(c("Phecode",nrow(GWASCatalog_KnownHitsDF)));
GWASCatalog_StudysetDF$MAPPED_TRAIT_URI <- sub(".*/", "", GWASCatalog_StudysetDF$MAPPED_TRAIT_URI);
print(c("raw1",nrow(GWASCatalog_StudysetDF)));
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[
  !grepl(" x ", GWASCatalog_StudysetDF$CHR_ID, ignore.case = TRUE),
]#remove SNP x SNP interaction
print(c("raw2",nrow(GWASCatalog_StudysetDF)));
GWASCatalog_StudysetDF1= GWASCatalog_StudysetDF %>%
  filter(!grepl(";", CHR_ID));
print(c("raw51",nrow(GWASCatalog_StudysetDF1)));
GWASCatalog_StudysetDF2= GWASCatalog_StudysetDF %>%
  filter(grepl(";", CHR_ID));
print(c("raw5",nrow(GWASCatalog_StudysetDF2)));
# Step 1: convert to character
GWASCatalog_StudysetDF2 <- GWASCatalog_StudysetDF2 %>%
  mutate(
    CHR_ID = as.character(CHR_ID),
    CHR_POS = as.character(CHR_POS)
  )

# Step 2: remove completely empty positions
GWASCatalog_StudysetDF2 <- GWASCatalog_StudysetDF2 %>%
  filter(CHR_ID != "" & CHR_POS != "")

# Step 3: split multi-locus
GWASCatalog_StudysetDF2 <- GWASCatalog_StudysetDF2 %>%
  mutate(
    CHR_ID = strsplit(CHR_ID, ";"),
    CHR_POS = strsplit(CHR_POS, ";")
  ) %>%
  filter(lengths(CHR_ID) == lengths(CHR_POS)) %>%  # keep only valid rows
  unnest(c(CHR_ID, CHR_POS)) %>%
  mutate(
    CHR_ID = trimws(CHR_ID),
    CHR_POS = as.numeric(trimws(CHR_POS))
  )
GWASCatalog_StudysetDF=rbind(GWASCatalog_StudysetDF1,GWASCatalog_StudysetDF2);
rownames(GWASCatalog_StudysetDF)=NULL;
GWASCatalog_StudysetDF=as.data.frame(GWASCatalog_StudysetDF);
print(c("raw3",nrow(GWASCatalog_StudysetDF)));
colnames(GWASCatalog_StudysetDF)[colnames(GWASCatalog_StudysetDF) == "DISEASE.TRAIT"] <- "traitFromSource";
colnames(GWASCatalog_StudysetDF)[colnames(GWASCatalog_StudysetDF) == "MAPPED_TRAIT_URI"] <- "diseaseIds";
colnames(GWASCatalog_StudysetDF)[colnames(GWASCatalog_StudysetDF) == "CHR_ID"] <- "chromosome";
colnames(GWASCatalog_StudysetDF)[colnames(GWASCatalog_StudysetDF) == "CHR_POS"] <- "position";
#GWASCatalog_KnownHitsDF$trait <- sub("\\s*\\(.*\\)", "", GWASCatalog_KnownHitsDF$traitFromSource)
#GWASCatalog_KnownHitsDF$phecode <- sub(".*\\(|\\)", "", GWASCatalog_KnownHitsDF$traitFromSource)
#GWASCatalog_KnownHitsDF$phecode <- sub("\\)", "", GWASCatalog_KnownHitsDF$phecode)
#GWASCatalog_KnownHitsDF$phecode <- gsub("^PheCode\\s*", "", GWASCatalog_KnownHitsDF$phecode);
#GWASCatalog_KnownHitsDF$phecode  =as.numeric(as.character(GWASCatalog_KnownHitsDF$phecode));
#GWASCatalog_KnownHitsDF=GWASCatalog_KnownHitsDF[,c("diseaseIds","phecode","traitFromSource","CHR_ID","CHR_POS")];
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[
  !apply(GWASCatalog_StudysetDF[ , c("chromosome","position")], 1,
         function(x) all(is.na(x) | x=="")),
]#remove the row with chr="" and pos=""
print(c("raw4",nrow(GWASCatalog_StudysetDF)));
GWASCatalog_StudysetDF=unique(GWASCatalog_StudysetDF);
#print(str(GWASCatalog_StudysetDF));
print(str(GWASCatalog_StudysetDF));
GWASCatalog_StudysetDF=GWASCatalog_StudysetDF[,c("diseaseIds",
                "chromosome","position","traitFromSource")];

print(colnames(GWASCatalog_StudysetDF));
print("By DiseaseID");
FinalMergedData_GWASCatalog=merge(OTG_GWASCata_Disease_PhecodeDF,
				  GWASCatalog_StudysetDF,
				  by="diseaseIds");
print(nrow(FinalMergedData_GWASCatalog));
GWASCatalog_StudysetDF <- na.omit(GWASCatalog_StudysetDF);
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[GWASCatalog_StudysetDF$chromosome != "", ]
print("remove missing row");
print(nrow(FinalMergedData_GWASCatalog));
FinalMergedData_GWASCatalog=FinalMergedData_GWASCatalog[,c("Phecode1.2","chromosome","position")];
fwrite(
   FinalMergedData_GWASCatalog,
   "AGD35K_Phecode122_0308_2026/Known_PheCode_GWASHits_GWASCatabyDiseaseIDs_0308_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );

print("By Trait");
GWASCatalog_StudysetDF$trait_Catalog <- sub("\\s*\\(.*\\)", "", GWASCatalog_StudysetDF$traitFromSource);
FinalMergedData_GWASCatalog=merge(OTG_GWASCata_Disease_PhecodeDF,GWASCatalog_StudysetDF,
		      by.x="trait_OTG",by.y="trait_Catalog");

print(nrow(FinalMergedData_GWASCatalog));
print("remove missing row");
GWASCatalog_StudysetDF <- na.omit(GWASCatalog_StudysetDF);
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[GWASCatalog_StudysetDF$chromosome != "", ]
print(nrow(FinalMergedData_GWASCatalog));

FinalMergedData_GWASCatalog=FinalMergedData_GWASCatalog[,c("Phecode1.2","chromosome","position")];
fwrite(
   FinalMergedData_GWASCatalog,
   "AGD35K_Phecode122_0308_2026/Known_PheCode_GWASHits_GWASCatabyTraitFromSource_0308_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );
