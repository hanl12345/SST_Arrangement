#install.packages("arrow")
#install.packages("R.utils")
library(arrow)
library(dplyr)
library(data.table)
library(tidyr)
library(R.utils);
library(readr);
library(stringr)
library(dplyr)
library(tidyr)
# -----------------------------
# 1. Check input path
# -----------------------------
GWASCatalogueDF <- fread("AGD35K_Phecode122_0308_2026/gwas-catalog-download-associations-alt-full.tsv.gz", sep = "\t",quote = "");
GWASCatalogueDF  =as.data.frame(GWASCatalogueDF);
GWASCatalogueDF=GWASCatalogueDF[,c(8,35,36)];
#DISEASE/TRAIT[8],MAPPED_TRAIT[35]    MAPPED_TRAIT_URI[36]
colnames(GWASCatalogueDF)=c("traitFromSource","name","EFO_URL");
#"name","description"
#print(str(GWASCatalogueDF));
print(nrow(GWASCatalogueDF));
# filter rows containing "Phecode"
GWASCatalogueDF=unique(GWASCatalogueDF);

GWASCatalogueDF$diseaseIds <- sub(".*/", "", GWASCatalogueDF$EFO_URL);
GWASCatalogueDF=GWASCatalogueDF[,c("traitFromSource","diseaseIds","name")];
GWASCatalogueDF$description ="NA";
print("GWASCatalogueDF");
print(str(GWASCatalogueDF));
GWASCatalogueDF=unique(GWASCatalogueDF);
print("GWASCatalogueDF dimension");
print(c("Number of GWASCatalogueDF unique",nrow(GWASCatalogueDF)));

print(colnames(GWASCatalogueDF));
#diseaseIds","traitFromSource","phecode"

OTGDVIDDF <- fread("AGD35K_Phecode122_0308_2026/disease_variant_study_Merged_OTGV25.12.tsv.gz",
				       	sep = "\t",quote = "");
#print(str(OTGDVIDDF));
OTGDVIDDF=OTGDVIDDF[,c("traitFromSource","diseaseIds","name","description")];
print("OTGDVIDDF dimension");
OTGDVIDDF=unique(OTGDVIDDF);
#OTGDVIDDF=OTGDVIDDF[,c("diseaseIds","trait","phecode")];
#OTGDVIDDF1=OTGDVIDDF[,c("diseaseIds","phecode")];
#Remove the duplicate line with same diseaseIds and phecode  that have existed in OTGDVIDDF 
#At first keep the lines with its phecode from agd369 phecode
print("OTGDVIDDF dimension");
print(c("Number of OTGDVIDDF unique",nrow(OTGDVIDDF)));
#print(str(OTGDVIDDF));
AGD35K_LabValueDF <- fread("AGD35K_Phecode122_0308_2026/InferredMatch_GWASTraitVsLabValue_BySameString.tsv.gz",sep = "\t",quote = "");
AGD35K_LabValueDF =as.data.frame(AGD35K_LabValueDF);
#print(colnames(AGD35K_LabValueDF));
AGD35K_OCTset=merge(AGD35K_LabValueDF,OTGDVIDDF,by="traitFromSource");
AGD35K_OCTset$MatchMethod="OTG_Exact";
#print(str(AGD35K_OCTset));
AGD35K_GWASCataset=merge(AGD35K_LabValueDF,GWASCatalogueDF,by="traitFromSource");
AGD35K_GWASCataset <- AGD35K_GWASCataset %>%
  anti_join(AGD35K_OCTset, by = c("diseaseIds", "traitFromSource"));
print("AGD35K_GWASCataset_after removing OTGset");
AGD35K_GWASCataset$MatchMethod="GWASCata_Exact";
#print(str(AGD35K_GWASCataset));
OTG_GWASCata_LabValueDF=rbind(AGD35K_OCTset,AGD35K_GWASCataset);
rownames(OTG_GWASCata_LabValueDF)=NULL;
OTG_GWASCata_LabValueDF=as.data.frame(OTG_GWASCata_LabValueDF);
print(str(OTG_GWASCata_LabValueDF));
colnames(OTG_GWASCata_LabValueDF)[colnames(OTG_GWASCata_LabValueDF) == "LabTrait"] <- "LabTrait_AGD35K";
colnames(OTG_GWASCata_LabValueDF)[colnames(OTG_GWASCata_LabValueDF) == "GWAStrait"] <- "name_AGD35K";
colnames(OTG_GWASCata_LabValueDF)[colnames(OTG_GWASCata_LabValueDF) == "name"] <- "name_EFO";
colnames(OTG_GWASCata_LabValueDF)[colnames(OTG_GWASCata_LabValueDF) == "description"] <- "description_EFO";

OTG_GWASCata_LabValueDF$trait_OTG=OTG_GWASCata_LabValueDF$traitFromSource;
#colnames(OTG_GWASCata_LabValueDF)[colnames(OTG_GWASCata_LabValueDF) == "traitFromSource"] <- "trait_OTG";
#colnames(OTG_GWASCata_Disease_PhecodeDF)[colnames(OTG_GWASCata_Disease_PhecodeDF) == "description"] <- "description_OTG";
#traitFromSource LabTrait_AGD35K name_AGD35K     diseaseIds      name    description     MatchMethod
#phecode diseaseIds      trait_OTG       MatchMethod     name_OTG        description_OTG phecode1.2_simpleLabel  phecode1.2_category
ColumnOrder=c("LabTrait_AGD35K","diseaseIds","trait_OTG","MatchMethod",
	      "name_EFO", "description_EFO","name_AGD35K");
OTG_GWASCata_LabValueDF1=OTG_GWASCata_LabValueDF[,ColumnOrder];
fwrite(
   OTG_GWASCata_LabValueDF1,
  "AGD35K_Phecode122_0308_2026/AGD35K_OTG_GWASCata_Disease_Labvalue.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );
print("creating file");

OTGKnownGWASSignalDF <- fread("AGD35K_Phecode122_0308_2026/disease_variant_study_Merged_OTGV25.12.tsv.gz");
OTGKnownGWASSignalDF =as.data.frame(OTGKnownGWASSignalDF);
OTGKnownGWASSignalDF=OTGKnownGWASSignalDF[,c("diseaseIds","variantId",
		"chromosome","position","traitFromSource")];
print(colnames(OTGKnownGWASSignalDF));
OTG_GWASCata_LabValueDF$trait_OTG=NULL;
print("OTG_GWASCata_LabValueDF");
print(colnames(OTG_GWASCata_LabValueDF));
OTG_GWASCata_LabValue_KnownHitsMatchBytraitFromSource=merge(OTGKnownGWASSignalDF,OTG_GWASCata_LabValueDF,by="traitFromSource");
print("here");
print(colnames(OTG_GWASCata_LabValue_KnownHitsMatchBytraitFromSource));
OTG_GWASCata_LabValue_KnownHitsMatchBytraitFromSource=OTG_GWASCata_LabValue_KnownHitsMatchBytraitFromSource[,c("LabTrait_AGD35K","chromosome","position")];

fwrite(
   OTG_GWASCata_LabValue_KnownHitsMatchBytraitFromSource,
  "AGD35K_Phecode122_0308_2026/Known_LabValue_GWASHits_OTG_GWASCatabyTraitFromSource_0204_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );

OTG_GWASCata_LabValue_KnownHitsMatchByDiseaseIds=merge(OTGKnownGWASSignalDF,OTG_GWASCata_LabValueDF,by="diseaseIds");
print(colnames(OTG_GWASCata_LabValue_KnownHitsMatchByDiseaseIds));
OTG_GWASCata_LabValue_KnownHitsMatchByDiseaseIds=OTG_GWASCata_LabValue_KnownHitsMatchByDiseaseIds[,c("LabTrait_AGD35K","chromosome","position")];
fwrite(
   OTG_GWASCata_LabValue_KnownHitsMatchByDiseaseIds,
   "AGD35K_Phecode122_0308_2026/Known_LabValue_GWASHits_OTG_GWASCatabyDiseaseIDs_0204_2026.tsv.gz",
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
print(c("raw4",nrow(GWASCatalog_StudysetDF)));
#print(unique(GWASCatalog_StudysetDF$chromosome));
#print(unique(GWASCatalog_StudysetDF$position));

GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[
  !apply(GWASCatalog_StudysetDF[ , c("chromosome","position")], 1,
         function(x) all(is.na(x) | x=="")),
]#remove the row with chr="" and pos=""
print(c("raw4",nrow(GWASCatalog_StudysetDF)));
GWASCatalog_StudysetDF=unique(GWASCatalog_StudysetDF);
#print(str(GWASCatalog_StudysetDF));
GWASCatalog_StudysetDF=GWASCatalog_StudysetDF[,c("diseaseIds",
                "chromosome","position","traitFromSource")];
GWASDiseaseIDs= unique(GWASCatalog_StudysetDF$diseaseIds);
print("By DiseaseID_1");
#print(str(OTG_GWASCata_LabValueDF));
print(str(OTG_GWASCata_LabValueDF))

FinalMergedData_GWASCatalog=merge(OTG_GWASCata_LabValueDF,GWASCatalog_StudysetDF,by="diseaseIds");
print(nrow(FinalMergedData_GWASCatalog));
GWASCatalog_StudysetDF <- na.omit(GWASCatalog_StudysetDF);
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[GWASCatalog_StudysetDF$chromosome != "", ]
print(nrow(FinalMergedData_GWASCatalog));

FinalMergedData_GWASCatalog=FinalMergedData_GWASCatalog[,c("LabTrait_AGD35K","chromosome","position")];
fwrite(
   FinalMergedData_GWASCatalog,
   "AGD35K_Phecode122_0308_2026/Known_LabValue_GWASHits_GWASCatabyDiseaseIDs_0206_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );

print("By Trait");
FinalMergedData_GWASCatalog=merge(OTG_GWASCata_LabValueDF,GWASCatalog_StudysetDF,by="traitFromSource");

print(nrow(FinalMergedData_GWASCatalog));
print("remove missing row");
GWASCatalog_StudysetDF <- na.omit(GWASCatalog_StudysetDF);
GWASCatalog_StudysetDF <- GWASCatalog_StudysetDF[GWASCatalog_StudysetDF$chromosome != "", ]
print(nrow(FinalMergedData_GWASCatalog));

FinalMergedData_GWASCatalog=FinalMergedData_GWASCatalog[,c("LabTrait_AGD35K","chromosome","position")];
fwrite(
   FinalMergedData_GWASCatalog,
   "AGD35K_Phecode122_0308_2026/Known_LabValue_GWASHits_GWASCatabyTraitFromSource_0206_2026.tsv.gz",
   sep = "\t",
   quote = FALSE,
   na = "NA"
 );
