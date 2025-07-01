library(readxl)
library(dplyr)
library(openxlsx)
library(googlesheets4)
library(readr)

# %%%%%%%%%%%%%%%%%%%%%%%%%% master qPCR for mya arenaria %%%%%%%%%%%%%%%%%%%%%%%%%%%%

qPCR <- read_excel('/Users/carissamayo/Metzger Lab Dropbox/Master_qPCR_datasheets/master_mya_qPCRdata_current.xlsx')

options(
  gargle_oauth_cache = ".secrets",
  gargle_oauth_email = TRUE,
  gargle_verbosity = "debug"
)
gs4_auth(cache = ".secrets", email = "cmayo@uw.edu", scopes = "https://www.googleapis.com/auth/drive")

meta <- read_sheet("https://docs.google.com/spreadsheets/d/1hmiQo29kPK_uOXzycjMn1C4htUoKEbiUBVRp3cRB9tQ/edit?usp=sharing",
                    sheet = "meta_master")

merged <- merge(qPCR, meta, by = "Clam_ID", all.x = "TRUE")
merged_df <- as.data.frame(merged)

write_csv(merged_df, "/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Raw/merged_mya_qpcr.csv")

# %%%%%%%%%%%%%%%%%%%%%%%% master eDNA for survey sites %%%%%%%%%%%%%%%%%%%%%%%%%%%%

eDNA <- read_excel('/Users/carissamayo/Metzger Lab Dropbox/Master_qPCR_datasheets/master_eDNA_qPCRdata_current.xlsx')

write_csv(eDNA, "/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Raw/edna_qpcr.csv")
