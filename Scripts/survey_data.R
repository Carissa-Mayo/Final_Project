library(tidyverse)
library(readr)
library(fuzzyjoin)

# %%%%%%%%%%%%%%%%%%%%%%%%%%% cleaning raw data %%%%%%%%%%%%%%%%%%%%%%%

###### population data ###### 

maya_qpcr_raw_dat <- read_csv("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Raw/merged_mya_qpcr.csv")

all_dat <- maya_qpcr_raw_dat %>% filter(Sample_type != "eDNA")
columns <- c("Clam_ID", "Sample_date", "Repeat_number", "CANCER_value1", "CANCER_value2", "CANCER_value3", 
             "N1N2_value1", "N1N2_value2", "N1N2_value3", "Date_received", "Source", "Location")
all_dat <-subset(maya_qpcr_raw_dat, select = columns)
names(all_dat) <- make.names(names(all_dat))

# calculating cancer fractions for each sample
LOD = 1
all_dat <- all_dat %>% 
  mutate_at(c('N1N2_value1', 'N1N2_value2', 'N1N2_value3', 'CANCER_value1', 'CANCER_value2', 'CANCER_value3'), as.numeric)
all_dat <- all_dat %>%
  rowwise() %>%
  mutate(N1N2_average = mean(c(N1N2_value1, N1N2_value2, N1N2_value3), na.rm=TRUE)) %>%
  mutate(CANCER_average = ifelse(mean(c(CANCER_value1, CANCER_value2, CANCER_value3), na.rm=TRUE) >= LOD, mean(c(CANCER_value1, CANCER_value2, CANCER_value3), na.rm=TRUE), 0)) %>%
  mutate(Amp = sum(!is.na(c(CANCER_value1, CANCER_value2, CANCER_value3)) & c(CANCER_value1, CANCER_value2, CANCER_value3) > 0)) %>% 
  mutate(Fraction = ((CANCER_average/N1N2_average)/(1-(CANCER_average/N1N2_average))))
all_dat <- all_dat[!is.na(all_dat$CANCER_average), ]
all_dat <- all_dat[is.na(all_dat$Repeat_number), ]

# only want the initial samples from bi-monthly collections
all_dat$Date_received <- as.Date(all_dat$Date_received, format="%m/%d/%Y")
all_dat$Sample_date <- as.Date(all_dat$Sample_date, format="%Y-%m-%d")
all_dat <- all_dat[which(all_dat$Date_received == all_dat$Sample_date), ]

# infectious threshold found from earlier analysis
infect_threshold <- 0.2412

###### eDNA data ###### 

edna_site_raw_dat <- read_csv("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Raw/edna_qpcr.csv")
all_edna <- edna_site_raw_dat %>% 
  mutate_at(c('CANCER_value1', 'CANCER_value2', 'CANCER_value3'), as.numeric)

LOD = 1
all_edna <- all_edna %>%
  rowwise() %>%
  mutate(CANCER_average = ifelse(mean(c(CANCER_value1, CANCER_value2, CANCER_value3), na.rm=TRUE) >= LOD, mean(c(CANCER_value1, CANCER_value2, CANCER_value3), na.rm=TRUE), LOD)) %>%
  mutate(Amp = sum(!is.na(c(CANCER_value1, CANCER_value2, CANCER_value3)) & 
                     c(CANCER_value1, CANCER_value2, CANCER_value3) > 0)) %>%
  mutate(CANCER_sd = sd(c(CANCER_value1, CANCER_value2, CANCER_value3), na.rm = TRUE),) %>% 
  mutate(CopyPerMl = case_when(CANCER_average/10 > 0.1 ~ CANCER_average/10,
                               CANCER_average/10 <= 0.1 ~ 0.1))

# %%%%%%%%%%%%%%%%%%%%%% Site-specific Prep Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_counts <- function(data) {
  df <- data %>% group_by(Date_received) %>% 
    summarise(
      S_num = sum(Fraction == 0), 
      E_num = sum((Fraction > 0) & (Fraction < infect_threshold)), 
      I_num = sum(Fraction >= infect_threshold),
      count = n()
    )
  return(df)
}

get_props <- function(data) {
  df <- data %>% 
    rowwise() %>% 
    mutate(S = (S_num / count)) %>% 
    mutate(E = (E_num / count)) %>% 
    mutate(I = (I_num / count))
  df <- df %>% select(-c("S_num", "E_num", "I_num"))
  return(df)
}

get_edna <- function(data) {
  df <- data %>% group_by(Sample_date) %>% 
    summarize(mean_CopyPerML = mean(CopyPerMl, na.rm = TRUE), sd_CopyPerML = sd(CopyPerMl, na.rm = TRUE)) 
  return(df)
}

get_stan_df <- function(data, edna, write, name) {
  df <- fuzzy_left_join(data, edna, by = c("Date_received" = "Sample_date"),
                        match_fun = list(function(x,y) abs(difftime(x, y, units = "days")) <= 5))
  df <- df %>% select(-c("Sample_date"))
  
  if (write == 1) {
    write_csv(df, file = paste("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Data/", name))
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%% Forest River Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gmgi_forest_dat <- all_dat %>% filter(Source == "Gloucester Marine Genomics Institute"| Source == "GMGI", Location == "Forest River Park, MA")
gmgi_forest_edna <- all_edna %>% filter(grepl("^MAFR-\\d$", Sample_ID))

#### count data ####
gmgi_forest_count <- get_counts(gmgi_forest_dat)

#### proportional data ####
gmgi_forest_prop <- get_props(gmgi_forest_count)

#### eDNA ####
gmgi_forest_edna <- get_edna(gmgi_forest_edna)

#### all together ####
gmgi_forest_result_eDNA <- get_stan_df(gmgi_forest_count, gmgi_forest_edna, 1, "count_FR_edna.csv")

# %%%%%%%%%%%%%%%%%%%%%%%%%%% Pavilion Beach Data %%%%%%%%%%%%%%%%%%%%%%%%%%% 

gmgi_pb_dat <- all_dat %>% filter(Source == "Gloucester Marine Genomics Institute"| Source == "GMGI", Location == "Pavilion Beach, MA")
gmgi_pb_edna <- all_edna %>% filter(grepl("^MAPB-\\d$", Sample_ID))

#### count data ####
gmgi_pb_count <- get_counts(gmgi_pb_dat)

#### proportional data ####
gmgi_pb_prop <- get_props(gmgi_pb_count)

#### eDNA ####
gmgi_pb_edna <- get_edna(gmgi_pb_edna)

#### all together ####
gmgi_pb_result_eDNA23 <- get_stan_df(gmgi_pb_count, gmgi_pb_edna, 1, "count_PB_edna.csv")

