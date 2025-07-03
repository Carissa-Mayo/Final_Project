library(tidyverse)
library(readr)
library(fuzzyjoin)
library(ggplot2)
library(patchwork)
library(reshape2)

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
    write_csv(df, file = paste0("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Data/", name))
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
gmgi_pb_result_eDNA23 <- get_stan_df(gmgi_pb_count, gmgi_pb_edna, 0, "count_PB_edna.csv")


# %%%%%%%%%%%%%%%%%%%%%%%%% Figure Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prop_plot <- function(data, write, name, title) {
  long <- data %>%
    mutate(
      S_prop = S * 100,
      E_prop = E * 100,
      I_prop = I * 100
    ) %>%
    select(Date_received, count, S_prop, E_prop, I_prop) %>%
    pivot_longer(
      cols = ends_with("_prop"),
      names_to = "Compartment",
      values_to = "Proportion"
    ) %>%
    mutate(
      Compartment = recode(Compartment, S_prop = "S", E_prop = "E", I_prop = "I"),
      Compartment = factor(Compartment, levels = c("S", "E", "I")),
      Date_received = as.Date(Date_received)
    )
  
  p <- ggplot(long, aes(x = Date_received, y = Proportion, fill = Compartment)) +
    geom_text(
      data = data,
      aes(x = as.Date(Date_received), y = 105, label = count),
      inherit.aes = FALSE, size = 3
    ) +
    geom_bar(stat = "identity", position = "stack", width = 30) +
    labs(
      title = paste("Susceptible, Exposed, and Infectious Proportions - 23/24", title),
      x = "Collection Date",
      y = "Proportion of Sample",
      fill = ""
    ) +
    scale_fill_manual(values = c("S" = "palegreen3", "E" = "darkorange", "I" = "firebrick3")) +
    scale_x_date(date_labels = "%m-%d-%y", breaks = data$Date_received[seq(1, nrow(data), by = 2)]) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(angle = 18, hjust = 1),
      plot.title = element_text(size = 24)
    )
  
  if (write == 1) {
    ggsave(
      filename = paste0("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Figures/SURVEY_prop_", name, ".png"),
      plot = p, width = 16, height = 5, dpi = 300
    )
  }
  
  return(p)
}

count_plot <- function(data, write, name, title) {
  colnames(data)[2] <- "S"
  colnames(data)[3] <- "E"
  colnames(data)[4] <- "I"
  df_melt <- melt(data, id.vars = "Date_received")
  df_melt <- df_melt[df_melt$variable != "count", ]
  
  p <- ggplot(df_melt, aes(x = Date_received, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Susceptible, Exposed, and Infectious Counts - 23/24 ", title),
         x = "Collection Date",
         y = "Count",
         fill = "") +
    scale_fill_manual(values = c("S" = "palegreen3", "E" = "darkorange", "I" = "firebrick3")) +
    scale_x_date(date_labels = "%m-%d-%y", breaks = data$Date_received[seq(1, nrow(data), by = 2)]) +
    theme_minimal(base_family = "Arial") +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 18, hjust = 1),
          plot.title = element_text(size = 24))
  
  if(write == 1) {
    ggsave(filename = paste0("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Figures/SURVEY_count_", name, ".png"),
           p, width = 16, height = 5, dpi = 300)
    message("Saved count")
  }
  
  return(p)
}

edna_plot <- function(data, write, name, title) {
  p <- ggplot(data, aes(x = Sample_date, y = mean_CopyPerML)) +
    geom_line(color = "#104862") +
    geom_point(size = 2, color = "#104862") +
    #geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    #annotate("text", x = data$Sample_date[2], y = 0.15, label = "LOD = 0.1", color = "red", hjust = -0.05, size = 3) +
    #scale_y_log10() +
    scale_x_date(date_labels = "%m-%d-%y", breaks = data$Sample_date[seq(1, nrow(data), by = 2)]) +
    #scale_y_continuous(breaks = c(0.00, 0.25, 0.5, 0.75, 1.00), limits = c(0, 1.05)) +
    labs(title = paste("eDNA Concentration", title),
         x = "Collection Date",
         y = "Copies per mL") +
    theme_minimal(base_family = "Arial") +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 18, hjust = 1),
          plot.title = element_text(size = 24))
  
  if(write == 1) {
    ggsave(filename = paste0("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Figures/SURVEY_edna_", name, ".png"),
           p, width = 16, height = 5, dpi = 300)
    message("Saved edna plot")
  }
  
  return(p)
}

both_plot <- function(p1, p2, write, name, title) {
  p1 <- p1 + ggtitle("Susceptible, Exposed, and Infectious Proportions")
  p2 <- p2 + ggtitle("eDNA Concentration")
  p <- p1 + p2 + plot_layout(ncol = 2) +
    plot_annotation(
      title = paste("23/24 Survey Results from ", title),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 24,
          family = "Arial"
        )
      )
    )
  
  if(write == 1) {
    ggsave(filename = paste0("/Users/carissamayo/Metzger Lab Dropbox/Carissa Mayo/Final_Project/Figures/SURVEY_both_", name, ".png"),
           p, width = 16, height = 5, dpi = 300)
  }
  
  return(p)
}
  
# %%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%

##### Forest River #####
gmgi_forest_prop$Date_received <- gmgi_forest_edna$Sample_date
gmgi_forest_count$Date_received <- gmgi_forest_edna$Sample_date

p_fr_prop <- prop_plot(gmgi_forest_prop, 1, "FR", "Forest River, MA")
p_fr_count <- count_plot(gmgi_forest_count, 1, "FR", "Forest River, MA")
p_fr_edna <- edna_plot(gmgi_forest_edna, 1, "FR", "Forest River, MA")
both_plot(p_fr_prop, p_fr_edna, 1, "FR", "Forest River, MA")

##### Pavilion Beach #####
gmgi_pb_prop$Date_received <- gmgi_pb_edna$Sample_date[1:length(gmgi_pb_prop$Date_received)]
gmgi_pb_count$Date_received <- gmgi_pb_edna$Sample_date[1:length(gmgi_pb_prop$Date_received)]

p_pb_prop <- prop_plot(gmgi_pb_prop[1:6, ], 1, "PB", "Pavilion Beach, MA")
p_fr_count <- count_plot(gmgi_pb_count[1:6, ], 1, "PB", "Pavilion Beach, MA")
p_pb_edna <- edna_plot(gmgi_pb_edna, 1, "PB", "Pavilion Beach, MA")
both_plot(p_pb_prop, p_pb_edna, 1, "PB", "Pavilion Beach, MA")
