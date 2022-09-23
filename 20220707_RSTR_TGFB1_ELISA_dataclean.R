#Project: RSTR the sequel - Luminex Data
#Nathan Kieswetter


# RSTR_TGFb1 ELISA analysis
#Aim: Assess cytokine secretion in the face of Mtbl, and PP! stimulation from the RSTR cohort

################################################
#NOTE: Well assignment for ELISA plate layout #
###############################################

# Load in the required libraries

library(tidyverse)
library(here)
library(readr)
library(tibble)
library(readxl) # Used to Read in excel document with multiple sheets

# First, specify the path of the folder containing the ELISA output files

path <- "G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/data/elisa/20220722_RSTR_TGFB1_batch7.xlsx"

batch <- str_sub(path,-6,-6)

TGFb1_batch <- readxl::read_excel(path = path, range = "A11:F104") %>%
  mutate(batch = batch)

# Remove standard and blank wells 

TGFb1_batch <- TGFb1_batch %>%
  filter(!str_detect(Well, 'A'),
         !str_detect(Well, 'H')
  )

#Change the headings to something more R-friendly

TGFb1_batch <- TGFb1_batch %>%
rename(blank_corr_OD = "Difference based on Blank corrected and Blank corrected (calculated)",
       conc_5P_ng.mL = "5-Parameter fit based on Difference in ng/mL (calculated)",
       blank_corr_450 = "Blank corrected based on Raw Data (450 1)",
       blank_corr_570 = "Blank corrected based on Raw Data (570 2)")

# Add PTID information based on row designation (from plate layout) - this will change on a per plate basis!
# Add PTID information

data <- TGFb1_batch

if(batch == 1){
  data <- data %>%
  mutate(PTID = case_when(
    startsWith(Well, "B") ~ "52",
    startsWith(Well, "C") ~ "56",
    startsWith(Well, "D") ~ "76",
    startsWith(Well, "E") ~ "95",
    startsWith(Well, "F") ~ "96",
    startsWith(Well, "G") ~ "97",
  )
)
} else if(batch == 2){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "181",
      startsWith(Well, "C") ~ "183",
      startsWith(Well, "D") ~ "323",
      startsWith(Well, "E") ~ "332",
      startsWith(Well, "F") ~ "346",
      startsWith(Well, "G") ~ "358",
    )
)
} else if(batch == 3){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "58",
      startsWith(Well, "C") ~ "88",
      startsWith(Well, "D") ~ "128",
      startsWith(Well, "E") ~ "133",
      startsWith(Well, "F") ~ "148",
      startsWith(Well, "G") ~ "180",
    )
    )
}else if(batch == 4){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "198",
      startsWith(Well, "C") ~ "209",
      startsWith(Well, "D") ~ "217",
      startsWith(Well, "E") ~ "235",
      startsWith(Well, "F") ~ "241",
      startsWith(Well, "G") ~ "248",
    )
    )
}else if(batch == 5){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "294",
      startsWith(Well, "C") ~ "297",
      startsWith(Well, "D") ~ "301",
      startsWith(Well, "E") ~ "310",
      startsWith(Well, "F") ~ "360",
      startsWith(Well, "G") ~ "380",
    )
    )
}else if(batch == 6){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "111",
      startsWith(Well, "C") ~ "150",
      startsWith(Well, "D") ~ "377",
      startsWith(Well, "E") ~ "187",
      startsWith(Well, "F") ~ "284",
      startsWith(Well, "G") ~ "410",
    )
    )
}else if(batch == 7){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "B") ~ "361",
      startsWith(Well, "C") ~ "254"
    )
    )
}

#Add stimulation condition based on well designation (from plate layout)

dmso_wellid <- c("B01", "C01", "D01", "E01", "F01", "G01",
                 "B02","C02", "D02", "E02", "F02", "G02",
                 "B03","C03", "D03", "E03", "F03", "G03",
                 "B04","C04", "D04", "E04", "F04", "G04"
)

pp1_wellid <- c("B05","C05", "D05", "E05", "F05", "G05",
                "B06","C06", "D06", "E06", "F06", "G06",
                "B07","C07", "D07", "E07", "F07", "G07",
                "B08","C08", "D08", "E08", "F08", "G08"
)

mtbl_wellid <- c("B09","C09", "D09", "E09", "F09", "G09",
                 "B10","C10", "D10", "E10", "F10", "G10",
                 "B11","C11", "D11", "E11", "F11", "G11",
                 "B12","C12", "D12", "E12", "F12", "G12"
)

data <- data %>%
  mutate(stimulation = case_when(Well %in% dmso_wellid ~ "DMSO",
                                 Well %in% pp1_wellid ~ "PP1",
                                 Well %in% mtbl_wellid ~ "MTBL"
  ))

#  Add time point information based on well designation (from plate layout)

timepoint_6h_wellid <- c("B01","C01","D01","E01", "F01", "G01",
                         "B05","C05","D05","E05", "F05", "G05",
                         "B09","C09","D09","E09", "F09", "G09"
)

timepoint_12h_wellid <- c("B02","C02","D02","E02", "F02", "G02",
                          "B06","C06","D06","E06", "F06", "G06",
                          "B10","C10","D10","E10", "F10", "G10"
)

timepoint_24h_wellid <- c("B03","C03","D03","E03", "F03", "G03",
                          "B07","C07","D07","E07", "F07", "G07",
                          "B11","C11","D11","E11", "F11", "G11"
)

timepoint_48h_wellid <- c("B04","C04","D04","E04", "F04", "G04",
                          "B08","C08","D08","E08", "F08", "G08", 
                          "B12","C12","D12","E12", "F12", "G12"
)

data <- data %>%
  mutate(timepoint = case_when(Well %in% timepoint_6h_wellid ~ "6",
                               Well %in% timepoint_12h_wellid ~ "12",
                               Well %in% timepoint_24h_wellid ~ "24",
                               Well %in% timepoint_48h_wellid ~ "48"
  ))

# Assign correct dilution factor (df_exp) based on experimental method:
# 6 h = 1
# 12 h = 1.25
# 24 h = 1.66
# 48 h = 2.5

data <- data %>%
  mutate("df_exp" = case_when(Well %in% timepoint_6h_wellid ~ "1",
                              Well %in% timepoint_12h_wellid ~ "1.25",
                              Well %in% timepoint_24h_wellid ~ "1.66",
                              Well %in% timepoint_48h_wellid ~ "2.5"
  ))

# Remove "<< Y range", replace with a 0 value coerce them to class numeric

data$conc_5P_ng.mL <- gsub("<< Y range",0, data$conc_5P_ng.mL)

# Convert to pg/mL and correct for dilution factor i.e. the real experimental concentration for both mfi and conc.

data2 <- data %>%
  mutate(conc_pg_df_cor = as.numeric(conc_5P_ng.mL) * 1000 / df_exp) %>%
  mutate(mfi_final =  as.numeric(blank_corr_OD) / df_exp)

# Remove background DMSO signal from corresponding conditions 

TGFb1_dmso <- data %>% 
  filter(stimulation == "DMSO") %>% 
  slice(rep(1:n(),3)) %>% 
  select(conc_pg_df_cor, mfi_final) %>%
  rename(DMSO_final_conc = conc_pg_df_cor,
         DMSO_final_MFI = mfi_final)

# Next subtract the DMSO from the cytokine signal and create a new column filled with this data

data <- bind_cols(data, TGFb1_dmso)

data <- data %>% 
  mutate(bkg_corr_pg.mL = conc_pg_df_cor - DMSO_final_conc) %>%
  mutate(bkg_corr_MFI = mfi_final - DMSO_final_MFI)
  
# Export data for further analysis

if(batch == 1){
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch1_clean.csv", row.names = FALSE) 
} else if (batch == 2) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch2_clean.csv", row.names = FALSE)
}else if (batch == 3) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch3_clean.csv", row.names = FALSE)
}else if (batch == 4) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch4_clean.csv", row.names = FALSE)
}else if (batch == 5) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch5_clean.csv", row.names = FALSE)
}else if (batch == 6) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch6_clean.csv", row.names = FALSE)
}else if (batch == 7) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\elisa\\elisa_batch7_clean.csv", row.names = FALSE)
}


