#Project: RSTR the sequel - Luminex Data
#Nathan Kieswetter


# RSTR_luminex analysis
#Aim: Assess cytokine secretion in the face of Mtbl, and PP! stimulation from the RSTR cohort

############################################################################
#NOTE: Well assignment for plate layout 1 (i.e. duplicate row of standards)#
############################################################################

#Load in the required libraries

library(tidyverse)
library(here)
library(readr)
library(tibble)
library(readxl) # Read in excel document with multiple sheets

#First, specify the path of the file you're interested in cleaning up

path <- "G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/data/luminex/Nathan Kieswetter Invitrogen Hu 27plex 071922_batch6.xlsx"

batch <- str_sub(path,-6,-6)

# Before reading data, we will return the names of the sheets for later use:

sheets <- excel_sheets(path)

# Next,we're going to extract all the sheet, define the data rectangle (A10:M94) and concatenate these into one data frame.

xl_list <- lapply(excel_sheets(path), read_excel, path = path, range = "A10:M94")

xl_list <- lapply(seq_along(sheets), function(i) {
  data.frame(sheet = I(sheets[i]), xl_list[[i]])
})

xl_list <- do.call(rbind, xl_list)

length(xl_list)

# Add PTID information based on row designation (from plate layout) - this will change on a per plate basis!
# Plate layout is saved in "Metadata" folder

data <- xl_list

if(batch == 1){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "C") ~ "52",
      startsWith(Well, "D") ~ "56",
      startsWith(Well, "E") ~ "76",
      startsWith(Well, "F") ~ "95",
      startsWith(Well, "G") ~ "96",
      startsWith(Well, "H") ~ "97"
    ))
} else if(batch == 3){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "C") ~ "58",
      startsWith(Well, "D") ~ "88",
      startsWith(Well, "E") ~ "128",
      startsWith(Well, "F") ~ "133",
      startsWith(Well, "G") ~ "148",
      startsWith(Well, "H") ~ "180"
    ))  
} else if(batch == 5){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "C") ~ "294",
      startsWith(Well, "D") ~ "297",
      startsWith(Well, "E") ~ "301",
      startsWith(Well, "F") ~ "310",
      startsWith(Well, "G") ~ "360",
      startsWith(Well, "H") ~ "380"
    ))
} else if(batch == 6){
  data <- data %>%
    mutate(PTID = case_when(
      startsWith(Well, "C") ~ "111",
      startsWith(Well, "D") ~ "150",
      startsWith(Well, "E") ~ "377",
      startsWith(Well, "F") ~ "187",
      startsWith(Well, "G") ~ "284",
      startsWith(Well, "H") ~ "410"
    ))
} else {
  print("Not an appropriate plate layout for this pipeline. Double check file path and bacth number")
}

#Add stimulation condition based on well designation (from plate layout)

dmso_wellid <- c("C1", "D1", "E1", "F1", "G1", "H1",
                 "C2", "D2", "E2", "F2", "G2", "H2",
                 "C3", "D3", "E3", "F3", "G3", "H3",
                 "C4", "D4", "E4", "F4", "G4", "H4"
)


pp1_wellid <- c("C5", "D5", "E5", "F5", "G5", "H5",
                "C6", "D6", "E6", "F6", "G6", "H6",
                "C7", "D7", "E7", "F7", "G7", "H7",
                "C8", "D8", "E8", "F8", "G8", "H8"
)

mtbl_wellid <- c("C9", "D9", "E9", "F9", "G9", "H9",
                 "C10", "D10", "E10", "F10", "G10", "H10",
                 "C11", "D11", "E11", "F11", "G11", "H11",
                 "C12", "D12", "E12", "F12", "G12", "H12"
)

data <- data %>%
  mutate(stimulation = case_when(Well %in% dmso_wellid ~ "DMSO",
                                 Well %in% pp1_wellid ~ "PP1",
                                 Well %in% mtbl_wellid ~ "MTBL"
                              
  ))

# Add time point information based on well designation (from plate layout)

timepoint_6h_wellid <- c("C1","D1","E1", "F1", "G1", "H1",
                         "C5","D5","E5", "F5", "G5", "H5",
                         "C9","D9","E9", "F9", "G9", "H9"
                         )

timepoint_12h_wellid <- c("C2","D2","E2", "F2", "G2", "H2",
                          "C6","D6","E6", "F6", "G6", "H6",
                          "C10","D10","E10", "F10", "G10", "H10"
)

timepoint_24h_wellid <- c("C3","D3","E3", "F3", "G3", "H3",
                         "C7","D7","E7", "F7", "G7", "H7",
                         "C11","D11","E11", "F11", "G11", "H11"
)

timepoint_48h_wellid <- c("C4","D4","E4", "F4", "G4", "H4",
                         "C8","D8","E8", "F8", "G8", "H8",
                         "C12","D12","E12", "F12", "G12", "H12"
)

data <- data %>%
  mutate(timepoint = case_when(Well %in% timepoint_6h_wellid ~ "6",
                               Well %in% timepoint_12h_wellid ~ "12",
                               Well %in% timepoint_24h_wellid ~ "24",
                               Well %in% timepoint_48h_wellid ~ "48"
  ))

# Assign correct dilution factor based on experimental method:
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

#Remove standard wells

std_wells <- c("A12,B12", "A1,B1", "A2,B2", "A3,B3", "A4,B4", "A5,B5", "A6,B6",
               "A7,B7", "A8,B8", "A9,B9", "A10,B10", "A11,B11")

data <- data[!(data$Well == std_wells),]

# Remove "*", "OOR <" and "OOR >" from the Obs.conc column of data and coerce them to class numeric

data$Obs.Conc <- gsub("OOR <", 0, data$Obs.Conc)
data$Obs.Conc <- gsub("\\*", "", data$Obs.Conc)
data$Obs.Conc <- gsub("OOR >",9999, data$Obs.Conc)
data$Obs.Conc <- as.numeric(data$Obs.Conc)


# Upper limits to use "OOR >"

# data$Obs.Conc = "OOR >"
# 'GM-CSF (44)' = 79413.894
# 'INF alpha (48)' = 3165.323
# 'INF gamma (43)' = 2515.398
# 'IL-1 alpha (62)' = 491.069
# 'IL-1 beta (18)' = 11443.806
# 'IL-10 (28)' = 7406.246
# 'IL-12p70 (34)' = 8016.500
# 'IL-13 (35)' = 17507.536
# 'IL-15 (65)' = 13669.643
# 'IL-17A (36)' = 9012.441
# 'IL-18 (66)' = 60374.196
# 'IL-1RA (38)' = 122339.068
# 'IL-2 (19)' = 32364.926
# 'IL-21 (72)' = 63516.448
# 'IL-22 (76)' = 106537.096
# 'IL-23 (63)' = 61155.405
# 'IL-27 (14)' = 3941.955
# 'IL-31 (37)' = 59699.378
# 'IL-4 (20)' = 52637.925
# 'IL-5 (21)' = 37429.828
# 'IL-6 (25)' = 11592.523
# 'IL-7 (26)' = 2889.529
# 'IL-8 (27)' = 9326.592
# 'IL-9 (52)' = 33890.535
# 'IP-10 (22)' = 5487.886
# 'TNF-alpha (45)' = 20859.008
# 'TNF-beta (54)' = 33150.281

# Calculate correct concentration/MF (i.e. taking into account he DF)

data <- mutate(data, final_conc =  as.numeric(data$Obs.Conc)/as.numeric(data$df_exp))
data <- mutate(data, final_MFI =  as.numeric(data$FI...Bkgd)/as.numeric(data$df_exp))

                                    
#Assign batch number to each plate

data <- data %>%
       mutate(batch = batch)  #Change this as you go!

#Remove uninformative columns

data <- select(data, -c(sheet, Type, Std.Dev, X.CV, Exp.Conc, X.Obs.Exp....100, Dilution))

data[is.na(data)] = 0 #Maybe don't use this

# Remove background DMSO signal from corresponding conditions - FIGURE THIS OUT For loop this for each cytokine

GMCSF_dmso <- data %>% 
  filter(Analyte == "GM-CSF (44)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))

## Create a separate data frame that can be used to calculate the background (DMSO) corrected cytokine values/MFI:

df_DMSO <- GMCSF_dmso

IFNA_dmso <- data %>% 
  filter(Analyte == "INF alpha (48)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IFNA_dmso)

IFNY_dmso <- data %>% 
  filter(Analyte == "INF gamma (43)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IFNY_dmso)

IL1a_dmso <- data %>% 
  filter(Analyte == "IL-1 alpha (62)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL1a_dmso)

IL1b_dmso <- data %>% 
  filter(Analyte == "IL-1 beta (18)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL1b_dmso)

IL10_dmso <- data %>% 
  filter(Analyte == "IL-10 (28)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL10_dmso)

IL12p70_dmso <- data %>% 
  filter(Analyte == "IL-12p70 (34)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL12p70_dmso)

IL13_dmso <- data %>% 
  filter(Analyte == "IL-13 (35)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL13_dmso)

IL15_dmso <- data %>% 
  filter(Analyte == "IL-15 (65)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL15_dmso)

IL17A_dmso <- data %>% 
  filter(Analyte == "IL-17A (36)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL17A_dmso)

IL18_dmso <- data %>% 
  filter(Analyte == "IL-18 (66)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL18_dmso)

IL1RA_dmso <- data %>% 
  filter(Analyte == "IL-1RA (38)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL1RA_dmso)

IL2_dmso <- data %>% 
  filter(Analyte == "IL-2 (19)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL2_dmso)

IL21_dmso <- data %>% 
  filter(Analyte == "IL-21 (72)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL21_dmso)

IL22_dmso <- data %>% 
  filter(Analyte == "IL-22 (76)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL22_dmso)

IL23_dmso <- data %>% 
  filter(Analyte == "IL-23 (63)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL23_dmso)

IL27_dmso <- data %>% 
  filter(Analyte == "IL-27 (14)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL27_dmso)

IL31_dmso <- data %>% 
  filter(Analyte == "IL-31 (37)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL31_dmso)

IL4_dmso <- data %>% 
  filter(Analyte == "IL-4 (20)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL4_dmso)

IL5_dmso <- data %>% 
  filter(Analyte == "IL-5 (21)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL5_dmso)

IL6_dmso <- data %>% 
  filter(Analyte == "IL-6 (25)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL6_dmso)

IL7_dmso <- data %>% 
  filter(Analyte == "IL-7 (26)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL7_dmso)

IL8_dmso <- data %>% 
  filter(Analyte == "IL-8 (27)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL8_dmso)

IL9_dmso <- data %>% 
  filter(Analyte == "IL-9 (52)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IL9_dmso)

IP10_dmso <- data %>% 
  filter(Analyte == "IP-10 (22)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, IP10_dmso)

TNFa_dmso <- data %>% 
  filter(Analyte == "TNF-alpha (45)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, TNFa_dmso)

TNFb_dmso <- data %>% 
  filter(Analyte == "TNF-beta (54)" & stimulation == "DMSO") %>% 
  slice(rep(1:n(),3))
df_DMSO <- bind_rows(df_DMSO, TNFb_dmso)

df_DMSO[is.na(df_DMSO)] = 0

df_DMSO <- df_DMSO %>% 
  select(final_conc, final_MFI) %>%
  rename("DMSO_final_conc" = "final_conc",
         "DMSO_final_MFI" = "final_MFI")

## Next subtract the DMSO from the cytokine signal and create a new column filled with this data

data <- bind_cols(data, df_DMSO)

data <- data %>% 
        mutate(conc_bkgd_corr = data$final_conc - data$DMSO_final_conc) %>%
        mutate(mfi_bkgd_corr = as.numeric(data$final_MFI) - as.numeric(data$DMSO_final_MFI)) 
  
#Export data for further analysis
if(batch == 1){
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\batch1_clean.csv", row.names = FALSE) 
} else if (batch == 3) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\batch3_clean.csv", row.names = FALSE)
}else if (batch == 5) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\batch5_clean.csv", row.names = FALSE)
}else if (batch == 6) {
  write.csv(data,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\batch6_clean.csv", row.names = FALSE)
}
