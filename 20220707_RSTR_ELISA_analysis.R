library(readr)
library(tidyverse)
library(ggpubr)
library(DescTools)
library(bayestestR)
library(readxl)

date <- 20220707

# Set WD

setwd('G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/out/elisa')

# Read in all batches from a single directory and row bind them.

df_RSTR_comb <- list.files(path = 'G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/out/elisa') %>% 
  lapply(read_csv) %>% 
  bind_rows

# Add the correct PTID information so that it matches the flow PTIDS

df_RSTR_comb <- df_RSTR_comb %>% 
                mutate(PTID = recode(PTID, '52' = 'RS102052',
                                           '56' = 'RS102056',
                                           '76' = 'RS102076',
                                           '95' = 'RS102095',
                                           '96' = 'RS102096',
                                           '97' = 'RS102097',
                                           '181' = 'RS102181',
                                           '183' = 'RS102183',
                                           '323' = 'RS102323',
                                           '332' = 'RS102332',
                                           '346' = 'RS102346',
                                           '358' = 'RS102358',
                                           '361' = 'RS102361',
                                           '58' = 'RS102058',
                                           '88' = 'RS102088',
                                           '133' = 'RS102133',
                                           '148' = 'RS102148',
                                           '180' = 'RS102180',
                                           '198' = 'RS102198',
                                           '209' = 'RS102209',
                                           '217' = 'RS102217',
                                           '235' = 'RS102235',
                                           '241' = 'RS102241',
                                           '248' = 'RS102248',
                                           '254' = 'RS102254',
                                           '297' = 'RS102297',
                                           '301' = 'RS102301',
                                           '310' = 'RS102310',
                                           '360' = 'RS102360',
                                           '380' = 'RS102380',
                                           '111' = 'RS102111',
                                           '150' = 'RS102150',
                                           '377' = 'RS102377',
                                           '187' = 'RS102187',
                                           '248' = 'RS102248',
                                           '410' = 'RS102410',
                                           '128' = 'RS102128',
                                           '294' = 'RS102294',
                                           '284' = 'RS102284',
))

# Assign Group (TST+ or P_neg) to all patients

p_neg_ptids <- c("RS102052", "RS102058", "RS102096","RS102097","RS102128",
                 "RS102133","RS102148","RS102150","RS102161","RS102181",
                 "RS102183","RS102198", "RS102209","RS102241","RS102294",
                 "RS102310","RS102346","RS102351","RS102360","RS102377")

tst_pos_ptids <- c("RS102056","RS102076","RS102088","RS102095","RS102111",
                   "RS102180","RS102187","RS102217","RS102235","RS102248",
                   "RS102254","RS102284","RS102297","RS102301","RS102323",
                   "RS102332","RS102358","RS102361","RS102380","RS102410")

df_RSTR_comb <- df_RSTR_comb %>%
  mutate(group = case_when(PTID %in% p_neg_ptids ~ "P_neg",
                             PTID %in% tst_pos_ptids ~ "TST_pos"
  ))

# Plot the data

# Assess background corrected conc. over time and stratified by stimulation

TGFb1_mDMSO <- df_RSTR_comb %>%
            filter(stimulation == c("PP1", "MTBL"))

# Plot the background corrected concentration data

ggplot(TGFb1_mDMSO, aes(x = factor(timepoint), 
                     y = bkg_corr_pg.mL,
                     color = group,
                     )) +
                    geom_boxplot() +
                    facet_grid(. ~ stimulation) +
                    scale_y_continuous(trans = 'sqrt') +
                    xlab("Timepoint") + 
                    ylab("TGFb1 Concentration (pg/mL)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  theme(text = element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns"))
  )

# Plot the non-background corrected concentration data

ggplot(df_RSTR_comb, aes(x = factor(timepoint), 
                  y = conc_pg_df_cor,
                  color = group)) +
  geom_boxplot() +
  facet_grid(. ~ stimulation) +
  scale_y_continuous(trans='sqrt') +
  xlab("Timepoint") + 
  ylab("TGFb1 Concentration (pg/mL)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(text = element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns"))
                     )

# Assess background corrected OD (MFI)

ggplot(TGFb1_mDMSO, 
       aes(x = factor(timepoint), 
           y = bkg_corr_MFI,
           color = group)) + 
  labs(fill = 'Group') +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  facet_grid(. ~ stimulation) +
#  scale_y_continuous(trans='sqrt') +
  xlab("Stimulation") + 
  ylab("Blank Corrected OD450") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(text = element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))

# Assess non-background corrected OD (MFI)

ggplot(df_RSTR_comb, 
       aes(x = factor(timepoint), 
           y = mfi_final,
           color = group)) + 
  labs(fill = 'Group') +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  facet_grid(. ~ stimulation) +
  #  scale_y_continuous(trans='sqrt') +
  xlab("Stimulation") + 
  ylab("Blank Corrected OD450") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(text = element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))

# Calculate the AUC

auc_split <- split(df_RSTR_comb, f = list(df_RSTR_comb$PTID, 
                                          df_RSTR_comb$stimulation,
                                          df_RSTR_comb$batch))

# Function to calculate the AUC per the factors defines above.

## Define the AUC (MFI) for the non-background corrected secreted data.

auc_fun1 <- function(x) { # Function to define the AUC 
  timept <- x$timepoint
  mfi <- x$mfi_final
  out <- AUC(x = timept, 
             y = mfi,
             method = "trapezoid",
             absolutearea = TRUE)
  return(out)
}

df_auc1 <- sapply(auc_split, auc_fun1) %>%
  as.data.frame() %>%
  na.omit() %>%
  rownames_to_column() %>%
  setNames(c('rowname', 'AUC_nbc')) # AUC_nbc = not background corrected

df_auc_meta <- str_split_fixed(df_auc1$rowname, "\\.", 3) %>%
  as.data.frame() %>%
  setNames(c('PTID', 'stimulation', 'batch'))

df_auc <- bind_cols(df_auc1, df_auc_meta)

df_auc <- df_auc %>%
  select(-rowname)

# The following function will calculate the MFI for the background corrected AUC of the MFI over time

auc_fun2 <- function(x) {
  timept <- x$timepoint
  mfi <- x$bkg_corr_MFI
  out <- AUC(x = timept, 
             y = mfi,
             method = "trapezoid",
             absolutearea = TRUE)
  return(out)
}

df_auc2 <- sapply(auc_split, auc_fun2) %>%
  as.data.frame() %>%
  na.omit() %>%
  rownames_to_column()%>%
  setNames(c('rowname', 'AUC_bc')) %>% # AUC_bc =  background corrected
  select(-rowname)

df_auc <- bind_cols(df_auc, df_auc2)

# Re-add group metadata

p_neg_ptids <- c("RS102052", "RS102058", "RS102096","RS102097","RS102128",
                 "RS102133","RS102148","RS102150","RS102161","RS102181",
                 "RS102183","RS102198", "RS102209","RS102241","RS102294",
                 "RS102310","RS102346","RS102351","RS102360","RS102377")

tst_pos_ptids <- c("RS102056","RS102076","RS102088","RS102095","RS102111",
                   "RS102180","RS102187","RS102217","RS102235","RS102248",
                   "RS102254","RS102284","RS102297","RS102301","RS102323",
                   "RS102332","RS102358","RS102361","RS102380","RS102410")

df_auc <- df_auc %>%
  mutate(group = case_when(PTID %in% p_neg_ptids ~ "P_neg",
                           PTID %in% tst_pos_ptids ~ "TST_pos"
  ))

#Add additional Metadata using mutating joins

RSTR_metadata <- read_excel('G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/data/metadata/RSTR_sequel_metadata.xlsx') %>%
  rename('PTID' = 'Sample ID')

df_auc <- inner_join(RSTR_metadata, df_auc, by = "PTID") %>%
  select(-c('Status', 'Position', 'Random Sequence'))

# Export the AUC_file

write.csv(df_auc,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\data\\elisa\\TGFb_all_batches_auc_mfi.csv", row.names = FALSE)

# Plot the AUC values for the PTIDS by stimulation

df_auc <- read.csv("G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\data\\elisa\\TGFb_all_batches_auc_mfi.csv")

df_auc_mDMSO <- df_auc %>%
  filter(stimulation != "DMSO")

## Plot the non-background controlled AUC data

ggplot(df_auc, 
       aes(x = factor(stimulation), 
             y = AUC_nbc,
             color = group)) + 
    labs(fill = 'Group') +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    #    facet_grid(. ~ group) +
    scale_y_continuous(trans='log10') +
    xlab("Stimulation") + 
    ylab("AUC (MFI)") +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw() +
    theme(text = element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
    stat_compare_means(method = "wilcox.test")  # This works really well and returns the same result as in the stats section below

## Plot the background corrected AUC based on the MFI data

ggplot(df_auc_mDMSO, 
         aes(x = factor(stimulation), 
             y = AUC_bc,
             color = group)) + 
    labs(fill = 'Group') +
    geom_boxplot(outlier.shape = NA, fill = "white") +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) + #adds the individual data points
    #    facet_grid(. ~ stimulation) +
    scale_y_continuous(trans='log10') +
    xlab("Stimulation") + 
    ylab("AUC (Bkgnd Corrected)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw() +
    theme(text = element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#984EA3", "#4DAF4A")) + 
    stat_compare_means(method = "wilcox.test")

# UMAP

## Clean the data for UMAP using the concentration data

set.seed(date)

df_auc_umap <- df_auc %>% 
  drop_na() %>%
  #  filter(Analyte == "IL-6") %>%
  mutate(ID = row_number()) 

df_auc_meta <- df_auc_umap %>%
  select(ID, batch, group, PTID) #ADULT_RISK_TBINF, M0_KCVSEX, KCHCA_AGE_YR_CURRENT)

## Run UMAP

umap_fit <- df_auc_umap %>%
  select(-c(batch)) %>%
  select(where(is.numeric)) %>%
  #  select(ID, AUC_bc) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  uwot::umap(n_neighbors = 100     #n_neighbors = nearest neighbors = perplexity. Default = 15
  )

## Generate the UMAP plot

umap_df <- umap_fit %>%   #if using umap::umap: umap_fit$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID = row_number()) %>%
  inner_join(df_auc_meta, by = "ID")

umap_df %>%
#  filter(Analyte == "IL-4") %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = factor(batch),
             shape = group)) +
  geom_point(size = 3, position = position_jitter(h = 0.15, w = 0.15)) +
  labs(shape = "Group",
       color = "Batch",
       x = "UMAP1",
       y = "UMAP2") +
  ggtitle("UMAP Plot")+
  theme(plot.title = element_text(hjust = 0.5))

## UMAP of bulk MFI data (i.e. not area under the curve (AUC)

set.seed(date)

df_mfi_umap <- df_RSTR_comb %>% 
  drop_na() %>%
  #  filter(Analyte == "IL-6") %>%
  mutate(ID = row_number()) 

df_mfi_meta <- df_mfi_umap %>%
  select(ID, batch, group, PTID) #ADULT_RISK_TBINF, M0_KCVSEX, KCHCA_AGE_YR_CURRENT)

## Run UMAP

umap_fit <- df_mfi_umap %>%
  select(ID, mfi_final) %>%
  select(where(is.numeric)) %>%
  #  select(ID, AUC_bc) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  uwot::umap(n_neighbors = 15   #n_neighbors = nearest neighbors = perplexity. Default = 15
  )

## Generate the UMAP plot

umap_df <- umap_fit %>%   #if using umap::umap: umap_fit$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID = row_number()) %>%
  inner_join(df_mfi_meta, by = "ID")

umap_df %>%
  #  filter(Analyte == "IL-4") %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = factor(batch),
             shape = group)) +
  geom_point(size = 3, position = position_jitter(h = 0.15, w = 0.15)) +
  labs(shape = "Group",
       color = "Batch",
       x = "UMAP1",
       y = "UMAP2") +
  ggtitle("UMAP Plot")+
  theme(plot.title = element_text(hjust = 0.5))

