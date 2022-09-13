#Project: RSTR the sequel - Luminex Data
#Nathan Kieswetter
# RSTR_TGFb1 Luminex analysis (AUC)
# Aim: Assess cytokine secretion in the face of Mtbl, and PP! stimulation from the RSTR cohort

# Load in the required libraries
library(tidyverse)
library(umap)
library(DescTools)
library(bayestestR)
library(ggpubr)
library(broom)
library(ComplexHeatmap)

date <- 20220809

# Aim: Calculate the AUC for the MFI values

#Read in the data

path <- 'G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/out/luminex/luminex_all_batches.csv'
  
df_RSTR_comb <- read_csv(path)

## Make multiple data frames per patient, per analyte to use for AUC calculations.

auc_split <- split(df_RSTR_comb, f = list(df_RSTR_comb$PTID, 
                                          df_RSTR_comb$Analyte, 
                                          df_RSTR_comb$stimulation,
                                          df_RSTR_comb$batch))

# Function to calculate the AUC per the factors defines above.

## Define the AUC (MFI) for the non-background corrected secreted data.

auc_fun1 <- function(x) { # Function to define the AUC 
  timept <- x$timepoint
  mfi <- x$final_MFI
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

df_auc_meta <- str_split_fixed(df_auc1$rowname, "\\.", 4) %>%
  as.data.frame() %>%
  setNames(c('PTID', 'Analyte', 'stimulation', 'batch'))

df_auc <- bind_cols(df_auc1, df_auc_meta)

df_auc <- df_auc %>%
  select(-rowname)

# The following function will calculate the MFI for the background corrected AUC of the MFI over time

auc_fun2 <- function(x) {
  timept <- x$timepoint
  mfi <- x$mfi_bkgd_corr
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

df_auc <- bind_cols(df_auc, df_auc2) # merge the two df's together for downstream analysis/graphing

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

write.csv(df_auc,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\luminex_all_batches_auc_mfi.csv", row.names = FALSE)

# Plot the AUC values for the PTIDS by stimulation

df_auc <- read.csv("G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\luminex_all_batches_auc_mfi.csv")

analyte_auc_list <- split(df_auc, f = df_auc$Analyte)

tmp <- df_auc %>%
  filter(stimulation != "DMSO")

analyte_auc_list_no_DMSO <- split(tmp, f = df_auc$Analyte)

## Plot the non-background controlled AUC data

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

plotdata3 <- function(x) { 
  ggplot(analyte_auc_list[[x]], 
         aes(x = factor(stimulation), 
             y = AUC_nbc,
             color = group)) + 
    ggtitle(names(analyte_auc_list)[[x]]) +
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
}

lapply(seq_len(length(analyte_auc_list)), plotdata3)

## Plot the background corrected AUC based on the MFI data

plotdata4 <- function(x) { 
  ggplot(analyte_auc_list_no_DMSO[[x]], 
         aes(x = factor(stimulation), 
             y = AUC_bc,
             color = group)) + 
    ggtitle(names(analyte_auc_list_no_DMSO)[[x]]) +
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
}

lapply(seq_len(length(analyte_auc_list_no_DMSO)), plotdata4)

# Plot all cytokines - non background corrected

tmp_pp1 <- df_auc %>%
  filter(stimulation != "DMSO") %>%
  filter(stimulation == "PP1")

tmp_mtbl <- df_auc %>%
  filter(stimulation != "DMSO") %>%
  filter(stimulation == "MTBL")

plot_pp1 <- ggplot(tmp_pp1, 
         aes(x = factor(Analyte), 
             y = AUC_nbc,
             color = group)) + 
    labs(color = 'Group') +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    #    facet_grid(. ~ group) +
    ggtitle("PP1") +
    scale_y_continuous(trans='log10') +
    xlab("Stimulation") + 
    ylab("AUC (MFI)") +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw() +
    theme(text = element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif", 
#                      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
#                                        symbols = c("****", "***", "**", "*", "ns"))
                       )

plot_mtbl <- ggplot(tmp_mtbl, 
                    aes(x = factor(Analyte), 
                        y = AUC_nbc,
                        color = group)) + 
  labs(color = 'Group') +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  #    facet_grid(. ~ group) +
  ggtitle("Mtb Lysate") +
  scale_y_continuous(trans='log10') +
  xlab("Stimulation") + 
  ylab("AUC (MFI)") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw() +
  theme(text = element_text(family="Arial"),
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
  stat_compare_means(method = "wilcox.test", 
#                     label = "p.signif", 
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
#                                        symbols = c("****", "***", "**", "*", "ns"))
                     )
plot_pp1
plot_mtbl

# UMAP

## Clean the data for UMAP using the concentration data

set.seed(date)

df_auc_umap <- df_auc %>% 
  drop_na() %>%
#  filter(Analyte == "IL-6") %>%
  mutate(ID = row_number()) 

df_auc_meta <- df_auc_umap %>%
  select(ID, Analyte, batch, group, PTID) #ADULT_RISK_TBINF, M0_KCVSEX, KCHCA_AGE_YR_CURRENT)

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
  select(ID, Analyte, batch, group, PTID) #ADULT_RISK_TBINF, M0_KCVSEX, KCHCA_AGE_YR_CURRENT)

## Run UMAP

umap_fit <- df_mfi_umap %>%
  select(ID, final_MFI, FI...Bkgd, FI) %>%
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



## Format the data for assessment of background corrected AUC values

df_auc_mdmso <- df_auc %>%
  filter(stimulation != "DMSO")

df_auc_mdmso$group <- as.factor(df_auc_mdmso$group)

df_auc_mdmso$stimulation <- as.factor(df_auc_mdmso$stimulation)

## Create a list of DF contain per analyte, per stimulation

##Background controlled

stim_auc_list <- df_auc_mdmso %>%
  split(f = list(df_auc_mdmso$stimulation,
                 df_auc_mdmso$Analyte)
  )

## Non-background controlled

### Format the data for analysis

df_auc$group <- as.factor(df_auc$group)

df_auc$stimulation <- as.factor(df_auc$stimulation)

stim_auc_list_nbc <- df_auc %>%
  split(f = list(df_auc$stimulation,
                 df_auc$Analyte)
)

# Wilcoxon test

## Background controlled data

wilcox_results <- lapply(stim_auc_list, function(x) with(x, wilcox.test(AUC_bc ~ group)))

names_list <- names(wilcox_results)

df_wilcox_results <- wilcox_results %>% 
  map_df(broom::tidy) %>%
  cbind(names_list) %>%
  separate(names_list, c('stimulation', 'analyte'),"[.]") %>%
  select(stimulation, analyte, statistic, p.value, method, alternative)

wilcox_sig_ttest <- df_wilcox_results %>% # Returns only significant observations
  filter(p.value <= 0.05)

## Non-background controlled data

wilcox_results_nbc <- lapply(stim_auc_list_nbc, function(x) with(x, wilcox.test(AUC_nbc ~ group)))

names_list <- names(wilcox_results_nbc)

df_wilcox_results_nbc <- wilcox_results_nbc %>%
  map_df(broom::tidy) %>%
  cbind(names_list) %>%
  separate(names_list, c('stimulation', 'analyte'),"[.]") %>%
  select(stimulation, analyte, statistic, p.value, method, alternative)

wilcox_sig_nbs <- df_wilcox_results_nbc %>%
  filter(p.value <= 0.05)
