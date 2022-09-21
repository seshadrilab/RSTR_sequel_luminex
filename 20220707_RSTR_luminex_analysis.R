library(dplyr)
library(readr)
library(pheatmap)
library(tidyverse)
library(umap)
library(DescTools)
library(bayestestR)
library(ggpubr)
library(broom)
library(purrr)

date <- 20220707

# Read in all batches from a single directory and row bind them.

#setwd('G:/Shared drives/Seshadri Lab/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/out/luminex')
setwd('~/GoogleDrive/Lab Members/Kieswetter_Nathan/data/here/RSTR_Luminex/out/luminex')
df_RSTR_comb <- list.files() %>% 
  lapply(read_csv) %>%
  map(~ mutate(.x, across("PTID", as.character))) %>% 
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
                                           '284' = 'RS102284'
))

# Check that all patients have been included i.e. there should be 18 and all batches (6) - should return true

length(unique(df_RSTR_comb$PTID)) == 38
length(unique(df_RSTR_comb$batch)) == 6

# Add a better cytokine name

df_RSTR_comb <- df_RSTR_comb %>% 
  mutate(Analyte = recode(Analyte, 
                       'GM-CSF (44)' = 'GM-CSF',
                       'INF alpha (48)' = 'IFNa',
                       'INF gamma (43)' = 'IFNy',
                       'IL-1 alpha (62)' = 'IL-1a',
                       'IL-1 beta (18)' = 'IL-1b',
                       'IL-10 (28)' = 'IL-10',
                       'IL-12p70 (34)' = 'IL-12p70',
                       'IL-13 (35)' = 'IL-13',
                       'IL-15 (65)' = 'IL-15',
                       'IL-17A (36)' = 'IL-17A',
                       'IL-18 (66)' = 'IL-18',
                       'IL-1RA (38)' = 'IL-1RA',
                       'IL-2 (19)' = 'IL-2',
                       'IL-21 (72)' = 'IL-21',
                       'IL-22 (76)' = 'IL-22',
                       'IL-23 (63)' = 'IL-23',
                       'IL-27 (14)' = 'IL-27',
                       'IL-31 (37)' = 'IL-31',
                       'IL-4 (20)' = 'IL-4',
                       'IL-5 (21)' = 'IL-5',
                       'IL-6 (25)' = 'IL-6',
                       'IL-7 (26)' = 'IL-7',
                       'IL-8 (27)' = 'IL-8',
                       'IL-9 (52)' = 'IL-9',
                       'IP-10 (22)' = 'IP-10',
                       'TNF-alpha (45)' = 'TNFa',
                       'TNF-beta (54)' = 'TNFb',
  ))

#Check that all 27 analytes are present - should return true.

length(unique(df_RSTR_comb$Analyte)) == 27

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

# Export the counts

write.csv(df_RSTR_comb,"G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\luminex_all_batches.csv", row.names = FALSE)


# Read in the counts

df_RSTR_comb <- read.csv("G:\\Shared drives\\Seshadri Lab\\Lab Members\\Kieswetter_Nathan\\data\\here\\RSTR_Luminex\\out\\luminex\\luminex_all_batches.csv")

# Plotting the Data - Bar Graphs (including DMSO)
### NOTE: Whilst the below code works, MFI was uses as the dependent variable. 

analyte_list <- split(df_RSTR_comb, f = df_RSTR_comb$Analyte) #splits df into multiple df's based on analyte

tmp <- df_RSTR_comb %>%
  filter(stimulation != "DMSO")
analyte_list_noDMSO <-  split(tmp, f = tmp$Analyte) #splits df into multiple df's based on analyte (but excludes DMSO)

## (NOT IN USE) Plot all cytokines as a measure  background corrected concentration values over time (i.e. DMSO subtracted)

plotdata1 <- function(x) { 
  ggplot(analyte_list_noDMSO[[x]], 
         aes(x = factor(timepoint), 
             y = bkg_corr_conc,
             fill = group)) + 
    ggtitle(names(analyte_list_noDMSO)[[x]]) + 
    geom_boxplot() +
    #geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    facet_grid(. ~ stimulation) +
    scale_y_continuous(trans='log10') +
    xlab("Timepoint") + 
    ylab("Concentration (pg/mL)") +
    theme_bw() +  #Jolies theme and color pallet 
    theme(text = element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_fill_manual(values = c("#984EA3", "#4DAF4A")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols = c("****", "***", "**", "*", "ns"))
    )
}

lapply(seq_len(length(analyte_list_noDMSO)), plotdata1)

## (NOT IN USE) Plot the raw data over time (i.e. no background correction)

plotdata2 <- function(x) { 
  ggplot(analyte_list[[x]], 
         aes(x = factor(timepoint), 
             y = final_conc,
             fill = group)) + 
    ggtitle(names(analyte_list)[[x]]) + 
    geom_boxplot()+
    facet_grid(. ~ stimulation) +
    scale_y_continuous(trans='log10') +
    xlab("Timepoint") + 
    ylab("Concentration (pg/mL)") +
    theme_bw() +
    theme(text = element_text(family="Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_fill_manual(values = c("#984EA3", "#4DAF4A"))
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))
  )

}

lapply(seq_len(length(analyte_list)), plotdata2)

## (IN USE) Plot all cytokines as a measure  background corrected MFI values over time (i.e. DMSO subtracted)

plotdata3 <- function(x) { 
  ggplot(analyte_list_noDMSO[[x]], 
         aes(x = factor(timepoint), 
             y = mfi_bkgd_corr,
             color = group)) + 
    ggtitle(names(analyte_list_noDMSO)[[x]]) +
    labs(fill = 'Group') +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    facet_grid(. ~ stimulation) +
    scale_y_continuous(trans='log10') +
    xlab("Stimulation") + 
    ylab("Background Subtracted MFI") +
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
}

lapply(seq_len(length(analyte_list_noDMSO)), plotdata3)

## (IN USE) Plot all cytokines as a measure raw MFI values over time

plotdata4 <- function(x) { 
  ggplot(analyte_list[[x]], 
         aes(x = factor(timepoint), 
             y = final_MFI,
             color = group)) + 
    ggtitle(names(analyte_list)[[x]]) +
    labs(fill = 'Group') +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
   facet_grid(. ~ stimulation) +
    scale_y_continuous(trans='log10') +
    xlab("Stimulation") + 
    ylab("MFI") +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw() +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#984EA3", "#4DAF4A")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols = c("****", "***", "**", "*", "ns"))
    )
}

lapply(seq_len(length(analyte_list)), plotdata4)

# UMAP Code 
### Looks weird due to low dimensionality (I think)

## Clean the data for UMAP using the concentration data

set.seed(date)

df_RSTR_umap <- df_RSTR_comb %>% 
  drop_na() %>%
  mutate(ID = row_number()) 

df_RSTR_meta <- df_RSTR_umap %>%
  select(ID, Analyte, batch, group, PTID, timepoint)

## Run UMAP

umap_fit <- df_RSTR_umap %>%
  select(where(is.numeric)) %>%
  select(ID, final_conc) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  uwot::umap(
    n_neighbors = 15,
    min_dist = 0.01
  )

## Generate the UMAP plot

umap_df <- umap_fit%>%   #if using umap::umap: umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(df_RSTR_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = batch,
             shape = group)) +
  geom_point() +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") +
  theme(plot.title = element_text(hjust = 0.5))

