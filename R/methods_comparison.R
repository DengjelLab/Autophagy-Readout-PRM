library(MSstats)
library(tidyverse)
library(MSstatsLOBD)
library(patchwork)
library(data.table)
library(ComplexHeatmap)
library(factoextra)
library(RColorBrewer)
library(ggpubr)
library(progress)

pwd <- c('')
setwd(pwd)
source("R/fn.R")

################################################################
###               Part 0: General parameters                 ###
################################################################

analysis_name = c("methods_comparison")
set.seed(21)
PRM_file = c("analyses/experiment4_rerun/protein_data_imputated.csv")
DDA_file = c("data/method_comparison_dda.csv")
DIA_file = c("data/method_comparison_dia.csv")
DIAn_file = c("data/method_comparison_dia_normalized.csv")
ccurve_file = c("data/cc_peptide_quantification.csv")

################################################################
###                    Part 1: Data input                    ###
################################################################

dir.create(paste("analyses", analysis_name, sep = "/"))
prm = read_csv(PRM_file)
dda = read_csv(DDA_file)
dia = read_csv(DIA_file)
dian = read_csv(DIAn_file)
target_list <- read.csv(ccurve_file) %>% 
  as_tibble() %>% 
  filter(!Protein.Gene %in% c("", "TUBB4B", "KRT18")) %>% #Tubulin is a normalizer not a real target
  pull(Protein.Gene) %>% 
  unique()

################################################################
###                    Part 2: Shared hits                   ###
################################################################

dia_targets = dia %>% 
  filter(PG.Genes %in% target_list) %>%
  select(PG.Genes, "[1] 20231022_AL_HBSS_pooled_DIA1.raw.PG.Quantity", "[2] 20231022_AL_HBSS_pooled_DIA2.raw.PG.Quantity", "[3] 20231022_AL_HBSS_pooled_DIA3.raw.PG.Quantity", "[4] 20231022_AL_HBSSbaf_pooled_DIAbaf1.raw.PG.Quantity",
         "[5] 20231022_AL_HBSSbaf_pooled_DIAbaf2.raw.PG.Quantity", "[6] 20231022_AL_HBSSbaf_pooled_DIAbaf3.raw.PG.Quantity") %>% 
  drop_na() %>% 
  pull(PG.Genes) %>% 
  unique()

dian_targets = dian %>% 
  filter(PG.Genes %in% target_list) %>%
  select(PG.Genes, "[1] 20231022_AL_HBSS_pooled_DIA1.raw.PG.Quantity", "[2] 20231022_AL_HBSS_pooled_DIA2.raw.PG.Quantity", "[3] 20231022_AL_HBSS_pooled_DIA3.raw.PG.Quantity", 
         "[4] 20231022_AL_HBSSbaf_pooled_DIAbaf1.raw.PG.Quantity", "[5] 20231022_AL_HBSSbaf_pooled_DIAbaf2.raw.PG.Quantity", "[6] 20231022_AL_HBSSbaf_pooled_DIAbaf3.raw.PG.Quantity") %>% 
  drop_na() %>% 
  pull(PG.Genes) %>% 
  unique()

ddam_targets = dda %>% 
  filter(Gene %in% target_list) %>% 
  select(Gene, `HBSS_1 MaxLFQ Intensity`, `HBSS_2 MaxLFQ Intensity`, `HBSS_3 MaxLFQ Intensity`,
         `HBSS_BafA1_1 MaxLFQ Intensity`, `HBSS_BafA1_2 MaxLFQ Intensity`, `HBSS_BafA1_3 MaxLFQ Intensity`) %>% 
  filter(if_any(c(`HBSS_1 MaxLFQ Intensity`, `HBSS_2 MaxLFQ Intensity`, `HBSS_3 MaxLFQ Intensity`,
                  `HBSS_BafA1_1 MaxLFQ Intensity`, `HBSS_BafA1_2 MaxLFQ Intensity`, `HBSS_BafA1_3 MaxLFQ Intensity`), function(x){x>0})) %>%
  pull(Gene) %>% 
  unique()

dda_targets = dda %>% 
  filter(Gene %in% target_list) %>% 
  select(Gene, `HBSS_1 Intensity`, `HBSS_2 Intensity`, `HBSS_3 Intensity`,
         `HBSS_BafA1_1 Intensity`, `HBSS_BafA1_2 Intensity`, `HBSS_BafA1_3 Intensity`) %>% 
  filter(if_any(c(`HBSS_1 Intensity`, `HBSS_2 Intensity`, `HBSS_3 Intensity`,
                  `HBSS_BafA1_1 Intensity`, `HBSS_BafA1_2 Intensity`, `HBSS_BafA1_3 Intensity`), function(x){x>0})) %>%
  pull(Gene) %>% 
  unique()

prm_targets = prm %>% 
  filter(Protein.Gene %in% target_list,
         ! is.na(Protein.Intensity)) %>%
  pull(Protein.Gene) %>% 
  unique()

length(prm_targets)
length(dian_targets)
length(dia_targets)
length(dda_targets)
length(ddam_targets)

dian_targets_cv = dian %>% 
  filter(PG.Genes %in% target_list) %>% 
  select(PG.Genes, "[1] 20231022_AL_HBSS_pooled_DIA1.raw.PG.Quantity", "[2] 20231022_AL_HBSS_pooled_DIA2.raw.PG.Quantity", "[3] 20231022_AL_HBSS_pooled_DIA3.raw.PG.Quantity", 
         "[4] 20231022_AL_HBSSbaf_pooled_DIAbaf1.raw.PG.Quantity", "[5] 20231022_AL_HBSSbaf_pooled_DIAbaf2.raw.PG.Quantity", "[6] 20231022_AL_HBSSbaf_pooled_DIAbaf3.raw.PG.Quantity") %>% 
  pivot_longer(!PG.Genes) %>% 
  mutate(Condition = ifelse(grepl("HBSSbaf", name), "HBSS + BafA1", "HBSS")) %>%
  group_by(PG.Genes, Condition) %>% 
  summarise(cv_cond = sd(value)/mean(value)) %>% 
  ungroup() %>%
  group_by(PG.Genes) %>% 
  summarise(cv = mean(cv_cond)) %>% 
  mutate(method = "DIA - Spectronaut \n Normalized") %>% 
  rename(Protein = PG.Genes)

dia_targets_cv = dia %>% 
  filter(PG.Genes %in% target_list) %>% 
  select(PG.Genes, "[1] 20231022_AL_HBSS_pooled_DIA1.raw.PG.Quantity", "[2] 20231022_AL_HBSS_pooled_DIA2.raw.PG.Quantity", "[3] 20231022_AL_HBSS_pooled_DIA3.raw.PG.Quantity", 
         "[4] 20231022_AL_HBSSbaf_pooled_DIAbaf1.raw.PG.Quantity", "[5] 20231022_AL_HBSSbaf_pooled_DIAbaf2.raw.PG.Quantity", "[6] 20231022_AL_HBSSbaf_pooled_DIAbaf3.raw.PG.Quantity") %>% 
  pivot_longer(!PG.Genes) %>% 
  mutate(Condition = ifelse(grepl("HBSSbaf", name), "HBSS + BafA1", "HBSS")) %>%
  group_by(PG.Genes, Condition) %>% 
  summarise(cv_cond = sd(value)/mean(value)) %>% 
  ungroup() %>%
  group_by(PG.Genes) %>% 
  summarise(cv = mean(cv_cond)) %>% 
  mutate(method = "DIA - Spectronaut") %>% 
  rename(Protein = PG.Genes)

ddam_targets_cv = dda %>% 
  filter(Gene %in% target_list) %>% 
  select(Gene, `HBSS_1 MaxLFQ Intensity`, `HBSS_2 MaxLFQ Intensity`, `HBSS_3 MaxLFQ Intensity`,
         `HBSS_BafA1_1 MaxLFQ Intensity`, `HBSS_BafA1_2 MaxLFQ Intensity`, `HBSS_BafA1_3 MaxLFQ Intensity`) %>% 
  pivot_longer(!Gene) %>%
  mutate(Condition = ifelse(grepl("BafA1", name), "HBSS + BafA1", "HBSS")) %>% 
  group_by(Gene, Condition) %>% 
  summarise(cv_cond = sd(value)/mean(value)) %>% 
  ungroup() %>%
  group_by(Gene) %>% 
  summarise(cv = mean(cv_cond)) %>% 
  mutate(method = "DDA \n MaxLFQ") %>% 
  rename(Protein = Gene)

dda_targets_cv = dda %>% 
  filter(Gene %in% target_list) %>% 
  select(Gene, `HBSS_1 Intensity`, `HBSS_2 Intensity`, `HBSS_3 Intensity`,
         `HBSS_BafA1_1 Intensity`, `HBSS_BafA1_2 Intensity`, `HBSS_BafA1_3 Intensity`) %>% 
  pivot_longer(!Gene) %>%
  mutate(Condition = ifelse(grepl("BafA1", name), "HBSS + BafA1", "HBSS")) %>% 
  group_by(Gene, Condition) %>% 
  summarise(cv_cond = sd(value)/mean(value)) %>% 
  ungroup() %>%
  group_by(Gene) %>% 
  summarise(cv = mean(cv_cond)) %>% 
  mutate(method = "DDA") %>% 
  rename(Protein = Gene)

prm_targets_cv = prm %>% 
  filter(Protein.Gene %in% target_list,
         ! is.na(Protein.Intensity)) %>% 
  group_by(Protein.Gene, Condition) %>% 
  summarise(cv_cond = sd(Protein.Intensity)/mean(Protein.Intensity)) %>% 
  ungroup() %>% 
  group_by(Protein.Gene) %>% 
  summarise(cv = mean(cv_cond)) %>% 
  mutate(method = "PRM") %>% 
  rename(Protein = Protein.Gene)

print(mean(dian_targets_cv$cv))
print(mean(dia_targets_cv$cv))
print(mean(ddam_targets_cv$cv, na.rm = TRUE))
print(mean(dda_targets_cv$cv, na.rm = TRUE))
print(mean(prm_targets_cv$cv))

cvs = rbind(dia_targets_cv, dian_targets_cv, dda_targets_cv ,ddam_targets_cv, prm_targets_cv)
write.csv(cvs ,paste("analyses", analysis_name, "CV_methods_comparisons.csv", sep = "/"))

