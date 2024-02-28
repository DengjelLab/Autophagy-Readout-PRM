library(MSstats)
library(tidyverse)
library(MSstatsLOBD)
library(data.table)
library(progress)

pwd <- c()
setwd(pwd)
source("R/fn.R")

################################################################
###               Part 0: General parameters                 ###
################################################################

analysis_name = c("experiment2")
set.seed(21)
ccurve_file = c("data/cc_peptide_quantification.csv")
experiment_file = c("data/fragment_quantification_exp12.csv")

#Ploting calibration curve results
plot_cc = FALSE

#Output files names:
pkg_linear_fits_file = c("MSStatsLOBD_linear_fits")
pkg_nonlinear_fits_file = c("MSStatsLOBD_nonlinear_fits")
FOM_file = c("nonlinear_fits_FOM")
exp_ccurve_match_file = c("experiment_calibrationcurve_matched")
peptide_data_file = c("peptide_data")
protein_data_file = c("protein_data")

################################################################
###           Part 1: Input data & filtering ions            ###
################################################################

#Cleanup files:
ccurve = read.csv(ccurve_file) %>% 
  as_tibble() %>% 
  filter(Peptide.Note %in% c("wrongblank", "")) %>% #the wrongblank label was used for peptides where the blank was incorrectly measured. Removed from analysis
  filter(Sample.Type != "Unknown")  %>% #Sample type was set to unknown for washes and RT test samples
  mutate(Total.Area = as.numeric(Total.Area)) %>% 
  filter(!is.na(Total.Area)) %>% 
  mutate(Fragment.Ion = paste0(Fragment.Ion, " ", Product.Charge, "+"))#Fragment Ion does not account for charge. Add it if 2 charge states were monitored.

runs_to_remove = c("PA15_D3_4","PA15_D3_5", "PA15_D4_4","PA15_D4_5", "PA15_D5_4", "PA15_D5_5", #These were specific to experiment 1
                   "PA15_DB1_2", "PA15_DB1_3", "PA15_DB2_2", "PA15_DB2_3", "PA15_DB3_2", "PA15_DB3_3")

experiment = read.csv(experiment_file) %>% 
  filter(Peptide.Note == "", Condition != "") %>% 
  as_tibble() %>% 
  mutate(Fragment.Ion = paste0(Fragment.Ion, " ", Product.Charge, "+")) %>% 
  filter(!Replicate.Name %in% runs_to_remove)

peptides <- unique(experiment$Peptide)

#Keeping only the ions present in the experiment in the calibration curve data
ccurve_filtered <- FilterIons(experiment, ccurve)

##Merge the fragments together to get peptide intensity and add controls back
summed_areas = ccurve_filtered %>%
  mutate(Area = as.numeric(Area)) %>% 
  group_by(Replicate.Name, Peptide) %>%
  summarise(Summed.Area = sum(Area)) %>% 
  ungroup()

standards_areas = ccurve[ccurve$Standard.Type == "Global Standard",] %>% 
  mutate(Area = as.numeric(Area)) %>% 
  group_by(Replicate.Name, Peptide) %>% 
  summarise(Summed.Area = sum(Area)) %>% 
  ungroup()

summed_areas = rbind(summed_areas, standards_areas)
ccurve = select(ccurve, Peptide, Replicate.Name, Analyte.Concentration, Standard.Type) %>% 
  distinct()
ccurve = left_join(summed_areas, ccurve, by = c("Replicate.Name", "Peptide"))

################################################################
###     Part 2: Normalize the signal of each ccurve run      ###
################################################################

std = ccurve %>% 
  filter(Standard.Type == "Global Standard", ! is.na(Summed.Area)) %>% 
  group_by(Replicate.Name) %>% 
  summarise(StdReplicate = sum(Summed.Area))

std$avg_standard = mean(std$StdReplicate)

ccurve = left_join(ccurve, std, by = "Replicate.Name")
ccurve = mutate(ccurve, NormArea = Summed.Area/(StdReplicate/avg_standard))

################################################################
###      Part 3: Run MSStats LOBD to get FOM and plots       ###
################################################################
#Preparing input for MSStats
df_out <- ccurve %>%
  filter(Standard.Type == "") %>% 
  select(NormArea, Analyte.Concentration, Peptide, Replicate.Name) %>% 
  rename(INTENSITY = NormArea,
         CONCENTRATION = Analyte.Concentration,
         NAME = Peptide,
         REPLICATE = Replicate.Name)

dir.create(paste("analyses", analysis_name, sep = "/"))
if(plot_cc){
  dir.create(paste("analyses", analysis_name ,"plots_cc", sep = "/"))
}

#Looping over all the peptides to fit regressions and generate plots
df_out = df_out %>% 
  mutate(INTENSITY = sapply(INTENSITY, FUN = ReplaceByNonNull))

#Function fitting both linear and non-linear fits to the different peptides
fits = FitAndPlotCalibrationCurves(df_out, plot_results = plot_cc, adress = paste("analyses", analysis_name, "plots_cc", sep = "/"))

#Function recalculating the slope for nonlinear fits. Linear regression of points above the noise.
peptides <- unique(df_out$NAME)
cor_fits = CorrectingSlope(fits, plot_results = plot_cc, adress = paste("analyses", analysis_name, "plots_cc", sep = "/"))

#Exporting fits results
linear_lod <- bind_rows(fits[[1]])
nonlinear_lod <- bind_rows(cor_fits)
fom <- WriteFom(linear_lod, nonlinear_lod)
fom_out <- fom %>%
  select(NAME, NONLIN_LOB, NONLIN_LOD, NONLIN_SLOPE, NONLIN_INTERCEPT) %>% 
  rename("Peptide" = "NAME", "LOD" = "NONLIN_LOB", "LOQ" = "NONLIN_LOD", "slope" = "NONLIN_SLOPE", "intercept" = "NONLIN_INTERCEPT")
write.csv(fom_out, file = paste("analyses", analysis_name, paste0(FOM_file, ".csv"), sep = "/"))
write.csv(linear_lod, file = paste("analyses", analysis_name, paste0(pkg_linear_fits_file, ".csv"), sep = "/"))
write.csv(nonlinear_lod, file = paste("analyses", analysis_name, paste0(pkg_nonlinear_fits_file, ".csv"), sep = "/"))

################################################################
###               Part 4: Overall cal. curve plots           ###
################################################################

if(plot_cc){
  #fom <- read.csv(paste("analyses", analysis_name,paste0(FOM_file, ".csv"), sep = "/"))
  linear_lod <- read.csv(paste("analyses", analysis_name, paste0(pkg_linear_fits_file, ".csv"), sep = "/"))
  nonlinear_lod <- read.csv(paste("analyses", analysis_name, paste0(pkg_nonlinear_fits_file, ".csv"), sep = "/"))
  fom <- WriteFom(linear_lod, nonlinear_lod)
  fits_comparison <- CompareFits(fom)
  ggsave(paste("analyses", analysis_name, "plots_cc", "fits_comparison.svg", sep = "/"),
         plot = fits_comparison,
         width = 28,
         height = 16,
         units = "cm")
  overview <- PlotFomOverview(fom)
  ggsave(paste("analyses", analysis_name, "plots_cc", "FOM_overview.svg", sep = "/"),
         plot = overview,
         width = 28,
         height = 16,
         units = "cm")
}

################################################################
###            Part 5: Experiment matching to curve          ###
################################################################
nonlin = read.csv(paste("analyses", analysis_name, paste0(FOM_file, ".csv"), sep = "/"))

experiment = experiment %>%
  mutate(Area = as.numeric(Area)) %>% 
  select(Protein.Gene, Peptide, Precursor, Isotope.Label.Type,Total.Area,
         Normalized.Area, Replicate.Name, Condition, BiolReplicates,
         Fragment.Ion, Area)

summed_fragments = experiment %>% 
  group_by(Peptide, Isotope.Label.Type, Replicate.Name) %>%
  summarise(Summed.Areas = sum(Area)) %>% 
  ungroup()

experiment = left_join(experiment, summed_fragments, by = c("Peptide", "Isotope.Label.Type", "Replicate.Name"))

normalized_areas = experiment %>%
  select(Peptide, Replicate.Name, Summed.Areas, Isotope.Label.Type) %>%
  distinct() %>% 
  group_by(Peptide, Replicate.Name) %>% 
  summarise(Normalized.Signal = Summed.Areas[Isotope.Label.Type == "light"]/Summed.Areas[Isotope.Label.Type == "heavy"]) %>% 
  ungroup()

df = left_join(experiment, normalized_areas, by = c("Peptide", "Replicate.Name")) %>% 
  left_join(nonlin, by = "Peptide") %>%
  mutate(Calculted.Concentration = (Summed.Areas - intercept)/slope) %>% 
  mutate(Above.LOD = Calculted.Concentration >= LOD, Above.LOQ = Calculted.Concentration >= LOQ)

peptide_fom <- PlotValidity(df)
ggsave(file = paste("analyses", analysis_name, "PeptidesFoM.svg",sep = "/"),
       plot = peptide_fom,
       width = 18,
       height = 12,
       units = "cm")

write.csv(df, file = paste("analyses", analysis_name ,paste0(exp_ccurve_match_file, ".csv"), sep = "/"))

################################################################
###         Part 6: Peptides & Protein quantification        ###
################################################################
df <- read.csv(paste("analyses", analysis_name, paste0(exp_ccurve_match_file, ".csv") , sep = "/"))

#Make a peptide-level dataframe but need to calculate protein amount first since this is based on fragment data
df = df %>%
  select(Protein.Gene, Peptide, Isotope.Label.Type, Normalized.Signal, Replicate.Name, Condition, LOQ, Above.LOQ, Fragment.Ion, Area, BiolReplicates)

df[df$Above.LOQ == FALSE, ]$Normalized.Signal = NA


protein_abundances = CalculateProteinAbundanceLOQimpute(df, Protein.Gene, Peptide, Replicate.Name, Isotope.Label.Type, Area, Above.LOQ)
df = left_join(df, protein_abundances, by = c("Protein.Gene", "Replicate.Name"))
prot_df = df %>%  
  select(Protein.Gene, Replicate.Name, Protein.Intensity,Condition, Imputed_data, BiolReplicates) %>% 
  distinct()
write.csv(prot_df, paste("analyses", analysis_name, paste0(protein_data_file, "_imputated.csv"), sep = "/"))

val_pep = df %>%
  group_by(Peptide, Replicate.Name) %>% 
  summarise(valid_pep = ifelse(mean(Above.LOQ) == 1, TRUE, FALSE)) %>% 
  ungroup()
df = left_join(df, val_pep, by = c("Peptide", "Replicate.Name"))
df[df$valid_pep == FALSE,]$Normalized.Signal = NA
pep_df = df %>% 
  filter(Isotope.Label.Type == "light") %>% 
  mutate(Peptide.Label = paste(Protein.Gene, Peptide)) %>% 
  select(Protein.Gene, Peptide.Label, Peptide, Normalized.Signal, Replicate.Name, Condition, valid_pep, BiolReplicates) %>% 
  distinct()
write.csv(pep_df, paste("analyses", analysis_name, paste0(peptide_data_file,"LOQ_mod_imputation", ".csv"), sep = "/"))
