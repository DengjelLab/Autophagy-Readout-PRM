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

analysis_name = c("cc_analysis")
set.seed(21)
ccurve_file = c("data/cc_peptide_quantification.csv")

#Ploting calibration curve results
plot_cc = TRUE

#Output files names:
pkg_linear_fits_file = c("MSStatsLOBD_linear_fits")
pkg_nonlinear_fits_file = c("MSStatsLOBD_nonlinear_fits")
FOM_file = c("nonlinear_fits_FOM")
exp_ccurve_match_file = c("experiment_calibrationcurve_matched")
peptide_data_file = c("peptide_data")
protein_data_file = c("protein_data")

################################################################
###                    Part 1: Input data                    ###
################################################################

#Cleanup files:
ccurve = read.csv(ccurve_file) %>% 
  as_tibble() %>% 
  filter(Peptide.Note %in% c("wrongblank", "")) %>% #the wrongblank label was used for peptides where the blank was incorrrectly measured. Removed from analysis
  filter(Sample.Type != "Unknown")  %>% #Sample type was set to unknown for washes and RT test samples
  mutate(Total.Area = as.numeric(Total.Area)) %>% 
  filter(!is.na(Total.Area)) %>% 
  mutate(Fragment.Ion = paste0(Fragment.Ion, " ", Product.Charge, "+"))#Fragment Ion does not mention charge. Add it if 2 charge state were monitored.

##Merge the fragments together to get peptide intensity and add controls back
summed_areas = ccurve %>%
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


#Looping over all the peptides to fit the regressions and generate the plots

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
protein_annotations <- ccurve %>% 
  select(Protein, Protein.Gene, Peptide) %>% 
  filter(! Protein.Gene %in% c("KRT18", "TUBB4B")) %>% 
  distinct() %>% 
  mutate(Protein = str_extract(Protein, "\\|(.*?)\\|"),
         Protein = str_replace(Protein, "\\|", ""),
         Protein = str_replace(Protein, "\\|", ""))
fom_out <- fom %>% 
  left_join(protein_annotations, by = "Peptide")

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
               width = 10,
               height = 10,
               units = "cm")
}

#Figure for calibration curve example:
# library(patchwork)
# i <- df_out[df_out$NAME == "IQLPSEK",]
# f <- fits[[2]][[52]]
# t <- plot_quantlim(i,f)
# t[[1]] <- t[[1]] +
#   theme(legend.position = "none",
#         plot.title = element_text(size = 10, hjust = 0.5),
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8),
#         axis.text.x = element_text(size = 6),
#         axis.text.y = element_text(size = 6))
# t[[2]] <- t[[2]] +
#   ggtitle("") +
#   theme(legend.position = "none",
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8),
#         axis.text.x = element_text(size = 6),
#         axis.text.y = element_text(size = 6),
#         plot.margin = margin(t = -5, r = 0.1, b = 0, l = 0),
#         plot.background = element_rect(fill = "transparent", color = NA))
#   
# 
# o <- t[[1]] + inset_element(t[[2]], left = 0.005, bottom = 0.45, right = 0.5, top = 1)
# ggsave("figures/gbrl_curve.svg",
#        plot = o,
#        width = 105,
#        height = 297*1/3,
#        unit = "mm",
#        dpi = 300)