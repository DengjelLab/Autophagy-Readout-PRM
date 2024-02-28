################################################################
###               Functions for PRM Analysis                 ###
################################################################

#Make a heatmap with ComplexHeatmap from a df containing the following columns:
# - A val column which contains the numeric values to be displayed as a color in the heatmap
# - A col column which contains the samples
# - A row column which contains the features
MakeHeatmap = function(df, row, col, val, ht_title, ht_rowtitle, ht_columntitle, ht_annotations){
  matrix = df %>% 
    select({{col}}, {{row}}, {{val}}) %>% 
    spread(key = {{col}}, value = {{val}})
  rownames(matrix) = matrix[,1]
  matrix = matrix[,-1]
  matrix = na.omit(matrix)
  matrix = t(scale(t(as.matrix(matrix))))
  ht = Heatmap(matrix,
               name = ht_title,
               width = ncol(matrix)*unit(3, "mm"),
               height = nrow(matrix)*unit(3, "mm"),
               border_gp = gpar(col = "white", lwd = 1),
               rect_gp = gpar(col = "white", lwd = 1),
               row_title = ht_rowtitle,
               column_title = ht_columntitle,
               row_names_gp = gpar(fontsize = 5.5),
               column_names_gp = gpar(fontsize = 5.5),
               top_annotation = ht_annotations)
  return(ht)
}

#Make a PCA plot using with FactoExtra from a df containing the following columns:
# - A val column which contains the numeric values to be displayed as a color in the heatmap
# - A col column which contains the samples
# - A row column which contains the features
PlotPCASamples = function(df, row, col, val, groups, plot_title){
  matrix = df %>% 
    select({{col}}, {{row}}, {{val}}) %>% 
    spread(key = {{col}}, value = {{val}})
  rownames(matrix) = matrix[,1]
  matrix = matrix[,-1]
  matrix = na.omit(matrix)
  matrix = t(scale(t(as.matrix(matrix))))
  pca = prcomp(t(matrix), scale = FALSE)
  pca_plot = fviz_pca_ind(pca, col.ind = groups, repel = TRUE, title = plot_title)
  return(pca_plot)
}


#A function calculating the protein abundance the same way Skyline does, i.e. mean(light ions)/mean(heavy ions) for
#all light ions with heavy counterpart.
#From following input data:
#df = dataframe with peptide and ions data
#protein_col = df column identifying the proteins
#replicate_col = df column identifying the replicates
#ion_id = df column identifying ions. CAUTION: For skylie documents use name + charge as just fragment names can be repeated
#ion_quant = df column with quantitative value associated to ion
#isotope_col = df column specifying isotope label type
CalculateProteinAbundance = function(df, protein_col, replicate_col, isotope_col,ion_quant_col){
  prot = df %>% 
    pull({{protein_col}}) %>% 
    unique()
  rep = df %>% 
    pull({{replicate_col}}) %>% 
    unique()
  prot_v = c()
  replicate_v = c()
  int_v = c()
  
  pBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(prot),
                           complete = "=",
                           incomplete = "-",
                           current = ">",
                           clear = FALSE,
                           width = 100)
  
  for(p in 1:length(prot)){
    pBar$tick()
    for(r in 1:length(rep)){
      s = df %>% 
        filter({{protein_col}} == prot[p] & {{replicate_col}} == rep[r]) %>% 
        pivot_wider(names_from = {{isotope_col}}, values_from = {{ion_quant_col}}) %>% 
        filter(heavy != 0)
      light = mean(s$light, na.rm = TRUE)
      heavy = mean(s$heavy)
      norm_prot = light/heavy
      prot_v = c(prot_v, prot[p])
      replicate_v = c(replicate_v, rep[r])
      int_v = c(int_v, norm_prot)
    }
  }
  protein_abundances = data.frame(prot_v, replicate_v, int_v) %>% 
    rename(Protein.Gene = prot_v, Replicate.Name = replicate_v, Protein.Intensity = int_v)
  return(protein_abundances)
}



#A function to fit both linear and non-linear regressions to the calibration curve results and plot the results
#Takes a MSStat formatted df as input
FitAndPlotCalibrationCurves = function(df, plot_results = TRUE, adress = ""){
  peptides = unique(df$NAME)
  lfs = list()
  nlfs = list()
  lps = list()
  nlps = list()

  pBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(peptides),
                           complete = "=",
                           incomplete = "-",
                           current = ">",
                           clear = FALSE,
                           width = 100)
  
  for(i in 1:length(peptides)){
    pBar$tick()
    peptide_name = peptides[i]
    #print(paste(i, "/", length(peptides), ": ", peptide_name))
    peptide_data = df %>% filter(NAME == peptide_name)
    lf = linear_quantlim(peptide_data)
    lfs[[i]] = lf
    nlf = nonlinear_quantlim(peptide_data)
    nlfs[[i]] = nlf
    lp = plot_quantlim(spikeindata = peptide_data, quantlim_out = lf, address = FALSE)
    nlp = plot_quantlim(spikeindata = peptide_data, quantlim_out = nlf, address = FALSE)
    #print("Data fitted")
    if(plot_results == TRUE){
      for(j in 1:length(lp)){
        lp[[j]] = lp[[j]] + theme(legend.position = "none")
      }
      lp[[1]] = lp[[1]] + theme(text = element_text(size = 5))
      lp[[2]] = lp[[2]] +
        theme(text = element_text(size = 3.5)) +
        ggtitle("")
      for(j in 1:length(nlp)){
        nlp[[j]] = nlp[[j]] + theme(legend.position = "none")
      }
      nlp[[1]] = nlp[[1]] + theme(text = element_text(size = 5))
      nlp[[2]] = nlp[[2]] +
        theme(text = element_text(size = 3.5)) +
        ggtitle("")
      olp = lp[[1]] + inset_element(lp[[2]], 0.05, 0.5 , 0.65, 1)
      onlp = nlp[[1]] + inset_element(nlp[[2]], 0.05, 0.5 , 0.65, 1)
      op = olp + onlp
      op = op + plot_annotation(title = peptide_name, subtitle = "Linear and Non-Linear fits")
      lps[[i]] = olp
      nlps[[i]] = onlp
      #print("Plots drawn")
      ggsave(file = paste(adress,paste0("Fits_", peptide_name,".svg"), sep = "/"),
             plot = op,
             width = 46,
             height = 30,
             units = "cm")
      #print("Plots exported")
    }
  }
  return(list(lfs, nlfs, lps, nlps))
}


#A functions that corrects the fits to get the right slope
CorrectingSlope = function(ls, plot_results = TRUE, adress = ""){
  nlfs = ls[[2]]
  nlps = ls[[4]]
  pBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(nlfs),
                           complete = "=",
                           incomplete = "-",
                           current = ">",
                           clear = FALSE,
                           width = 100)
  for(i in 1:length(nlfs)){
    pBar$tick()
    #Correction
    peptide_name = nlfs[[i]]$NAME[1]
    lod = max(nlfs[[i]]$LOB)
    linear_part = nlfs[[i]][nlfs[[i]]$CONCENTRATION > lod, ]
    lf = lm(linear_part$MEAN ~ linear_part$CONCENTRATION)
    nlfs[[i]]$intercept = lf$coefficients[1]
    nlfs[[i]]$slope = lf$coefficients[2]
    nlfs[[i]]$linreg = nlfs[[i]]$slope * nlfs[[i]]$CONCENTRATION + nlfs[[i]]$intercept
    #Plots
    if (plot_results == TRUE){
      mean3lod = nlfs[[i]][which.min(abs(nlfs[[i]]$CONCENTRATION - 3 * lod)), ]$MEAN
      fl = ggplot(data = nlfs[[i]], aes(x = CONCENTRATION)) +
        geom_point(aes(y = MEAN, col = "package_fit"), alpha = 0.5) +
        geom_line(aes(y = linreg, col = "correction"),linetype = "dotted", linewidth = 1, alpha = 0.5) +
        scale_color_manual(name = "slope", breaks = c("package_fit", "alex_correction"),
                           values = c("package_fit" = "black", "correction" = "blue")) +
        theme_bw() +
        ggtitle(paste(peptides[i], "whole range"))
      fi = fl +
        xlim(0, 3*lod) +
        ylim(0, mean3lod) +
        ggtitle(paste(peptides[i], "inset"))
      f = fl + inset_element(fi, 0.05, 0.5 , 0.65, 1)
      nlp = nlps[[i]]
      of = nlp + f
      ggsave(file = paste(adress, paste0("Correction_fit_", peptide_name,".svg"), sep = "/"),
             plot = of,
             width = 46,
             height = 30,
             units = "cm")
      
    }
  }
  return(nlfs)
}

#Function to replace the 0 by small random values to enable standard deviation calculation in the calibration curve
ReplaceByNonNull = function(x){
  if(x == 0){
    abs(rnorm(1,0,1))
  } else{
    x
  }
}

#Write FOM from the linear and non-linear fit files
WriteFom <- function(lf, nlf){
  lin <- lf %>% 
    select(NAME, LOD, LOB, SLOPE, INTERCEPT) %>% 
    rename(LIN_LOD = LOD, LIN_LOB=LOB, LIN_SLOPE = SLOPE, LIN_INTERCEPT = INTERCEPT) %>% 
    distinct() %>% 
    group_by(NAME) %>% 
    arrange(desc(LIN_LOD)) %>% 
    slice(1)
  nonlin <- nlf %>% 
    select(NAME, LOD, LOB, slope, intercept) %>% 
    rename(NONLIN_LOD = LOD, NONLIN_LOB=LOB, NONLIN_SLOPE = slope, NONLIN_INTERCEPT = intercept) %>% 
    distinct() %>% 
    group_by(NAME) %>% 
    arrange(desc(NONLIN_LOD)) %>% 
    slice(1)
  fom = lin %>% 
    left_join(nonlin, by = "NAME")
}

#Comparing the LOB and LOD values obtained by both fitting methods
CompareFits <- function(fom){
  lob_graph = ggplot(fom, aes(x = LIN_LOB, y = NONLIN_LOB)) +
    geom_point()+
    annotate("segment", x = 1, y = 1, xend = 200000, yend = 200000) +
    theme_bw() +
    scale_x_continuous(name = "Linear LOB", breaks = c(10e0, 10e1, 10e2, 10e3, 10e4, 10e5), trans = "log") +
    scale_y_continuous(name = "Non-linear LOB", breaks = c(10e0, 10e1, 10e2, 10e3, 10e4, 10e5), trans = "log") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5))+
    ggtitle("LOB comparisons")
  lod_graph = ggplot(fom, aes(x = LIN_LOD, y = NONLIN_LOD)) +
    geom_point()+
    annotate("segment", x = 1, y = 1, xend = 200000, yend = 200000) +
    theme_bw()+
    scale_x_continuous(name = "Linear LOD", breaks = c(10e0, 10e1, 10e2, 10e3, 10e4, 10e5), trans = "log") +
    scale_y_continuous(name = "Non-linear LOD", breaks = c(10e0, 10e1, 10e2, 10e3, 10e4, 10e5), trans = "log") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5))+
    ggtitle("LOD comparisons")
  plt = lob_graph + lod_graph
  return(plt)
}

#Overview of LOD and LOQ for a given fit
PlotFomOverview <- function(fom){
  lob_lod_comparisons = lob_lod_comparisons <- ggplot(fom, aes(x = NONLIN_LOB, y = NONLIN_LOD)) +
    geom_abline(slope = 1, lty = 2, color = "#BDD7E7") +
    geom_point(shape = 21, fill = "#08519C") +
    theme_bw() +
    ylim(0, 1400) +
    xlim(0, 1400) +
    labs(x = "Limit of detection (amol)", y = "Limit of quantification (amol)") +
    ggtitle("Peptides analytical values") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "none")
  return(lob_lod_comparisons)
}


#Function to filter the calibration curve ions, so that they only contain ions matching the ones measured in the experiment. So a calibration curve can be calculated based
#on the experimental data only. Takes experiment dataframe, extract valid ions for each peptides and filters the calibration curve for these ions only.
FilterIons <- function(exp_df, cc_df){
  experiment_ions <- exp_df %>% 
    group_by(Peptide) %>% 
    select(Peptide, Fragment.Ion) %>% 
    unique()
  output_cc <- data.frame(matrix(nrow = 0, ncol = ncol(cc_df)))
  colnames(output_cc) <- colnames(cc_df)
  for(i in 1:length(unique(exp_df$Peptide))){
    p <- peptides[i]
    valid_ions <- experiment_ions[experiment_ions$Peptide == p, ]$Fragment.Ion
    filtered_ions <- cc_df[(cc_df$Peptide == p) & (cc_df$Fragment.Ion %in% valid_ions), ]
    output_cc <- rbind(output_cc, filtered_ions)
  }
  return(output_cc)
}

#Plot showing peptides abundances and how they relate to the calculated LOD an LOQ.
#There is most likely a more elegant way of doing it but I did not find it fast enough
PlotValidity <- function(df){
  lad <- length(unique(df[df$Above.LOD == TRUE & df$Isotope.Label.Type == "light", ]$Peptide))
  had <- length(unique(df[df$Above.LOD == TRUE & df$Isotope.Label.Type == "heavy", ]$Peptide))
  laq <- length(unique(df[df$Above.LOQ == TRUE & df$Isotope.Label.Type == "light", ]$Peptide))
  haq <- length(unique(df[df$Above.LOQ == TRUE & df$Isotope.Label.Type == "heavy", ]$Peptide))
  lbd <- length(unique(df[df$Above.LOD == FALSE & df$Isotope.Label.Type == "light", ]$Peptide))
  hbd <- length(unique(df[df$Above.LOD == FALSE & df$Isotope.Label.Type == "heavy", ]$Peptide))
  lbq <- length(unique(df[df$Above.LOQ == FALSE & df$Isotope.Label.Type == "light", ]$Peptide))
  hbq <- length(unique(df[df$Above.LOQ == FALSE & df$Isotope.Label.Type == "heavy", ]$Peptide))
  counts <- c(lad, had, laq, haq, lbd, hbd, lbq, hbq)
  groups <- c("light > LOD", "heavy > LOD", "light > LOQ", "heavy > LOQ", "light > LOD", "heavy > LOD", "light > LOQ", "heavy > LOQ")
  values <- c(T,T,T,T,F,F,F,F)
  plot_df <- bind_cols(counts = counts, groups = groups, values = values)
  peptide_fom <- ggplot(plot_df, aes(x = groups, y = counts, fill = values)) +
    geom_bar(stat = "identity", width = 0.4, position = position_dodge()) +
    theme_bw() +
    ggtitle("Peptide measurements analytical validity")
  return(peptide_fom)
}

#Quantifying proteins based on all ions à la Skyline. To use within the Quantification function because it needs some filtering of the df first
QuantifyProtein <- function(df, isotope_col, ion_quant_col){
  df <- df %>% 
    pivot_wider(names_from = {{isotope_col}}, values_from = {{ion_quant_col}}) %>% 
    filter(heavy != 0)
  L <- mean(df$light, na.rm = TRUE)
  H <- mean(df$heavy)
  return(L/H)
}



#A function that calculates the protein abundance the same way Skyline does, i.e. mean(light ions)/mean(heavy ions) for
#all light ions with heavy counterpart. This one taking into account that some peptides are below LOQ, in this case they are filtered out
#for a given prot/replicate combination but the protein abundance is calculated with other peptides when possible.It will also impute the missing
#value at the protein level with a downshift strategy
#From following input data:
#df = dataframe with peptide and ions data
#protein_col = df column identifying the proteins
#replicate_col = df column identifying the replicates
#ion_id = df column identifying ions. CAUTION: For skylie documents use name + charge as just fragment names can be repeated
#ion_quant = df column with quantitative value associated to ion
#isotope_col = df column specifying isotope label type
#loq_col is a boolean column telling if loq is reached or not
CalculateProteinAbundanceLOQimpute = function(df, protein_col, peptide_col, replicate_col, isotope_col,ion_quant_col, loq_col){
  prot = df %>% 
    pull({{protein_col}}) %>% 
    unique()
  replicate = df %>% 
    pull({{replicate_col}}) %>% 
    unique()
  #Checks if there is at least 1 pep > LOQ in the experiment
  val_pep_across = df %>%
    group_by({{peptide_col}}) %>% 
    summarise(validity_across_experiment = mean({{loq_col}})) %>% 
    ungroup()
  df = left_join(df, val_pep_across, by = c("Peptide"))
  #Checks whether all measurements within one replicate (light + heavy) are above LOQ
  val_pep_in = df %>%
    group_by({{peptide_col}}, {{replicate_col}}) %>% 
    summarise(validity_light_heavy = mean({{loq_col}})) %>% 
    ungroup()
  df = left_join(df, val_pep_in, by = c("Peptide", "Replicate.Name"))
  
  #Initialize empty vector for results output
  prot_v = c()
  replicate_v = c()
  int_v = c()
  
  #Iterate through proteins, if all peptides are valid -> calculate abundance in each replicate.
  #If not, calculates when possible and sets the rest as NA
  for(p in 1:length(prot)){
    current_prot = df %>% 
      filter({{protein_col}} == prot[p])
    valid_peptides = current_prot %>% filter(validity_across_experiment == 1) %>% pull({{peptide_col}}) %>% unique()
    if(length(valid_peptides)>0){
      current_prot = current_prot %>% 
        filter({{peptide_col}} %in% valid_peptides)
      for(r in 1:length(replicate)){
        ions_values = current_prot %>%
          filter({{replicate_col}} == replicate[r])
        prot_v = c(prot_v, prot[p])
        replicate_v = c(replicate_v, replicate[r])
        int_v = c(int_v, QuantifyProtein(ions_values, isotope_col = {{isotope_col}}, ion_quant_col = {{ion_quant_col}}))
      }
    } else {
      for(r in 1:length(replicate)){
        valid_peptides_subset = current_prot %>% filter({{replicate_col}} == replicate[r] & validity_light_heavy == 1)
        if (nrow(valid_peptides_subset > 0)){
          prot_v = c(prot_v, prot[p])
          replicate_v = c(replicate_v, replicate[r])
          int_v = c(int_v, QuantifyProtein(valid_peptides_subset, isotope_col = {{isotope_col}}, ion_quant_col = {{ion_quant_col}}))
        } else{
          prot_v = c(prot_v, prot[p])
          replicate_v = c(replicate_v, replicate[r])
          int_v = c(int_v, NA)
        }
      }
    }
  }
  protein_abundances = data.frame(prot_v, replicate_v, int_v) %>% 
    rename(Protein.Gene = prot_v, Replicate.Name = replicate_v, Protein.Intensity = int_v) %>% 
    mutate(Imputed_data = FALSE)
  
  #Now imputation part
  #Only calculate if missing values
  condition = df %>% 
    select(Protein.Gene, Condition, Replicate.Name)
  
  for(p in 1:length(prot)){
    if(sum(is.na(protein_abundances[protein_abundances$Protein.Gene == prot[p], ]$Protein.Intensity)) == 0){
      next
    } else{ #Make sure the all NA results are filtered out
      print(paste(prot[p], "-Imputation"))
      #Group by conditions
      current_prot = protein_abundances[protein_abundances$Protein.Gene == prot[p], ]
      n_missing = sum(is.na(current_prot$Protein.Intensity))
      if(n_missing == length(unique(current_prot$Replicate.Name))){
        print(paste(prot[p], "has no valid measurements at all"))
      } else {
        current_prot = left_join(current_prot, condition, by = c("Protein.Gene", "Replicate.Name"))
        grouped_sd = current_prot %>%
          dplyr::group_by(Condition) %>% 
          summarise(grouped_sd = sd(Protein.Intensity, na.rm = TRUE))
        print(grouped_sd)
        average_sd = mean(grouped_sd$grouped_sd, na.rm = TRUE)
        print(average_sd)
        min_val = min(current_prot$Protein.Intensity, na.rm = TRUE)
        print(min_val)
        mean_val = min_val * 0.8
        filling_values = rnorm(n_missing, mean_val, average_sd)
        print(filling_values)
        protein_abundances[protein_abundances$Protein.Gene == prot[p] & is.na(protein_abundances$Protein.Intensity), ]$Imputed_data = TRUE
        protein_abundances[protein_abundances$Protein.Gene == prot[p] & is.na(protein_abundances$Protein.Intensity), ]$Protein.Intensity = filling_values
        
      }
    }
  }
  return(protein_abundances)
}
