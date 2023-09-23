#V2

library(tidyverse)
library(summarytools) #Not sure still using this?
library(broom)  # devtools::install_github("tidymodels/broom")
library(vegan)
library(RRPP)
require(pvclust)
library(Hmisc)

#graphics
library(cowplot)
library(ggrepel)
library(ggh4x)
library(GGally)
library(patchwork)

library(knitr)
library(kableExtra)

set.seed(28738)

setwd(
  "C:/Users/prade/OneDrive - The University of Liverpool/Obsidian/Liverpool/ARDAT/WP1 Mechanistic understanding of innate immune response to AAV/PH009/A1/Manually Debarcoded_PHOnly"
)

# setwd(
#   "C:/Users/prade/OneDrive - The University of Liverpool/Obsidian/Liverpool/ARDAT/WP1 Mechanistic understanding of innate immune response to AAV/PH009/A1/Manually Debarcoded_LHOnly/Boots/"
# )

# setwd("C:/Users/prade/OneDrive - The University of Liverpool/Obsidian/Liverpool/ARDAT/WP1 Mechanistic understanding of innate immune response to AAV/PH009/Adam Lister/manually debarcoded")


#Read Input====
chrCSVList <- list.files(getwd(), pattern = "\\.csv$")
for (i in 1:length(chrCSVList))
{
  assign(substr(chrCSVList[i], 1, nchar(chrCSVList[i]) - 4), read_csv(chrCSVList[i]))
}

rm(i, chrCSVList)

#Helper functions====
fVectorChangeArrowPos <- function (data, vecChangeArrowPos)
{
  data$vc_x_new = data$vc_xend + vecChangeArrowPos * (data$vc_x - data$vc_xend)
  data$vc_y_new = data$vc_yend + vecChangeArrowPos * (data$vc_y - data$vc_yend)
  data$vc_z_new = data$vc_zend + vecChangeArrowPos * (data$vc_z - data$vc_zend)
  return(data)
}

plotPCA <-
  function(Data = pca_fit,
           DataAug = ggPCAData,
           X,
           Y,
           Title = graphTitle,
           Shape = shapeBy,
           Colour = colorBy) {
    pcX = paste0(".fittedPC", X)
    pcY = paste0(".fittedPC", Y)
    
    ggPlot <- DataAug %>%
      ggplot(., aes(.data[[pcX]], .data[[pcY]])) +
      ggtitle(Title) +
      scale_shape_manual(values = c(18, 20)) +
      theme_half_open(12) +
      theme(legend.title = element_blank()) +
      geom_point(
        data = DataAug %>% filter(status == "summary") %>% arrange(patient_id),
        size = 4,
        aes(shape = .data[[Shape]],  color = .data[[Colour]])
      )
    
    if (class(Data) == "prcomp") {
      print("4")
      ggPlot <- ggPlot +
        ylab(paste("PC", Y, paste0(round(
          tidy(Data, matrix = "eigenvalues")[Y, 3] * 100
        ), "%"))) +
        xlab(paste("PC", X, paste0(round(
          tidy(Data, matrix = "eigenvalues")[X, 3] * 100
        ), "%")))
    }
    
    else if (class(Data) == "princmp")
    {
      ggPlot <- ggPlot +
        ylab(paste("PC", Y, paste0(round(
          Data$vars[Y] / sum(Data$vars) * 100
        ), "%"))) +
        xlab(paste("PC", X, paste0(round(
          Data$vars[X] / sum(Data$vars) * 100
        ), "%")))
    }
    
    
    return(ggPlot)
  }

savePlot <-
  function(plot,
           title = ggPCA$labels$title,
           subtitle = "_5_Components",
           h = 5,
           w = 6,
           background = "white") {
    ggsave(
      plot,
      filename = paste0(title, subtitle, ".png"),
      device = "png",
      dpi = 800,
      width = w,
      height = h,
      units = "in",
      bg = background
    )
  }

plotVecChange <-
  function(plot,
           Data = vecChange,
           X,
           Y,
           Colour = colorBy,
           arrowStyle = arrowStyleVecChange) {
    ggVec <-  plot + geom_segment(
      data = Data,
      aes(
        x = .data[[paste0(X, "_new")]],
        y = .data[[paste0(Y, "_new")]],
        xend = .data[[paste0(X, "end")]],
        yend = .data[[paste0(Y, "end")]],
        color = .data[[Colour]],
        alpha = 0.4
      ),
      inherit.aes = FALSE,
      arrow = arrowStyle,
      linetype = "solid",
      show.legend = FALSE
    ) +
      geom_segment(
        data = Data,
        aes(
          xend = .data[[paste0(X, "_new")]],
          yend = .data[[paste0(Y, "_new")]],
          x = .data[[X]],
          y = .data[[Y]],
          color = .data[[Colour]],
          alpha = 0.4
        ),
        inherit.aes = FALSE,
        linetype = "solid",
        show.legend = FALSE
      )
    return(ggVec)
  }

arrowStyleLoadings <- arrow(
  angle = 40,
  ends = "first",
  type = "closed",
  length = grid::unit(8, "pt")
)

arrowStyleVecChange <- arrow(
  angle = 20,
  ends = "first",
  type = "closed",
  length = grid::unit(6, "pt")
)

dataCleanup <- function(tblColumn) {
  #Replace all 0s with minimum value / 3
  minValue <- min(tblColumn[tblColumn > 0], na.rm = TRUE)
  
  if (is.finite(minValue) & !is.na(minValue)) {
    tblColumn[tblColumn == 0 | NA] <- minValue / 3
  }
  else {
    #Case where all values in column are 0 or something else, add in very small number
    tblColumn[tblColumn == 0 | NA] <- 0.000000000001
  }
  
  #If variance of column is 0, replace the first element with a smaller number
  if (is.numeric(tblColumn) & var(tblColumn) == 0) {
    tblColumn[1] <- tblColumn[1] - tblColumn[1] / 1000
  }
  return(tblColumn)
}



# groupData <- allData %>%
#   group_by(condition, metacluster_id, patient_id, sero, sample_id) %>%
#   add_column(grp_id = group_indices(.), .before = ".cell")
#
#Bootstrapping====
{
  # splitData <- groupData %>% ungroup() %>% nest(.by = grp_id)
  #
  # {
  #   summaryStats = NULL
  #   cellDistribution = NULL
  #   clusterDistribution = NULL
  #
  #   nGroups = length(splitData[[1]])
  #   J = 500 # no of bootstrap repeats
  #   I <- nGroups
  # }
  #
  # rm(allData, groupData)
  # gc(full = TRUE)
  #
  # performBootstrap <- function(I) {
  #   for (i in I)
  #   {
  #     N = length(splitData[[2]][[i]][[1]])
  #
  #     temp <-
  #       splitData[[2]][[i]] %>%
  #       group_by(condition, metacluster_id, patient_id, sero, sample_id)
  #
  #     groupTime <- Sys.time()
  #     print("========================")
  #     print(paste("Group ", i, "of ", nGroups))
  #     print(paste("Sampling/Replacing ", N, "cells"))
  #     print("========================")
  #
  #     summaryStats = NULL
  #     cellDistribution = NULL
  #     clusterDistribution = NULL
  #
  #     for (j in 1:J) {
  #       temp <- temp %>%
  #         slice_sample(n = N, replace = TRUE) %>%
  #         ungroup()
  #
  #       #temp %>%
  #       #select(1:7) %>% mutate(bootstrap = j) %>% bind_rows(cellDistribution, .) -> cellDistribution
  #
  #       temp %>%
  #         select(!c(.cell, cluster_id)) %>%
  #         summarise(across(where(is.character), unique),
  #                   #across(where(is.numeric), ~ dataCleanup(.)  %>% median(na.rm = TRUE))) %>%
  #                   across(where(is.numeric), median)) %>%
  #         mutate(bootstrap = j) %>%
  #         bind_rows(summaryStats, .) -> summaryStats
  #
  #       print(paste("Bootstrap ", j , "of ", J, "; of group ", i, "of ", nGroups))
  #
  #     }
  #
  #     # write_csv(
  #     #   cellDistribution,
  #     #   file = paste("CellDistributions_", temp[[2]][[1]], ".csv"),
  #     #   append = TRUE
  #     # )
  #     write_csv(
  #       summaryStats,
  #       file = paste("Medians_", temp[[2]][[1]], "_", temp[[7]][[1]], ".csv"),
  #       append = FALSE
  #     )
  #
  #     rm(summaryStats)
  #     rm(cellDistribution)
  #     gc(verbose = TRUE, full = TRUE)
  #
  #   }
  # }
  #
  
  # hist(table(
  #   cellDistribution %>% filter(sample_id == "LV172_UT_HS", metacluster_id == "pDCs") %>% .$.cell
  # ))
  
  # bootData = tibble()
  # chrCSVList <-
  #   list.files(getwd(), pattern = "^(Medians).*\\.csv$")
  #
  # for (i in 1:length(chrCSVList))
  # {
  #   read_csv(chrCSVList[i], col_names = FALSE) %>%
  #     # mutate(across(where(is.character), unique),
  #     #        across(where(is.numeric), dataCleanup)) %>%
  #     bind_rows(bootData, .) -> bootData
  #
  #   print(paste("Reading ....", chrCSVList[i]))
}

# colnames(bootData) <- read_csv("allData.csv") %>%
#   select(-c(".cell", "cluster_id")) %>%
#   colnames() %>%
#   c(., "bootstrap_no")

# colnames(bootData) <- c("sample_id", "condition", "sero", "patient_id", "metacluster_id", "bootstrap_no",
#                         "CD11c", "CD123", "CD127 (IL-7R)", "CD14_AND_TNFA", "CD16", "CD183 (CXCR3)",
#                         "CD19", "CD25 (IL-2R)", "CD279 (PD-1)", "CD28", "CD3", "CD4",
#                         "CD45RA", "CD56 (NCAM)", "CD69", "CD8a", "CXCR5", "FoxP3", "Granzyme B",
#                         "HLA-DR", "IFNbeta", "IFNGamma", "IL-10", "IL-17A", "IL-2", "Ki 67",
#                         "Perforin")
# write_csv(bootData, "Medians_Compiled.csv")


# rm(i, chrCSVList)
#}

#Read files====
# {
groupData <- Medians_Compiled %>%
  select(-c(bootstrap_no)) %>% #bootstrap_no for Bootstraped data; "cluster_id", ".cell" for all data
  group_by(condition, metacluster_id, sample_id)

#
summaryData <-
  groupData %>%
  summarise(across(where(is.character), unique),
            across(where(is.numeric), median)) %>%
  add_column(status = "summary", .before = "sample_id") %>%
  write_csv(file = "summaryData.csv")

groupData %>%
  add_column(status = "raw", .before = "sample_id") %>%
  bind_rows(summaryData) -> bindData
#
allData_Clean <- bindData %>%
  ungroup() %>%
  group_by(metacluster_id, condition) %>%
  mutate(across(where(is.numeric), dataCleanup)) %>%
  write_csv("allData_Clean.csv")

allData_CleanR <- allData_Clean %>%
  mutate(across(where(is.numeric), ~ decostand(.x, method = "range")))
#
#
# }
#
# {
#Testing vars====
dataVar <-
  summaryData_Clean %>% group_by(status,
                                 sample_id,
                                 condition,
                                 sero,
                                 patient_id,
                                 metacluster_id) %>% filter(sero != "AAV2" &
                                                              metacluster_id == "Monocytes")

loadingsKey = StainingKey %>% filter(marker_class == "state") %>%
  select(-c(fcs_colname, marker_class)) %>%
  distinct()
graphTitle = paste(" NK Cells")
loadingsNumber = 3
loadingsMultiple = 1
labelIndividuals = FALSE
drawVectorChange = TRUE
vecChangeArrowPos = 2 / 3
vecChangeArrowLength = 8
graphTitle = "Title"
statEllipseBy = 'sero'
shapeBy = 'condition'
colorBy = 'sero'
imageBG = "white"
X = 1
Y = 2
Z = 3


#Main Function ====
PCA <-
  function(dataVar,
           X = 1,
           Y = 2,
           Z = 3,
           genComponentPlots = TRUE,
           genMultiANOVA = TRUE,
           performCluster = TRUE,
           startIndex = 9,
           loadingsKey,
           loadingsNumber = 3,
           loadingsMultiple = 1,
           labelIndividuals = FALSE,
           perIndividualEllipse = TRUE,
           drawVectorChange = TRUE,
           vecChangeArrowPos = 2 / 3,
           graphTitle = "Title",
           pairsPlot = FALSE,
           sparse = FALSE,
           statEllipseBy = 'sero',
           shapeBy = 'condition',
           colorBy = 'sero',
           imageBG = "white")
  {
    loadingsKey <- loadingsKey[[1]]
    pcX = paste0(".fittedPC", X)
    pcY = paste0(".fittedPC", Y)
    pcZ = paste0(".fittedPC", Z)
    
    pca_fit = NULL
    
    #sparse PCA fit ====
    if (sparse) {
      dfDataVar <- dataVar %>%
        ungroup() %>%
        select(where(is.numeric)) %>% # retain only numeric columns
        scale() %>%
        as.data.frame()
      
      pca_fit <-
        princmp(~ .,
                data = dfDataVar,
                method = "sparse",
                sw = TRUE)
      
      
      #Data Prep ====
      Rotations <- pca_fit$scoef %>%
        as_tibble(rownames = "column") %>%
        rename_with(~ (gsub("Comp.", ".fittedPC", .x, fixed = TRUE))) %>% #Renaming for better integration
        mutate(column = gsub("`", "", column)) %>% #Removing the `` added by princmp to names
        select(all_of(c("column", pcX, pcY))) %>%
        filter(column %in% loadingsKey) %>%
        distinct()
      
      allLoadings <- bind_rows(
        Rotations %>%
          slice_min(order_by = .[[pcX]], n = loadingsNumber),
        Rotations %>%
          slice_max(order_by = .[[pcX]], n = loadingsNumber),
        Rotations %>%
          slice_min(order_by = .[[pcY]], n = loadingsNumber),
        Rotations %>%
          slice_max(order_by = .[[pcY]], n = loadingsNumber)
      ) %>% distinct()
      
      ggPCAData <-
        dataVar %>%
        bind_cols(pca_fit$scores %>% as_tibble()) %>%  # add original dataset back in
        group_by(across(c(all_of(
          dataVar %>% group_vars()
        )))) %>%
        rename_with(~ (gsub("Comp.", ".fittedPC", .x, fixed = TRUE)))
    }
    else{
      #PCA fit ====
      pca_fit <- dataVar %>%
        ungroup() %>%
        select(where(is.numeric)) %>% # retain only numeric columns
        prcomp(scale = TRUE) # do PCA on scaled data; data is centered by default.
      
      #Data Prep ====
      Rotations <- pca_fit %>%
        tidy(matrix = "rotation") %>%
        pivot_wider(
          names_from = "PC",
          names_prefix = ".fittedPC",
          values_from = "value"
        ) %>%
        select(c(column, pcX, pcY)) %>%
        filter(column %in% loadingsKey) %>%
        distinct()
      
      allLoadings <- bind_rows(
        Rotations %>%
          slice_min(order_by = .[[pcX]], n = loadingsNumber),
        Rotations %>%
          slice_max(order_by = .[[pcX]], n = loadingsNumber),
        Rotations %>%
          slice_min(order_by = .[[pcY]], n = loadingsNumber),
        Rotations %>%
          slice_max(order_by = .[[pcY]], n = loadingsNumber)
      ) %>% distinct()
      
      #Data augmentation
      ggPCAData <-
        pca_fit %>%
        augment(dataVar) %>%  # add original dataset back in
        group_by(across(c(
          all_of(dataVar %>% group_vars()), .rownames
        )))
      
    }
    
    
    
    
    #P6 PCA plotting====
    {
      ggPCA <-
        plotPCA(pca_fit, ggPCAData, X, Y, graphTitle, shapeBy, colorBy)
      
      #P7 Components plotting====
      if (genComponentPlots)
      {
        ggPCA0 <-
          ggPCA + guides(shape = guide_legend(nrow = 1)) + guides(fill = guide_legend(nrow = 1))
        
        ggPCA1 <-
          plotPCA(pca_fit, ggPCAData, X, Z, NULL, shapeBy, colorBy) + guides(color = guide_legend(nrow = 1))
        
        ggPCA2 <-
          plotPCA(pca_fit, ggPCAData, Y, Z, NULL, shapeBy, colorBy)
        
        ggComponents <-
          ggPCA0 + (ggPCA1 / ggPCA2) + plot_annotation(tag_levels = c('A', '1'), tag_sep = '.') + plot_layout(widths = c(2, 1), guides = 'collect') &
          theme(legend.position = "bottom", legend.box = "horizontal")
        
        savePlot(
          ggComponents,
          ggPCA$labels$title,
          background = imageBG,
          subtitle = "_7_Components"
        )
      }
      
      if (labelIndividuals) {
        ggPCA <-
          ggPCA +
          geom_label_repel(
            data = ggPCAData %>% filter(status == "summary"),
            aes(label = sample_id),
            min.segment.length = 0,
            box.padding = 2,
            segment.alpha = 0.5,
            segment.linetype = 3,
            segment.curvature = -1e-30,
            alpha = 0.8,
            show.legend = FALSE
          )
      }
      #Per donor ellipse====
      if (perIndividualEllipse) {
        ggPCA <-
          ggPCA + stat_ellipse(
            data = ggPCAData %>% filter(status != "summary"),
            aes(.data[[paste0(".fittedPC", X)]], .data[[paste0(".fittedPC", Y)]], group = sample_id, fill = .data[[colorBy]]),
            inherit.aes = FALSE,
            type = "norm",
            geom = "polygon",
            alpha = 0.1,
            show.legend = FALSE
          )
      }
      
      else{
        ggPCA <-
          ggPCA + stat_ellipse(
            data = ggPCAData %>% filter(status == "summary"),
            aes(.data[[paste0(".fittedPC", X)]], .data[[paste0(".fittedPC", Y)]], fill = .data[[colorBy]]),
            inherit.aes = FALSE,
            type = "norm",
            geom = "polygon",
            alpha = 0.1,
            show.legend = FALSE
          )
      }
      
      savePlot(ggPCA,
               ggPCA$labels$title,
               background = imageBG,
               subtitle = "_6_Raw")
      
      #P8 Biplot====
      internalMultiple = min(min(abs(layer_scales(ggPCA)$x$range$range)) / max(abs(allLoadings[[pcX]])),
                             min(abs(layer_scales(ggPCA)$y$range$range)) / max(abs(allLoadings[[pcY]])))
      
      ggPCAL <- ggPCA + geom_segment(
        data = allLoadings,
        aes(
          .data[[pcX]] * internalMultiple * loadingsMultiple,
          .data[[pcY]] * internalMultiple * loadingsMultiple,
          xend = 0,
          yend = 0
        ),
        arrow = arrowStyleLoadings,
        inherit.aes = FALSE
      ) +
        geom_label_repel(
          data = allLoadings,
          aes(
            .data[[pcX]] * internalMultiple * loadingsMultiple,
            .data[[pcY]] * internalMultiple * loadingsMultiple,
            label = column
          ),
          inherit.aes = FALSE,
          min.segment.length = 0,
          box.padding = 1,
          segment.alpha = 0.5,
          segment.linetype = 5,
          segment.curvature = -1e-30,
          alpha = 0.8
        )
      
      savePlot(ggPCAL,
               ggPCA$labels$title,
               background = imageBG,
               subtitle = "_8_Load")
      
      #P5 Variances====
      if (class(pca_fit) == "prcomp") {
        featurePercentage <- pca_fit$rotation %>%
          abs() %>%
          sweep(2, colSums(.), "/") %>%
          as_tibble(rownames = "variable") %>%
          pivot_longer(cols = 2:(length(.)), names_to = "PC") %>%
          filter(!grepl("^CD|CXCR5|HLA-DR|\\.cell|cluster_id", variable)) %>%
          filter(PC %in% c(paste0("PC", X),
                           paste0("PC", Y),
                           paste0("PC", Z))) %>%
          #arrange(desc(value)) %>% mutate(variable = factor(variable, levels = variable)) %>%   # This trick update the factor levels
          ggplot(aes(value * 100, variable)) +
          geom_col() +
          xlab(paste0("% Variance explained")) +
          ylab("") +
          theme_half_open(12) +
          facet_wrap2( ~ PC)
        
        loadingsPercentage <- pca_fit %>%
          tidy(matrix = "pcs") %>%
          mutate(percent = percent * 100) %>%
          mutate(cumulative = cumulative * 100) %>%
          ggplot(aes(PC, percent)) +
          geom_col() +
          #geom_line(aes(PC, cumulative)) +
          scale_x_continuous("PCs",
                             breaks = c(1:10),
                             limits = c(0, 10)) +
          ylab("% Variance explained") +
          theme_half_open(12)
        
        savePlot(
          featurePercentage / loadingsPercentage + plot_annotation(title = paste(
            ggPCA$labels$title, "variances explained"
          )),
          #Possible bug
          ggPCA$labels$title,
          background = imageBG,
          subtitle = "_5_Var"
        )
      }
      #Calculate vector change coordinates====
      if (drawVectorChange) {
        {
          temp <- ggPCAData %>%
            filter(status == "summary") %>%
            arrange(sero, sample_id)
          vecChange <- data.frame(
            patient_id = temp$patient_id[c(TRUE, FALSE)],
            sample_id = temp$sample_id[c(TRUE, FALSE)],
            sero = temp$sero[c(TRUE, FALSE)],
            metacluster_id = temp$metacluster_id[c(TRUE, FALSE)],
            status = temp$status[c(TRUE, FALSE)],
            condition = temp$condition[c(TRUE, FALSE)],
            
            # vc_xend = temp$.fittedPC1[c(TRUE, FALSE)],
            # vc_x = temp$.fittedPC1[c(FALSE, TRUE)],
            # vc_yend = temp$.fittedPC2[c(TRUE, FALSE)],
            # vc_y = temp$.fittedPC2[c(FALSE, TRUE)]
            ##Uncomment to flip arrow direction
            
            vc_x = temp[[pcX]][c(TRUE, FALSE)],
            vc_xend = temp[[pcX]][c(FALSE, TRUE)],
            vc_y = temp[[pcY]][c(TRUE, FALSE)],
            vc_yend = temp[[pcY]][c(FALSE, TRUE)],
            vc_z = temp[[pcZ]][c(TRUE, FALSE)],
            vc_zend = temp[[pcZ]][c(FALSE, TRUE)]
          ) %>% as_tibble()
          
          rm(temp)
        }
        
        vecChange <-
          vecChange %>% fVectorChangeArrowPos(vecChangeArrowPos)
        
        #P1 Vector change====
        ggPCAV <- plotVecChange(ggPCA,
                                vecChange,
                                "vc_x",
                                "vc_y",
                                Colour = colorBy)
        
        savePlot(ggPCAV,
                 ggPCA$labels$title,
                 background = imageBG,
                 subtitle = "_1_Vec")
        
        ggPCALV <- plotVecChange(ggPCAL,
                                 vecChange,
                                 "vc_x",
                                 "vc_y",
                                 Colour = colorBy)
        
        #P2 Vector Load change====
        savePlot(ggPCALV,
                 ggPCA$labels$title,
                 background = imageBG,
                 subtitle = "_2_LoadVec")
        
        if (genComponentPlots)
        {
          ggPCA0Vec <- plotVecChange(ggPCA0,
                                     vecChange,
                                     "vc_x",
                                     "vc_y",
                                     Colour = colorBy)
          
          ggPCA1Vec <- plotVecChange(ggPCA1,
                                     vecChange,
                                     "vc_x",
                                     "vc_z",
                                     Colour = colorBy)
          
          ggPCA2Vec <- plotVecChange(ggPCA2,
                                     vecChange,
                                     "vc_y",
                                     "vc_z",
                                     Colour = colorBy)
          
          ggComponentsVec <-
            ggPCA0Vec + (ggPCA1Vec / ggPCA2Vec) + plot_annotation(tag_levels = c('A', '1'), tag_sep = '.') + plot_layout(widths = c(2, 1), guides = 'collect') &
            theme(legend.position = "bottom", legend.box = "horizontal")
          
          #P3 Components Vec====
          savePlot(
            ggComponentsVec,
            ggPCA$labels$title,
            background = imageBG,
            subtitle = "_3_ComponentsVec"
          )
          
        }
        
      }
      
      #P9,10 Pairs plot====
      if (pairsPlot) {
        #Pairs plot
        ggPairsSero <- ggpairs(ggPCAData,
                               columns = 35:39,
                               aes(color = condition, shape = sero))
        savePlot(ggPairsSero,
                 ggPCA$labels$title,
                 background = imageBG,
                 subtitle = "_9_PairsPlot_Tmt")
        
        ggPairsTmt <- ggpairs(ggPCAData,
                              columns = 35:39,
                              aes(color = sero, shape = patient_id))
        savePlot(ggPairsTmt,
                 ggPCA$labels$title,
                 background = imageBG,
                 subtitle = "_10_PairsPlot_Sero")
      }
    }
    
    #P4 Univariate graphs====
    {
      tblDataLong <-
        pivot_longer(
          dataVar %>% filter(status == "summary") %>% ungroup(),
          cols = -all_of(dataVar %>% group_vars()),
          names_to = "variable"
        ) %>%
        filter(!grepl("^CD|CXCR5|HLA-DR", variable))
      
      ggPairwise <-
        ggplot(tblDataLong ,
               aes(x = interaction(condition, sero, sep = "!"), y = value)) +
        theme_half_open(12) +
        scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Seropositivity") +
        stat_summary(fun.y = mean, geom = "point") +
        stat_summary(fun.data = mean_se, geom = "errorbar") +
        geom_point(position = "jitter", aes(color = patient_id)) +
        ggtitle(graphTitle) +
        facet_wrap2( ~ variable)
      
      savePlot(
        ggPairwise,
        ggPCA$labels$title,
        background = imageBG,
        w = 12,
        subtitle = "_4_Pairwise"
      )
    }
    
    #P11 Hierarchical clustering ====
    if (performCluster) {
      png(
        filename = paste0(ggPCA$labels$title, "_10_HeirClus.png"),
        width = 6,
        height = 5,
        units = "in",
        res = 800
      )
      ggHClust <- ggPCAData %>%
        ungroup() %>%
        mutate(.rownames = paste(sample_id, " ", sero)) %>%
        filter(status == "summary") %>%
        select(.rownames | all_of(starts_with(".fitted"))) %>%
        column_to_rownames(var = ".rownames") %>%
        t() %>%
        as.matrix() %>%
        pvclust::pvclust(
          method.hclust = "complete",
          method.dist = "correlation",
          nboot = 1000,
          parallel = TRUE
        ) %>% plot()
      dev.off()
    }
    
  }

# #Multivariate ANOVA====
dataVar <- summaryData_Clean %>%
  filter(sero != "AAV2") %>%   #select(which(!grepl("^CD|CXCR5|HLA-DR", names(.)))) %>%
  ungroup()

performStats <- function(dataVar, permutations = 1, title) {
  print("....Calculating distance matrix")
  
  brayDataVar <- dataVar %>%
    ungroup() %>%
    select(where(is.numeric)) %>%
    vegdist()
  
  rrppSummaryData <- rrpp.data.frame(
    distances = brayDataVar ,
    sero = dataVar$sero,
    condition = dataVar$condition,
    patient_id = dataVar$patient_id,
    sample_id = dataVar$sample_id,
    metacluster_id = dataVar$metacluster_id
  )
  
  print("....Fitting model")
  
  rrppModelFitTypeIII <-
    lm.rrpp(
      distances ~ condition * metacluster_id,
      data = rrppSummaryData,
      iter = 1,
      SS.type = "III"
    )
  
  print("....Analysing model")
  rrppModelFitTypeIII %>%
    anova.lm.rrpp() %>%
    .$table %>%
    kable(label = "Q", caption = "R") %>%
    kable_classic() %>%
    save_kable(file = paste0(title, "_RRPPModel_ANOVATable.png"),
               density = 800)
  
  print("....Generating pairwise")
  pairwise(rrppModelFitTypeIII, groups = rrppSummaryData$sample_id) %>%
    summary.pairwise() %>%
    kable() %>%
    kable_classic() %>%
    save_kable(file = paste0(title, "_RRPPModel_Pairwise.png"),
               density = 800) #Pairwise comparision of LS means
  
  getANOVAStats(rrppModelFitTypeIII, stat = "cohenf")[["cohenf"]][, 1] ^
    2  %>%
    kable() %>%
    kable_classic()
  
  #u = df of predictors
  getANOVAStats(rrppModelFitTypeIII, stat = "cohenf")[["cohenf"]][, 1] ^
    2  %>%
    pwr.f2.test(
      u = 1,
      f2 = .,
      power = 0.8,
      sig.level = 0.05
    )
  
  print("....Calculating dispersion matrix")
  summaryDispersal <-
    betadisper(brayDataVar, group = dataVar$sample_id, type = "centroid")
  
  print("Analysing distance matrix")
  permutest(summaryDispersal)$tab %>%
    kable() %>%
    kable_classic() %>%
    save_kable(file = paste0(title, "_BetaDis_ANOVATable.png"),
               density = 800)
  
  print("Generating pairwise")
  TukeyHSD(summaryDispersal) %>%
    kable() %>%
    kable_classic() %>%
    save_kable(file = paste0(title, "_BetaDis_Pairwise.png"),
               density = 800)
  
}

#Wrapper function====
performPCA <-
  function (dataVar,
            PC_X = 2,
            PC_Y = 3,
            loadingKey = StainingKey,
            prefix = "Subset",
            loadingsNumber = 3,
            drawVecChange = TRUE,
            labIndividuals = FALSE,
            elipIndividuals = TRUE,
            componentPlot = TRUE,
            sparsePCA = TRUE) {
    print("Preprocessing....")
    
    StainingKeyType <-
      loadingKey %>% filter(marker_class == "type") %>%
      select(-c(fcs_colname, marker_class)) %>%
      distinct()
    
    StainingKeyState <-
      loadingKey %>% filter(marker_class == "state") %>%
      filter(!row_number() %in% c(1:9)) %>%
      select(-c(fcs_colname, marker_class)) %>%
      distinct()
    
    print("Lineage....")
    PCA(
      X = PC_X,
      Y = PC_Y,
      dataVar,
      loadingsKey = StainingKeyType,
      graphTitle = paste(prefix, "Lineages"),
      statEllipseBy = "metacluster_id",
      colorBy = "metacluster_id",
      drawVectorChange = FALSE,
      pairsPlot = FALSE,
      genMultiANOVA = FALSE,
      perIndividualEllipse = FALSE,
      performCluster = FALSE,
      sparse = sparsePCA
    )
    
    print("CD8s....")
    PCA(
      subset(dataVar, metacluster_id == 'CD8 T cell'),
      X = PC_X,
      Y = PC_Y,
      loadingsKey = StainingKeyState,
      graphTitle = paste(prefix, "CD8 T Cell"),
      drawVectorChange = drawVecChange,
      labelIndividuals = labIndividuals,
      pairsPlot = componentPlot,
      perIndividualEllipse = elipIndividuals,
      sparse = sparsePCA
    )
    
    # print("CD4s....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'CD4 T cell'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "CD4 T Cell"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    #
    # print("Bs....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'B cells'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "B Cell"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    #
    # print("Monocytes....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'Monocytes'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "Monocytes"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    #
    # print("pDCs....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'pDCs'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "pDCs"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    #
    # print("NKs....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'NK Cells'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "NK Cells"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    
    #Luisa's vars
    # print("Misc....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'Misc'),
    #   X = PC_X,
    #   Y = PC_Y,
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, "Misc"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals,
    #   pairsPlot = componentPlot,
    #   perIndividualEllipse = elipIndividuals,
    #   sparse = sparsePCA
    # )
    
    #Adam's vars
    # print("NK Ts....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'NK T Cells'),
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, " NK T Cells"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals
    # )
    #
    # print("CD3 only....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'CD3 Only'),
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, " CD3 Only"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals
    # )
    #
    # print("CD3 CD4 CD8 only....")
    # PCA(
    #   subset(dataVar, metacluster_id == 'CD3 CD4 CD8'),
    #   loadingsKey = StainingKeyState,
    #   graphTitle = paste(prefix, " CD3 CD4 CD8"),
    #   drawVectorChange = drawVecChange,
    #   labelIndividuals = labIndividuals
    # )
    
  }

#Run wrapper====
performPCA(
  allData_Clean %>%
    group_by(status, sample_id, condition, sero, patient_id, metacluster_id) %>%
    filter(patient_id != "LV081"),
  #%>% select(which(!grepl("^CD|CXCR5|HLA-DR", names(.)))),
  prefix = "Sparse",
  PC_X = 1,
  PC_Y = 2,
  drawVecChange = TRUE,
  labIndividuals = FALSE,
  componentPlot = FALSE,
  elipIndividuals = TRUE,
  sparsePCA = TRUE
)
