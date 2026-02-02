# count_cells_over_median.R
# R script to use mean intensity data from CellProfiler object results files
# to count cells exceeding a threshold intensity value defined as the median of the control group
# and summarize by treatment and sample

# Theresa Swayne, Columbia University, 2025-26
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# To use: Run the script.
# Will be prompted for a file.

# ---- Setup ----
require(tidyverse)
require(readr)
require(stringr)
require(reshape2)

# ---- Prompt for an object file ----
# no message will be displayed
objectFile <- file.choose()

# Read the data from the file
objectData <- read_csv(objectFile,
                    locale = locale())

# template provided for calculating vimentin positivity -- but this measurement is not appropriate if cells are detected using vimentin

# determine threshold for positivity by median of control group
# medians <- objectData %>% 
# group_by(Metadata_Treatment) %>% 
# summarise(nCells = n(),
#           medianSMA = median(Intensity_MeanIntensity_SMA),
#           medianPDGFR = median(Intensity_MeanIntensity_PDGFRalpha),
#           medianVimentin = median(Intensity_MeanIntensity_Vimentin))

medians <- objectData %>% 
  group_by(Metadata_Treatment) %>% 
  summarise(nCells = n(),
            medianSMA = median(Intensity_MeanIntensity_SMA),
            medianPDGFR = median(Intensity_MeanIntensity_PDGFRalpha))

# get the control group median from the first row of the summary
threshSMA <- medians$medianSMA[1]
threshPDGFR <- medians$medianPDGFR[1]
#threshVim <- medians$medianVimentin[1]
            
# ---- Count positive cells within each sample and treatment ---
# counts <- objectData %>% 
#   group_by(Metadata_Treatment, Metadata_Sample) %>% 
#   summarise(nCells = n(),
#             PositiveSMA = sum(Intensity_MeanIntensity_SMA > threshSMA),
#             PercentPositiveSMA = PositiveSMA/nCells,
#             PositivePDGFR = sum(Intensity_MeanIntensity_PDGFRalpha > threshPDGFR),
#             PercentPositivePDGFR = PositivePDGFR/nCells,
#             PositiveVimentin = sum(Intensity_MeanIntensity_Vimentin > threshVim),
#             PercentPositiveVimentin = PositiveVimentin/nCells)

counts <- objectData %>% 
  group_by(Metadata_Treatment, Metadata_Sample) %>% 
  summarise(nCells = n(),
            PositiveSMA = sum(Intensity_MeanIntensity_SMA > threshSMA),
            PercentPositiveSMA = PositiveSMA/nCells,
            PositivePDGFR = sum(Intensity_MeanIntensity_PDGFRalpha > threshPDGFR),
            PercentPositivePDGFR = PositivePDGFR/nCells)

# make the column headers show the threshold used
SMA_name <- paste0("SMA Mean >",round(threshSMA, digits = 5))
PDGFR_name <- paste0("PDGFR Mean >", round(threshPDGFR, digits = 5))
#Vim_name <- paste0("Vimentin Mean >", round(threshVim, digits = 5))

# prepare a file for saving with the informative column headers
counts_save <- counts
# counts_save <- counts_save %>% 
#   rename(!!SMA_name := PositiveSMA) %>%
#   rename(!!PDGFR_name := PositivePDGFR) %>%
#   rename(!!Vim_name := PositiveVimentin)

counts_save <- counts_save %>% 
  rename(!!SMA_name := PositiveSMA) %>%
  rename(!!PDGFR_name := PositivePDGFR)

# ---- Count positive cells by treatment group, pooling all samples ---

# pooled_counts <- objectData %>% 
#   group_by(Metadata_Treatment) %>% 
#   summarise(nCells = n(),
#             PositiveSMA = sum(Intensity_MeanIntensity_SMA > threshSMA),
#             PercentPositiveSMA = PositiveSMA/nCells,
#             PositivePDGFR = sum(Intensity_MeanIntensity_PDGFRalpha > threshPDGFR),
#             PercentPositivePDGFR = PositivePDGFR/nCells,
#             PositiveVimentin = sum(Intensity_MeanIntensity_Vimentin > threshVim),
#             PercentPositiveVimentin = PositiveVimentin/nCells)

pooled_counts <- objectData %>% 
  group_by(Metadata_Treatment) %>% 
  summarise(nCells = n(),
            PositiveSMA = sum(Intensity_MeanIntensity_SMA > threshSMA),
            PercentPositiveSMA = PositiveSMA/nCells,
            PositivePDGFR = sum(Intensity_MeanIntensity_PDGFRalpha > threshPDGFR),
            PercentPositivePDGFR = PositivePDGFR/nCells)
            
pooled_counts_save <- pooled_counts
# pooled_counts_save <- pooled_counts_save %>% 
#   rename(!!SMA_name := PositiveSMA) %>%
#   rename(!!PDGFR_name := PositivePDGFR) %>%
#   rename(!!Vim_name := PositiveVimentin)

pooled_counts_save <- pooled_counts_save %>% 
  rename(!!SMA_name := PositiveSMA) %>%
  rename(!!PDGFR_name := PositivePDGFR)

# ---- Save new file ----
objectName <- str_sub(basename(objectFile), 1, -5) # name of the file without higher levels or extension
parentDir <- dirname(objectFile) # parent of the file
outputFile = paste(objectName, "_count_summary.csv") # spaces will be inserted
write_csv(counts_save,file.path(parentDir, outputFile))
outputPooledFile = paste(objectName, "_pooled_count_summary.csv") # spaces will be inserted
write_csv(pooled_counts_save,file.path(parentDir, outputPooledFile))

# ---- Plots ----
p_countSMA <- ggplot(counts,
                 aes(Metadata_Treatment, PercentPositiveSMA)) + 
  geom_col() +
  facet_wrap(~Metadata_Sample) +
  theme(text=element_text(size=20))


#  theme(text=element_text(size=20))

outputPlotSMA = paste(objectName, "SMA count plot.pdf")

# plot size can be changed according to the number of panels
ggsave(file.path(parentDir, outputPlotSMA), width=20, height = 14)

p_countPDGFR <- ggplot(counts,
                       aes(Metadata_Treatment, PercentPositivePDGFR)) + 
  geom_col() +
  facet_wrap(~Metadata_Sample)+
  theme(text=element_text(size=20))

outputPlotPDGFR = paste(objectName, "PDGFR count plot.pdf")
ggsave(file.path(parentDir, outputPlotPDGFR), width=20, height = 14)

# template for plotting vimentin positivity -- but this measurement is not appropriate if cells are detected using vimentin
# 
# p_countVim <- ggplot(counts,
#                        aes(Metadata_Treatment, PercentPositiveVimentin)) + 
#   geom_col() +
#   facet_wrap(~Metadata_Sample)+
#   theme(text=element_text(size=20))
# 
# outputPlotVim = paste(objectName, "Vimentin count plot.pdf")
# ggsave(file.path(parentDir, outputPlotVim), width=20, height = 14)

# to get both on the same plot we have to reshape the data
# pct_pos <- counts %>% 
#   select(-c(PositiveSMA, PositivePDGFR, PositiveVimentin)) %>%
#   rename(SMA = PercentPositiveSMA, PDGFR = PercentPositivePDGFR, Vim = PercentPositiveVimentin)

pct_pos <- counts %>% 
  select(-c(PositiveSMA, PositivePDGFR)) %>%
  rename(SMA = PercentPositiveSMA, PDGFR = PercentPositivePDGFR)

pct_pos_pivot <- pct_pos %>%
  pivot_longer(-c(Metadata_Sample,Metadata_Treatment, nCells), names_to="Marker", values_to="PercentPositive")

# p_counts <- ggplot(pct_pos_pivot,
#                        aes(x = Metadata_Treatment, y = PercentPositive, fill = Marker)) + 
#   geom_col(position = "dodge") +
#   facet_wrap(~Metadata_Sample)+
#   theme(text=element_text(size=20)) +
#   scale_fill_manual(values = c("#1ABC9C", "#FF5733"))

p_counts <- ggplot(pct_pos_pivot,
                   aes(x = Marker, y = PercentPositive, fill = Metadata_Treatment)) + 
  geom_col(position = "dodge") +
  facet_wrap(~Metadata_Sample)+
  theme(text=element_text(size=20))  +
  scale_fill_manual(values = c("#000000", "#FF0000"))

outputCountPlot = paste(objectName, "combined count plot.pdf")
ggsave(file.path(parentDir, outputCountPlot), width=20, height = 14)

# ---- Plots of pooled data ----
p_pooledSMA <- ggplot(pooled_counts,
                      aes(Metadata_Treatment, PercentPositiveSMA)) + 
  geom_col()+
  theme(text=element_text(size=20))

outputPooledSMA = paste(objectName, "SMA pooled count plot.pdf")

# plot size can be changed according to the number of panels
ggsave(file.path(parentDir, outputPooledSMA), width=20, height = 14)

p_pooledPDGFR <- ggplot(pooled_counts,
                        aes(Metadata_Treatment, PercentPositivePDGFR)) + 
  geom_col()+
  theme(text=element_text(size=20))

outputPooledPDGFR = paste(objectName, "PDGFR pooled count plot.pdf")
ggsave(file.path(parentDir, outputPooledPDGFR), width=20, height = 14)

# template for plotting vimentin positivity -- but this measurement is not appropriate if cells are detected using vimentin
# 
# p_pooledVim <- ggplot(pooled_counts,
#                         aes(Metadata_Treatment, PercentPositiveVim)) + 
#   geom_col()+
#   theme(text=element_text(size=20))
# 
# outputPooledVim = paste(objectName, "Vimentin pooled count plot.pdf")
# ggsave(file.path(parentDir, outputPooledVim), width=20, height = 14)
