# summarize_CP_tables_SD_plot.R
# R script to select mean intensity data from CellProfiler object results files
# and summarize by treatment and sample
# To use: Run the script.
# Will be prompted for a file

# ---- Setup ----
require(tidyverse)
require(readr)
require(stringr)


# ---- Prompt for an object file ----
# no message will be displayed
objectFile <- file.choose()

# Read the data from the file
objectData <- read_csv(objectFile,
                    locale = locale())


mean_intens <- objectData %>% 
  group_by(Metadata_Treatment, Metadata_Sample) %>% 
  summarise(nCells = n(),
            MeanSMA = mean(Intensity_MeanIntensity_SMA),
            MeanSMA_SD = sd(Intensity_MeanIntensity_SMA),
            MeanPDGFR = mean(Intensity_MeanIntensity_PDGFRalpha),
            MeanPDGFR_SD = sd(Intensity_MeanIntensity_PDGFRalpha))


mean_pooled_intens <- objectData %>% 
  group_by(Metadata_Treatment) %>% 
  summarise(nCells = n(),
            MeanSMA = mean(Intensity_MeanIntensity_SMA),
            MeanSMA_SD = sd(Intensity_MeanIntensity_SMA),
            MeanPDGFR = mean(Intensity_MeanIntensity_PDGFRalpha),
            MeanPDGFR_SD = sd(Intensity_MeanIntensity_PDGFRalpha))

# ---- Save new file ----
objectName <- str_sub(basename(objectFile), 1, -5) # name of the file without higher levels or extension
parentDir <- dirname(objectFile) # parent of the file
outputFile = paste(objectName, "_summary.csv") # spaces will be inserted
write_csv(mean_intens,file.path(parentDir, outputFile))
outputPooledFile = paste(objectName, "_pooledsummary.csv") # spaces will be inserted
write_csv(mean_pooled_intens,file.path(parentDir, outputPooledFile))

# ---- Plots ----
p_meanSMA <- ggplot(objectData,
                 aes(x=Metadata_Treatment,
                     y=Intensity_MeanIntensity_SMA)) + 
  geom_violin(trim=FALSE) +
  facet_wrap(~Metadata_Sample) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red")

outputPlotSMA = paste(objectName, "SMA mean plot.pdf")

# plot size can be changed according to the number of panels
ggsave(file.path(parentDir, outputPlotSMA), width=20, height = 14)

p_meanPDGFR <- ggplot(objectData,
                    aes(x=Metadata_Treatment,
                        y=Intensity_MeanIntensity_PDGFRalpha)) + 
  geom_violin(trim=FALSE) +
  facet_wrap(~Metadata_Sample) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red")

outputPlotPDGFR = paste(objectName, "PDGFR mean plot.pdf")
ggsave(file.path(parentDir, outputPlotPDGFR), width=20, height = 14)


# ---- Plots of pooled data ----
p_pooledSMA <- ggplot(objectData,
                    aes(x=Metadata_Treatment,
                        y=Intensity_MeanIntensity_SMA)) + 
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red")

outputPooledSMA = paste(objectName, "SMA mean pooled plot.pdf")

# plot size can be changed according to the number of panels
ggsave(file.path(parentDir, outputPooledSMA), width=20, height = 14)

p_pooledPDGFR <- ggplot(objectData,
                      aes(x=Metadata_Treatment,
                          y=Intensity_MeanIntensity_PDGFRalpha)) + 
  geom_violin(trim=FALSE) +
  stat_summary(fun.data=mean_sdl,
               geom="pointrange", color="red")

outputPooledPDGFR = paste(objectName, "PDGFR mean pooled plot.pdf")
ggsave(file.path(parentDir, outputPooledPDGFR), width=20, height = 14)
