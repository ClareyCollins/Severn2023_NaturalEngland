# ---
# Title: "Severn Phase III - Species Diversity Bubble Plot"
# Author: "Clare E. Collins, Graham Sellers and Nathan P. Griffiths"
# Date: "31 January 2024"
# ---
# 
# ## Prepare working environment
# 
# Clear R memory <<rm(list=ls())>>, set working directory (check <<getwd()>>; Set <<setwd()>>) and load required packages.
# 

##Packages------

library(plyr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(patchwork)
library(cowplot)


# To ensure reproducibility, print details about the version of R being used for analysis.
sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)

#### Start with 12S verts data ####

## Import raw Tapirs outputs NE data
NE2023 = read.csv('data/NE23_mifish_blast98_denoise.csv', sep = ',',  header=TRUE)

## Make '.' a '-'
colnames(NE2023) = gsub('\\.', '-', colnames(NE2023))

## Rename first column
colnames(NE2023)[1] = "assignment"

## Remove underscore from species names
NE2023$assignment = gsub("_", " ", NE2023$assignment)

##### Controls #####

## Create separate dataframe of controls
NEcontrols2023 = NE2023[grep('assignment|NEG|POS|CSNeg', colnames(NE2023))]

## Remove species that do not feature in any of the controls
NEcontrols2023_species = NEcontrols2023[rowSums(NEcontrols2023[(2:10)]) > 0,]

## Turn table into long format
NEcontrols2023_species = cbind(NEcontrols2023_species[1], stack(NEcontrols2023_species[-1]))

## Change column names
colnames(NEcontrols2023_species)[2:3] = c("Reads","Controls") 
## Add a column labelling all values in this run as "Severn Verts III"
NEcontrols2023_species$Run = "Severn Verts III"

## Order data by taxonomic assignment
NEcontrols2023_species = NEcontrols2023_species[order(NEcontrols2023_species$assignment),]

## Create column specifying type of negative control
NEcontrols2023_species$Type = ifelse(grepl("CSNeg", NEcontrols2023_species$Controls), "Extraction",
                                     ifelse(grepl("NEG", NEcontrols2023_species$Controls), "Neg",
                                            "Pos"))

## Create factor to order heatmap by
NEcontrols2023_species$fType = factor(NEcontrols2023_species$Type, 
                               levels=c("Extraction",
                                        "Neg", "Pos"))

## Plot contamination found in negative controls
options(scipen=999)
plot_control = ggplot(NEcontrols2023_species, aes(x=Controls, 
                                    y=fct_rev(as_factor(assignment)), 
                                    fill=Reads))
plot_control = plot_control + geom_tile(colour="black")
plot_control = plot_control + scale_fill_gradientn(name="Total read counts", 
                                limits=c(0,100000),
                                breaks=c(0,50000, 100000),
                                colours=c("white","grey","black"), 
                                values=c(0,0.1,1))
plot_control = plot_control + labs(x="Process controls", y="Taxonomic assignment")
plot_control = plot_control + theme_bw()
plot_control = plot_control + theme(panel.grid.major = element_line(colour="white"),
                 panel.grid.minor = element_line(colour="white"), 
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(colour="black", face = "italic"),
                 text = element_text(size=12),
                 legend.key.size = unit(1, 'lines'))
plot_control = plot_control + facet_grid(Run ~ fType, scale = "free", space = "free")
plot_control

png("output/Severn_III_NE2023_Unfiltered_Contamination_Verts.png", width = 5000, height = 2500, units = "px", res = 345)
plot_control
dev.off()

##### Samples #####

## Make assignment column the row names for these transposed data
rownames(NE2023) = NE2023$assignment
NE2023 = NE2023[, -grep('assignment', colnames(NE2023))]

## Sort by taxonomy
NE2023 = NE2023[order(NE2023$taxonomy),]

## isolate fish and others to continue working on fish in NE2023
NE2023_nonfish = NE2023[-grep('Actinopteri|Petromyzontiformes', NE2023$taxonomy),]
NE2023 = NE2023[grep('Actinopteri|Petromyzontiformes', NE2023$taxonomy),]

## Grab taxonomy, before removing column below
NE2023_tax = NE2023$taxonomy

## Remove negatives, positives, filtration and extraction blanks, and taxonomy (leaves 163 variables)
NE2023 = NE2023[, -grep('taxonomy|NEG|POS|CSNeg', colnames(NE2023))]

## Check positive doesn't appear in samples and if 0 remove positive control (Maylandia zebra) if does appear, export this data to warning file
NE2023_POSITIVEWARNING = NULL
if (all(NE2023['Maylandia zebra',] <= 0)) { 
  ## Remove the positive control row
  NE2023 = NE2023[!(rownames(NE2023) == 'Maylandia zebra'),]
} else {
  NE2023_POSITIVEWARNING = subset(NE2023, rownames(NE2023) == 'Maylandia zebra')
}

write.csv(NE2023_POSITIVEWARNING, "output/POSITIVE_WARNING.csv")

## Now remove positive row if still there
NE2023 = NE2023[!(rownames(NE2023) == 'Maylandia zebra'),]

## Noise filter data at 0.1% (but record data removed in a df NE2023_noisefiltered):
## This function transposes (t) df (swapping columns and rownames) and divides each element transposed in NE2023 by the corresponding column sum, normalising the columns of the matrix and then checks whether each normalised element is less (or more) than 0.001 returning logical (TRUE/FALSE), the true values are returned as 0 and the false values are left as before in the corresponding df
NE2023_noisefiltered = NE2023
NE2023_noisefiltered[t(t(NE2023)/colSums(NE2023)) > 0.001] = 0
NE2023[t(t(NE2023)/colSums(NE2023)) < 0.001] = 0
write.csv(NE2023_noisefiltered, "output/NE2023_noisefiltered.csv")

## Drop (but record - species in these dropped samples only) samples with less than 500 total (i.e. sum across all species) reads:
NE2023_samplesdropped = NE2023[!colSums(NE2023) > 500]
NE2023_samplesdropped = NE2023_samplesdropped[rowSums(NE2023_samplesdropped[])>0,]
write.csv(NE2023_samplesdropped,"output/NE2023_droppedless500reads.csv")

NE2023 = NE2023[!colSums(NE2023) < 500]

## Remove assignments to Cyprinidae, Percidae and unassigned
NE2023 = NE2023[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(NE2023)),]

## Change Alosa_alosa to Alosa spp.
rownames(NE2023) = sub('Alosa alosa', 'Alosa spp.', rownames(NE2023))

## Merge Leuciscus to Leuciscus leuciscus and change to Leuciscus spp.
Leuciscus_leuciscus = NE2023['Leuciscus leuciscus',]
Leuciscus = NE2023['Leuciscus',]
NE2023['Leuciscus leuciscus',] = Leuciscus_leuciscus + Leuciscus
NE2023 = NE2023[!(rownames(NE2023) == 'Leuciscus'),]
rownames(NE2023) = sub('Leuciscus leuciscus', 'Leuciscus spp.', rownames(NE2023))

## Change Lampetra fluviatilis to Lampetra_spp.
rownames(NE2023) = sub('Lampetra fluviatilis', 'Lampetra spp.', rownames(NE2023))

## Remove (but record) species with no assignments
NE2023_speciesnoassignments = NE2023[rowSums(NE2023) == 0,]
write.csv(NE2023_speciesnoassignments, "output/NE2023_speciesnoassignments.csv")
NE2023 = NE2023[rowSums(NE2023) > 0,]

## Change column names from 'tube ids' to 'sample ids' from dictionary file
sample_names = read.csv('data/NE2023_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube ## Gives the sample the name of the tube so that this can be matched and replaced with the sample name below
colnames(NE2023) = samples[colnames(NE2023)]

## Sort by Genus
NE2023 = NE2023[order(row.names(NE2023)),]

##### Bubble Plot #####

## Transpose and sort (create new df, copy rownames to a column, move that column to the start and remove rownames, then transpose df)
bubble_df = NE2023
bubble_df$assignment = rownames(bubble_df)
bubble_df = bubble_df[,c(58,1:57)]
rownames(bubble_df) = NULL
bubble_df = data.frame(t(bubble_df))

## Assign assignment row as column names and delete assignment row
colnames(bubble_df) = bubble_df['assignment',]
bubble_df = bubble_df[!(rownames(bubble_df) == 'assignment'),]

## Assign rownames column as site column, move it to the start and delete the rownames
bubble_df$site = rownames(bubble_df)
bubble_df = bubble_df [,c(26,1:25)]
rownames(bubble_df) = NULL

## Remove hyphen and sample number from site ID
bubble_df$site = gsub(".{2}$", "", bubble_df$site)

## Convert characters to numeric for all columns except the first
bubble_df[,-1] = lapply(bubble_df[,-1], function(x) as.numeric(as.character(x)))

## Merge (sum) reads by site (will be site per date e.g. Site01-A) PACKAGE plyr
bubble_df = ddply(bubble_df, .(site), numcolwise(sum))

## Add assigned reads column (get dimensions (rows and columns) of the df; print the number of columns and use this to set the larger number needed to add reads after rowSums)
r=dim(bubble_df)
r[2]
bubble_df$assigned = rowSums(bubble_df[2:26])

## Add site as rownames (and then remove site column)
rownames(bubble_df) = bubble_df$site
bubble_df = bubble_df[, -grep('site', colnames(bubble_df))]

## Make proportional reads eDNA dataset 
site = rownames(bubble_df)
prop_reads_df = bubble_df
prop_reads_df = prop_reads_df/prop_reads_df$assigned
r = dim(prop_reads_df)
prop_reads_df = prop_reads_df[1:(r[2]-1)]
prop_reads_df = cbind(site, prop_reads_df) 
dim(prop_reads_df)
prop_reads_df = melt(prop_reads_df, id=(1)) ## melt from reshape PACKAGE

## Change variable to taxonomic
colnames(prop_reads_df) = sub('variable', 'taxonomic', colnames(prop_reads_df))

## Make bubble plot
### Make Bubble Plot Figure 
set1 = subset(prop_reads_df, prop_reads_df$value>0)
set1$site = factor(set1$site, levels = site)

## Add common names, importing them from dictionary file, matching to taxonomic names and then inserting into this df
common_names_df = read.csv('data/NE2023_CommonNames.txt', sep = '\t', header = F, col.names = c('taxonomic', 'common'))
taxonomic_common = match(set1$taxonomic, common_names_df$taxonomic)
set1$common = common_names_df$common[taxonomic_common]

## Organise sampling events 
set1$event = set1$site
set1 = set1[,c(1, 5, 4, 3, 2)]

## Remove sample event from site ID 
set1$site = gsub(".{2}$", "", set1$site)

## Check species richness 
fr = table(set1$event)
fr = as.data.frame(fr)
## re-order fr
x = fr[c(1,3,5,7,9,11,13),]
z = fr[c(2,4,6,8,10,12,14),]
y = rbind(x, z)


## Remove sample ID from event
set1$event = gsub("^.....", "", set1$event)

## Ggplot bubble plot

event.labs = c("9-10 May 2023", "22 May 2023", "6 June 2023", "19 June 2023")
names(event.labs) = c("A", "B", "C", "D")

## Taxonomic name plot
bubble_plot = ggplot (set1, aes(x=site, y=taxonomic, size=value*100)) +
  geom_point(alpha=0.9, shape=21, fill="#70cf57", color="black") +
  scale_size(range = c(1, 10), name="Relative fish reads (%)") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Sampling Site") +
  ylab("") +
  theme(axis.text.y = element_text(face = "italic", size = 12))+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_x_discrete(position = "top", labels = c("00", "01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(set1$taxonomic))) + 
  facet_wrap(~event, ncol=4, labeller = labeller(event = event.labs))


png("output/Fish2023_Summary_taxonomic.png", width = 4500, height = 2500, units = "px", res = 345)
bubble_plot
dev.off()

## Common name bubbleplot with 00 site not showing for sampling days A, C and D

# Specify the order of levels for the 'common' variable
common_order <- c(
  "Atlantic salmon", "Brown trout", "Chub", "Common barbel", "Common bleak",
  "Common bream", "Common minnow", "Common roach", "Dace spp.", "European bullhead", "European eel", "European flounder",
  "European perch", "Grass carp", "Grayling", "Gudgeon",
  "Northern pike", "Rudd", "Ruffe", "Shad spp.", "Silver bream",
  "Stone loach", "Tench", "Three-spined stickleback", "Zander"
)
set1$common <- factor(set1$common, levels = common_order)

# Create plots for each event
bubble_plot_A = ggplot(subset(set1, event == "A"), aes(x = site, y = common, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#70cf57", color = "black") +
  scale_size(range = c(1, 10), name = "Relative fish reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(set1$common))
    ) +
  labs (x = "Sample sites")


bubble_plot_B = ggplot(subset(set1, event == "B"), aes(x = site, y = common, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#70cf57", color = "black") +
  scale_size(range = c(1, 10), name = "Relative fish reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("00", "01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(set1$common))
  ) + 
  labs (x = "Sample sites")

bubble_plot_C = ggplot(subset(set1, event == "C"), aes(x = site, y = common, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#70cf57", color = "black") +
  scale_size(range = c(1, 10), name = "Relative fish reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(set1$common))
  ) +
  labs (x = "Sample sites")

bubble_plot_D = ggplot(subset(set1, event == "D"), aes(x = site, y = common, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#70cf57", color = "black") +
  scale_size(range = c(1, 10), name = "Relative fish reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(set1$common))
  ) +
  labs (x = "Sample sites")

# Create a combined legend
combined_legend = cowplot::get_legend(bubble_plot_A + theme(legend.position="right"))  # Using plot A as a reference

# Arrange plots using patchwork
combined_plots = ((bubble_plot_A +
                     facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
  (bubble_plot_B +
     facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
  (bubble_plot_C +
  facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs)))  /
  (bubble_plot_D +
  facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
  combined_legend) +
  plot_layout(ncol = 5, widths = c(1, 1.2, 1, 1, 0.8))

combined_plots
# Save the combined plot
png("output/Fish2023_Summary_common.png", width = 4500, height = 2500, units = "px", res = 345)
combined_plots
dev.off()


### Invertebrate Data ####
## Import raw Tapirs outputs
NEInverts2023 = read.csv("data/NE23_leese_blast95_denoise.csv", header=TRUE)

## Grab taxonomy, before removing column below
NEInverts2023_tax = NEInverts2023$taxonomy

## Make '.' a '-'
colnames(NEInverts2023) = gsub('\\.', '-', colnames(NEInverts2023))

## Rename first column
colnames(NEInverts2023)[1] = "assignment"

## Remove underscore from species names
NEInverts2023$assignment = gsub("_", " ", NEInverts2023$assignment)


##### Inverts controls #####

## Create separate dataframe of controls
NEInvertscontrols2023 = NEInverts2023[grep('assignment|NEG|POS|CSNeg', colnames(NEInverts2023))]

## Remove species that do not feature in any of the controls
NEInvertscontrols2023_species = NEInvertscontrols2023[rowSums(NEInvertscontrols2023[(2:10)]) > 0,]

## Turn table into long format
NEInvertscontrols2023_species = cbind(NEInvertscontrols2023_species[1], stack(NEInvertscontrols2023_species[-1]))

## Change column names
colnames(NEInvertscontrols2023_species)[2:3] = c("reads","controls") 
## Add a column labelling all values in this run as "Severn Inverts III"
NEInvertscontrols2023_species$run = "Severn Inverts III"

## Order data by taxonomic assignment
# NEInvertscontrols2023_species = NEInvertscontrols2023_species[order(NEInvertscontrols2023_species$assignment),]

## Create column specifying type of negative control
NEInvertscontrols2023_species$type = ifelse(grepl("CSNeg", NEInvertscontrols2023_species$controls), "Extraction",
                                     ifelse(grepl("NEG", NEInvertscontrols2023_species$controls), "Neg",
                                            "Pos"))

## Create factor to order heatmap by
NEInvertscontrols2023_species$ftype = factor(NEInvertscontrols2023_species$type, 
                                      levels=c("Extraction",
                                               "Neg", "Pos"))

## Plot contamination found in negative controls
options(scipen=999)
plot_contr_inverts = ggplot(NEInvertscontrols2023_species, aes(x=controls, 
                                                  y=fct_rev(as_factor(assignment)), 
                                                  fill=reads))
plot_contr_inverts = plot_contr_inverts + geom_tile(colour="black")
plot_contr_inverts = plot_contr_inverts + scale_fill_gradientn(name="Total read counts", 
                                                   limits=c(0,100000),
                                                   breaks=c(0,50000, 100000),
                                                   colours=c("white","grey","black"), 
                                                   values=c(0,0.1,1))
plot_contr_inverts = plot_contr_inverts + labs(x="Process controls", y="Taxonomic assignment")
plot_contr_inverts = plot_contr_inverts + theme_bw()
plot_contr_inverts = plot_contr_inverts + theme(panel.grid.major = element_line(colour="white"),
                                    panel.grid.minor = element_line(colour="white"), 
                                    axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                                    axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                                    axis.ticks.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.text.y = element_text(colour="black", face = "italic"),
                                    text = element_text(size=12),
                                    legend.key.size = unit(1, 'lines'))
plot_contr_inverts = plot_contr_inverts + facet_grid(run ~ ftype, scale = "free", space = "free")
plot_contr_inverts

png("output/Severn_III_NE2023_Unfiltered_Contamination_Inverts.png", width = 5000, height = 2500, units = "px", res = 345)
plot_contr_inverts
dev.off()

###Inverts samples####
## Make assignment column the row names for these transposed data
rownames(NEInverts2023) = NEInverts2023$assignment
NEInverts2023 = NEInverts2023[, -grep('assignment', colnames(NEInverts2023))]

## Remove negatives, positives, filtration and extraction blanks, and taxonomy (leaves 163 variables)
NEInverts2023 = NEInverts2023[, -grep('taxonomy|NEG|POS|CSNeg', colnames(NEInverts2023))]

## Check positive doesn't appear in samples and if 0 remove positive control if does appear, export this data to warning file
NEInverts2023_POSITIVEWARNING = NULL
if (all(NEInverts2023['Calliphora vicina',] <= 0)) { 
  ## Remove the positive control row
  NEInverts2023 = NEInverts2023[!(rownames(NEInverts2023) == 'Calliphora vicina'),]
} else {
  NEInverts2023_POSITIVEWARNING = subset(NEInverts2023, rownames(NEInverts2023) == 'Calliphora vicina')
}

write.csv(NEInverts2023_POSITIVEWARNING, "output/POSITIVE_WARNING.csv")

## Now remove positive row if still there
NEInverts2023 = NEInverts2023[!(rownames(NEInverts2023) == 'Calliphora vicina'),]

## Remove any reads under 20 in any cell 
NEInverts2023[NEInverts2023 < 20] = 0

## Remove (but record) species with no assignments
NEInverts2023_speciesnoassignments = NEInverts2023[rowSums(NEInverts2023) == 0,]
write.csv(NEInverts2023_speciesnoassignments, "output/NEInverts2023_speciesnoassignments.csv")
NEInverts2023 = NEInverts2023[rowSums(NEInverts2023) > 0,]

## Drop (but record - species in these dropped samples only) samples with less than 500 total (i.e. sum across all species) reads:
NEInverts2023_samplesdropped = NEInverts2023[!colSums(NEInverts2023) > 500]
NEInverts2023_samplesdropped = NEInverts2023_samplesdropped[rowSums(NEInverts2023_samplesdropped[])>0,]
write.csv(NEInverts2023_samplesdropped,"output/NEInverts2023_droppedless500reads.csv")
NEInverts2023 = NEInverts2023[!colSums(NEInverts2023) < 500]


## Change column names from 'tube ids' to 'sample ids' from dictionary file
sample_names = read.csv('data/NE2023_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube ## Gives the sample the name of the tube so that this can be matched and replaced with the sample name below
colnames(NEInverts2023) = samples[colnames(NEInverts2023)]

## Add taxonomic information in to get order
## Split the strings in NEInverts2023_tax
NEInverts2023_taxonomicdf = data.frame(taxonomic = c(NEInverts2023_tax))
NEInverts2023_taxonomicdf$order = NEInverts2023_taxonomicdf
## Extract the "o__text" part from each element of NEInverts2023_taxonomic_df
NEInverts2023_taxonomicdf$order = sapply(NEInverts2023_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^o__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^o__", "", order_part)
  } else {
    order_text = NA  # If "o__" part not found, return NA
  }
  return(order_text)
})

## Extract the s__ part to get species column
NEInverts2023_taxonomicdf$species = sapply(NEInverts2023_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^s__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^s__", "", order_part)
  } else {
    order_text = NA  # If "s__" part not found, return NA
  }
  return(order_text)
})

## Remove underscore from species column
NEInverts2023_taxonomicdf$species = gsub("_", " ", NEInverts2023_taxonomicdf$species)

## Match species in NEInverts2023 to those in taxonomicdf and return the correct order in an order column
NEInverts2023 = merge(NEInverts2023, NEInverts2023_taxonomicdf, by.x = "row.names", by.y = "species", all.x = TRUE)

## Change rownames column heading to assignment and move order to the front and remove the taxonomic column
colnames(NEInverts2023) [1] = c("assignment")
NEInverts2023 = NEInverts2023[,c(66,1:65)]
NEInverts2023 = NEInverts2023[!(colnames(NEInverts2023) == 'taxonomic')]

## Add order to Ampullaceana balthica from NBN Atlas and those recorded to genus level only
NEInverts2023[4, "order"] = "Hygrophila"
NEInverts2023[36, "order"] = "Trichoptera"
NEInverts2023[39, "order"] = "Diptera"
NEInverts2023[71, "order"] = "Diptera"

## Remove all rows where order is unknown or NA (includes family assignments and unassigned)
NEInverts2023 = NEInverts2023[!is.na(NEInverts2023$order) & NEInverts2023$order != "" & NEInverts2023$order != "unknown", ]

## Barplot Species detected per order
invertcount = table(NEInverts2023$order) 
invertcountdf = as.data.frame(invertcount)
as.factor(invertcountdf$Var1)
sum(invertcountdf$Freq)

NEInvert2023_barplot = ggplot(data = invertcountdf, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black", fill="#78004F") + coord_flip() +
  theme_bw() + geom_text(aes(label=Freq), hjust=1.6, color="white", size=4.5) +
  xlab("") + ylab("Species detected per order") +  
  theme(axis.text.y = element_text(face = "italic", size = 12))+
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(face = "bold", size = 12)) + 
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
  scale_x_discrete(limits=rev(levels(invertcountdf$Var1))) + ggtitle("River Severn 2023")


png("output/Severn2023_invert_barplot.png", width = 2500, height = 2500, units = "px", res = 345)
NEInvert2023_barplot
dev.off()

############# Bubble Plot inverts 2023 ################

## Transpose and sort (create new df, then remove assignment column, combine order rows into the same row (summing data) and then transpose df (no rownames so removed lines from NG code))
inverts_bp_df = NEInverts2023
inverts_bp_df = inverts_bp_df[!(colnames(inverts_bp_df) == 'assignment')]
inverts_bp_df = aggregate(. ~ order, data = inverts_bp_df, sum)
inverts_bp_df = data.frame(t(inverts_bp_df))

## Assign order row as column names and delete order row
colnames(inverts_bp_df) = inverts_bp_df['order',]
inverts_bp_df = inverts_bp_df[!(rownames(inverts_bp_df) == 'order'),]

## Assign rownames column as site column, move it to the start and delete the rownames
inverts_bp_df$site = rownames(inverts_bp_df)
inverts_bp_df = inverts_bp_df [,c(18,1:17)]
rownames(inverts_bp_df) = NULL

## Convert characters to numeric for all columns except the first
inverts_bp_df[,-1] = lapply(inverts_bp_df[,-1], function(x) as.numeric(as.character(x)))

## Remove sample number (1,2,3) and merge (sum) reads by site (will be site per date e.g. Site01-A) PACKAGE plyr
inverts_bp_df$site = gsub(".{2}$", "", inverts_bp_df$site)
inverts_bp_df = ddply(inverts_bp_df, .(site), numcolwise(sum))

## Add assigned reads column (get dimensions (rows and columns) of the df; print the number of columns and use this to set the larger number needed to add reads after rowSums)
r=dim(inverts_bp_df)
r[2]
inverts_bp_df$assigned = rowSums(inverts_bp_df[2:18])

## Make proportional reads eDNA dataset 
site = inverts_bp_df$site
rownames(inverts_bp_df) = inverts_bp_df$site
inverts_bp_df = inverts_bp_df[, -grep('site', colnames(inverts_bp_df))]
prop_reads_df = inverts_bp_df
prop_reads_df = prop_reads_df/prop_reads_df$assigned
r = dim(prop_reads_df)
prop_reads_df = prop_reads_df[1:(r[2]-1)]
prop_reads_df = cbind(site, prop_reads_df) 
dim(prop_reads_df)
prop_reads_df = melt(prop_reads_df, id=(1)) ## melt from reshape PACKAGE

## Change variable to order
colnames(prop_reads_df) = sub('variable', 'order', colnames(prop_reads_df))

## Make bubble plot
### Make Bubble Plot Figure 
inverts_set1 = subset(prop_reads_df, prop_reads_df$value>0)
inverts_set1$site = factor(inverts_set1$site, levels = site)

## Organise sampling events 
inverts_set1$event = inverts_set1$site
inverts_set1 = inverts_set1[,c(1, 4, 2, 3)]

## Remove sample event from site ID 
inverts_set1$site = gsub(".{2}$", "", inverts_set1$site)

## Check species richness 
inverts_fr = table(inverts_set1$event)
inverts_fr = as.data.frame(inverts_fr)
## re-order fr WHY?
x = inverts_fr[c(1,3,5,7,9,11,13),]
z = inverts_fr[c(2,4,6,8,10,12,14),]
y = rbind(x, z)

## Remove sample ID from event
inverts_set1$event = gsub("^.....", "", inverts_set1$event)

## Taxonomic name plot
inverts_bubble_plot = ggplot (inverts_set1, aes(x=site, y=order, size=value*100)) +
  geom_point(alpha=0.9, shape=21, fill="#78004F", color="black") +
  scale_size(range = c(1, 10), name="Relative invert reads (%)") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Sampling Site Code") +
  ylab("") +
  theme(axis.text.y = element_text(face = "italic", size = 12))+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 10))+
  scale_x_discrete(position = "top", labels = c("00", "01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(inverts_set1$order))) + 
  facet_wrap(~event, ncol=4, labeller = labeller(event = event.labs))


png("output/Inverts2023_Summary_taxonomic.png", width = 4500, height = 2500, units = "px", res = 345)
inverts_bubble_plot
dev.off()

## Bubbleplot with 00 site not showing for sampling days A, C and D

# Create plots for each event
inverts_bubble_plot_A = ggplot(subset(inverts_set1, event == "A"), aes(x = site, y = order, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#78004F", color = "black") +
  scale_size(range = c(1, 10), name = "Relative invert reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(inverts_set1$order))) +
  labs(x = "Sample Sites")  # Modify x-axis title

inverts_bubble_plot_B = ggplot(subset(inverts_set1, event == "B"), aes(x = site, y =order, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#78004F", color = "black") +
  scale_size(range = c(1, 10), name = "Relative invert reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("00", "01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(inverts_set1$order))) +
  labs(x = "Sample Sites")  # Modify x-axis title

inverts_bubble_plot_C = ggplot(subset(inverts_set1, event == "C"), aes(x = site, y = order, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#78004f", color = "black") +
  scale_size(range = c(1, 10), name = "Relative invert reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(inverts_set1$order))) +
  labs(x = "Sample Sites")  # Modify x-axis title

inverts_bubble_plot_D = ggplot(subset(inverts_set1, event == "D"), aes(x = site, y = order, size = value * 100)) +
  geom_point(alpha = 0.9, shape = 21, fill = "#78004f", color = "black") +
  scale_size(range = c(1, 10), name = "Relative invert reads (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_x_discrete(position = "top", labels = c("01", "02", "03", "04", "05")) +
  scale_y_discrete(limits = rev(levels(inverts_set1$order))) +
  labs(x = "Sample Sites")  # Modify x-axis title

# Create a combined legend
inverts_combined_legend = cowplot::get_legend(inverts_bubble_plot_A + theme(legend.position="right"))  # Using plot A as a reference

# Arrange plots using patchwork
inverts_combined_plots = ((inverts_bubble_plot_A +
                     facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
                    (inverts_bubble_plot_B +
                       facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
                    (inverts_bubble_plot_C +
                       facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs)))  /
                    (inverts_bubble_plot_D +
                       facet_wrap(~event, ncol = 4, labeller = labeller(event = event.labs))) /
                    inverts_combined_legend) +
  plot_layout(ncol = 5, widths = c(1, 1.2, 1, 1, 0.8))

# Save the combined plot
png("output/Inverts2023_Summary_taxonomic.png", width = 4500, height = 2500, units = "px", res = 345)
inverts_combined_plots
dev.off()
