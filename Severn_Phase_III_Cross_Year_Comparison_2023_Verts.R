# ---
# Title: "Severn Phase III - Cross-year Comparison Species Diversity Bubble Plot - VERTS"
# Author: "Clare E. Collins, Graham Sellers and Nathan P. Griffiths"
# Date: "31 January 2024"
# ---
# 
# ## Prepare working environment
# 
# Clear R memory <<rm(list=ls())>>, set working directory (check <<getwd()>>; Set <<setwd()>>) and load required packages.
# 

##Packages------

library('reshape2')
library('vegan')
library(ggplot2)
library(stringr)

# 2018 VERTS DATA BIT --------------------------------------------------------------------------------

# read in the data
dat_2018 = read.csv('data/comparison/2018_reanalysis_MIFISH_verts_blast98_denoise.tsv', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

# grab taxonomy, remove from data
dat_2018_tax = dat_2018$taxonomy
dat_2018 = dat_2018[-grep('taxonomy', colnames(dat_2018))]

# noise filter data at 0.1%
dat_2018[t(t(dat_2018)/colSums(dat_2018)) < 0.001] = 0

# drop samples with less than 500 reads:
dat_2018 = dat_2018[,!colSums(dat_2018) < 500]

# add back taxonomy
dat_2018$taxonomy = dat_2018_tax

# remove assignments to Cyprinidae, Percidae and unassigned
dat_2018 = dat_2018[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(dat_2018)),]

# isolate fish and remove taxonomy (again)
dat_2018 = dat_2018[grep('Actinopteri|Petromyzontiformes', dat_2018$taxonomy),]
dat_2018 = dat_2018[-grep('taxonomy', colnames(dat_2018))]

# remove species with no assignments
dat_2018 = dat_2018[rowSums(dat_2018) > 0,]

# change Alosa_alosa to Alosa_spp.
rownames(dat_2018) = sub('Alosa_alosa', 'Alosa_spp.', rownames(dat_2018))

# change Leuciscus to Leuciscus_spp.
rownames(dat_2018) = sub('Leuciscus', 'Leuciscus_spp.', rownames(dat_2018))

# remove filtration and extraction blanks
dat_2018 = dat_2018[-grep('^N|N', colnames(dat_2018))]

# remove samples with no assignments
dat_2018 = dat_2018[colSums(dat_2018) > 0]

# select only River Severn sites
dat_2018 = dat_2018[-grep('TB|PB', colnames(dat_2018))]

# isolate site repeated between all 3 years (Diglis)
Diglis_2018 = dat_2018[grep('DW', colnames(dat_2018))]
Sabrina_2018 = dat_2018[grep('WB', colnames(dat_2018))]
Stourport_2018 = dat_2018[grep('L', colnames(dat_2018))]

# pool sites
pooled_2018 = data.frame(rowSums(Stourport_2018), rowSums(Sabrina_2018), rowSums(Diglis_2018))
pooled_2018 = pooled_2018[rowSums(pooled_2018) > 0,]
pooled_2018 = pooled_2018[sort(rownames(pooled_2018)),]
names(pooled_2018) = c('Stourport_2018', 'Sabrina_2018', 'Diglis_2018')


# 2021 VERTS DATA BIT--------------------------------------------------------------------------------

# read in the data
dat_2021 = read.csv('data/comparison/NE2021_MIFISH_fish_blast98_denoise.tsv', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

# remove negatives, positives, filtration and extraction blanks, and taxonomy
dat_2021 = dat_2021[, -grep('taxonomy|NEG|POS|EB|FB', colnames(dat_2021))]

# noise filter data at 0.1%
dat_2021[t(t(dat_2021)/colSums(dat_2021)) < 0.001] = 0

# drop samples with less than 500 reads:
dat_2021 = dat_2021[!colSums(dat_2021) < 500]

# remove positive control (Maylandia zebra)
dat_2021 = dat_2021[-grep('Maylandia_zebra', rownames(dat_2021)),]

# remove assignments to Cyprinidae, Percidae and unassigned
dat_2021 = dat_2021[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(dat_2021)),]

# change Alosa_alosa to Alosa_spp.
rownames(dat_2021) = sub('Alosa', 'Alosa_spp.', rownames(dat_2021))

# change Leuciscus to Leuciscus_spp.
rownames(dat_2021) = sub('Leuciscus', 'Leuciscus_spp.', rownames(dat_2021))

# change Lampetra fluviatilis to Lampetra_spp.
rownames(dat_2021) = sub('Lampetra_fluviatilis', 'Lampetra_spp.', rownames(dat_2021))

# make '.' a '-' as intended (why does read.csv do this to the original input file???)
colnames(dat_2021) = gsub('\\.', '-', colnames(dat_2021))

# correct samples with improved 1:10 dilution results
dilutes = c('NE23', 'NE26', 'NE89')

for(i in dilutes){dat_2021[i] = dat_2021[paste(i, '_2', sep = '')]}

# remove duplicated samples (NE63-2-1 and NE63-2-2)
dat_2021 = dat_2021[-grep('NE63-2-1|NE63-2-2', colnames(dat_2021))]

colnames(dat_2021)[which(names(dat_2021) == 'NE63-1-1')] = 'NE63'
colnames(dat_2021)[which(names(dat_2021) == 'NE63-1-2')] = 'NE65'

# change column names from 'tube ids' to 'sample ids'
sample_names = read.csv('data/comparison/NE2021_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube
colnames(dat_2021) = samples[colnames(dat_2021)]

# remove species with no assignments
dat_2021 = dat_2021[rowSums(dat_2021) > 0,]

# remove 1:10 dilutions and lab filtration
dat_2021 = dat_2021[, -grep('_2|_lab', colnames(dat_2021))]

# isolate sites repeated between all years (Stourport (S10), Sabrina (S11), Diglis (S13))
Stourport_2021 = dat_2021[grep('S10', colnames(dat_2021))]
Sabrina_2021 = dat_2021[grep('S11', colnames(dat_2021))]
Diglis_2021 = dat_2021[grep('S13', colnames(dat_2021))]

# pool sites
pooled_2021 = data.frame(rowSums(Stourport_2021), rowSums(Sabrina_2021), rowSums(Diglis_2021))
pooled_2021 = pooled_2021[rowSums(pooled_2021) > 0,]
pooled_2021 = pooled_2021[sort(rownames(pooled_2021)),]
names(pooled_2021) = c('Stourport_2021', 'Sabrina_2021', 'Diglis_2021')


# 2022 VERTS DATA BIT--------------------------------------------------------------------------------

# read in the data
dat_2022 = read.csv('data/comparison/FISH_98_Final.csv')

# make '.' a '-' as intended
colnames(dat_2022) = gsub('\\.', '-', colnames(dat_2022))

# grab taxonomy, remove from data
dat_2022_tax = dat_2022$taxonomy
dat_2022 = dat_2022[-grep('taxonomy', colnames(dat_2022))]

# sort row names
rownames(dat_2022) = dat_2022[,1]
dat_2022 <- dat_2022[, -1]

# noise filter data at 0.1%
dat_2022[t(t(dat_2022)/colSums(dat_2022)) < 0.001] = 0

# drop samples with less than 500 reads:
dat_2022 = dat_2022[,!colSums(dat_2022) < 500]

# add back taxonomy
dat_2022$taxonomy = dat_2022_tax

# remove assignments to Cyprinidae, Percidae and unassigned
dat_2022 = dat_2022[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(dat_2022)),]

# isolate fish and remove taxonomy (again)
dat_2022 = dat_2022[grep('Actinopteri|Petromyzontiformes', dat_2022$taxonomy),]
dat_2022 = dat_2022[-grep('taxonomy', colnames(dat_2022))]

# remove species with no assignments
dat_2022 = dat_2022[rowSums(dat_2022) > 0,]

# change Alosa_alosa to Alosa_spp.
rownames(dat_2022) = sub('Alosa', 'Alosa_spp.', rownames(dat_2022))

# change Leuciscus to Leuciscus_spp.
rownames(dat_2022) = sub('Leuciscus', 'Leuciscus_spp.', rownames(dat_2022))

# change Lampetra fluviatilis to Lampetra_spp.
rownames(dat_2022) = sub('Lampetra_fluviatilis', 'Lampetra_spp.', rownames(dat_2022))

# select only River Severn sites
dat_2022 <- dat_2022[grep('S.-', colnames(dat_2022))]

# isolate sites repeated between all years (Stourport (S5), Sabrina (S2), Diglis (S1))
Stourport_2022 = dat_2022[grep('S5', colnames(dat_2022))]
Sabrina_2022 = dat_2022[grep('S2', colnames(dat_2022))]
Diglis_2022 = dat_2022[grep('S1', colnames(dat_2022))]

# pool sites
pooled_2022 = data.frame(rowSums(Stourport_2022), rowSums(Sabrina_2022), rowSums(Diglis_2022))
pooled_2022 = pooled_2022[rowSums(pooled_2022) > 0,]
pooled_2022 = pooled_2022[sort(rownames(pooled_2022)),]
names(pooled_2022) = c('Stourport_2022', 'Sabrina_2022', 'Diglis_2022')



# 2023 VERTS DATA BIT--------------------------------------------------------------------------------

# read in the data
dat_2023 = read.csv('data/comparison/NE23_mifish_blast98_denoise.csv')

# make '.' a '-' as intended
colnames(dat_2023) = gsub('\\.', '-', colnames(dat_2023))

# grab taxonomy, remove from data
dat_2023_tax = dat_2023$taxonomy
dat_2023 = dat_2023[-grep('taxonomy', colnames(dat_2023))]

# sort row names
rownames(dat_2023) = dat_2023[,1]
dat_2023 <- dat_2023[, -1]

# noise filter data at 0.1%
dat_2023[t(t(dat_2023)/colSums(dat_2023)) < 0.001] = 0

# drop samples with less than 500 reads:
dat_2023 = dat_2023[,!colSums(dat_2023) < 500]

# add back taxonomy
dat_2023$taxonomy = dat_2023_tax

# remove assignments to Cyprinidae, Percidae and unassigned
dat_2023 = dat_2023[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(dat_2023)),]

# isolate fish and remove taxonomy (again)
dat_2023 = dat_2023[grep('Actinopteri|Petromyzontiformes', dat_2023$taxonomy),]
dat_2023 = dat_2023[-grep('taxonomy', colnames(dat_2023))]

# change column names from 'tube ids' to 'sample ids'
sample_names = read.csv('data/comparison/NE2023_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube
colnames(dat_2023) = samples[colnames(dat_2023)]

# remove species with no assignments
dat_2023 = dat_2023[rowSums(dat_2023) > 0,]

## Change Alosa_alosa to Alosa spp.
rownames(dat_2023) = sub('Alosa_alosa', 'Alosa spp.', rownames(dat_2023))

## Merge Leuciscus to Leuciscus leuciscus and change to Leuciscus spp.
Leuciscus_leuciscus = dat_2023['Leuciscus_leuciscus',]
Leuciscus = dat_2023['Leuciscus',]
dat_2023['Leuciscus_leuciscus',] = Leuciscus_leuciscus + Leuciscus
dat_2023 = dat_2023[!(rownames(dat_2023) == 'Leuciscus'),]
rownames(dat_2023) = sub('Leuciscus_leuciscus', 'Leuciscus spp.', rownames(dat_2023))

# isolate sites repeated between all years (Stourport (S5), Sabrina (S2), Diglis (S1))
Stourport_2023 = dat_2023[grep('CS03', colnames(dat_2023))]
Sabrina_2023 = dat_2023[grep('CS02', colnames(dat_2023))]
Diglis_2023 = dat_2023[grep('CS01', colnames(dat_2023))]

# pool sites
pooled_2023 = data.frame(rowSums(Stourport_2023), rowSums(Sabrina_2023), rowSums(Diglis_2023))
pooled_2023 = pooled_2023[rowSums(pooled_2023) > 0,]
pooled_2023 = pooled_2023[sort(rownames(pooled_2023)),]
names(pooled_2023) = c('Stourport_2023', 'Sabrina_2023', 'Diglis_2023')

# Verts Data Sorting --------------------------------------------------

## Make em' proportional so they are comparable between different years

## Make proportional reads eDNA dataset 
#2018
det1 = as.data.frame(t(pooled_2018))
det1$assigned=rowSums(det1)
det1=det1/det1$assigned
r=dim(det1)
det1=det1[1:(r[2]-1)]
det1$site = rownames(det1)
r=dim(det1)
det1 = melt(det1, id=(r[2]))

#2021
det2 = as.data.frame(t(pooled_2021))
det2$assigned=rowSums(det2)
det2=det2/det2$assigned
r=dim(det2)
det2=det2[1:(r[2]-1)]
det2$site = rownames(det2)
r=dim(det2)
det2 = melt(det2, id=(r[2]))

#2022
det3 = as.data.frame(t(pooled_2022))
det3$assigned=rowSums(det3)
det3=det3/det3$assigned
r=dim(det3)
det3=det3[1:(r[2]-1)]
det3$site = rownames(det3)
r=dim(det3)
det3 = melt(det3, id=(r[2]))

#2023
det4 = as.data.frame(t(pooled_2023))
det4$assigned=rowSums(det4)
det4=det4/det4$assigned
r=dim(det4)
det4=det4[1:(r[2]-1)]
det4$site = rownames(det4)
r=dim(det4)
det4 = melt(det4, id=(r[2]))


## Merge 
fin <- rbind(det1, det2, det3, det4)

## Add year
fin$year <- str_extract(fin$site, ".{4}$")

## Add Site 
fin$loc <- gsub(".{5}$", "", fin$site)
## Change Site to d/s Diglis, u/s Diglis, u/s Lincomb
fin$loc <- gsub("Stourport", "u/s Lincomb weir (03)", fin$loc)
fin$loc <- gsub("Diglis", "d/s Diglis weir (01)", fin$loc)
fin$loc <- gsub("Sabrina", "u/s Diglis weir (02)", fin$loc)


# Remove Underscore
fin$variable = gsub('_', ' ', fin$variable)

# remove zeros
fin = fin[fin$value > 0,]

## Add common names, importing them from dictionary file, matching to taxonomic names and then inserting into this df
common_names_df = read.csv('data/NE2023_CommonNames.txt', sep = '\t', header = F, col.names = c('taxonomic', 'common'))
taxonomic_common_fin = match(fin$variable, common_names_df$taxonomic)
fin$common = common_names_df$common[taxonomic_common_fin]

#Bubble plot Verts Common Names --------------------------------
comparison_common_bp = ggplot(fin, aes(x=year, y=common, size=value*100)) +
  geom_point(alpha=0.9, shape=21, fill="#70cf57", color="black") +
  scale_size(range = c(1, 10), name="Relative fish reads (%)") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Sampling year") +
  ylab("") +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))+
  scale_x_discrete(position = "top") +
  facet_wrap(~loc, scales = 'free_x') + scale_y_discrete(limits=rev)

png("output/Fish_comparison_2018_2021_2022_2023_commonnames.png", width = 3500, height = 2500, units = "px", res = 345)
comparison_common_bp
dev.off()

#Bubble plot Verts Taxonomic Names --------------------------------
comparison_bp = ggplot(fin, aes(x=year, y=variable, size=value*100)) +
  geom_point(alpha=0.9, shape=21, fill="#70cf57", color="black") +
  scale_size(range = c(1, 10), name="Relative fish reads (%)") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Sampling year") +
  ylab("") +
  theme(axis.text.y = element_text(face = "italic", size = 12))+
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))+
  scale_x_discrete(position = "top") +
  facet_wrap(~loc, scales = 'free_x') + scale_y_discrete(limits=rev)

png("output/Fish_comparison_2018_2021_2022_2023_taxonomic.png", width = 3500, height = 2500, units = "px", res = 345)
comparison_bp
dev.off()