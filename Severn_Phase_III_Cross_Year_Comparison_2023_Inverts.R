# ---
# Title: "Severn Phase III - Cross-year Comparison Species Diversity Bubble Plot - INVERTS"
# Author: "Clare E. Collins, Graham Sellers and Nathan P. Griffiths"
# Date: "19 April 2024"
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

# 2021 INVERTS DATA BIT-----------------

# ==============================================================================

# invert_dat_2021A CLEANUP - from Graham Sellers 2021 Script

# ==============================================================================

# read in the invert_dat_2021a:
invert_dat_2021 = read.csv('data/comparison/NE2021_leese_ncbi_blast95_denoise.tsv', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

# make '.' a '-' as intended:
colnames(invert_dat_2021) = gsub('\\.', '-', colnames(invert_dat_2021))

# remove duplicated samples (NE63-2-1 and NE63-2-2):
invert_dat_2021 = invert_dat_2021[-grep('NE63-2-1|NE63-2-2', colnames(invert_dat_2021))]

# rename the duplicated samples:
colnames(invert_dat_2021)[which(names(invert_dat_2021) == 'NE63-1-1')] = 'NE63'
colnames(invert_dat_2021)[which(names(invert_dat_2021) == 'NE63-1-2')] = 'NE65'

# isolate taxonomy:
invert_dat_2021_tax = invert_dat_2021$taxonomy

# remove negatives, positives, filtration and extraction blanks, and taxonomy:
invert_dat_2021 = invert_dat_2021[, -grep('taxonomy|NEG|POS', colnames(invert_dat_2021))]

# change column names from 'tube ids' to 'sample ids':
sample_names = read.csv('data/comparison/NE2021_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube
colnames(invert_dat_2021) = samples[colnames(invert_dat_2021)]

# remove 1:10 dilutions and lab filtration:
if(length(grep('_2|_lab', colnames(invert_dat_2021)))){
  invert_dat_2021 = invert_dat_2021[, -grep('_2|_lab', colnames(invert_dat_2021))]
}

OTUs_pre_threshold = nrow(invert_dat_2021) - 1

# noise filter invert_dat_2021a for < 20 reads:
invert_dat_2021[invert_dat_2021 < 20] = 0

# drop samples with less than 500 reads:
invert_dat_2021 = invert_dat_2021[!colSums(invert_dat_2021) < 500]

# remove assignments to unassigned:
invert_dat_2021 = invert_dat_2021[-grep('unassigned', rownames(invert_dat_2021)),]

# # remove unassigned from taxonomy:
# invert_dat_2021_tax = invert_dat_2021_tax[-grep('unassigned', invert_dat_2021_tax)]
# 
# # update taxonomy:
# invert_dat_2021_tax = invert_dat_2021_tax[rowSums(invert_dat_2021) > 0]

# remove species with no assignments:
invert_dat_2021 = invert_dat_2021[rowSums(invert_dat_2021) > 0,]

# order by name and season:
ordered = c(sort(colnames(invert_dat_2021[grep('Summer', colnames(invert_dat_2021))])),
            sort(colnames(invert_dat_2021[grep('Autumn', colnames(invert_dat_2021))])))

invert_dat_2021 = invert_dat_2021[,ordered]

OTUs_post_threshold = nrow(invert_dat_2021)

## Now working with clean data invert_dat_2021 -----

## add taxonomic information in to get order
## split the strings in invert_dat_2021_tax
invert_dat_2021_taxonomicdf = data.frame(taxonomic = c(invert_dat_2021_tax))
invert_dat_2021_taxonomicdf$order = invert_dat_2021_taxonomicdf
## extract the "o__text" part from each element of invert_dat_2021_taxonomic_df
invert_dat_2021_taxonomicdf$order = sapply(invert_dat_2021_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^o__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^o__", "", order_part)
  } else {
    order_text = NA  # If "o__" part not found, return NA
  }
  return(order_text)
})

## extract the s__ part to get species column
invert_dat_2021_taxonomicdf$species = sapply(invert_dat_2021_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^s__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^s__", "", order_part)
  } else {
    order_text = NA  # If "s__" part not found, return NA
  }
  return(order_text)
})

## match species in invert_dat_2021 to those in invert_dat_2021_taxonomicdf and return the correct order in an order column
invert_dat_2021 = merge(invert_dat_2021, invert_dat_2021_taxonomicdf, by.x = "row.names", by.y = "species", all.x = TRUE)

## Change rownames column heading to assignment and move order to the front and remove the taxonomic column
colnames(invert_dat_2021) [1] = c("assignment")
invert_dat_2021 = invert_dat_2021[,c(89,1:88)]
invert_dat_2021 = invert_dat_2021[!(colnames(invert_dat_2021) == 'taxonomic')]

## remove positive across all samples despite it probably being correct data [CC]
invert_dat_2021 = invert_dat_2021[invert_dat_2021$assignment != 'Calliphora_vicina',]

## show columns where order is unknown
rows_with_unknown_order21 = invert_dat_2021[invert_dat_2021$order == 'unknown', ]

## add order for above species from NBN Atlas
invert_dat_2021$order[invert_dat_2021$assignment == "Ampullaceana_balthica"] = "Hygrophila"
invert_dat_2021$order[invert_dat_2021$assignment == "Physella_acuta"] = "Hygrophila"
invert_dat_2021$order[invert_dat_2021$assignment == "Radix_auricularia"] = "Hygrophila"

## show columns where order is NA
rows_with_na_order21 = invert_dat_2021[is.na(invert_dat_2021$order), ]

## add order for rows recorded to genus level only
invert_dat_2021$order[invert_dat_2021$assignment == "Lucilia"] = "Diptera"
invert_dat_2021$order[invert_dat_2021$assignment == "Simulium"] = "Diptera"

## remove all rows where order is unknown or NA (includes family assignments and unassigned)
invert_dat_2021 = invert_dat_2021[!is.na(invert_dat_2021$order) & invert_dat_2021$order != "" & invert_dat_2021$order != "unknown", ]

## remove assignment column and combine order rows into the same row (summing data)
invert_dat_2021 = invert_dat_2021[!(colnames(invert_dat_2021) == 'assignment')]
invert_dat_2021 <- aggregate(. ~ order, data = invert_dat_2021, sum)

## make order column rownames
rownames(invert_dat_2021) = invert_dat_2021$order
invert_dat_2021 = invert_dat_2021[, -grep("order", colnames(invert_dat_2021))]

# isolate sites repeated between all years (Stourport (S5), Sabrina (S2), Diglis (S1))
Stourport_2021_inverts = invert_dat_2021[grep('S13', colnames(invert_dat_2021))]
Sabrina_2021_inverts = invert_dat_2021[grep('S11', colnames(invert_dat_2021))]
Diglis_2021_inverts = invert_dat_2021[grep('S10', colnames(invert_dat_2021))]

# pool sites
pooled_inverts_2021 = data.frame(rowSums(Stourport_2021_inverts), rowSums(Sabrina_2021_inverts), rowSums(Diglis_2021_inverts))
pooled_inverts_2021 = pooled_inverts_2021[rowSums(pooled_inverts_2021) > 0,]
pooled_inverts_2021 = pooled_inverts_2021[sort(rownames(pooled_inverts_2021)),]
names(pooled_inverts_2021) = c('Stourport_2021_inverts', 'Sabrina_2021_inverts', 'Diglis_2021_inverts')



# 2022 INVERTS DATA BIT-----------------

# read in the data
invert_dat_2022 = read.csv('data/comparison/NE22_INVERTS_95_Final.csv')

# make '.' a '-'
colnames(invert_dat_2022) = gsub('\\.', '-', colnames(invert_dat_2022))

## remove positive across all samples despite it probably being correct data [CC]
invert_dat_2022 = invert_dat_2022[invert_dat_2022$`X-OTU_ID` != 'Eudasyphora_cyanicolor',]

## remove negatives, positives, filtration and extraction blanks and Thames samples and OTU, taxonomic
invert_dat_2022 = invert_dat_2022[, -grep('X.OTU_ID|taxonomy|EB|NC|PC|T', colnames(invert_dat_2022))]

## remove any reads under 20 in any cell 
invert_dat_2022[invert_dat_2022 < 20] = 0

## combine order rows into the same row (summing data)
invert_dat_2022 <- aggregate(. ~ Order, data = invert_dat_2022, sum)

## make order column the row names
rownames(invert_dat_2022) = invert_dat_2022$Order
invert_dat_2022 = invert_dat_2022[, -grep("Order", colnames(invert_dat_2022))]

## drop samples with less than 500 reads (exclude first non-numeric column):
invert_dat_2022 = invert_dat_2022[!colSums(invert_dat_2022) < 500]

## remove orders with no assignments
invert_dat_2022 = invert_dat_2022[rowSums(invert_dat_2022) > 0, ]

## remove unassigned
invert_dat_2022 = invert_dat_2022[row.names(invert_dat_2022) != 'u__unassigned',]

# isolate sites repeated between all years (Stourport (S5), Sabrina (S2), Diglis (S1))
Stourport_2022_inverts = invert_dat_2022[grep('S5', colnames(invert_dat_2022))]
Sabrina_2022_inverts = invert_dat_2022[grep('S2', colnames(invert_dat_2022))]
Diglis_2022_inverts = invert_dat_2022[grep('S1', colnames(invert_dat_2022))]

# pool sites
pooled_inverts_2022 = data.frame(rowSums(Stourport_2022_inverts), rowSums(Sabrina_2022_inverts), rowSums(Diglis_2022_inverts))
pooled_inverts_2022 = pooled_inverts_2022[rowSums(pooled_inverts_2022) > 0,]
pooled_inverts_2022 = pooled_inverts_2022[sort(rownames(pooled_inverts_2022)),]
names(pooled_inverts_2022) = c('Stourport_2022_inverts', 'Sabrina_2022_inverts', 'Diglis_2022_inverts')


# 2023 INVERTS DATA BIT----------------

# read in the data
invert_dat_2023 = read.csv('data/comparison/NE23_leese_blast95_denoise.csv')

# make '.' a '-'
colnames(invert_dat_2023) = gsub('\\.', '_', colnames(invert_dat_2023))

# grab taxonomy, remove from data
invert_dat_2023_tax = invert_dat_2023$taxonomy
invert_dat_2023 = invert_dat_2023[-grep('taxonomy', colnames(invert_dat_2023))]

## make assignment column the row names
rownames(invert_dat_2023) = invert_dat_2023$X_OTU_ID
invert_dat_2023 = invert_dat_2023[, -grep("X_OTU_ID", colnames(invert_dat_2023))]

## add taxonomic information in to get order
## split the strings in invert_dat_2023_tax
invert_dat_2023_taxonomicdf = data.frame(taxonomic = c(invert_dat_2023_tax))
invert_dat_2023_taxonomicdf$order = invert_dat_2023_taxonomicdf
## extract the "o__text" part from each element of invert_dat_2023_taxonomic_df
invert_dat_2023_taxonomicdf$order = sapply(invert_dat_2023_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^o__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^o__", "", order_part)
  } else {
    order_text = NA  # If "o__" part not found, return NA
  }
  return(order_text)
})

## extract the s__ part to get species column
invert_dat_2023_taxonomicdf$species = sapply(invert_dat_2023_taxonomicdf$taxonomic, function(x) {
  parts = unlist(strsplit(x, "; "))
  order_part = parts[grep("^s__", parts)]
  if (length(order_part) > 0) {
    order_text = gsub("^s__", "", order_part)
  } else {
    order_text = NA  # If "s__" part not found, return NA
  }
  return(order_text)
})


## remove negatives, positives, filtration and extraction blanks and taxonomic
invert_dat_2023 = invert_dat_2023[, -grep('taxonomy|NEG|POS|CSNeg', colnames(invert_dat_2023))]

# remove positive across all samples
invert_dat_2023 = invert_dat_2023[row.names(invert_dat_2023) != 'Calliphora_vicina',]

## remove any reads under 20 in any cell 
invert_dat_2023[invert_dat_2023 < 20] = 0

# drop samples with less than 500 reads:
invert_dat_2023 = invert_dat_2023[!colSums(invert_dat_2023) < 500]

## remove species with no assignments
invert_dat_2023 = invert_dat_2023[rowSums(invert_dat_2023) > 0, ]

# change column names from 'tube ids' to 'sample ids' but first change column names to have a - instead of _
colnames(invert_dat_2023) = gsub('_', '-', colnames(invert_dat_2023))
sample_names = read.csv('data/comparison/NE2023_sample_names.txt', sep = '\t', header = F, col.names = c('tube', 'sample'))
samples = sample_names$sample
names(samples) = sample_names$tube
colnames(invert_dat_2023) = samples[colnames(invert_dat_2023)]

## match species in invert_dat_2023 to those in invert_dat_2023_taxonomicdf and return the correct order in an order column
invert_dat_2023 = merge(invert_dat_2023, invert_dat_2023_taxonomicdf, by.x = "row.names", by.y = "species", all.x = TRUE)

## Change rownames column heading to assignment and move order to the front and remove the taxonomic column
colnames(invert_dat_2023) [1] = c("assignment")
invert_dat_2023 = invert_dat_2023[,c(66,1:65)]
invert_dat_2023 = invert_dat_2023[!(colnames(invert_dat_2023) == 'taxonomic')]


## show columns where order is unknown
rows_with_unknown_order23 = invert_dat_2023[invert_dat_2023$order == 'unknown', ]

## add order for above species from NBN Atlas
invert_dat_2023$order[invert_dat_2023$assignment == "Ampullaceana_balthica"] = "Hygrophila"
invert_dat_2023$order[invert_dat_2023$assignment == "Radix_auricularia"] = "Hygrophila"

## show columns where order is NA
rows_with_na_order23 = invert_dat_2023[is.na(invert_dat_2023$order), ]

## add order for rows recorded to genus level only
invert_dat_2023$order[invert_dat_2023$assignment == "Limnephilus"] = "Trichoptera"
invert_dat_2023$order[invert_dat_2023$assignment == "Lucilia"] = "Diptera"
invert_dat_2023$order[invert_dat_2023$assignment == "Simulium"] = "Diptera"

## remove all rows where order is unknown or NA (includes family assignments and unassigned)
invert_dat_2023 = invert_dat_2023[!is.na(invert_dat_2023$order) & invert_dat_2023$order != "" & invert_dat_2023$order != "unknown", ]

## remove assignment column and combine order rows into the same row (summing data)
invert_dat_2023 = invert_dat_2023[!(colnames(invert_dat_2023) == 'assignment')]
invert_dat_2023 <- aggregate(. ~ order, data = invert_dat_2023, sum)

## make order column rownames
rownames(invert_dat_2023) = invert_dat_2023$order
invert_dat_2023 = invert_dat_2023[, -grep("order", colnames(invert_dat_2023))]

# isolate sites repeated between all years (Stourport (S5), Sabrina (S2), Diglis (S1))
Stourport_2023_inverts = invert_dat_2023[grep('CS03', colnames(invert_dat_2023))]
Sabrina_2023_inverts = invert_dat_2023[grep('CS02', colnames(invert_dat_2023))]
Diglis_2023_inverts = invert_dat_2023[grep('CS01', colnames(invert_dat_2023))]

# pool sites
pooled_inverts_2023 = data.frame(rowSums(Stourport_2023_inverts), rowSums(Sabrina_2023_inverts), rowSums(Diglis_2023_inverts))
pooled_inverts_2023 = pooled_inverts_2023[rowSums(pooled_inverts_2023) > 0,]
pooled_inverts_2023 = pooled_inverts_2023[sort(rownames(pooled_inverts_2023)),]
names(pooled_inverts_2023) = c('Stourport_2023_inverts', 'Sabrina_2023_inverts', 'Diglis_2023_inverts')



# Inverts Data Sorting --------------------------------------------------

## Make em' proportional so they are comparable between different years

## Make proportional reads eDNA dataset 

#2021
det2 = as.data.frame(t(pooled_inverts_2021))
det2$assigned=rowSums(det2)
det2=det2/det2$assigned
r=dim(det2)
det2=det2[1:(r[2]-1)]
det2$site = rownames(det2)
r=dim(det2)
det2 = melt(det2, id=(r[2]))

#2022
det3 = as.data.frame(t(pooled_inverts_2022))
det3$assigned=rowSums(det3)
det3=det3/det3$assigned
r=dim(det3)
det3=det3[1:(r[2]-1)]
det3$site = rownames(det3)
r=dim(det3)
det3 = melt(det3, id=(r[2]))

#2023
det4 = as.data.frame(t(pooled_inverts_2023))
det4$assigned=rowSums(det4)
det4=det4/det4$assigned
r=dim(det4)
det4=det4[1:(r[2]-1)]
det4$site = rownames(det4)
r=dim(det4)
det4 = melt(det4, id=(r[2]))


## Merge 
fin = rbind(det2, det3, det4)

## Add year
fin$year_inverts = str_extract(fin$site, ".{12}$")
fin$year = substr(fin$year_inverts, 1,4)
fin = fin[, -grep('year_inverts', colnames(fin))]

## Add Site by taking off the 13 letters after the site name 
fin$loc <- gsub(".{13}$", "", fin$site)
## Change Site to d/s Diglis, u/s Diglis, u/s Lincomb
fin$loc <- gsub("Stourport", "u/s Lincomb weir (03)", fin$loc)
fin$loc <- gsub("Diglis", "d/s Diglis weir (01)", fin$loc)
fin$loc <- gsub("Sabrina", "u/s Diglis weir (02)", fin$loc)


# remove zeros
fin = fin[fin$value > 0,]

#Bubble plot inVerts Orders --------------------------------
inverts_comparison_bp = ggplot(fin, aes(x=year, y=variable, size=value*100)) +
  geom_point(alpha=0.9, shape=21, fill="#78004F", color="black") +
  scale_size(range = c(1, 10), name="Relative invert reads (%)") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Sampling year") +
  ylab("") +
  theme(axis.text.y = element_text(face = "italic", size = 12))+
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))+
  scale_x_discrete(position = "top") +
  facet_wrap(~loc, scales = 'free_x') + scale_y_discrete(limits=rev)

png("output/Inverts_comparison_2021_2022_2023.png", width = 3500, height = 2500, units = "px", res = 345)
inverts_comparison_bp
dev.off()