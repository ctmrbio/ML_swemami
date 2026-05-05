## SweMaMi PTB pre-ML analysis

### Note, this code is was repeated for both vaginal and fecal sequences, separately. To make the code easier to decipher, here we present one example.

#### 1) cleaning up the files
       Will separate after decontamination.
       Also let's try this with the newly annotated samples. See how many we get.
   
#### 3) Alpha diversity -- will use this for the ML trials
#### 4) Finding the significant metadata variables using Adonis.
#### 5) Choosing the metadata variables to use.
#### 6) Running ancombc to find the differentially abundant organisms.
#### 7) The functional analysis of the samples.
#### 6) preparing for machine learning analysis.

#First step is loading all the required packages.


```R
library(ANCOMBC)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(DT)
library(mia)
library(caret)
library(microbiome)
library(microViz)
library(patchwork)
library(ggplot2)
library(mixOmics)
library(devtools)
library(vegan)
library(scutr)
library(BiocManager)
library(decontam)
library(mixOmics)
library(compositions)
library(MASS)
library(MLeval)
```

 
    


## Reading in the files and paths and making the dataframe with all the samples


```R
## reading in the lists for the resequenced samples and the previously sequenced metaphlan samples

sample_resequenced_samples <- list.files(path="/PATH", pattern="metaphlan.txt", all.files=T, full.names=T)

metaphlan_paths_sample <- list.files(path="/PATH", pattern="metaphlan.txt", all.files=T, full.names=T)

metaphlan_input_files <- c(sample_resequenced_samples,metaphlan_paths_sample)
## reading in the case control

case_control_original_file <- read.csv("/PATH/PTB_sETB_case-contrl_apr24.csv",sep=";")

## reading in the metadata information

metadata <- read.csv("/PATH/metadata_SweMaMi.csv",sep=",")
metadata_Q1_cleaned <- read.csv("/PATH/metadata_Q1_for_ML.csv",sep=",")
metadata_Q2_cleaned <- read.csv("/PATH/metadata_Q2_for_ML.csv",sep=",")
metadata_ab_kit_studienummer <- read.csv("/PATH/kit_antibiotic_use_information.csv",sep=",")
new_case_control <- read.csv("/PATH/240627_SweMaMi_PTB_case-control.csv",sep=";")

# rename colnames
colnames(new_case_control) <- c("type","Case","Control1","Control2")

```

# File preparation


```R
## combining everything and then removing all the first instances of duplicates
# and then get rid of those in this file and the combine the resequenced removed with the ones that we
# currently have in here.

## files are read in order of flowcell ID

n <- 0


### So, the goal here is to clean up the files and later separate the T1 and T2 sample analysis files.
for (species_file in metaphlan_input_files){
  
   species_input <- read.csv(species_file,sep = "\t")
   if (n == 0) {combined_species_files <- species_input} else {combined_species_files <- merge(combined_species_files,species_input,by="clade_name",all=T)}
    n <- n + 1  
    print(species_file)
  
}

### getting rid of NAs which occur when certain samples do not have a species in them. We simply set them to zero.
combined_species_files[is.na(combined_species_files)] <- 0

### for breaking down the long names for the samples to only the sample ID
   
sample_full_names <- colnames(combined_species_files)
sample_full_names <- data.frame(sample_full_names)

### this leads to a dataframe of the names split into 3 rows. The second rows is the sample ID and the one we need.
sample_IDs <- lapply(sample_full_names, function(x) str_split(x,"__"))
sample_IDs <- data.frame(sample_IDs)
                                 
sample_IDs_2 <- t(sample_IDs)
sample_IDs_2 <- data.frame(sample_IDs_2)
                    

### Now get rid of the subscripts of the name
                     
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_1","",x))
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_2","",x))
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_3","",x))

                     
### Now here I give the new cleaned column names to my columns
combined_species_with_sample_IDs <- combined_species_files
colnames(combined_species_with_sample_IDs) <- sample_IDs_2$X2
                          
                          
combined_species_with_sample_IDs_2 <- combined_species_with_sample_IDs                           
                          
rownames(combined_species_with_sample_IDs_2) <- combined_species_with_sample_IDs_2$clade_name
combined_species_with_sample_IDs_2 <- combined_species_with_sample_IDs_2[,-c(1)] 
                          
### since the files were read in flowcell order, the final repeat is the most recently sequenced one
## so going to keep that, so will remove the first occurrences of the repeated files.
                          
      
                          
reversed_df <- combined_species_with_sample_IDs_2[, rev(colnames(combined_species_with_sample_IDs_2))]
colnames(reversed_df) <- sub("\\.\\d+$", "", colnames(reversed_df))
reversed_df_filtered <- reversed_df[, !duplicated(colnames(reversed_df))]
combined_samples_unique <- reversed_df_filtered
                          
```




#### Decontamination and renormalization (given the decontam does not remove any species when used with metaphlan we skipped it in this case)

#### fixing names


```R
### let's narrow them down to species in the rownames
metaphlan_input <- combined_samples_unique
species_metaphlan_input <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("s__", rownames(metaphlan_input)[i])== TRUE && grepl("t__",rownames(metaphlan_input)[i]) == FALSE) {

        species_metaphlan_input <- rbind(species_metaphlan_input,metaphlan_input[i,])
    } 

}
```


```R
### let's narrow them down to species in the rownames
metaphlan_input <- combined_samples_unique
species_metaphlan_input <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("s__", rownames(metaphlan_input)[i])== TRUE && grepl("t__",rownames(metaphlan_input)[i]) == FALSE) {

        species_metaphlan_input <- rbind(species_metaphlan_input,metaphlan_input[i,])
    } 

}
```

#### decontamination and normalization


```R



metaphlan_input_filtered <- species_metaphlan_input
metaphlan_input_filtered$species_sum <- rowSums(metaphlan_input_filtered[,1:ncol(metaphlan_input_filtered)])
metaphlan_input_filtered <- metaphlan_input_filtered %>% filter(species_sum > 0.1)

###

metaphlan_species_all_samples <- metaphlan_input_filtered[ , -which(names(metaphlan_input_filtered) %in% c("species_sum"))]


### Here, I'm going to remove any column that has nonezero elements in only 10 or less of its entries.

metaphlan_nonzero_rows <- rowSums(species_metaphlan_input != 0) > 10
metaphlan_species_all_samples <- metaphlan_species_all_samples[metaphlan_nonzero_rows, ]


metaphlan_species_all_samples <- na.omit(metaphlan_species_all_samples)



### Adding the unclassified row back in

### adding the unclassified row (the UNCLASSIFED row is present in the initial file)
metaphlan_species_all_samples <- rbind(metaphlan_species_all_samples,metaphlan_input[nrow(metaphlan_input),])

## removing samples were the sum is equal to zero

metaphlan_species_all_samples <- metaphlan_species_all_samples[, colSums(metaphlan_species_all_samples) != 0]

### renormalizing and changing sums to 1

metaphlan_species_all_samples_normalized <- metaphlan_species_all_samples
for (i in 1:ncol(metaphlan_species_all_samples_normalized)) {
    metaphlan_species_all_samples_normalized[,i] <- 
    metaphlan_species_all_samples_normalized[i]/colSums(metaphlan_species_all_samples_normalized[i])
}

## sanity check of the column sums
print(colSums(metaphlan_species_all_samples_normalized))
```



#### CLRing: keep in mind that this is for each column as sample so I can use it on the above.


```R
metaphlan_samples_sample <- metaphlan_species_all_samples_normalized + 0.00000001
metaphlan_samples_sample_CLR <- clr(metaphlan_samples_sample)
metaphlan_samples_sample_CLR <- data.frame(metaphlan_samples_sample_CLR)
# fixing colnames since they've been renamed and were in number format initially
colnames(metaphlan_samples_sample_CLR) <- gsub("^X", "", colnames(metaphlan_samples_sample_CLR))

```





## Separating timepoints and making the case and control files

#### Putting them in one column along with their keys


```R
case_control_PT <- gather(new_case_control,"onset","studienummer",2:4)
colnames(case_control_PT) <- c("type","key","Studienummer")
case_control_PT$key <- gsub("1","",case_control_PT$key)
case_control_PT$key <- gsub("2","",case_control_PT$key)
#table(case_control_PT$key)
```

First we need transpose the metaphlan file above in order to be able to use if for separation at timepoints 1 and 2.


```R

### First step is transposing the data so I can use the merge function
metaphlan_sample_output_t <- t(metaphlan_samples_sample_CLR)
metaphlan_sample_output_t <- data.frame(metaphlan_sample_output_t)

### now CLR

metaphlan_species_non_CLR <- t(metaphlan_species_all_samples)
metaphlan_species_non_CLR <- data.frame(metaphlan_species_non_CLR)
```

### non CLRed


```R
## separating all case and controls for TP1 and TP2 (using Studienummer)

### TP1
TP1_all_samples_non_clr <- merge(studienummer_fec_kits[,c(1,2)],metaphlan_species_non_CLR,by.y=0,by.x="kit1.faecal_sample.barcode",all.x=TRUE,all.y=TRUE)


TP1_case_control_non_clr <- merge(case_control_PT,TP1_all_samples_non_clr,by="Studienummer",all.x=TRUE,all.y=FALSE)
TP1_case_control_non_clr <- na.omit(TP1_case_control_non_clr)

### TP2

TP2_all_samples_mon_clr <- merge(studienummer_fec_kits[,c(1,3)],metaphlan_species_non_CLR,by.y=0,by.x="kit2.faecal_sample.barcode",all.x=TRUE,all.y=TRUE)

TP2_case_control_non_clr <- merge(case_control_PT,TP2_all_samples_mon_clr,by="Studienummer",all.x=TRUE,all.y=TRUE)
TP2_case_control_non_clr <- na.omit(TP2_case_control_non_clr)


```

### CLRed


```R
## separating all case and controls for TP1 and TP2 (using Studienummer)

### TP1
TP1_all_samples <- merge(studienummer_fec_kits[,c(1,2)],metaphlan_sample_output_t,by.y=0,by.x="kit1.faecal_sample.barcode",all.x=TRUE,all.y=TRUE)


TP1_case_control_0 <- merge(case_control_PT,TP1_all_samples,by="Studienummer",all.x=TRUE,all.y=FALSE)
TP1_case_control_0 <- na.omit(TP1_case_control_0)

### TP2

TP2_all_samples <- merge(studienummer_fec_kits[,c(1,3)],metaphlan_sample_output_t,by.y=0,by.x="kit2.faecal_sample.barcode",all.x=TRUE,all.y=TRUE)

TP2_case_control_0 <- merge(case_control_PT,TP2_all_samples,by="Studienummer",all.x=TRUE,all.y=TRUE)
TP2_case_control_0 <- na.omit(TP2_case_control_0)

#write.csv(TP1_case_control_0,"/ceph/projects/010_SweMaMi/analyses/nicole/input_files_PTB/TP1_case_control_CLR.csv",row.names=FALSE)

```

### non CLRed


```R

## type of preterm TP1
TP1_case_control_late_non_clr <- filter(TP1_case_control_non_clr,type %in% c("late_w34-36"))
TP1_case_control_mod_non_clr <- filter(TP1_case_control_non_clr,type %in% c("moderate_w32-33"))
TP1_case_control_very_non_clr <- filter(TP1_case_control_non_clr,type %in% c("very_w28_31"))
TP1_case_control_ext_non_clr <- filter(TP1_case_control_non_clr,type %in% c("extremely_lt28"))

## type of preterm TP2
TP2_case_control_late_non_clr <- filter(TP2_case_control_non_clr,type %in% c("late_w34-36"))
TP2_case_control_mod_non_clr <- filter(TP2_case_control_non_clr,type %in% c("moderate_w32-33"))
TP2_case_control_very_non_clr <- filter(TP2_case_control_non_clr,type %in% c("very_w28_31"))
TP2_case_control_ext_non_clr <- filter(TP2_case_control_non_clr,type %in% c("extremely_lt28"))

## too few extremes (only 1) so we combine them with the very.
TP2_case_control_very_non_clr <- rbind(TP2_case_control_very_non_clr,TP2_case_control_ext_non_clr)
write.csv(TP2_case_control_very_non_clr,"/PATH/F1_V1_species.csv")
```

### CLRed


```R

## type of preterm TP1
TP1_case_control_late_0 <- filter(TP1_case_control_0,type %in% c("late_w34-36"))
TP1_case_control_mod_0 <- filter(TP1_case_control_0,type %in% c("moderate_w32-33"))
TP1_case_control_very_0 <- filter(TP1_case_control_0,type %in% c("very_w28_31"))
TP1_case_control_ext_0 <- filter(TP1_case_control_0,type %in% c("extremely_lt28"))

## type of preterm TP2
TP2_case_control_late_0 <- filter(TP2_case_control_0,type %in% c("late_w34-36"))
TP2_case_control_mod_0 <- filter(TP2_case_control_0,type %in% c("moderate_w32-33"))
TP2_case_control_very_0 <- filter(TP2_case_control_0,type %in% c("very_w28_31"))
TP2_case_control_ext_0 <- filter(TP2_case_control_0,type %in% c("extremely_lt28"))

## too few extremes (only 1) so we combine them with the very.
TP2_case_control_very_0 <- rbind(TP2_case_control_very_0,TP2_case_control_ext_0)

```

### TOP SPECIES

#### Going to find the top 50 species for PTB CC TP1 and TP2. Separate them by case and control though.



```R

################################################# All samples ##################################################


TP1_sample_composition <- TP1_case_control_non_clr[,-c(1:4)]
rownames(TP1_sample_composition) <- TP1_case_control_non_clr$Studienummer
TP1_composition_t <- t(TP1_sample_composition)

TP1_composition_t <- data.frame(TP1_composition_t)

####
TP1_composition_t$sum_species <- rowSums(TP1_composition_t[,1:ncol(TP1_composition_t)])
TP1_composition_t <- TP1_composition_t[order(TP1_composition_t$sum_species, decreasing = TRUE),]
TP1_composition_top_50 <- TP1_composition_t[1:50,]

## removing the sum_species column

TP1_composition_top_50 <- TP1_composition_top_50[ , -which(names(TP1_composition_top_50) %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP1_other_species_value <- NULL

for (i in 1:ncol(TP1_composition_top_50)) {
  TP1_other_species_value <-cbind(TP1_other_species_value,(1-sum(TP1_composition_top_50[,i])))
}

TP1_other_species_value <- data.frame(TP1_other_species_value)
colnames(TP1_other_species_value) <- colnames(TP1_composition_top_50)
TP1_composition_top_50 <- rbind(TP1_composition_top_50,TP1_other_species_value)
row.names(TP1_composition_top_50)[row.names(TP1_composition_top_50) == "1"] <- "other.species"


############################################ Only cases ########################################################

TP1_case <- filter(TP1_case_control_non_clr, key %in% c("Case"))

TP1_sample_composition_case <- TP1_case[,-c(1:4)]
rownames(TP1_sample_composition_case) <- TP1_case$Studienummer
TP1_composition_case_t <- t(TP1_sample_composition_case)

TP1_composition_case_t <- data.frame(TP1_composition_case_t)

####
TP1_composition_case_t$sum_species <- rowSums(TP1_composition_case_t[,1:ncol(TP1_composition_case_t)])
TP1_composition_case_t <- TP1_composition_case_t[order(TP1_composition_case_t$sum_species, decreasing = TRUE),]
TP1_composition_case_top_50 <- TP1_composition_case_t[1:50,]

## removing the sum_species column

TP1_composition_case_top_50 <- TP1_composition_case_top_50[ , -which(names(TP1_composition_case_top_50)
                                                                   %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP1_other_species_value_case <- NULL

for (i in 1:ncol(TP1_composition_case_top_50)) {
  TP1_other_species_value_case <-cbind(TP1_other_species_value_case,(1-sum(TP1_composition_case_top_50[,i])))
}

TP1_other_species_value_case <- data.frame(TP1_other_species_value_case)
colnames(TP1_other_species_value_case) <- colnames(TP1_composition_case_top_50)
TP1_composition_case_top_50 <- rbind(TP1_composition_case_top_50,TP1_other_species_value_case)
row.names(TP1_composition_case_top_50)[row.names(TP1_composition_case_top_50) == "1"] <- "other.species"


############################################ Only controls ########################################################

TP1_control <- filter(TP1_case_control_non_clr, key %in% c("Control"))

TP1_sample_composition_control <- TP1_control[,-c(1:4)]
rownames(TP1_sample_composition_control) <- TP1_control$Studienummer
TP1_composition_control_t <- t(TP1_sample_composition_control)

TP1_composition_control_t <- data.frame(TP1_composition_control_t)

####
TP1_composition_control_t$sum_species <- rowSums(TP1_composition_control_t[,1:ncol(TP1_composition_control_t)])
TP1_composition_control_t <- TP1_composition_control_t[order(TP1_composition_control_t$sum_species, decreasing = TRUE),]
TP1_composition_control_top_50 <- TP1_composition_control_t[1:50,]

## removing the sum_species column

TP1_composition_control_top_50 <- TP1_composition_control_top_50[ , -which(names(TP1_composition_control_top_50)
                                                                   %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP1_other_species_value_control <- NULL

for (i in 1:ncol(TP1_composition_control_top_50)) {
  TP1_other_species_value_control <-cbind(TP1_other_species_value_control,(1-sum(TP1_composition_control_top_50[,i])))
}

TP1_other_species_value_control <- data.frame(TP1_other_species_value_control)
colnames(TP1_other_species_value_control) <- colnames(TP1_composition_control_top_50)
TP1_composition_control_top_50 <- rbind(TP1_composition_control_top_50,TP1_other_species_value_control)
row.names(TP1_composition_control_top_50)[row.names(TP1_composition_control_top_50) == "1"] <- "other.species"


```

#### Finding top 50 for TP2 all case and control and then case and control separately

Next, I'm going to look at all the top 50 for the TP2 samples.


```R

################################################# All samples ##################################################


TP2_sample_composition <- TP2_case_control_non_clr[,-c(1:4)]
rownames(TP2_sample_composition) <- TP2_case_control_non_clr$Studienummer
TP2_composition_t <- t(TP2_sample_composition)

TP2_composition_t <- data.frame(TP2_composition_t)

####
TP2_composition_t$sum_species <- rowSums(TP2_composition_t[,1:ncol(TP2_composition_t)])
TP2_composition_t <- TP2_composition_t[order(TP2_composition_t$sum_species, decreasing = TRUE),]
TP2_composition_top_50 <- TP2_composition_t[1:50,]

## removing the sum_species column

TP2_composition_top_50 <- TP2_composition_top_50[ , -which(names(TP2_composition_top_50) %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP2_other_species_value <- NULL

for (i in 1:ncol(TP2_composition_top_50)) {
  TP2_other_species_value <-cbind(TP2_other_species_value,(1-sum(TP2_composition_top_50[,i])))
}

TP2_other_species_value <- data.frame(TP2_other_species_value)
colnames(TP2_other_species_value) <- colnames(TP2_composition_top_50)
TP2_composition_top_50 <- rbind(TP2_composition_top_50,TP2_other_species_value)
row.names(TP2_composition_top_50)[row.names(TP2_composition_top_50) == "1"] <- "other.species"


############################################ Only cases ########################################################

TP2_case <- filter(TP2_case_control_non_clr, key %in% c("Case"))

TP2_sample_composition_case <- TP2_case[,-c(1:4)]
rownames(TP2_sample_composition_case) <- TP2_case$Studienummer
TP2_composition_case_t <- t(TP2_sample_composition_case)

TP2_composition_case_t <- data.frame(TP2_composition_case_t)

####
TP2_composition_case_t$sum_species <- rowSums(TP2_composition_case_t[,1:ncol(TP2_composition_case_t)])
TP2_composition_case_t <- TP2_composition_case_t[order(TP2_composition_case_t$sum_species, decreasing = TRUE),]
TP2_composition_case_top_50 <- TP2_composition_case_t[1:50,]

## removing the sum_species column

TP2_composition_case_top_50 <- TP2_composition_case_top_50[ , -which(names(TP2_composition_case_top_50)
                                                                   %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP2_other_species_value_case <- NULL

for (i in 1:ncol(TP2_composition_case_top_50)) {
  TP2_other_species_value_case <-cbind(TP2_other_species_value_case,(1-sum(TP2_composition_case_top_50[,i])))
}

TP2_other_species_value_case <- data.frame(TP2_other_species_value_case)
colnames(TP2_other_species_value_case) <- colnames(TP2_composition_case_top_50)
TP2_composition_case_top_50 <- rbind(TP2_composition_case_top_50,TP2_other_species_value_case)
row.names(TP2_composition_case_top_50)[row.names(TP2_composition_case_top_50) == "1"] <- "other.species"


############################################ Only controls ########################################################

TP2_control <- filter(TP2_case_control_non_clr, key %in% c("Control"))

TP2_sample_composition_control <- TP2_control[,-c(1:4)]
rownames(TP2_sample_composition_control) <- TP2_control$Studienummer
TP2_composition_control_t <- t(TP2_sample_composition_control)

TP2_composition_control_t <- data.frame(TP2_composition_control_t)

####
TP2_composition_control_t$sum_species <- rowSums(TP2_composition_control_t[,1:ncol(TP2_composition_control_t)])
TP2_composition_control_t <- TP2_composition_control_t[order(TP2_composition_control_t$sum_species, decreasing = TRUE),]
TP2_composition_control_top_50 <- TP2_composition_control_t[1:50,]

## removing the sum_species column

TP2_composition_control_top_50 <- TP2_composition_control_top_50[ , -which(names(TP2_composition_control_top_50)
                                                                   %in% c("sum_species"))]

### adding the line for other species cause we're getting rid of a ton of species

TP2_other_species_value_control <- NULL

for (i in 1:ncol(TP2_composition_control_top_50)) {
  TP2_other_species_value_control <-cbind(TP2_other_species_value_control,(1-sum(TP2_composition_control_top_50[,i])))
}

TP2_other_species_value_control <- data.frame(TP2_other_species_value_control)
colnames(TP2_other_species_value_control) <- colnames(TP2_composition_control_top_50)
TP2_composition_control_top_50 <- rbind(TP2_composition_control_top_50,TP2_other_species_value_control)
row.names(TP2_composition_control_top_50)[row.names(TP2_composition_control_top_50) == "1"] <- "other.species"


```

## Alpha Diversity 

#### TP1 alpha diversity

For TP1 case and control alpha diversity. We will later add the case and control to them. 


```R
### all case and controls finding the species only now

TP1_case_control_species <- TP1_case_control_non_clr[,-c(1:4)]
rownames(TP1_case_control_species) <- TP1_case_control_non_clr[,1]

##### observed species all case controls
TP1_cc_richness <- specnumber(TP1_case_control_species)

TP1_cc_richness <- data.frame(TP1_cc_richness)
colnames(TP1_cc_richness) <- c("richness")


#### shannon all case controls

TP1_cc_shannon <- diversity(TP1_case_control_species)

TP1_cc_shannon <- data.frame(TP1_cc_shannon)
colnames(TP1_cc_shannon) <- c("shannon")


#### Pielou's evenness for all case controls

TP1_cc_pielou <- TP1_cc_shannon/log(TP1_cc_richness)

TP1_cc_pielou <- data.frame(TP1_cc_pielou)

colnames(TP1_cc_pielou) <- c("pielou")

#### inverse simpson all case controls

TP1_cc_invsimp <- diversity(TP1_case_control_species, index = "invsimpson")

TP1_cc_invsimp <- data.frame(TP1_cc_invsimp)
colnames(TP1_cc_invsimp) <- c("invsimp")


###### Now combining all alpha diversity indices

TP1_cc_alpha_intermediate <- merge(TP1_cc_richness,TP1_cc_shannon,by=0)
TP1_cc_alpha <- merge(TP1_cc_alpha_intermediate,TP1_cc_pielou,by.x="Row.names",by.y=0)
TP1_cc_alpha <- merge(TP1_cc_alpha,TP1_cc_invsimp,,by.x="Row.names",by.y=0)
names(TP1_cc_alpha)[names(TP1_cc_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP1_cc_alpha_key <- merge(TP1_case_control_non_clr[,c("Studienummer","key")],TP1_cc_alpha,by="Studienummer")


```


```R
write.csv(TP1_cc_alpha_key,"/PATH/V1_F1_all_species.csv",sep=",")
```

##### alpha for late preterm


```R
### all late preterm case and controls finding the species only now

TP1_case_control_late_species <- TP1_case_control_late_non_clr[,-c(1:4)]
rownames(TP1_case_control_late_species) <- TP1_case_control_late_non_clr[,1]

##### observed species all case controls
TP1_cc_late_richness <- specnumber(TP1_case_control_late_species)

TP1_cc_late_richness <- data.frame(TP1_cc_late_richness)

#### shannon all case controls

TP1_cc_late_shannon <- diversity(TP1_case_control_late_species)

TP1_cc_late_shannon <- data.frame(TP1_cc_late_shannon)

#### Pielou's evenness for all case controls

TP1_cc_late_pielou <- TP1_cc_late_shannon/log(TP1_cc_late_richness)

TP1_cc_late_pielou <- data.frame(TP1_cc_late_pielou)

colnames(TP1_cc_late_pielou) <- c("TP1_cc_late_pielou")

#### inverse simpson all case controls

TP1_cc_late_invsimp <- diversity(TP1_case_control_late_species, index = "invsimpson")

TP1_cc_late_invsimp <- data.frame(TP1_cc_late_invsimp)

###### Now combining all alpha diversity indices

TP1_cc_late_alpha_intermediate <- merge(TP1_cc_late_richness,TP1_cc_late_shannon,by=0)
TP1_cc_late_alpha <- merge(TP1_cc_late_alpha_intermediate,TP1_cc_late_pielou,by.x="Row.names",by.y=0)
TP1_cc_late_alpha <- merge(TP1_cc_late_alpha,TP1_cc_late_invsimp,,by.x="Row.names",by.y=0)
names(TP1_cc_late_alpha)[names(TP1_cc_late_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP1_cc_late_alpha_key <- merge(TP1_case_control_late_non_clr[,c("Studienummer","key")],TP1_cc_late_alpha,by="Studienummer")

```

##### Alpha for moderate


```R
### all mod preterm case and controls finding the species only now

TP1_case_control_mod_species <- TP1_case_control_mod_non_clr[,-c(1:4)]
rownames(TP1_case_control_mod_species) <- TP1_case_control_mod_non_clr[,1]

##### observed species all case controls
TP1_cc_mod_richness <- specnumber(TP1_case_control_mod_species)

TP1_cc_mod_richness <- data.frame(TP1_cc_mod_richness)

#### shannon all case controls

TP1_cc_mod_shannon <- diversity(TP1_case_control_mod_species)

TP1_cc_mod_shannon <- data.frame(TP1_cc_mod_shannon)

#### Pielou's evenness for all case controls

TP1_cc_mod_pielou <- TP1_cc_mod_shannon/log(TP1_cc_mod_richness)

TP1_cc_mod_pielou <- data.frame(TP1_cc_mod_pielou)

colnames(TP1_cc_mod_pielou) <- c("TP1_cc_mod_pielou")

#### inverse simpson all case controls

TP1_cc_mod_invsimp <- diversity(TP1_case_control_mod_species, index = "invsimpson")

TP1_cc_mod_invsimp <- data.frame(TP1_cc_mod_invsimp)

###### Now combining all alpha diversity indices

TP1_cc_mod_alpha_intermediate <- merge(TP1_cc_mod_richness,TP1_cc_mod_shannon,by=0)
TP1_cc_mod_alpha <- merge(TP1_cc_mod_alpha_intermediate,TP1_cc_mod_pielou,by.x="Row.names",by.y=0)
TP1_cc_mod_alpha <- merge(TP1_cc_mod_alpha,TP1_cc_mod_invsimp,,by.x="Row.names",by.y=0)
names(TP1_cc_mod_alpha)[names(TP1_cc_mod_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP1_cc_mod_alpha_key <- merge(TP1_case_control_mod_non_clr[,c("Studienummer","key")],TP1_cc_mod_alpha,by="Studienummer")


```

##### Alpha for very preterm


```R
### all very preterm case and controls finding the species only now

TP1_case_control_very_species <- TP1_case_control_very_non_clr[,-c(1:4)]
rownames(TP1_case_control_very_species) <- TP1_case_control_very_non_clr[,1]

##### observed species all case controls
TP1_cc_very_richness <- specnumber(TP1_case_control_very_species)

TP1_cc_very_richness <- data.frame(TP1_cc_very_richness)

#### shannon all case controls

TP1_cc_very_shannon <- diversity(TP1_case_control_very_species)

TP1_cc_very_shannon <- data.frame(TP1_cc_very_shannon)

#### Pielou's evenness for all case controls

TP1_cc_very_pielou <- TP1_cc_very_shannon/log(TP1_cc_very_richness)

TP1_cc_very_pielou <- data.frame(TP1_cc_very_pielou)

colnames(TP1_cc_very_pielou) <- c("TP1_cc_very_pielou")

#### inverse simpson all case controls

TP1_cc_very_invsimp <- diversity(TP1_case_control_very_species, index = "invsimpson")

TP1_cc_very_invsimp <- data.frame(TP1_cc_very_invsimp)

###### Now combining all alpha diversity indices

TP1_cc_very_alpha_intermediate <- merge(TP1_cc_very_richness,TP1_cc_very_shannon,by=0)
TP1_cc_very_alpha <- merge(TP1_cc_very_alpha_intermediate,TP1_cc_very_pielou,by.x="Row.names",by.y=0)
TP1_cc_very_alpha <- merge(TP1_cc_very_alpha,TP1_cc_very_invsimp,,by.x="Row.names",by.y=0)
names(TP1_cc_very_alpha)[names(TP1_cc_very_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP1_cc_very_alpha_key <- merge(TP1_case_control_very_non_clr[,c("Studienummer","key")],TP1_cc_very_alpha,by="Studienummer")



```

###### alpha for ext preterm


```R
### all ext preterm case and controls finding the species only now

TP1_case_control_ext_species <- TP1_case_control_ext_non_clr[,-c(1:4)]
rownames(TP1_case_control_ext_species) <- TP1_case_control_ext_non_clr[,1]

##### observed species all case controls
TP1_cc_ext_richness <- specnumber(TP1_case_control_ext_species)

TP1_cc_ext_richness <- data.frame(TP1_cc_ext_richness)

#### shannon all case controls

TP1_cc_ext_shannon <- diversity(TP1_case_control_ext_species)

TP1_cc_ext_shannon <- data.frame(TP1_cc_ext_shannon)

#### Pielou's evenness for all case controls

TP1_cc_ext_pielou <- TP1_cc_ext_shannon/log(TP1_cc_ext_richness)

TP1_cc_ext_pielou <- data.frame(TP1_cc_ext_pielou)

colnames(TP1_cc_ext_pielou) <- c("TP1_cc_ext_pielou")

#### inverse simpson all case controls

TP1_cc_ext_invsimp <- diversity(TP1_case_control_ext_species, index = "invsimpson")

TP1_cc_ext_invsimp <- data.frame(TP1_cc_ext_invsimp)

###### Now combining all alpha diversity indices

TP1_cc_ext_alpha_intermediate <- merge(TP1_cc_ext_richness,TP1_cc_ext_shannon,by=0)
TP1_cc_ext_alpha <- merge(TP1_cc_ext_alpha_intermediate,TP1_cc_ext_pielou,by.x="Row.names",by.y=0)
TP1_cc_ext_alpha <- merge(TP1_cc_ext_alpha,TP1_cc_ext_invsimp,,by.x="Row.names",by.y=0)
names(TP1_cc_ext_alpha)[names(TP1_cc_ext_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP1_cc_ext_alpha_key <- merge(TP1_case_control_ext_non_clr[,c("Studienummer","key")],TP1_cc_ext_alpha,by="Studienummer")

```

#### wilcoxon test now


```R

TP1_wilcoxon_test <- list()

############# all case and control

TP1_alpha_case <- filter(TP1_cc_alpha_key, key %in% c("Case"))
TP1_alpha_control <- filter(TP1_cc_alpha_key, key %in% c("Control"))

### shannon

TP1_case_shannon <- TP1_alpha_case$TP1_cc_shannon
TP1_control_shannon <- TP1_alpha_control$TP1_cc_shannon
TP1_shannon_wilcoxon <- wilcox.test(TP1_case_shannon, TP1_control_shannon, alternative = "two.sided")
TP1_wilcoxon_test[[1]] <- TP1_shannon_wilcoxon
### richness

TP1_case_richness <- TP1_alpha_case$TP1_cc_richness
TP1_control_richness <- TP1_alpha_control$TP1_cc_richness
TP1_richness_wilcoxon <- wilcox.test(TP1_case_richness, TP1_control_richness, alternative = "two.sided")
TP1_wilcoxon_test[[2]] <- TP1_richness_wilcoxon

### Pielou

TP1_case_pielou <- TP1_alpha_case$TP1_cc_pielou
TP1_control_pielou <- TP1_alpha_control$TP1_cc_pielou
TP1_pielou_wilcoxon <- wilcox.test(TP1_case_pielou, TP1_control_pielou, alternative = "two.sided")
TP1_wilcoxon_test[[3]] <- TP1_pielou_wilcoxon

### Inverse Simpson

TP1_case_invsimp <- TP1_alpha_case$TP1_cc_invsimp
TP1_control_invsimp <- TP1_alpha_control$TP1_cc_invsimp
TP1_invsimp_wilcoxon <- wilcox.test(TP1_case_invsimp, TP1_control_invsimp, alternative = "two.sided")
TP1_wilcoxon_test[[4]] <- TP1_invsimp_wilcoxon

```

##### for late preterm


```R

TP1_wilcoxon_late_test <- list()

############# all case_late and control_late

TP1_alpha_case_late <- filter(TP1_cc_late_alpha_key, key %in% c("Case"))
TP1_alpha_control_late <- filter(TP1_cc_late_alpha_key, key %in% c("Control"))

### shannon

TP1_case_late_shannon <- TP1_alpha_case_late$TP1_cc_late_shannon
TP1_control_late_shannon <- TP1_alpha_control_late$TP1_cc_late_shannon
TP1_shannon_wilcoxon_late <- wilcox.test(TP1_case_late_shannon, TP1_control_late_shannon, alternative = "two.sided")
TP1_wilcoxon_late_test[[1]] <- TP1_shannon_wilcoxon_late
### richness

TP1_case_late_richness <- TP1_alpha_case_late$TP1_cc_late_richness
TP1_control_late_richness <- TP1_alpha_control_late$TP1_cc_late_richness
TP1_richness_wilcoxon_late <- wilcox.test(TP1_case_late_richness, TP1_control_late_richness, alternative = "two.sided")
TP1_wilcoxon_late_test[[2]] <- TP1_richness_wilcoxon_late

### Pielou

TP1_case_late_pielou <- TP1_alpha_case_late$TP1_cc_late_pielou
TP1_control_late_pielou <- TP1_alpha_control_late$TP1_cc_late_pielou
TP1_pielou_wilcoxon_late <- wilcox.test(TP1_case_late_pielou, TP1_control_late_pielou, alternative = "two.sided")
TP1_wilcoxon_late_test[[3]] <- TP1_pielou_wilcoxon_late

### Inverse Simpson

TP1_case_late_invsimp <- TP1_alpha_case_late$TP1_cc_late_invsimp
TP1_control_late_invsimp <- TP1_alpha_control_late$TP1_cc_late_invsimp
TP1_invsimp_wilcoxon_late <- wilcox.test(TP1_case_late_invsimp, TP1_control_late_invsimp, alternative = "two.sided")
TP1_wilcoxon_late_test[[4]] <- TP1_invsimp_wilcoxon_late

```

##### for moderate preterm


```R

TP1_wilcoxon_mod_test <- list()

############# all case_mod and control_mod

TP1_alpha_case_mod <- filter(TP1_cc_mod_alpha_key, key %in% c("Case"))
TP1_alpha_control_mod <- filter(TP1_cc_mod_alpha_key, key %in% c("Control"))

### shannon

TP1_case_mod_shannon <- TP1_alpha_case_mod$TP1_cc_mod_shannon
TP1_control_mod_shannon <- TP1_alpha_control_mod$TP1_cc_mod_shannon
TP1_shannon_wilcoxon_mod <- wilcox.test(TP1_case_mod_shannon, TP1_control_mod_shannon, alternative = "two.sided")
TP1_wilcoxon_mod_test[[1]] <- TP1_shannon_wilcoxon_mod
### richness

TP1_case_mod_richness <- TP1_alpha_case_mod$TP1_cc_mod_richness
TP1_control_mod_richness <- TP1_alpha_control_mod$TP1_cc_mod_richness
TP1_richness_wilcoxon_mod <- wilcox.test(TP1_case_mod_richness, TP1_control_mod_richness, alternative = "two.sided")
TP1_wilcoxon_mod_test[[2]] <- TP1_richness_wilcoxon_mod

### Pielou

TP1_case_mod_pielou <- TP1_alpha_case_mod$TP1_cc_mod_pielou
TP1_control_mod_pielou <- TP1_alpha_control_mod$TP1_cc_mod_pielou
TP1_pielou_wilcoxon_mod <- wilcox.test(TP1_case_mod_pielou, TP1_control_mod_pielou, alternative = "two.sided")
TP1_wilcoxon_mod_test[[3]] <- TP1_pielou_wilcoxon_mod

### Inverse Simpson

TP1_case_mod_invsimp <- TP1_alpha_case_mod$TP1_cc_mod_invsimp
TP1_control_mod_invsimp <- TP1_alpha_control_mod$TP1_cc_mod_invsimp
TP1_invsimp_wilcoxon_mod <- wilcox.test(TP1_case_mod_invsimp, TP1_control_mod_invsimp, alternative = "two.sided")
TP1_wilcoxon_mod_test[[4]] <- TP1_invsimp_wilcoxon_mod


```


###### very preterm


```R

TP1_wilcoxon_very_test <- list()

############# all case_very and control_very

TP1_alpha_case_very <- filter(TP1_cc_very_alpha_key, key %in% c("Case"))
TP1_alpha_control_very <- filter(TP1_cc_very_alpha_key, key %in% c("Control"))

### shannon

TP1_case_very_shannon <- TP1_alpha_case_very$TP1_cc_very_shannon
TP1_control_very_shannon <- TP1_alpha_control_very$TP1_cc_very_shannon
TP1_shannon_wilcoxon_very <- wilcox.test(TP1_case_very_shannon, TP1_control_very_shannon, alternative = "two.sided")
TP1_wilcoxon_very_test[[1]] <- TP1_shannon_wilcoxon_very
### richness

TP1_case_very_richness <- TP1_alpha_case_very$TP1_cc_very_richness
TP1_control_very_richness <- TP1_alpha_control_very$TP1_cc_very_richness
TP1_richness_wilcoxon_very <- wilcox.test(TP1_case_very_richness, TP1_control_very_richness, alternative = "two.sided")
TP1_wilcoxon_very_test[[2]] <- TP1_richness_wilcoxon_very

### Pielou

TP1_case_very_pielou <- TP1_alpha_case_very$TP1_cc_very_pielou
TP1_control_very_pielou <- TP1_alpha_control_very$TP1_cc_very_pielou
TP1_pielou_wilcoxon_very <- wilcox.test(TP1_case_very_pielou, TP1_control_very_pielou, alternative = "two.sided")
TP1_wilcoxon_very_test[[3]] <- TP1_pielou_wilcoxon_very

### Inverse Simpson

TP1_case_very_invsimp <- TP1_alpha_case_very$TP1_cc_very_invsimp
TP1_control_very_invsimp <- TP1_alpha_control_very$TP1_cc_very_invsimp
TP1_invsimp_wilcoxon_very <- wilcox.test(TP1_case_very_invsimp, TP1_control_very_invsimp, alternative = "two.sided")
TP1_wilcoxon_very_test[[4]] <- TP1_invsimp_wilcoxon_very


```


##### for extreme


```R

TP1_wilcoxon_ext_test <- list()

############# all case_ext and control_ext

TP1_alpha_case_ext <- filter(TP1_cc_ext_alpha_key, key %in% c("Case"))
TP1_alpha_control_ext <- filter(TP1_cc_ext_alpha_key, key %in% c("Control"))

### shannon

TP1_case_ext_shannon <- TP1_alpha_case_ext$TP1_cc_ext_shannon
TP1_control_ext_shannon <- TP1_alpha_control_ext$TP1_cc_ext_shannon
TP1_shannon_wilcoxon_ext <- wilcox.test(TP1_case_ext_shannon, TP1_control_ext_shannon, alternative = "two.sided")
TP1_wilcoxon_ext_test[[1]] <- TP1_shannon_wilcoxon_ext
### richness

TP1_case_ext_richness <- TP1_alpha_case_ext$TP1_cc_ext_richness
TP1_control_ext_richness <- TP1_alpha_control_ext$TP1_cc_ext_richness
TP1_richness_wilcoxon_ext <- wilcox.test(TP1_case_ext_richness, TP1_control_ext_richness, alternative = "two.sided")
TP1_wilcoxon_ext_test[[2]] <- TP1_richness_wilcoxon_ext

### Pielou

TP1_case_ext_pielou <- TP1_alpha_case_ext$TP1_cc_ext_pielou
TP1_control_ext_pielou <- TP1_alpha_control_ext$TP1_cc_ext_pielou
TP1_pielou_wilcoxon_ext <- wilcox.test(TP1_case_ext_pielou, TP1_control_ext_pielou, alternative = "two.sided")
TP1_wilcoxon_ext_test[[3]] <- TP1_pielou_wilcoxon_ext

### Inverse Simpson

TP1_case_ext_invsimp <- TP1_alpha_case_ext$TP1_cc_ext_invsimp
TP1_control_ext_invsimp <- TP1_alpha_control_ext$TP1_cc_ext_invsimp
TP1_invsimp_wilcoxon_ext <- wilcox.test(TP1_case_ext_invsimp, TP1_control_ext_invsimp, alternative = "two.sided")
TP1_wilcoxon_ext_test[[4]] <- TP1_invsimp_wilcoxon_ext

```

#### violing plots for all


```R
TP1_alpha_plot<- ggplot(TP1_cc_alpha_key, aes(x = key, y = TP1_cc_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_shannon.pdf",plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_alpha_key, aes(x = key, y = TP1_cc_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_invsimp.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_alpha_key, aes(x = key, y = TP1_cc_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_pielou.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_alpha_key, aes(x = key, y = TP1_cc_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_richness.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 


```

###### violin for late


```R
TP1_alpha_plot<- ggplot(TP1_cc_late_alpha_key, aes(x = key, y = TP1_cc_late_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_late_shannon.pdf",plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_late_alpha_key, aes(x = key, y = TP1_cc_late_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_late_invsimp.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_late_alpha_key, aes(x = key, y = TP1_cc_late_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_late_pielou.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_late_alpha_key, aes(x = key, y = TP1_cc_late_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_late_richness.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 


```

##### violin for moderate


```R
TP1_alpha_plot<- ggplot(TP1_cc_mod_alpha_key, aes(x = key, y = TP1_cc_mod_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_mod_shannon.pdf",plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_mod_alpha_key, aes(x = key, y = TP1_cc_mod_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_mod_invsimp.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_mod_alpha_key, aes(x = key, y = TP1_cc_mod_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_mod_pielou.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_mod_alpha_key, aes(x = key, y = TP1_cc_mod_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_mod_richness.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 


```

##### violin for very


```R
TP1_alpha_plot<- ggplot(TP1_cc_very_alpha_key, aes(x = key, y = TP1_cc_very_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_very_shannon.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_very_alpha_key, aes(x = key, y = TP1_cc_very_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_very_invsimp.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_very_alpha_key, aes(x = key, y = TP1_cc_very_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_very_pielou.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_very_alpha_key, aes(x = key, y = TP1_cc_very_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_very_richness.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

```

##### violin for ext


```R
TP1_alpha_plot<- ggplot(TP1_cc_ext_alpha_key, aes(x = key, y = TP1_cc_ext_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_ext_shannon.pdf",plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_ext_alpha_key, aes(x = key, y = TP1_cc_ext_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_ext_invsimp.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_ext_alpha_key, aes(x = key, y = TP1_cc_ext_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_ext_pielou.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 

TP1_alpha_plot<- ggplot(TP1_cc_ext_alpha_key, aes(x = key, y = TP1_cc_ext_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP1_cc_ext_richness.pdf", plot = TP1_alpha_plot, width = 12, height = 8, dpi = 600) 


```

## TP2 alpha diversity


```R
### all case and controls finding the species only now

TP2_case_control_species <- TP2_case_control_non_clr[,-c(1:4)]
rownames(TP2_case_control_species) <- TP2_case_control_non_clr[,1]

##### observed species all case controls
TP2_cc_richness <- specnumber(TP2_case_control_species)

TP2_cc_richness <- data.frame(TP2_cc_richness)
colnames(TP2_cc_richness) <- c("richness")


#### shannon all case controls

TP2_cc_shannon <- diversity(TP2_case_control_species)

TP2_cc_shannon <- data.frame(TP2_cc_shannon)
colnames(TP2_cc_shannon) <- c("shannon")


#### Pielou's evenness for all case controls

TP2_cc_pielou <- TP2_cc_shannon/log(TP2_cc_richness)

TP2_cc_pielou <- data.frame(TP2_cc_pielou)

colnames(TP2_cc_pielou) <- c("pielou")

#### inverse simpson all case controls

TP2_cc_invsimp <- diversity(TP2_case_control_species, index = "invsimpson")

TP2_cc_invsimp <- data.frame(TP2_cc_invsimp)
colnames(TP2_cc_invsimp) <- c("invsimp")


###### Now combining all alpha diversity indices

TP2_cc_alpha_intermediate <- merge(TP2_cc_richness,TP2_cc_shannon,by=0)
TP2_cc_alpha <- merge(TP2_cc_alpha_intermediate,TP2_cc_pielou,by.x="Row.names",by.y=0)
TP2_cc_alpha <- merge(TP2_cc_alpha,TP2_cc_invsimp,,by.x="Row.names",by.y=0)
names(TP2_cc_alpha)[names(TP2_cc_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP2_cc_alpha_key <- merge(TP2_case_control_non_clr[,c("Studienummer","key")],TP2_cc_alpha,by="Studienummer")


```


```R
write.csv(TP2_cc_alpha_key,"/PATH/V2_F2_all_species.csv",sep=",")
```

##### later preterm


```R
### all late preterm case and controls finding the species only now

TP2_case_control_late_species <- TP2_case_control_late_non_clr[,-c(1:4)]
rownames(TP2_case_control_late_species) <- TP2_case_control_late_non_clr[,1]

##### observed species all case controls
TP2_cc_late_richness <- specnumber(TP2_case_control_late_species)

TP2_cc_late_richness <- data.frame(TP2_cc_late_richness)

#### shannon all case controls

TP2_cc_late_shannon <- diversity(TP2_case_control_late_species)

TP2_cc_late_shannon <- data.frame(TP2_cc_late_shannon)

#### Pielou's evenness for all case controls

TP2_cc_late_pielou <- TP2_cc_late_shannon/log(TP2_cc_late_richness)

TP2_cc_late_pielou <- data.frame(TP2_cc_late_pielou)

colnames(TP2_cc_late_pielou) <- c("TP2_cc_late_pielou")

#### inverse simpson all case controls

TP2_cc_late_invsimp <- diversity(TP2_case_control_late_species, index = "invsimpson")

TP2_cc_late_invsimp <- data.frame(TP2_cc_late_invsimp)

###### Now combining all alpha diversity indices

TP2_cc_late_alpha_intermediate <- merge(TP2_cc_late_richness,TP2_cc_late_shannon,by=0)
TP2_cc_late_alpha <- merge(TP2_cc_late_alpha_intermediate,TP2_cc_late_pielou,by.x="Row.names",by.y=0)
TP2_cc_late_alpha <- merge(TP2_cc_late_alpha,TP2_cc_late_invsimp,,by.x="Row.names",by.y=0)
names(TP2_cc_late_alpha)[names(TP2_cc_late_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP2_cc_late_alpha_key <- merge(TP2_case_control_late_non_clr[,c("Studienummer","key")],TP2_cc_late_alpha,by="Studienummer")

#write.csv(TP2_cc_late_alpha_key,"/PATH/TP2_cc_late_alpha_key.csv")

```

###### moderate preterm


```R
### all mod preterm case and controls finding the species only now

TP2_case_control_mod_species <- TP2_case_control_mod_non_clr[,-c(1:4)]
rownames(TP2_case_control_mod_species) <- TP2_case_control_mod_non_clr[,1]

##### observed species all case controls
TP2_cc_mod_richness <- specnumber(TP2_case_control_mod_species)

TP2_cc_mod_richness <- data.frame(TP2_cc_mod_richness)

#### shannon all case controls

TP2_cc_mod_shannon <- diversity(TP2_case_control_mod_species)

TP2_cc_mod_shannon <- data.frame(TP2_cc_mod_shannon)

#### Pielou's evenness for all case controls

TP2_cc_mod_pielou <- TP2_cc_mod_shannon/log(TP2_cc_mod_richness)

TP2_cc_mod_pielou <- data.frame(TP2_cc_mod_pielou)

colnames(TP2_cc_mod_pielou) <- c("TP2_cc_mod_pielou")

#### inverse simpson all case controls

TP2_cc_mod_invsimp <- diversity(TP2_case_control_mod_species, index = "invsimpson")

TP2_cc_mod_invsimp <- data.frame(TP2_cc_mod_invsimp)

###### Now combining all alpha diversity indices

TP2_cc_mod_alpha_intermediate <- merge(TP2_cc_mod_richness,TP2_cc_mod_shannon,by=0)
TP2_cc_mod_alpha <- merge(TP2_cc_mod_alpha_intermediate,TP2_cc_mod_pielou,by.x="Row.names",by.y=0)
TP2_cc_mod_alpha <- merge(TP2_cc_mod_alpha,TP2_cc_mod_invsimp,,by.x="Row.names",by.y=0)
names(TP2_cc_mod_alpha)[names(TP2_cc_mod_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP2_cc_mod_alpha_key <- merge(TP2_case_control_mod_non_clr[,c("Studienummer","key")],TP2_cc_mod_alpha,by="Studienummer")

#write.csv(TP2_cc_mod_alpha_key,"/PATH/TP2_cc_mod_alpha_key.csv")


```

##### very preterm


```R
### all very preterm case and controls finding the species only now

TP2_case_control_very_species <- TP2_case_control_very_non_clr[,-c(1:4)]
rownames(TP2_case_control_very_species) <- TP2_case_control_very_non_clr[,1]

##### observed species all case controls
TP2_cc_very_richness <- specnumber(TP2_case_control_very_species)

TP2_cc_very_richness <- data.frame(TP2_cc_very_richness)

#### shannon all case controls

TP2_cc_very_shannon <- diversity(TP2_case_control_very_species)

TP2_cc_very_shannon <- data.frame(TP2_cc_very_shannon)

#### Pielou's evenness for all case controls

TP2_cc_very_pielou <- TP2_cc_very_shannon/log(TP2_cc_very_richness)

TP2_cc_very_pielou <- data.frame(TP2_cc_very_pielou)

colnames(TP2_cc_very_pielou) <- c("TP2_cc_very_pielou")

#### inverse simpson all case controls

TP2_cc_very_invsimp <- diversity(TP2_case_control_very_species, index = "invsimpson")

TP2_cc_very_invsimp <- data.frame(TP2_cc_very_invsimp)

###### Now combining all alpha diversity indices

TP2_cc_very_alpha_intermediate <- merge(TP2_cc_very_richness,TP2_cc_very_shannon,by=0)
TP2_cc_very_alpha <- merge(TP2_cc_very_alpha_intermediate,TP2_cc_very_pielou,by.x="Row.names",by.y=0)
TP2_cc_very_alpha <- merge(TP2_cc_very_alpha,TP2_cc_very_invsimp,,by.x="Row.names",by.y=0)
names(TP2_cc_very_alpha)[names(TP2_cc_very_alpha) == "Row.names"] <- "Studienummer"

## Merging with the key column so we know which studienummer corresponds to case and which to control

## all case and controls
TP2_cc_very_alpha_key <- merge(TP2_case_control_very_non_clr[,c("Studienummer","key")],TP2_cc_very_alpha,by="Studienummer")

#write.csv(TP2_cc_very_alpha_key,"/PATH/TP2_cc_very_alpha_key.csv")


```

#### wilcoxon test now


```R
TP2_wilcoxon_test <- list()

############# all case and control

TP2_alpha_case <- filter(TP2_cc_alpha_key, key %in% c("Case"))
TP2_alpha_control <- filter(TP2_cc_alpha_key, key %in% c("Control"))

### shannon

TP2_case_shannon <- TP2_alpha_case$TP2_cc_shannon
TP2_control_shannon <- TP2_alpha_control$TP2_cc_shannon
TP2_shannon_wilcoxon <- wilcox.test(TP2_case_shannon, TP2_control_shannon, alternative = "two.sided")
TP2_wilcoxon_test[[1]] <- TP2_shannon_wilcoxon
### richness

TP2_case_richness <- TP2_alpha_case$TP2_cc_richness
TP2_control_richness <- TP2_alpha_control$TP2_cc_richness
TP2_richness_wilcoxon <- wilcox.test(TP2_case_richness, TP2_control_richness, alternative = "two.sided")
TP2_wilcoxon_test[[2]] <- TP2_richness_wilcoxon

### Pielou

TP2_case_pielou <- TP2_alpha_case$TP2_cc_pielou
TP2_control_pielou <- TP2_alpha_control$TP2_cc_pielou
TP2_pielou_wilcoxon <- wilcox.test(TP2_case_pielou, TP2_control_pielou, alternative = "two.sided")
TP2_wilcoxon_test[[3]] <- TP2_pielou_wilcoxon

### Inverse Simpson

TP2_case_invsimp <- TP2_alpha_case$TP2_cc_invsimp
TP2_control_invsimp <- TP2_alpha_control$TP2_cc_invsimp
TP2_invsimp_wilcoxon <- wilcox.test(TP2_case_invsimp, TP2_control_invsimp, alternative = "two.sided")
TP2_wilcoxon_test[[4]] <- TP2_invsimp_wilcoxon

#capture.output(TP2_wilcoxon_test,file="/PATH/TP2_wilcoxon_test_PTB.txt")

```

##### late preterm


```R

TP2_wilcoxon_late_test <- list()

############# all case_late and control_late

TP2_alpha_case_late <- filter(TP2_cc_late_alpha_key, key %in% c("Case"))
TP2_alpha_control_late <- filter(TP2_cc_late_alpha_key, key %in% c("Control"))

### shannon

TP2_case_late_shannon <- TP2_alpha_case_late$TP2_cc_late_shannon
TP2_control_late_shannon <- TP2_alpha_control_late$TP2_cc_late_shannon
TP2_shannon_wilcoxon_late <- wilcox.test(TP2_case_late_shannon, TP2_control_late_shannon, alternative = "two.sided")
TP2_wilcoxon_late_test[[1]] <- TP2_shannon_wilcoxon_late
### richness

TP2_case_late_richness <- TP2_alpha_case_late$TP2_cc_late_richness
TP2_control_late_richness <- TP2_alpha_control_late$TP2_cc_late_richness
TP2_richness_wilcoxon_late <- wilcox.test(TP2_case_late_richness, TP2_control_late_richness, alternative = "two.sided")
TP2_wilcoxon_late_test[[2]] <- TP2_richness_wilcoxon_late

### Pielou

TP2_case_late_pielou <- TP2_alpha_case_late$TP2_cc_late_pielou
TP2_control_late_pielou <- TP2_alpha_control_late$TP2_cc_late_pielou
TP2_pielou_wilcoxon_late <- wilcox.test(TP2_case_late_pielou, TP2_control_late_pielou, alternative = "two.sided")
TP2_wilcoxon_late_test[[3]] <- TP2_pielou_wilcoxon_late

### Inverse Simpson

TP2_case_late_invsimp <- TP2_alpha_case_late$TP2_cc_late_invsimp
TP2_control_late_invsimp <- TP2_alpha_control_late$TP2_cc_late_invsimp
TP2_invsimp_wilcoxon_late <- wilcox.test(TP2_case_late_invsimp, TP2_control_late_invsimp, alternative = "two.sided")
TP2_wilcoxon_late_test[[4]] <- TP2_invsimp_wilcoxon_late

#capture.output(TP2_wilcoxon_late_test,file="/PATH/TP2_wilcoxon_late_test_PTB.txt")

```

##### moderate preterm


```R

TP2_wilcoxon_mod_test <- list()

############# all case_mod and control_mod

TP2_alpha_case_mod <- filter(TP2_cc_mod_alpha_key, key %in% c("Case"))
TP2_alpha_control_mod <- filter(TP2_cc_mod_alpha_key, key %in% c("Control"))

### shannon

TP2_case_mod_shannon <- TP2_alpha_case_mod$TP2_cc_mod_shannon
TP2_control_mod_shannon <- TP2_alpha_control_mod$TP2_cc_mod_shannon
TP2_shannon_wilcoxon_mod <- wilcox.test(TP2_case_mod_shannon, TP2_control_mod_shannon, alternative = "two.sided")
TP2_wilcoxon_mod_test[[1]] <- TP2_shannon_wilcoxon_mod
### richness

TP2_case_mod_richness <- TP2_alpha_case_mod$TP2_cc_mod_richness
TP2_control_mod_richness <- TP2_alpha_control_mod$TP2_cc_mod_richness
TP2_richness_wilcoxon_mod <- wilcox.test(TP2_case_mod_richness, TP2_control_mod_richness, alternative = "two.sided")
TP2_wilcoxon_mod_test[[2]] <- TP2_richness_wilcoxon_mod

### Pielou

TP2_case_mod_pielou <- TP2_alpha_case_mod$TP2_cc_mod_pielou
TP2_control_mod_pielou <- TP2_alpha_control_mod$TP2_cc_mod_pielou
TP2_pielou_wilcoxon_mod <- wilcox.test(TP2_case_mod_pielou, TP2_control_mod_pielou, alternative = "two.sided")
TP2_wilcoxon_mod_test[[3]] <- TP2_pielou_wilcoxon_mod

### Inverse Simpson

TP2_case_mod_invsimp <- TP2_alpha_case_mod$TP2_cc_mod_invsimp
TP2_control_mod_invsimp <- TP2_alpha_control_mod$TP2_cc_mod_invsimp
TP2_invsimp_wilcoxon_mod <- wilcox.test(TP2_case_mod_invsimp, TP2_control_mod_invsimp, alternative = "two.sided")
TP2_wilcoxon_mod_test[[4]] <- TP2_invsimp_wilcoxon_mod

#capture.output(TP2_wilcoxon_mod_test,file="/PATH/TP2_wilcoxon_mod_test_PTB.txt")

```


##### very preterm


```R

TP2_wilcoxon_very_test <- list()

############# all case_very and control_very

TP2_alpha_case_very <- filter(TP2_cc_very_alpha_key, key %in% c("Case"))
TP2_alpha_control_very <- filter(TP2_cc_very_alpha_key, key %in% c("Control"))

### shannon

TP2_case_very_shannon <- TP2_alpha_case_very$TP2_cc_very_shannon
TP2_control_very_shannon <- TP2_alpha_control_very$TP2_cc_very_shannon
TP2_shannon_wilcoxon_very <- wilcox.test(TP2_case_very_shannon, TP2_control_very_shannon, alternative = "two.sided")
TP2_wilcoxon_very_test[[1]] <- TP2_shannon_wilcoxon_very
### richness

TP2_case_very_richness <- TP2_alpha_case_very$TP2_cc_very_richness
TP2_control_very_richness <- TP2_alpha_control_very$TP2_cc_very_richness
TP2_richness_wilcoxon_very <- wilcox.test(TP2_case_very_richness, TP2_control_very_richness, alternative = "two.sided")
TP2_wilcoxon_very_test[[2]] <- TP2_richness_wilcoxon_very

### Pielou

TP2_case_very_pielou <- TP2_alpha_case_very$TP2_cc_very_pielou
TP2_control_very_pielou <- TP2_alpha_control_very$TP2_cc_very_pielou
TP2_pielou_wilcoxon_very <- wilcox.test(TP2_case_very_pielou, TP2_control_very_pielou, alternative = "two.sided")
TP2_wilcoxon_very_test[[3]] <- TP2_pielou_wilcoxon_very

### Inverse Simpson

TP2_case_very_invsimp <- TP2_alpha_case_very$TP2_cc_very_invsimp
TP2_control_very_invsimp <- TP2_alpha_control_very$TP2_cc_very_invsimp
TP2_invsimp_wilcoxon_very <- wilcox.test(TP2_case_very_invsimp, TP2_control_very_invsimp, alternative = "two.sided")
TP2_wilcoxon_very_test[[4]] <- TP2_invsimp_wilcoxon_very

#capture.output(TP2_wilcoxon_very_test,file="/PATH/TP2_wilcoxon_very_test_PTB.txt")

```


##### Violin plots for all


```R
TP2_alpha_plot<- ggplot(TP2_cc_alpha_key, aes(x = key, y = TP2_cc_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_shannon.pdf",plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_alpha_key, aes(x = key, y = TP2_cc_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_invsimp.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_alpha_key, aes(x = key, y = TP2_cc_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_pielou.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_alpha_key, aes(x = key, y = TP2_cc_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_richness.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 



```

##### for late preterm


```R
TP2_alpha_plot<- ggplot(TP2_cc_late_alpha_key, aes(x = key, y = TP2_cc_late_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_late_shannon.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_late_alpha_key, aes(x = key, y = TP2_cc_late_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_late_invsimp.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_late_alpha_key, aes(x = key, y = TP2_cc_late_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_late_pielou.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_late_alpha_key, aes(x = key, y = TP2_cc_late_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_late_richness.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

     
```

##### for moderate


```R
TP2_alpha_plot<- ggplot(TP2_cc_mod_alpha_key, aes(x = key, y = TP2_cc_mod_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_mod_shannon.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_mod_alpha_key, aes(x = key, y = TP2_cc_mod_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_mod_invsimp.pdf",plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600)  

TP2_alpha_plot<- ggplot(TP2_cc_mod_alpha_key, aes(x = key, y = TP2_cc_mod_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_mod_pielou.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_mod_alpha_key, aes(x = key, y = TP2_cc_mod_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_mod_richness.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

```

##### for very preterm


```R
TP2_alpha_plot<- ggplot(TP2_cc_very_alpha_key, aes(x = key, y = TP2_cc_very_shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_very_shannon.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_very_alpha_key, aes(x = key, y = TP2_cc_very_invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_very_invsimp.pdf",plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600)  

TP2_alpha_plot<- ggplot(TP2_cc_very_alpha_key, aes(x = key, y = TP2_cc_very_pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_very_pielou.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

TP2_alpha_plot<- ggplot(TP2_cc_very_alpha_key, aes(x = key, y = TP2_cc_very_richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="red")

ggsave(filename = "/PATH/TP2_cc_very_richness.pdf", plot = TP2_alpha_plot, width = 12, height = 8, dpi = 600) 

```

### comparing case and control alpha diversity


```R

### add the timpoint column 

TP2_cases_alpha_2 <- TP2_alpha_case %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_cases_alpha_2 <- TP1_alpha_case %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_cases_alpha_2) <- alpha_column_names
colnames(TP1_cases_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

cases_alpha_species <- rbind(TP1_cases_alpha_2,TP2_cases_alpha_2)

### removing any values that are unique

cases_alpha_species_filtered <- cases_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

cases_alpha_species_filtered <- cases_alpha_species_filtered[order(cases_alpha_species_filtered$Studienummer),]

### In this example, the group_by function is used to group the alphata by the "Studienummer" column, 
### and the cur_group_id() function is used to assign a unique number to each distinct pair of 
### values in the "Studienummer" column.

cases_alpha_species_filtered <- cases_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### late preterm


```R


### add the timpoint column 

TP2_cases_late_alpha_2 <- TP2_alpha_case_late %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_cases_late_alpha_2 <- TP1_alpha_case_late %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_cases_late_alpha_2) <- alpha_column_names
colnames(TP1_cases_late_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

cases_late_alpha_species <- rbind(TP1_cases_late_alpha_2,TP2_cases_late_alpha_2)

### removing any values that are unique

cases_late_alpha_species_filtered <- cases_late_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

cases_late_alpha_species_filtered <- cases_late_alpha_species_filtered[order(cases_late_alpha_species_filtered$Studienummer),]

### In this example, the group_by function is used to group the alphata by the "Studienummer" column, 
### and the cur_group_id() function is used to assign a unique number to each distinct pair of 
### values in the "Studienummer" column.

cases_late_alpha_species_filtered <- cases_late_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### moderate preterm


```R


### add the timpoint column 

TP2_cases_mod_alpha_2 <- TP2_alpha_case_mod %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_cases_mod_alpha_2 <- TP1_alpha_case_mod %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_cases_mod_alpha_2) <- alpha_column_names
colnames(TP1_cases_mod_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

cases_mod_alpha_species <- rbind(TP1_cases_mod_alpha_2,TP2_cases_mod_alpha_2)

### removing any values that are unique

cases_mod_alpha_species_filtered <- cases_mod_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

cases_mod_alpha_species_filtered <- cases_mod_alpha_species_filtered[order(cases_mod_alpha_species_filtered$Studienummer),]



cases_mod_alpha_species_filtered <- cases_mod_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### very preterm


```R


### add the timpoint column 

TP2_cases_very_alpha_2 <- TP2_alpha_case_very %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_cases_very_alpha_2 <- TP1_alpha_case_very %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_cases_very_alpha_2) <- alpha_column_names
colnames(TP1_cases_very_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

cases_very_alpha_species <- rbind(TP1_cases_very_alpha_2,TP2_cases_very_alpha_2)

### removing any values that are unique

cases_very_alpha_species_filtered <- cases_very_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

cases_very_alpha_species_filtered <- cases_very_alpha_species_filtered[order(cases_very_alpha_species_filtered$Studienummer),]


cases_very_alpha_species_filtered <- cases_very_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### wilcoxon for time point comparison


```R
### all

cases_wilcoxon_test <- list()

############# all case and control

cases_TP1 <- filter(cases_alpha_species_filtered, timepoint %in% c("TP1"))
cases_alpha_TP2 <- filter(cases_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

cases_TP1_shannon <- cases_TP1$shannon
cases_TP2_shannon <- cases_alpha_TP2$shannon
cases_shannon_wilcoxon <- wilcox.test(cases_TP1_shannon, cases_TP2_shannon, alternative = "two.sided")
cases_wilcoxon_test[[1]] <- cases_shannon_wilcoxon
### richness

cases_TP1_richness <- cases_TP1$richness
cases_TP2_richness <- cases_alpha_TP2$richness
cases_richness_wilcoxon <- wilcox.test(cases_TP1_richness, cases_TP2_richness, alternative = "two.sided")
cases_wilcoxon_test[[2]] <- cases_richness_wilcoxon

### Pielou

cases_TP1_pielou <- cases_TP1$pielou
cases_TP2_pielou <- cases_alpha_TP2$pielou
cases_pielou_wilcoxon <- wilcox.test(cases_TP1_pielou, cases_TP2_pielou, alternative = "two.sided")
cases_wilcoxon_test[[3]] <- cases_pielou_wilcoxon

### Inverse Simpson

cases_TP1_invsimp <- cases_TP1$invsimp
cases_TP2_invsimp <- cases_alpha_TP2$invsimp
cases_invsimp_wilcoxon <- wilcox.test(cases_TP1_invsimp, cases_TP2_invsimp, alternative = "two.sided")
cases_wilcoxon_test[[4]] <- cases_invsimp_wilcoxon

#capture.output(cases_wilcoxon_test,file="/PATH/cases_wilcoxon_test_PTB.txt")

```

#### late preterm


```R
##### late

cases_wilcoxon_late_test <- list()

############# all TP2_late and control_late

cases_TP1_late <- filter(cases_late_alpha_species_filtered, timepoint %in% c("TP1"))
cases_alpha_TP2_late <- filter(cases_late_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

cases_TP1_late_shannon <- cases_TP1_late$shannon
cases_TP2_late_shannon <- cases_alpha_TP2_late$shannon
cases_shannon_wilcoxon_late <- wilcox.test(cases_TP1_late_shannon, cases_TP2_late_shannon, alternative = "two.sided")
cases_wilcoxon_late_test[[1]] <- cases_shannon_wilcoxon_late
### richness

cases_TP1_late_richness <- cases_TP1_late$richness
cases_TP2_late_richness <- cases_alpha_TP2_late$richness
cases_richness_wilcoxon_late <- wilcox.test(cases_TP1_late_richness, cases_TP2_late_richness, alternative = "two.sided")
cases_wilcoxon_late_test[[2]] <- cases_richness_wilcoxon_late

### Pielou

cases_TP1_late_pielou <- cases_TP1_late$pielou
cases_TP2_late_pielou <- cases_alpha_TP2_late$pielou
cases_pielou_wilcoxon_late <- wilcox.test(cases_TP1_late_pielou, cases_TP2_late_pielou, alternative = "two.sided")
cases_wilcoxon_late_test[[3]] <- cases_pielou_wilcoxon_late

### Inverse Simpson

cases_TP1_late_invsimp <- cases_TP1_late$invsimp
cases_TP2_late_invsimp <- cases_alpha_TP2_late$invsimp
cases_invsimp_wilcoxon_late <- wilcox.test(cases_TP1_late_invsimp, cases_TP2_late_invsimp, alternative = "two.sided")
cases_wilcoxon_late_test[[4]] <- cases_invsimp_wilcoxon_late

#capture.output(cases_wilcoxon_late_test,file="/PATH/cases_wilcoxon_late_test_PTB.txt")

```

#### moderate preterm


```R
##### mod

cases_wilcoxon_mod_test <- list()

############# all TP2_mod and control_mod

cases_TP1_mod <- filter(cases_mod_alpha_species_filtered, timepoint %in% c("TP1"))
cases_alpha_TP2_mod <- filter(cases_mod_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

cases_TP1_mod_shannon <- cases_TP1_mod$shannon
cases_TP2_mod_shannon <- cases_alpha_TP2_mod$shannon
cases_shannon_wilcoxon_mod <- wilcox.test(cases_TP1_mod_shannon, cases_TP2_mod_shannon, alternative = "two.sided")
cases_wilcoxon_mod_test[[1]] <- cases_shannon_wilcoxon_mod
### richness

cases_TP1_mod_richness <- cases_TP1_mod$richness
cases_TP2_mod_richness <- cases_alpha_TP2_mod$richness
cases_richness_wilcoxon_mod <- wilcox.test(cases_TP1_mod_richness, cases_TP2_mod_richness, alternative = "two.sided")
cases_wilcoxon_mod_test[[2]] <- cases_richness_wilcoxon_mod

### Pielou

cases_TP1_mod_pielou <- cases_TP1_mod$pielou
cases_TP2_mod_pielou <- cases_alpha_TP2_mod$pielou
cases_pielou_wilcoxon_mod <- wilcox.test(cases_TP1_mod_pielou, cases_TP2_mod_pielou, alternative = "two.sided")
cases_wilcoxon_mod_test[[3]] <- cases_pielou_wilcoxon_mod

### Inverse Simpson

cases_TP1_mod_invsimp <- cases_TP1_mod$invsimp
cases_TP2_mod_invsimp <- cases_alpha_TP2_mod$invsimp
cases_invsimp_wilcoxon_mod <- wilcox.test(cases_TP1_mod_invsimp, cases_TP2_mod_invsimp, alternative = "two.sided")
cases_wilcoxon_mod_test[[4]] <- cases_invsimp_wilcoxon_mod

#capture.output(cases_wilcoxon_mod_test,file="/PATH/cases_wilcoxon_mod_test_PTB.txt")

```

#### very preterm


```R
##### very

cases_wilcoxon_very_test <- list()

############# all TP2_very and control_very

cases_TP1_very <- filter(cases_very_alpha_species_filtered, timepoint %in% c("TP1"))
cases_alpha_TP2_very <- filter(cases_very_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

cases_TP1_very_shannon <- cases_TP1_very$shannon
cases_TP2_very_shannon <- cases_alpha_TP2_very$shannon
cases_shannon_wilcoxon_very <- wilcox.test(cases_TP1_very_shannon, cases_TP2_very_shannon, alternative = "two.sided")
cases_wilcoxon_very_test[[1]] <- cases_shannon_wilcoxon_very
### richness

cases_TP1_very_richness <- cases_TP1_very$richness
cases_TP2_very_richness <- cases_alpha_TP2_very$richness
cases_richness_wilcoxon_very <- wilcox.test(cases_TP1_very_richness, cases_TP2_very_richness, alternative = "two.sided")
cases_wilcoxon_very_test[[2]] <- cases_richness_wilcoxon_very

### Pielou

cases_TP1_very_pielou <- cases_TP1_very$pielou
cases_TP2_very_pielou <- cases_alpha_TP2_very$pielou
cases_pielou_wilcoxon_very <- wilcox.test(cases_TP1_very_pielou, cases_TP2_very_pielou, alternative = "two.sided")
cases_wilcoxon_very_test[[3]] <- cases_pielou_wilcoxon_very

### Inverse Simpson

cases_TP1_very_invsimp <- cases_TP1_very$invsimp
cases_TP2_very_invsimp <- cases_alpha_TP2_very$invsimp
cases_invsimp_wilcoxon_very <- wilcox.test(cases_TP1_very_invsimp, cases_TP2_very_invsimp, alternative = "two.sided")
cases_wilcoxon_very_test[[4]] <- cases_invsimp_wilcoxon_very

#capture.output(cases_wilcoxon_very_test,file="/PATH/cases_wilcoxon_very_test_PTB.txt")

```


##### Plots


```R
cases_alpha_shannon_plot<- ggplot(cases_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample Cases TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/TP_case_alpha_plot_shannon.pdf", plot = cases_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

cases_alpha_pielou_plot<- ggplot(cases_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample Cases TP1 and TP2 for Pielou")



ggsave(filename = "/PATH/TP_case_alpha_plot_pielou.pdf", plot = cases_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

cases_alpha_richness_plot<- ggplot(cases_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample Cases TP1 and TP2 for Richness")


ggsave(filename = "/PATH/TP_case_alpha_plot_richness.pdf", plot = cases_alpha_richness_plot, width = 12, height = 8, dpi = 600)

cases_alpha_invsimp_plot<- ggplot(cases_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample Cases TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/TP_case_alpha_plot_invsimp.pdf", plot = cases_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### late preterm


```R
cases_alpha_shannon_plot<- ggplot(cases_late_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample late Cases TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_case_late_alpha_plot_shannon.pdf", plot = cases_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

cases_alpha_pielou_plot<- ggplot(cases_late_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late Cases TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_case_late_alpha_plot_pielou.pdf", plot = cases_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

cases_alpha_richness_plot<- ggplot(cases_late_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late Cases TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_case_late_alpha_plot_richness.pdf", plot = cases_alpha_richness_plot, width = 12, height = 8, dpi = 600)

cases_alpha_invsimp_plot<- ggplot(cases_late_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late Cases TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_case_late_alpha_plot_invsimp.pdf", plot = cases_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### moderate preterm


```R
cases_alpha_shannon_plot<- ggplot(cases_mod_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample mod Cases TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_case_mod_alpha_plot_shannon.pdf", plot = cases_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

cases_alpha_pielou_plot<- ggplot(cases_mod_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod Cases TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_case_mod_alpha_plot_pielou.pdf", plot = cases_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

cases_alpha_richness_plot<- ggplot(cases_mod_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod Cases TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_case_mod_alpha_plot_richness.pdf", plot = cases_alpha_richness_plot, width = 12, height = 8, dpi = 600)

cases_alpha_invsimp_plot<- ggplot(cases_mod_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod Cases TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_case_mod_alpha_plot_invsimp.pdf", plot = cases_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### very preterm


```R
cases_alpha_shannon_plot<- ggplot(cases_very_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample very Cases TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_case_very_alpha_plot_shannon.pdf", plot = cases_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

cases_alpha_pielou_plot<- ggplot(cases_very_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very Cases TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_case_very_alpha_plot_pielou.pdf", plot = cases_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

cases_alpha_richness_plot<- ggplot(cases_very_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very Cases TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_case_very_alpha_plot_richness.pdf", plot = cases_alpha_richness_plot, width = 12, height = 8, dpi = 600)

cases_alpha_invsimp_plot<- ggplot(cases_very_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very Cases TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_case_very_alpha_plot_invsimp.pdf", plot = cases_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

### controls time point comparison


```R

### add the timpoint column 

TP2_controls_alpha_2 <- TP2_alpha_control %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_controls_alpha_2 <- TP1_alpha_control %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_controls_alpha_2) <- alpha_column_names
colnames(TP1_controls_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

controls_alpha_species <- rbind(TP1_controls_alpha_2,TP2_controls_alpha_2)

### removing any values that are unique

controls_alpha_species_filtered <- controls_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

controls_alpha_species_filtered <- controls_alpha_species_filtered[order(controls_alpha_species_filtered$Studienummer),]


controls_alpha_species_filtered <- controls_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### late preterm


```R
### add the timpoint column 

TP2_controls_late_alpha_2 <- TP2_alpha_control_late %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_controls_late_alpha_2 <- TP1_alpha_control_late %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_controls_late_alpha_2) <- alpha_column_names
colnames(TP1_controls_late_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

controls_late_alpha_species <- rbind(TP1_controls_late_alpha_2,TP2_controls_late_alpha_2)

### removing any values that are unique

controls_late_alpha_species_filtered <- controls_late_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

controls_late_alpha_species_filtered <- controls_late_alpha_species_filtered[order(controls_late_alpha_species_filtered$Studienummer),]


controls_late_alpha_species_filtered <- controls_late_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### moderate preterm


```R
### add the timpoint column 

TP2_controls_mod_alpha_2 <- TP2_alpha_control_mod %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_controls_mod_alpha_2 <- TP1_alpha_control_mod %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_controls_mod_alpha_2) <- alpha_column_names
colnames(TP1_controls_mod_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

controls_mod_alpha_species <- rbind(TP1_controls_mod_alpha_2,TP2_controls_mod_alpha_2)

### removing any values that are unique

controls_mod_alpha_species_filtered <- controls_mod_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

controls_mod_alpha_species_filtered <- controls_mod_alpha_species_filtered[order(controls_mod_alpha_species_filtered$Studienummer),]


controls_mod_alpha_species_filtered <- controls_mod_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### very preterm


```R

### add the timpoint column 

TP2_controls_very_alpha_2 <- TP2_alpha_control_very %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_controls_very_alpha_2 <- TP1_alpha_control_very %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

alpha_column_names <- c("timepoint","Studienummer","key","richness","shannon","pielou","invsimp")
colnames(TP2_controls_very_alpha_2) <- alpha_column_names
colnames(TP1_controls_very_alpha_2) <- alpha_column_names

### combining the two alphataframes TP1 and TP2

controls_very_alpha_species <- rbind(TP1_controls_very_alpha_2,TP2_controls_very_alpha_2)

### removing any values that are unique

controls_very_alpha_species_filtered <- controls_very_alpha_species %>%
  group_by(Studienummer) %>%
  filter(n() > 1) %>%
  ungroup()

### ordering the alphataframe by Studienummer

controls_very_alpha_species_filtered <- controls_very_alpha_species_filtered[order(controls_very_alpha_species_filtered$Studienummer),]


controls_very_alpha_species_filtered <- controls_very_alpha_species_filtered %>%
  group_by(Studienummer) %>%
  mutate(paired = cur_group_id())
```

#### Wilcoxon test


```R
### all


controls_wilcoxon_test <- list()

############# all case and control

controls_TP1 <- filter(controls_alpha_species_filtered, timepoint %in% c("TP1"))
controls_alpha_TP2 <- filter(controls_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

controls_TP1_shannon <- controls_TP1$shannon
controls_TP2_shannon <- controls_alpha_TP2$shannon
controls_shannon_wilcoxon <- wilcox.test(controls_TP1_shannon, controls_TP2_shannon, alternative = "two.sided")
controls_wilcoxon_test[[1]] <- controls_shannon_wilcoxon
### richness

controls_TP1_richness <- controls_TP1$richness
controls_TP2_richness <- controls_alpha_TP2$richness
controls_richness_wilcoxon <- wilcox.test(controls_TP1_richness, controls_TP2_richness, alternative = "two.sided")
controls_wilcoxon_test[[2]] <- controls_richness_wilcoxon

### Pielou

controls_TP1_pielou <- controls_TP1$pielou
controls_TP2_pielou <- controls_alpha_TP2$pielou
controls_pielou_wilcoxon <- wilcox.test(controls_TP1_pielou, controls_TP2_pielou, alternative = "two.sided")
controls_wilcoxon_test[[3]] <- controls_pielou_wilcoxon

### Inverse Simpson

controls_TP1_invsimp <- controls_TP1$invsimp
controls_TP2_invsimp <- controls_alpha_TP2$invsimp
controls_invsimp_wilcoxon <- wilcox.test(controls_TP1_invsimp, controls_TP2_invsimp, alternative = "two.sided")
controls_wilcoxon_test[[4]] <- controls_invsimp_wilcoxon

#capture.output(controls_wilcoxon_test,file="/PATH/controls_wilcoxon_test_PTB.txt")

```

#### later preterm


```R
##### late

controls_wilcoxon_late_test <- list()



controls_TP1_late <- filter(controls_late_alpha_species_filtered, timepoint %in% c("TP1"))
controls_alpha_TP2_late <- filter(controls_late_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

controls_TP1_late_shannon <- controls_TP1_late$shannon
controls_TP2_late_shannon <- controls_alpha_TP2_late$shannon
controls_shannon_wilcoxon_late <- wilcox.test(controls_TP1_late_shannon, controls_TP2_late_shannon, alternative = "two.sided")
controls_wilcoxon_late_test[[1]] <- controls_shannon_wilcoxon_late
### richness

controls_TP1_late_richness <- controls_TP1_late$richness
controls_TP2_late_richness <- controls_alpha_TP2_late$richness
controls_richness_wilcoxon_late <- wilcox.test(controls_TP1_late_richness, controls_TP2_late_richness, alternative = "two.sided")
controls_wilcoxon_late_test[[2]] <- controls_richness_wilcoxon_late

### Pielou

controls_TP1_late_pielou <- controls_TP1_late$pielou
controls_TP2_late_pielou <- controls_alpha_TP2_late$pielou
controls_pielou_wilcoxon_late <- wilcox.test(controls_TP1_late_pielou, controls_TP2_late_pielou, alternative = "two.sided")
controls_wilcoxon_late_test[[3]] <- controls_pielou_wilcoxon_late

### Inverse Simpson

controls_TP1_late_invsimp <- controls_TP1_late$invsimp
controls_TP2_late_invsimp <- controls_alpha_TP2_late$invsimp
controls_invsimp_wilcoxon_late <- wilcox.test(controls_TP1_late_invsimp, controls_TP2_late_invsimp, alternative = "two.sided")
controls_wilcoxon_late_test[[4]] <- controls_invsimp_wilcoxon_late

#capture.output(controls_wilcoxon_late_test,file="/PATH/controls_wilcoxon_late_test_PTB.txt")

```

#### moderate preterm


```R
 ##### mod

controls_wilcoxon_mod_test <- list()



controls_TP1_mod <- filter(controls_mod_alpha_species_filtered, timepoint %in% c("TP1"))
controls_alpha_TP2_mod <- filter(controls_mod_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

controls_TP1_mod_shannon <- controls_TP1_mod$shannon
controls_TP2_mod_shannon <- controls_alpha_TP2_mod$shannon
controls_shannon_wilcoxon_mod <- wilcox.test(controls_TP1_mod_shannon, controls_TP2_mod_shannon, alternative = "two.sided")
controls_wilcoxon_mod_test[[1]] <- controls_shannon_wilcoxon_mod
### richness

controls_TP1_mod_richness <- controls_TP1_mod$richness
controls_TP2_mod_richness <- controls_alpha_TP2_mod$richness
controls_richness_wilcoxon_mod <- wilcox.test(controls_TP1_mod_richness, controls_TP2_mod_richness, alternative = "two.sided")
controls_wilcoxon_mod_test[[2]] <- controls_richness_wilcoxon_mod

### Pielou

controls_TP1_mod_pielou <- controls_TP1_mod$pielou
controls_TP2_mod_pielou <- controls_alpha_TP2_mod$pielou
controls_pielou_wilcoxon_mod <- wilcox.test(controls_TP1_mod_pielou, controls_TP2_mod_pielou, alternative = "two.sided")
controls_wilcoxon_mod_test[[3]] <- controls_pielou_wilcoxon_mod

### Inverse Simpson

controls_TP1_mod_invsimp <- controls_TP1_mod$invsimp
controls_TP2_mod_invsimp <- controls_alpha_TP2_mod$invsimp
controls_invsimp_wilcoxon_mod <- wilcox.test(controls_TP1_mod_invsimp, controls_TP2_mod_invsimp, alternative = "two.sided")
controls_wilcoxon_mod_test[[4]] <- controls_invsimp_wilcoxon_mod

#capture.output(controls_wilcoxon_mod_test,file="/PATH/controls_wilcoxon_mod_test_PTB.txt")

```



#### very preterm


```R
##### very

controls_wilcoxon_very_test <- list()



controls_TP1_very <- filter(controls_very_alpha_species_filtered, timepoint %in% c("TP1"))
controls_alpha_TP2_very <- filter(controls_very_alpha_species_filtered, timepoint %in% c("TP2"))

### shannon

controls_TP1_very_shannon <- controls_TP1_very$shannon
controls_TP2_very_shannon <- controls_alpha_TP2_very$shannon
controls_shannon_wilcoxon_very <- wilcox.test(controls_TP1_very_shannon, controls_TP2_very_shannon, alternative = "two.sided")
controls_wilcoxon_very_test[[1]] <- controls_shannon_wilcoxon_very
### richness

controls_TP1_very_richness <- controls_TP1_very$richness
controls_TP2_very_richness <- controls_alpha_TP2_very$richness
controls_richness_wilcoxon_very <- wilcox.test(controls_TP1_very_richness, controls_TP2_very_richness, alternative = "two.sided")
controls_wilcoxon_very_test[[2]] <- controls_richness_wilcoxon_very

### Pielou

controls_TP1_very_pielou <- controls_TP1_very$pielou
controls_TP2_very_pielou <- controls_alpha_TP2_very$pielou
controls_pielou_wilcoxon_very <- wilcox.test(controls_TP1_very_pielou, controls_TP2_very_pielou, alternative = "two.sided")
controls_wilcoxon_very_test[[3]] <- controls_pielou_wilcoxon_very

### Inverse Simpson

controls_TP1_very_invsimp <- controls_TP1_very$invsimp
controls_TP2_very_invsimp <- controls_alpha_TP2_very$invsimp
controls_invsimp_wilcoxon_very <- wilcox.test(controls_TP1_very_invsimp, controls_TP2_very_invsimp, alternative = "two.sided")
controls_wilcoxon_very_test[[4]] <- controls_invsimp_wilcoxon_very

#capture.output(controls_wilcoxon_very_test,file="/PATH/controls_wilcoxon_very_test_PTB.txt")

```


#### violin plots


```R
controls_alpha_shannon_plot<- ggplot(controls_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample controls TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/TP_control_alpha_plot_shannon.pdf", plot = controls_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

controls_alpha_pielou_plot<- ggplot(controls_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample controls TP1 and TP2 for Pielou")



ggsave(filename = "/PATH/TP_control_alpha_plot_pielou.pdf", plot = controls_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

controls_alpha_richness_plot<- ggplot(controls_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample controls TP1 and TP2 for Richness")


ggsave(filename = "/PATH/TP_control_alpha_plot_richness.pdf", plot = controls_alpha_richness_plot, width = 12, height = 8, dpi = 600)

controls_alpha_invsimp_plot<- ggplot(controls_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample controls TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/TP_control_alpha_plot_invsimp.pdf", plot = controls_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### late preterm


```R
controls_alpha_shannon_plot<- ggplot(controls_late_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample late controls TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_control_late_alpha_plot_shannon.pdf", plot = controls_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

controls_alpha_pielou_plot<- ggplot(controls_late_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late controls TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_control_late_alpha_plot_pielou.pdf", plot = controls_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

controls_alpha_richness_plot<- ggplot(controls_late_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late controls TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_control_late_alpha_plot_richness.pdf", plot = controls_alpha_richness_plot, width = 12, height = 8, dpi = 600)

controls_alpha_invsimp_plot<- ggplot(controls_late_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample late controls TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_control_late_alpha_plot_invsimp.pdf", plot = controls_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### moderate preterm


```R
controls_alpha_shannon_plot<- ggplot(controls_mod_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample mod controls TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_control_mod_alpha_plot_shannon.pdf", plot = controls_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

controls_alpha_pielou_plot<- ggplot(controls_mod_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod controls TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_control_mod_alpha_plot_pielou.pdf", plot = controls_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

controls_alpha_richness_plot<- ggplot(controls_mod_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod controls TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_control_mod_alpha_plot_richness.pdf", plot = controls_alpha_richness_plot, width = 12, height = 8, dpi = 600)

controls_alpha_invsimp_plot<- ggplot(controls_mod_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample mod controls TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_control_mod_alpha_plot_invsimp.pdf", plot = controls_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

#### very preterm


```R
controls_alpha_shannon_plot<- ggplot(controls_very_alpha_species_filtered, aes(x = timepoint, y = shannon)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue") +
  labs(title = "Corresponding sample very controls TP1 and TP2 for Shannon diversity")


ggsave(filename = "/PATH/F_TP_control_very_alpha_plot_shannon.pdf", plot = controls_alpha_shannon_plot, width = 12, height = 8, dpi = 600)

controls_alpha_pielou_plot<- ggplot(controls_very_alpha_species_filtered, aes(x = timepoint, y = pielou)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very controls TP1 and TP2 for Pielou")


ggsave(filename = "/PATH/F_TP_control_very_alpha_plot_pielou.pdf", plot = controls_alpha_pielou_plot, width = 12, height = 8, dpi = 600)

controls_alpha_richness_plot<- ggplot(controls_very_alpha_species_filtered, aes(x = timepoint, y = richness)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very controls TP1 and TP2 for Richness")


ggsave(filename = "/PATH/F_TP_control_very_alpha_plot_richness.pdf", plot = controls_alpha_richness_plot, width = 12, height = 8, dpi = 600)

controls_alpha_invsimp_plot<- ggplot(controls_very_alpha_species_filtered, aes(x = timepoint, y = invsimp)) +
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="blue")+
  labs(title = "Corresponding sample very controls TP1 and TP2 for Inverse Simpson")


ggsave(filename = "/PATH/F_TP_control_very_alpha_plot_invsimp.pdf", plot = controls_alpha_invsimp_plot, width = 12, height = 8, dpi = 600)


```

# PERMANOVA/ Adonis

### For TP1. We're going to start by looping through it.


```R
## removing the troublesome variables from metadata which include compl_grops

metadata_Q1_cleaned <- metadata_Q1_cleaned[, !(names(metadata_Q1_cleaned) %in% c("stress_sum_score_Q1","Compl_grops"))]
metadata_Q1_cleaned[,c(4,140)] <- sapply(metadata_Q1_cleaned[,c(4,140)],function(x) as.character(x))
metadata_Q1_cleaned <- metadata_Q1_cleaned[,-c(1)]
TP1_adonis_clr <- merge(TP1_case_control_0,metadata_Q1_cleaned,by="Studienummer",all.x=TRUE,all.y=FALSE)

adonis_file_2 <- TP1_adonis_clr[,c(5:1869,1)]
adonis_species_test <- adonis_file_2[,1:1865]

```


```R
## ALready ran on interactive

adonis_TP1_list_clr <- list()
m <- 0
for (i in 1870:ncol(TP1_adonis_clr)) {
    m <- m+1
    adonis_file <- TP1_adonis_clr[,c(5:1869,i)] ### the section that only contains the species + the metadata we want to look into
    adonis_file <- na.omit(adonis_file)
    adonis_species <- adonis_file[,1:1865] ## the species only
    name_of_factor <- colnames(adonis_file[1866]) ## the metadatavariable we want to look into
    form2 <- as.formula(paste("adonis_species",name_of_factor,sep="~"))
    adonis_TP1_list_clr[[m]] <- adonis2(form2, data = adonis_file, na.action = na.omit, permutations = 999,method="euclidean")
    print(paste(m,name_of_factor,sep=", "))
}

#capture.output(adonis_TP1_list_clr,file="/PATH/adonis_TP1_list_clr_eclidean.txt")
#saveRDS(adonis_TP1_list_clr,file="/PATH/adonis_TP1_list_clr_eclidean.rds")


```

 
```R
TP1_adonis_R2_P <- NULL

for (i in 1:138) {
 df <- data.frame(adonis_TP1_list_clr[i])
 additional_row <- c(rownames(df)[1],df[1,3],df[1,5])  
 TP1_adonis_R2_P <- rbind(TP1_adonis_R2_P,additional_row)
}
TP1_adonis_R2_P <- data.frame(TP1_adonis_R2_P)
colnames(TP1_adonis_R2_P) <- c("variable","R2","p_val")
TP1_adonis_R2_P[,c(2,3)] <- sapply(TP1_adonis_R2_P[,c(2,3)], function(x) as.numeric(x))
TP1_adonis_R2_P <- data.frame(TP1_adonis_R2_P)
TP1_adonis_R2_P <- na.omit(TP1_adonis_R2_P)                                  
TP1_adonis_R2_P <- TP1_adonis_R2_P[order(TP1_adonis_R2_P$p_val,decreasing=FALSE),]

write.csv(TP1_adonis_R2_P,"/PATH/TP1_adonis_R2_P_2.csv")
                                  

```


### For TP2. Will loop through the metadata


```R
metadata_Q2_cleaned_2 <- metadata_Q2_cleaned[,-c(1,7)]

TP2_adonis_clr <- merge(TP2_case_control_0,metadata_Q2_cleaned_2,by="Studienummer",all.x=TRUE,all.y=FALSE)

 
## For the second questionnaire
adonis_TP2_list_clr <- list()
m <- 0
for (i in 1870:ncol(TP2_adonis_clr)) {
    m <- m+1
    adonis_file <- TP2_adonis_clr[,c(5:1869,i)] ### the section that only contains the species + the metadata we want to look into
    adonis_file <- na.omit(adonis_file)
    adonis_species <- adonis_file[,1:1865] ## the species only
    name_of_factor <- colnames(adonis_file[1866]) ## the metadatavariable we want to look into
    form2 <- as.formula(paste("adonis_species",name_of_factor,sep="~"))
    adonis_TP2_list_clr[[m]] <- adonis2(form2, data = adonis_file, na.action = na.omit, permutations = 999,method="euclidean")
    print(paste(m,name_of_factor,sep=", "))
}

#capture.output(adonis_TP2_list_clr,file="/PATH/adonis_TP2_list_clr_eclidean.txt")
#saveRDS(adonis_TP2_list_clr,file="/PATH/adonis_TP2_list_clr_eclidean.rds")

```


```R
TP2_adonis_R2_P <- NULL

for (i in 1:63) {
 df <- data.frame(adonis_TP2_list_clr[i])
 additional_row <- c(rownames(df)[1],df[1,3],df[1,5])  
 TP2_adonis_R2_P <- rbind(TP2_adonis_R2_P,additional_row)
}
TP2_adonis_R2_P <- data.frame(TP2_adonis_R2_P)
colnames(TP2_adonis_R2_P) <- c("variable","R2","p_val")
TP2_adonis_R2_P[,c(2,3)] <- sapply(TP2_adonis_R2_P[,c(2,3)], function(x) as.numeric(x))
TP2_adonis_R2_P <- data.frame(TP2_adonis_R2_P)
TP2_adonis_R2_P <- na.omit(TP2_adonis_R2_P)                                  
TP2_adonis_R2_P <- TP2_adonis_R2_P[order(TP2_adonis_R2_P$p_val,decreasing=FALSE),]

write.csv(TP2_adonis_R2_P,"/PATH/TP2_adonis_R2_P_2.csv")
                                  
TP2_adonis_R2_P

```



# ANCOMBC 

### Building the phyloseq object and ANCOMBC for TP1


```R
# for the otu matrix the row names are going to be the studienummer
TP1_species_for_ancombc <- TP1_case_control_non_clr[,c(5:ncol(TP1_case_control_non_clr))]
rownames(TP1_species_for_ancombc) <- TP1_case_control_non_clr[,1]

## prepping the OTU section
TP1_OTU_matrix <- as.matrix(TP1_species_for_ancombc)


### for the metadata file the rows names are the studienummers and the other columns are everything you 
## want to came in 
TP1_SMP <- TP1_case_control_non_clr[,c(1,3)] ## this is studienummer and keys
rownames(TP1_SMP) <- TP1_SMP[,c(1)] ## row is going to be studienummer

### now let's add the rows that we're interested in from Permanova
TP1_SMP <- merge(TP1_SMP,metadata_Q1_cleaned[,c("Studienummer","PUQUE_rating_Q1","Prev_PTB",
                                       "Prev_RPL","age_3_groups","BMI_3_groups")],by.x=0,by.y="Studienummer"
                ,all.x = TRUE, all.y = FALSE)

rownames(TP1_SMP) <- TP1_SMP$Row.names
TP1_SMP <- TP1_SMP[,-c(1,2)]


####### Now, putting them in the correct format for the phyloseq object
SMP <- sample_data(TP1_SMP)
OTU <- otu_table(TP1_OTU_matrix, taxa_are_rows = FALSE) 
      
TP1_phylo_object <- phyloseq(SMP,OTU)

### Now, building the tree summarized experiment from the phyloseq 
## object that I made in the previous sections of this file.

TP1_tse = mia::makeTreeSummarizedExperimentFromPhyloseq(TP1_phylo_object)

## running ancombc
output_ancombc2_TP1 = ancombc2(data = TP1_tse, assay_name = "counts",
                  fix_formula = "key", verbose = TRUE, p_adj_method="BH",group = NULL, struc_zero = FALSE, tax_level="Species",neg_lb = FALSE,
                  alpha = 0.1)

#saveRDS(output_ancombc2_TP1$res,"/PATH/output_ancombc2_TP1_key.rds")

```

### Building the phyloseq object and ANCOMBC for TP2


```R
# for the otu matrix the row names are going to be the studienummer
TP2_species_for_ancombc <- TP2_case_control_non_clr[,c(5:ncol(TP2_case_control_non_clr))]
rownames(TP2_species_for_ancombc) <- TP2_case_control_non_clr[,1]

## prepping the OTU section
TP2_OTU_matrix <- as.matrix(TP2_species_for_ancombc)


### for the metadata file the rows names are the studienummers and the other columns are everything you 
## want to came in 
TP2_SMP <- TP2_case_control_non_clr[,c(1,3)] ## this is studienummer and keys
rownames(TP2_SMP) <- TP2_SMP[,c(1)] ## row is going to be studienummer

### now let's add the rows that we're interested in from Permanova
TP2_SMP <- merge(TP2_SMP,metadata_Q2_cleaned[,c("Studienummer","Q2_neuro_medication"
                                       ,"BMI_prior")],by.x=0,by.y="Studienummer"
                ,all.x = TRUE, all.y = FALSE)

rownames(TP2_SMP) <- TP2_SMP$Row.names
TP2_SMP <- TP2_SMP[,-c(1,2)]


####### Now, putting them in the correct format for the phyloseq object
SMP <- sample_data(TP2_SMP)
OTU <- otu_table(TP2_OTU_matrix, taxa_are_rows = FALSE) 
      
TP2_phylo_object <- phyloseq(SMP,OTU)

### Now, building the tree summarized experiment from the phyloseq 
## object that I made in the previous sections of this file.

TP2_tse = mia::makeTreeSummarizedExperimentFromPhyloseq(TP2_phylo_object)

## running ancombc
output_ancombc2_TP2 = ancombc2(data = TP2_tse, assay_name = "counts",
                  fix_formula = "key", verbose = TRUE, p_adj_method="BH",group = NULL, struc_zero = FALSE, tax_level="Species",neg_lb = FALSE,
                  alpha = 0.1)

saveRDS(output_ancombc2_TP2$res,"/PATH/output_ancombc2_TP2_key.rds")


```

## time point comparison


```R
### Here separating into cases and then combining by timepoints

TP2_cases <- filter(TP2_case_control_non_clr, key %in% c("Case"))
TP1_cases <- filter(TP1_case_control_non_clr, key %in% c("Case"))

TP1_cases_2 <- TP1_cases[,-c(1,2,3,4)]
rownames(TP1_cases_2) <- TP1_cases[,c(1)]

TP2_cases_2 <- TP2_cases[,-c(1,2,3,4)]
rownames(TP2_cases_2) <- TP2_cases[,c(1)]


TP2_cases_2 <- TP2_cases_2 %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_cases_2 <- TP1_cases_2 %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

TP1_cases_t <- t(TP1_cases_2)
TP1_cases_t <- data.frame(TP1_cases_t)

TP2_cases_t <- t(TP2_cases_2)
TP2_cases_t <- data.frame(TP2_cases_t)

TP1_TP2_cases_t <- merge(TP1_cases_t,TP2_cases_t,by=0,all=TRUE)
TP1_TP2_cases_t[is.na(TP1_TP2_cases_t)] <- 0

TP1_TP2_cases <- t(TP1_TP2_cases_t)
TP1_TP2_cases <- data.frame(TP1_TP2_cases)

colnames(TP1_TP2_cases) <- TP1_TP2_cases[1,]
TP1_TP2_cases <- TP1_TP2_cases[-c(1),]

## prepping the phyloseq
cases_timpoints <- TP1_TP2_cases[,c("timepoint","UNCLASSIFIED")]
cases_species <- TP1_TP2_cases[ , -which(names(TP1_TP2_cases) %in% c("timepoint"))]
cases_species[,c(1:ncol(cases_species))] <- sapply(cases_species[,c(1:ncol(cases_species))], 
                                                   function(x) as.numeric(x))
cases_OTU <- as.matrix(cases_species)

## building phyloseq                                                   
                                                   
SMP <- sample_data(cases_timpoints)
OTU <- otu_table(cases_OTU, taxa_are_rows = FALSE) 
      
cases_object <- phyloseq(SMP,OTU)

### Now, building the tree summarized experiment from the phyloseq 
## object that I made in the previous sections of this file.

cases_tse = mia::makeTreeSummarizedExperimentFromPhyloseq(cases_object)
                                                   
## Running the ANCOMBC2
                                                   
output_ancombc2_cases = ancombc2(data = cases_tse, assay_name = "counts",
                  fix_formula = "timepoint", verbose = TRUE, p_adj_method="BH",
                                 group = NULL, struc_zero = FALSE, tax_level="Species",neg_lb = FALSE,
                  alpha = 0.1)  
                                                   
                                                   
saveRDS(output_ancombc2_cases$res,"/PATH/output_ancombc2_cases.rds")
                                                   
```



#### for controls TP1 and TP2 ANCOMBC2 comparison


```R
### Here separating into controls and then combining by timepoints

TP2_controls <- filter(TP2_case_control_non_clr, key %in% c("Control"))
TP1_controls <- filter(TP1_case_control_non_clr, key %in% c("Control"))

TP1_controls_2 <- TP1_controls[,-c(1,2,3,4)]
rownames(TP1_controls_2) <- TP1_controls[,c(1)]

TP2_controls_2 <- TP2_controls[,-c(1,2,3,4)]
rownames(TP2_controls_2) <- TP2_controls[,c(1)]


TP2_controls_2 <- TP2_controls_2 %>%
  mutate(timepoint = 'TP2') %>%
    relocate(timepoint)
TP1_controls_2 <- TP1_controls_2 %>%
  mutate(timepoint = 'TP1') %>%
    relocate(timepoint)

TP1_controls_t <- t(TP1_controls_2)
TP1_controls_t <- data.frame(TP1_controls_t)

TP2_controls_t <- t(TP2_controls_2)
TP2_controls_t <- data.frame(TP2_controls_t)

TP1_TP2_controls_t <- merge(TP1_controls_t,TP2_controls_t,by=0,all=TRUE)
TP1_TP2_controls_t[is.na(TP1_TP2_controls_t)] <- 0

TP1_TP2_controls <- t(TP1_TP2_controls_t)
TP1_TP2_controls <- data.frame(TP1_TP2_controls)

colnames(TP1_TP2_controls) <- TP1_TP2_controls[1,]
TP1_TP2_controls <- TP1_TP2_controls[-c(1),]

## prepping the phyloseq
controls_timpoints <- TP1_TP2_controls[,c("timepoint","UNCLASSIFIED")]
controls_species <- TP1_TP2_controls[ , -which(names(TP1_TP2_controls) %in% c("timepoint"))]
controls_species[,c(1:ncol(controls_species))] <- sapply(controls_species[,c(1:ncol(controls_species))], 
                                                   function(x) as.numeric(x))
controls_OTU <- as.matrix(controls_species)

## building phyloseq                                                   
                                                   
SMP <- sample_data(controls_timpoints)
OTU <- otu_table(controls_OTU, taxa_are_rows = FALSE) 
      
controls_object <- phyloseq(SMP,OTU)

### Now, building the tree summarized experiment from the phyloseq 
## object that I made in the previous sections of this file.

controls_tse = mia::makeTreeSummarizedExperimentFromPhyloseq(controls_object)
                                                   
## Running the ANCOMBC2
                                                   
output_ancombc2_controls = ancombc2(data = controls_tse, assay_name = "counts",
                  fix_formula = "timepoint", verbose = TRUE, p_adj_method="BH",
                                 group = NULL, struc_zero = FALSE, tax_level="Species",neg_lb = FALSE,
                  alpha = 0.1)         
                                                         
                                                         
saveRDS(output_ancombc2_controls$res,"/PATH/output_ancombc2_controls.rds")
                                                         
```


### making volcano plots for case and control ancombc2


```R
### making the volcano plot for time difference
F_control_ancombc_DA <- output_ancombc2_controls$res
F_control_ancombc_DA <- F_control_ancombc_DA[-c(1,2,5,6),]
taxon_names <- F_control_ancombc_DA[,c(1)]
taxon_names <- data.frame(taxon_names)
tax_IDs <- lapply(taxon_names, function(x) str_split(x,".s__"))
tax_IDs <- data.frame(tax_IDs)
tax_IDs_t <- t(tax_IDs)
tax_IDs_t <- data.frame(tax_IDs_t)   
F_control_ancombc_DA$species <- tax_IDs_t$X2

volcano_plot_control_time_comparison <- ggplot(F_control_ancombc_DA, aes(x = lfc_timepointTP2, y = -log10(p_timepointTP2), color = p_timepointTP2 < 0.05)) +
     geom_point() +
     scale_color_manual(values = c("red", "black")) +
     labs(title = "Volcano Plot sample controls T2 to T1 comparison",
          x = "Log Fold Change",
          y = "-log10(p-value)") +
  #   theme(legend.position = "none") +  
     geom_text(data = subset(F_control_ancombc_DA, p_timepointTP2 < 0.05),
                   aes(label = species), nudge_y = 0, nudge_x = -1, color = "blue",size=4)
      

ggsave(filename = "/PATH/volcano_plot_control_time_comparison_sample.pdf", plot = volcano_plot_control_time_comparison, width = 12, height = 8, dpi = 600) 

```




```R
### making the volcano plot for time difference
F_case_ancombc_DA <- output_ancombc2_cases$res
F_case_ancombc_DA <- F_case_ancombc_DA[-c(1,2,5,6),]
taxon_names <- F_case_ancombc_DA[,c(1)]
taxon_names <- data.frame(taxon_names)
tax_IDs <- lapply(taxon_names, function(x) str_split(x,".s__"))
tax_IDs <- data.frame(tax_IDs)
tax_IDs_t <- t(tax_IDs)
tax_IDs_t <- data.frame(tax_IDs_t)   
F_case_ancombc_DA$species <- tax_IDs_t$X2

volcano_plot_case_time_comparison <- ggplot(F_case_ancombc_DA, aes(x = lfc_timepointTP2, y = -log10(p_timepointTP2), color = p_timepointTP2 < 0.05)) +
     geom_point() +
     scale_color_manual(values = c("red", "black")) +
     labs(title = "Volcano Plot sample cases T2 to T1 comparison",
          x = "Log Fold Change",
          y = "-log10(p-value)") +
  #   theme(legend.position = "none") +  
     geom_text(data = subset(F_case_ancombc_DA, p_timepointTP2 < 0.05),
                   aes(label = species), nudge_y = 0, nudge_x = -1, color = "blue",size=4)
      

ggsave(filename = "/PATH/volcano_plot_case_time_comparison_sample.pdf", plot = volcano_plot_case_time_comparison, width = 12, height = 8, dpi = 600) 

```


## TP1


```R
### let's prepare the input file and write it here. using: TP1_case_control_clr, metadata_Q1_cleaned, 
# output_ancombc2_TP1$res, adonis_TP1_sig

### first fixing up the taxa so we only keep the ones that ancombc2 deemed as important. This cuts it down by 1/3.
TP1_cc_clr_da_taxa <- TP1_case_control_clr[output_ancombc2_TP1$res$taxon]
TP1_cc_clr_da_taxa <- cbind(TP1_case_control_clr[,c(1,2)],TP1_cc_clr_da_taxa)

### fixing the metadata so it only includes the studienummer and the significant varibales
Q1_metadata_sig_TP1 <- metadata_Q1_cleaned[adonis_TP1_sig$variable]
Q1_metadata_sig_TP1 <- cbind(metadata_Q1_cleaned[,c(2,3)],Q1_metadata_sig_TP1)
Q1_metadata_sig_TP1 <- Q1_metadata_sig_TP1[,-c(2)] ## removing the Key variable so we don't have two of it

### making sure all the variables are characters
Q1_metadata_sig_TP1[,c(1:ncol(Q1_metadata_sig_TP1))] <- sapply(Q1_metadata_sig_TP1[,c(1:ncol(Q1_metadata_sig_TP1))],
                                                             function(x) as.character(x))
Q1_metadata_sig_TP1 <- data.frame(Q1_metadata_sig_TP1)                                                             

# merging the metadata with the ancombc variables with the permanova determined metadata variables
                                                             
TP1_cc_clr_da_taxa_metadata <- merge(TP1_cc_clr_da_taxa,Q1_metadata_sig_TP1,
                                    by="Studienummer",all.x=TRUE,all.y=FALSE)
TP1_cc_clr_da_taxa_metadata[is.na(TP1_cc_clr_da_taxa_metadata)] <- 0

## removing the studienummer column since we don't need that for the ML variables                                                             
                                                             
TP1_cc_clr_da_taxa <- TP1_cc_clr_da_taxa[,-c(1)]                                                             
TP1_cc_clr_da_taxa_metadata <- TP1_cc_clr_da_taxa_metadata[,-c(1)]
                                                             
#write.csv(TP1_cc_clr_da_taxa_metadata,"/ceph/projects/010_SweMaMi/analyses/nicole/sample_outputs/prediction/Oct_2023/input_files/TP1_cc_clr_da_taxa_metadata.csv",row.names=FALSE)
#write.csv(TP1_cc_clr_da_taxa,"/ceph/projects/010_SweMaMi/analyses/nicole/sample_outputs/prediction/Oct_2023/input_files/TP1_cc_clr_da_taxa.csv",row.names=FALSE)
                                                             
                                                             
                                                             
```

## For TP2


```R
### let's prepare the input file and write it here. using: TP2_case_control_clr, metadata_Q2_cleaned, 
# output_ancombc2_TP2$res, adonis_TP2_sig

### first fixing up the taxa so we only keep the ones that ancombc2 deemed as important. This cuts it down by 1/3.
TP2_cc_clr_da_taxa <- TP2_case_control_clr[output_ancombc2_TP2$res$taxon]
TP2_cc_clr_da_taxa <- cbind(TP2_case_control_clr[,c(1,2)],TP2_cc_clr_da_taxa)

### fixing the metadata so it only includes the studienummer and the significant varibales
Q2_metadata_sig_TP2 <- metadata_Q2_cleaned[adonis_TP2_sig$variable]
Q2_metadata_sig_TP2 <- cbind(metadata_Q2_cleaned[,c(2,4)],Q2_metadata_sig_TP2) ## 2 is studienummer, 4 key, will remove key. putting two so it doesn't change the column name.
Q2_metadata_sig_TP2 <- Q2_metadata_sig_TP2[,-c(2)] ## removing the Key variable so we don't have two of it

### making sure all the variables are characters
Q2_metadata_sig_TP2[,c(1:ncol(Q2_metadata_sig_TP2))] <- sapply(Q2_metadata_sig_TP2[,c(1:ncol(Q2_metadata_sig_TP2))],
                                                             function(x) as.character(x))
Q2_metadata_sig_TP2 <- data.frame(Q2_metadata_sig_TP2)                                                             

# merging the metadata with the ancombc variables with the permanova determined metadata variables
                                                             
TP2_cc_clr_da_taxa_metadata <- merge(TP2_cc_clr_da_taxa,Q2_metadata_sig_TP2,
                                    by="Studienummer",all.x=TRUE,all.y=FALSE)
TP2_cc_clr_da_taxa_metadata[is.na(TP2_cc_clr_da_taxa_metadata)] <- 0

## removing the studienummer column since we don't need that for the ML variables                                                             
                                                             
TP2_cc_clr_da_taxa <- TP2_cc_clr_da_taxa[,-c(1)]                                                             
TP2_cc_clr_da_taxa_metadata <- TP2_cc_clr_da_taxa_metadata[,-c(1)]
                                                             
#write.csv(TP2_cc_clr_da_taxa_metadata,"/ceph/projects/010_SweMaMi/analyses/nicole/sample_outputs/prediction/Oct_2023/input_files/TP2_cc_clr_da_taxa_metadata.csv",row.names=FALSE)
#write.csv(TP2_cc_clr_da_taxa,"/ceph/projects/010_SweMaMi/analyses/nicole/sample_outputs/prediction/Oct_2023/input_files/TP2_cc_clr_da_taxa.csv",row.names=FALSE)
                                                             
                                                             

```

### Plots

#### Let's start with adonis bar plots for R2/p-val

##### TP1


```R
###### fixing the rds file of the list that we have so that we can get the necessary bits out as a df

## combinding TP1 adonis first line of each row
first_rows <- lapply(adonis_TP1_list_clr, function(tbl) tbl[1, ])
adonis_TP1_df <- do.call(rbind, first_rows)
adonis_TP1_df <- data.frame(adonis_TP1_df)                     
adonis_TP1_df <- na.omit(adonis_TP1_df)
adonis_TP1_df_ordered <- adonis_TP1_df[order(adonis_TP1_df$Pr..F., decreasing = FALSE),]

adonis_TP1_df_ordered$variable <- rownames(adonis_TP1_df_ordered)                     
colnames(adonis_TP1_df_ordered) <- c("Df","SumOfSqs","R2","F","p.val","variable")   
                     
## combinding TP2 adonis first line of each row
first_rows <- lapply(adonis_TP2_list_clr, function(tbl) tbl[1, ])
adonis_TP2_df <- do.call(rbind, first_rows)
adonis_TP2_df <- data.frame(adonis_TP2_df)                     
adonis_TP2_df <- na.omit(adonis_TP2_df)              
adonis_TP2_df_ordered <- adonis_TP2_df[order(adonis_TP2_df$Pr..F., decreasing = FALSE),]
                     
adonis_TP2_df_ordered$variable <- rownames(adonis_TP2_df_ordered)
                     
colnames(adonis_TP1_df_ordered) <- c("Df","SumOfSqs","R2","F","p.val","variable")   
colnames(adonis_TP2_df_ordered) <- c("Df","SumOfSqs","R2","F","p.val","variable")   
                     
```


```R
## now let's make the barplot for it
TP1_adonis_bar <- ggplot(adonis_TP1_df_ordered,aes(x=reorder(variable,-R2),y=R2,fill=p.val)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="blue",high="pink") +
  coord_flip() +
  xlab("F TP1 variables") +
    theme(axis.text.y = element_text(size = 5))

```


```R
#### now here I'm going to remove some unnecessary columns (not sig, and repeated and remake the plots)

adonis_TP1_sig_variables <- adonis_TP1_df_ordered %>% filter(p.val < 0.04)
TP1_adonis_fewer_bar <- ggplot(adonis_TP1_sig_variables,aes(x=reorder(variable,-R2),y=R2,fill=p.val)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="blue",high="pink") +
  coord_flip() +
  xlab("F TP1 variables") +
    theme(axis.text.y = element_text(size = 12))

```


### TP2


```R
## now let's make the barplot for it

TP2_adonis_bar <- ggplot(adonis_TP2_df_ordered,aes(x=reorder(variable,-R2),y=R2,fill=p.val)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="blue",high="pink") +
  coord_flip() +
  xlab("F TP2 variables") +
    theme(axis.text.y = element_text(size = 8))

```


    

    



```R
## now let's make the barplot for it for sig
adonis_TP2_sig_variables <- adonis_TP2_df_ordered %>% filter(p.val < 0.06)

TP2_adonis_bar <- ggplot(adonis_TP2_sig_variables,aes(x=reorder(variable,-R2),y=R2,fill=p.val)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low="blue",high="pink") +
  coord_flip() +
  xlab("F TP2 variables") +
    theme(axis.text.y = element_text(size = 15))

```


    



### PCOA and NMDS plots

#### PCOA for TP1


```R
TP1_adonis_clr_species <-TP1_adonis_clr[,c(3:1424)]
#adonis_TP1_significant_filtered <- filter(adonis_TP1, P.val < 0.051) ### this is only used to find the significant ones
TP1_adonis_significant <- TP1_adonis_clr[,adonis_TP1_sig$variable]
rownames(TP1_adonis_significant) <- TP1_adonis_clr[,1]
rownames(TP1_adonis_clr_species) <- TP1_adonis_clr[,1]
TP1_adonis_clr_species_mat <- as.matrix(TP1_adonis_clr_species)

### make sure all the metadata variables are categorical
TP1_adonis_significant[,c(1:ncol(TP1_adonis_significant))] <- sapply(TP1_adonis_significant[,c(1:ncol(TP1_adonis_significant))],
                                                                   function(x) as.character(x))

                                                                     TP1_diss_mat <- vegdist(TP1_adonis_clr_species_mat,method="euclidean")
TP1_pcoa_results <- cmdscale(TP1_diss_mat)
TP1_pcoa_results_df <- data.frame(TP1_pcoa_results)
colnames(TP1_pcoa_results_df) <- c("PCoA1","PCoA2")
TP1_pcoa_metadata <- merge(TP1_pcoa_results_df,TP1_adonis_significant,by=0)
                                                                     
for (i in 4:ncol(TP1_pcoa_metadata)) {
  name <- colnames(TP1_pcoa_metadata[i])
  pcoa_plot <- ggplot(TP1_pcoa_metadata,aes(PCoA1,PCoA2,color=TP1_pcoa_metadata[,i])) +
  geom_point() +
  labs(color = name)
  file_name <- paste("/PATH/TP1_",name,"_pcoa.png",sep="")
  ggsave(filename = file_name, plot = pcoa_plot, width = 12, height = 8, dpi = 600) 
}                                                                     
```

#### NMDS for TP1


```R
TP1_nmds_results <- metaMDS(TP1_diss_mat, distance="euclidean",k=2, try=500)
NMDS1 <- data.frame(TP1_nmds_results$points[,1])
NMDS2 <- data.frame(TP1_nmds_results$points[,2])
TP1_nmds_df <- cbind(NMDS1,NMDS2)
colnames(TP1_nmds_df) <- c("NMDS1","NMDS2")
TP1_NMDS_metadata <- merge(TP1_nmds_df,TP1_adonis_significant,by=0)

for (i in 4:ncol(TP1_NMDS_metadata)) {
  name <- colnames(TP1_NMDS_metadata[i])
  pcoa_plot <- ggplot(TP1_NMDS_metadata,aes(NMDS1,NMDS2,color=TP1_NMDS_metadata[,i])) +
  geom_point() +
  labs(color = name)
  file_name <- paste("/PATH/TP1_",name,"_nmds.png",sep="")
  ggsave(filename = file_name, plot = pcoa_plot, width = 12, height = 8, dpi = 600) 
}

```


## PCOA for TP2


```R
TP2_adonis_clr_species <-TP2_adonis_clr[,c(3:1424)]
#adonis_TP1_significant_filtered <- filter(adonis_TP1, P.val < 0.051) ### this is only used to find the significant ones
TP2_adonis_significant <- TP2_adonis_clr[,adonis_TP2_sig$variable]
rownames(TP2_adonis_significant) <- TP2_adonis_clr[,1]
rownames(TP2_adonis_clr_species) <- TP2_adonis_clr[,1]
TP2_adonis_clr_species_mat <- as.matrix(TP2_adonis_clr_species)

### make sure all the metadata variables are categorical
TP2_adonis_significant[,c(1:ncol(TP2_adonis_significant))] <- sapply(TP2_adonis_significant[,c(1:ncol(TP2_adonis_significant))],
                                                                   function(x) as.character(x))

                                                                     
TP2_diss_mat <- vegdist(TP2_adonis_clr_species_mat,method="euclidean")
TP2_pcoa_results <- cmdscale(TP2_diss_mat)
TP2_pcoa_results_df <- data.frame(TP2_pcoa_results)
colnames(TP2_pcoa_results_df) <- c("PCoA1","PCoA2")
TP2_pcoa_metadata <- merge(TP2_pcoa_results_df,TP2_adonis_significant,by=0)     
                                                                     
for (i in 4:ncol(TP2_pcoa_metadata)) {
  name <- colnames(TP2_pcoa_metadata[i])
  pcoa_plot <- ggplot(TP2_pcoa_metadata,aes(PCoA1,PCoA2,color=TP2_pcoa_metadata[,i])) +
  geom_point() +
  labs(color = name)
  file_name <- paste("/PATH/TP2_",name,"_pcoa.png",sep="")
  ggsave(filename = file_name, plot = pcoa_plot, width = 12, height = 8, dpi = 600) 
}                                                                     
```

#### NMDS for TP2


```R
TP2_nmds_results <- metaMDS(TP2_diss_mat, distance="euclidean",k=2, try=500)
NMDS1 <- data.frame(TP2_nmds_results$points[,1])
NMDS2 <- data.frame(TP2_nmds_results$points[,2])
TP2_nmds_df <- cbind(NMDS1,NMDS2)
colnames(TP2_nmds_df) <- c("NMDS1","NMDS2")
TP2_NMDS_metadata <- merge(TP2_nmds_df,TP2_adonis_significant,by=0)
```



```R
for (i in 4:ncol(TP2_NMDS_metadata)) {
  name <- colnames(TP2_NMDS_metadata[i])
  pcoa_plot <- ggplot(TP2_NMDS_metadata,aes(NMDS1,NMDS2,color=TP2_NMDS_metadata[,i])) +
  geom_point() +
  labs(color = name)
  file_name <- paste("/PATH/TP2_",name,"_nmds.png",sep="")
  ggsave(filename = file_name, plot = pcoa_plot, width = 12, height = 8, dpi = 600) 
}

```

### plots for ANCOMBC results

#### let's start by reading in the file that we made previously


```R
##### keeping W in here
TP1_ancombc_read_in_summary_DA <- TP1_ancombc_read_in %>% filter(q_keyControl < 0.05)
TP1_ancombc_read_in_summary_DA_2 <- TP1_ancombc_read_in_summary_DA %>% filter(lfc_keyControl > 1.5)
TP1_ancombc_read_in_summary_DA_3 <- TP1_ancombc_read_in_summary_DA %>% filter(lfc_keyControl < -1.5)
TP1_ancombc_DA <- rbind(TP1_ancombc_read_in_summary_DA_2,TP1_ancombc_read_in_summary_DA_3)
TP1_ancombc_DA <- TP1_ancombc_DA[!duplicated(TP1_ancombc_DA$taxon),]

### adding a taxon column for only the species
taxon_names <- TP1_ancombc_DA[,c(1)]
taxon_names <- data.frame(taxon_names)
tax_IDs <- lapply(taxon_names, function(x) str_split(x,".s__"))
tax_IDs <- data.frame(tax_IDs)
tax_IDs_t <- t(tax_IDs)
tax_IDs_t <- data.frame(tax_IDs_t)   
TP1_ancombc_DA$species <- tax_IDs_t$X2

```


```R
### making the volcano plot for TP1

volcano_plot_TP1 <- ggplot(TP1_ancombc_DA, aes(x = lfc_keyControl, y = -log10(q_keyControl), color = q_keyControl < 0.05)) +
     geom_point() +
     scale_color_manual(values = c("red", "black")) +
     labs(title = "Volcano Plot for sample TP1",
          x = "Log Fold Change",
          y = "-log10(p-value)") +
     theme(legend.position = "none") +  
     geom_text(data = subset(TP1_ancombc_DA, q_keyControl < 0.05),
                   aes(label = species), nudge_y = -0.02, nudge_x = 0.25, color = "blue",size=4)
      

ggsave(filename = "/PATH/volcano_plot_TP1.pdf", plot = volcano_plot_TP1, width = 12, height = 8, dpi = 600) 

```


    

    



```R
volcano_plot_TP1_w <- ggplot(TP1_ancombc_DA, aes(x = W_keyControl, y = -log10(q_keyControl), color = q_keyControl < 0.05)) +
     geom_point() +
     scale_color_manual(values = c("red", "black")) +
     labs(title = "effect size Volcano Plot for sample TP1",
          x = "Log Fold Change",
          y = "-log10(p-value)") +
     theme(legend.position = "none") +  
     geom_text(data = subset(TP1_ancombc_DA, q_keyControl < 0.05),
                   aes(label = species), nudge_y = -0.02, nudge_x = 0.25, color = "blue",size=2)
      

ggsave(filename = "/PATH/volcano_plot_TP1_effect_size.png", plot = volcano_plot_TP1_w, width = 12, height = 8, dpi = 600) 

```


    

    


## ANCOMBC for TP2


```R
### Taking the variables that change it most
TP2_ancombc_read_in_summary_DA <- TP2_ancombc_read_in %>% filter(q_keyControl < 0.05)
TP2_ancombc_read_in_summary_DA_2 <- TP2_ancombc_read_in_summary_DA %>% filter(lfc_keyControl > 1.5)
TP2_ancombc_read_in_summary_DA_3 <- TP2_ancombc_read_in_summary_DA %>% filter(lfc_keyControl < -1.5)
TP2_ancombc_DA <- rbind(TP2_ancombc_read_in_summary_DA_2,TP2_ancombc_read_in_summary_DA_3)
TP2_ancombc_DA <- TP2_ancombc_DA[!duplicated(TP2_ancombc_DA$taxon),]

### adding a taxon column for only the species

taxon_names <- TP2_ancombc_DA[,c(1)]
taxon_names <- data.frame(taxon_names)
tax_IDs <- lapply(taxon_names, function(x) str_split(x,".s__"))
tax_IDs <- data.frame(tax_IDs)
tax_IDs_t <- t(tax_IDs)
tax_IDs_t <- data.frame(tax_IDs_t)   
TP2_ancombc_DA$species <- tax_IDs_t$X2
```


```R
### making the volcano plot for TP2

volcano_plot_TP2 <- ggplot(TP2_ancombc_DA, aes(x = lfc_keyControl, y = -log10(q_keyControl), color = q_keyControl < 0.05)) +
     geom_point() +
     scale_color_manual(values = c("red", "black")) +
     labs(title = "Volcano Plot for sample TP2",
          x = "Log Fold Change",
          y = "-log10(q-value)") +
     theme(legend.position = "none") +  
     geom_text(data = subset(TP2_ancombc_DA, q_keyControl < 0.05),
                   aes(label = species), nudge_y = -0.05, nudge_x = -0.1, color = "blue",size=4)
      
ggsave(filename = "/PATH/volcano_plot_TP2.pdf", plot = volcano_plot_TP2, width = 12, height = 8, dpi = 600) 

```


    
    

