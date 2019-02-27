

##Code used to thin out the evidence file and find max CVs
##Could be made more efficient possibly by using the peptides file
library(dplyr)
library(tidyr)

##Read in necessary files
myEvidence <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/evidence.txt")
myPeptides <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/peptides.txt")


##Get rid of half of experiments so they can be used for validation, and remove modified and bad scoring peptides
myEvidenceHalved <- myEvidence[-(grep("_2$", myEvidence$Raw.file)), ]
workingEvidence <- myEvidenceHalved %>% filter(Modifications == "Unmodified", MS.MS.IDs != "-1")

##Remove missing values then find max intensity experiment for each Sequence charge combination
maxCVs <- workingEvidence %>% filter(is.na(Intensity) == FALSE) %>% group_by(Sequence, Charge) %>% summarise(myMaxCV = max(Intensity))
names(maxCVs) <- c(names(maxCVs[1:2], "Intensity"))


##merge Max CVs with other data to get raw file and amino acid counts
AllFeaturesMaxCV <- merge(maxCVs, workingEvidence, by = c("Sequence", "Charge", "Intensity"))
AllFeaturesMaxCV <- select(AllFeaturesMaxCV, Sequence:Length, Raw.file)
AllFeaturesMaxCV$Raw.file <- factor(AllFeaturesMaxCV$Raw.file)
AllFeaturesMaxCV <- merge(AllFeaturesMaxCV, myPeptides[, 1:36], by = "Sequence")
AllFeaturesMaxCV <- select(allfeatures, Sequence:Raw.file, A.count:O.count)

##Change names and write file containing max intensity CV, Charge, Sequence, AA counts and Length to a csv
names(AllFeaturesMaxCV) <- gsub("Raw\\.file", "Max Intensity CV", names(AllFeaturesMaxCV))
write.csv(AllFeaturesMaxCV, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv")
