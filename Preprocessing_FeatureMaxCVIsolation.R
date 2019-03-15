

##Code used to thin out the evidence file and find max CVs
##Could be made more efficient possibly by using the peptides file
library(dplyr)
library(tidyr)

####Preprocessing just to pull out info we want####

##Read in necessary files
myEvidence <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/evidence.txt")
myPeptides <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/peptides.txt")


##Get rid of half of experiments so they can be used for validation, and remove modified and bad scoring peptides
#myEvidenceHalved <- myEvidence[-(grep("_2$", myEvidence$Raw.file)), ]
workingEvidence <- myEvidence %>% filter(Modifications == "Unmodified", MS.MS.IDs != "-1") %>% filter(is.na(Intensity) == FALSE)

##Remove missing values then find max intensity experiment for each Sequence charge combination
maxCVs <- workingEvidence %>% group_by(Sequence, Charge) %>% summarise(myMaxCV = max(Intensity))
names(maxCVs) <- c(names(maxCVs[1:2]), "Intensity")


##merge Max CVs with other data to get raw file and amino acid counts
AllFeaturesMaxCV <- merge(maxCVs, workingEvidence, by = c("Sequence", "Charge", "Intensity"))
AllFeaturesMaxCV <- select(AllFeaturesMaxCV, Sequence:Length, Raw.file)
AllFeaturesMaxCV$Raw.file <- gsub(".*_CV_", "", AllFeaturesMaxCV$Raw.file)
AllFeaturesMaxCV$Raw.file <- factor(AllFeaturesMaxCV$Raw.file)
AllFeaturesMaxCV <- merge(AllFeaturesMaxCV, myPeptides[, 1:36], by = "Sequence")
AllFeaturesMaxCV <- select(AllFeaturesMaxCV, Sequence:Raw.file, A.Count:O.Count)

##Change names and write file containing max intensity CV, Charge, Sequence, AA counts and Length to a csv
AllFeaturesMaxCV$Raw.file <- gsub("_2$", "", AllFeaturesMaxCV$Raw.file)
names(AllFeaturesMaxCV) <- gsub("Raw\\.file", "Max Intensity CV", names(AllFeaturesMaxCV))
write.csv(AllFeaturesMaxCV, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv", row.names = F)



#### Looking at the distributions####
##Shows distribution of the peptides by charge and Max Intensity CV

library(RColorBrewer)

AllFeaturesMaxCV <- read.csv("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv")
names(AllFeaturesMaxCV) <- gsub("Max\\.Intensity\\.CV", "Peak_CV", names(AllFeaturesMaxCV))

##Code to make barplot of charge and Peak CV
barplot(table(AllFeaturesMaxCV$Charge, AllFeaturesMaxCV$Peak_CV), col = brewer.pal(7, "Spectral"))
legend("topleft", pch =15, col = brewer.pal(7, "Spectral"), legend = levels(factor(AllFeaturesMaxCV$Charge)))


#Make table showing numeric Charge and Peak CV distribution
CrossComparison <- as.data.frame(table(AllFeaturesMaxCV$Charge, AllFeaturesMaxCV$Peak_CV))
View(CrossComparison)
names(CrossComparison) <- c("Charge", "Peak_CV", "Freq")


##Function that allows for extraction of peptide intensities for peptide identified in the maximum number of experiments

FAIMSPeakBuilder <- function(myCharge, myCV, myFeatureFrame = AllFeaturesMaxCV){
  charge1CV25 <- filter(myFeatureFrame, Charge == myCharge, Peak_CV == myCV)
  charge1CV25 <- merge(charge1CV25, workingEvidence, by = c("Sequence", "Charge"))
  charge1CV25$Sequence <- factor(charge1CV25$Sequence)
  mymax <- max(table(charge1CV25$Sequence))
  mySeq <- labels(table(charge1CV25$Sequence)[(grep(mymax, table(charge1CV25$Sequence)))[1]])
  
  myPep <- charge1CV25 %>% filter(Sequence == mySeq) %>% select(Sequence:Peak_CV, PEP, Score, Raw.file, Intensity.x, Intensity.y, Retention.time, Retention.length) %>% 
    mutate(LogIntensity = log(Intensity.y, 2))
  myPep$Raw.file <- gsub("^.*_CV_|_2$", "", myPep$Raw.file)
  print(mymax)
  print(length(grep(mymax, table(charge1CV25$Sequence))))
  
  return(myPep)

}

Charge2by9 <- rbind(FAIMSPeakBuilder(2, 30),FAIMSPeakBuilder(2, 40),FAIMSPeakBuilder(2, 50),
                    FAIMSPeakBuilder(2, 60),FAIMSPeakBuilder(2, 70) ,FAIMSPeakBuilder(2, 80),
                    FAIMSPeakBuilder(2, 90),FAIMSPeakBuilder(2, 100),FAIMSPeakBuilder(2, 110))

ggplot(Charge2by9, aes(Retention.time, LogIntensity, color = Raw.file)) + geom_point() + geom_line() +
  facet_wrap(~Peak_CV, ncol = 3, nrow = 3, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90))


#Can do this over the different charge states probably would be helpful to add a smoother
Charge3by4 <- rbind(FAIMSPeakBuilder(3, 40),FAIMSPeakBuilder(3, 60),FAIMSPeakBuilder(3, 80),FAIMSPeakBuilder(3, 100))

ggplot(Charge3by4, aes(Retention.time, LogIntensity, color = Raw.file)) + geom_point() + geom_line() +
  facet_wrap(~Peak_CV, ncol = 2, nrow = 2, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90))


##Produces multicolor points and smoother for each facet

ggplot(Charge2by9, aes(Raw.file, LogIntensity)) + geom_point(color = factor(Charge2by9$Raw.file)) + 
  geom_smooth(method = "loess") + facet_wrap(~Sequence, ncol = 3, nrow = 3, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90))



Charge1ManyCV <- rbind(FAIMSPeakBuilder(1,25), FAIMSPeakBuilder(1,30), FAIMSPeakBuilder(1, 35))
ggplot(Charge1ManyCV, aes(Raw.file, LogIntensity, color = Sequence)) + geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(~Sequence, ncol = 1, nrow = 3)


#Histogram but puts all on same absolute scale. Not proportional
ggplot(AllFeaturesMaxCV, aes(Peak_CV, colour = Charge)) + geom_histogram(binwidth = 5) + 
  facet_wrap(~Charge, ncol = 2, nrow = 4)

#Reads in, isolates and makes density plots showing the error between the model and 
FeatureswError <- read.csv("FeaturesMaxCV_withErr.csv")
AllChargeby9 <- FeatureswError[grep("30|40|50|60|70|80|90|110|100", FeatureswError$y2), ]
ggplot(AllChargeby9, aes(x=Signd_Model_Error, group = Charge, fill = factor(Charge))) + geom_density(alpha = 0.4) + 
  facet_wrap(~y2, ncol = 3, nrow = 3)




####To normalize peptide counts####


#Read file in 
AllFeaturesMaxCV <- read.csv("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv")


NormPepCounts <- function(mydataframe, mydatacols){
  
  datacols <- mydataframe[, mydatacols]
  for(i in 1:nrow(datacols)){
    pepLength <- mydataframe$Length.x[i]
    for(j in 1:ncol(datacols)){
      datacols[i,j] <- datacols[i,j] / pepLength
      
    }
  }
  names(datacols) <- gsub("Count$", "NormCounts", names(datacols))
  newdataframe <- cbind(mydataframe, datacols)
  return(newdataframe)
}

##That took way too long so here we are

normCountsAllFeatures <- AllFeaturesMaxCV

normCountsAllFeatures <- normCountsAllFeatures %>% mutate(A.Norm = A.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(R.Norm = R.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(N.Norm = N.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(D.Norm = D.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(C.Norm = C.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(Q.Norm = Q.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(E.Norm = E.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(G.Norm = G.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(H.Norm = H.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(I.Norm = I.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(L.Norm = L.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(K.Norm = K.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(M.Norm = M.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(F.Norm = F.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(S.Norm = S.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(T.Norm = T.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(W.Norm = W.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(Y.Norm = Y.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(V.Norm = V.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(U.Norm = U.Count / Length.x)
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(O.Norm = O.Count / Length.x)

write.csv(normCountsAllFeatures, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/NormCountsFeatures_MaxCV.csv", row.names = F)

