
###Import necessary libraries####
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)


#########
##################         NEW      STRATEGY    FOR     PREPROCESSING      JGM(Fragpipe & Skyline)  DATA      ###################


#Make binary all values for CV settings

#Define binarizing function, applies 0.5 default
binarizer <- function(x, threshold = 0.5){
  
  if(x >= threshold){
    x <- 1
  }
  else{
    x <- 0
  }
}

####Reading in old dataset#####
#CleanReport <- read.delim("D:/Projects/FAIMS_MachineLearning/2020/January/InitialDataProcessingFromJesse/JMMdata_maxCVvalues.txt", stringsAsFactors = F)
#####GENERATE CLEAN DATASET WITH BINARY LABELS AND FEATURES####

#Read scaled intensity from CSV
CleanReport <- read.delim("P:/JGM_FAIMS_CVprediction/JMM_new_data/NEW_JMMdata_maxCVvalues.txt")

#Replace some column names to allow for use of some old code
colnames(CleanReport) <- gsub("z_modseq", "SeqCharge", colnames(CleanReport))

#removing acetylation notation from the modified sequences and replace with 'a' in modified sequence
CleanReport$SeqCharge <- gsub("\\[\\+42\\]", "a", CleanReport$SeqCharge)

#Parse out features such as Charge, Modified Sequence, and True Sequence, Length
CleanReport <- CleanReport %>% mutate(Charge = str_extract(SeqCharge, "[0-9]")) %>% 
  mutate(ModSequence = str_extract(SeqCharge, "[aA-zZ]+")) %>% mutate(Length = nchar(ModSequence))

#Make substitute canonical amino acids with modified amino acids
CleanReport$Sequence <- gsub("a|m", "M", CleanReport$ModSequence)


##Need to generate threshold table for cutoffs##
Rawjgmdata <- read.csv("D:/Projects/FAIMS_MachineLearning/2020/January/InitialDataProcessingFromJesse/Peptide_area_q_idotp.csv", 
                       stringsAsFactors = F)
CVindex <- grep("Total\\.Area", colnames(Rawjgmdata))
cvsetting <- c()
medianArea <- c()

for(CVindex in grep("Total\\.Area", colnames(Rawjgmdata))){
  
  cvsetting <- c(cvsetting, gsub("\\.Total\\.Area$", "", colnames(Rawjgmdata[CVindex])))
  medianArea <- c(medianArea, median(as.numeric(Rawjgmdata[Rawjgmdata[,CVindex] > 0, CVindex]), na.rm = T))
  
  
}

medianArea_df <- as.data.frame(cbind(cvsetting, medianArea)) %>% 
  mutate(medianArea = as.numeric(as.character(medianArea)), cvsetting = as.character(cvsetting)) %>%
  group_by(cvsetting) %>% mutate(Threshold = (mean(medianArea) - min(medianArea))/ (max(medianArea) - min(medianArea))) #%>% 


###Three options for the cutoff###

## 1)   Original, just cuts off at 50%####
##CleanReport[1:nrow(CleanReport),2:max(grep("X[0-9]+", colnames(CleanReport)))] <- 
##  as.data.frame(apply(CleanReport[1:nrow(CleanReport),2:max(grep("X[0-9]+", colnames(CleanReport)))], c(1,2), binarizer))

##Numbers of each label
## X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
## 778  5374 15498 25470 35397 37203 35591 31146 26324 15904 17986 11517  7343  4350  2475  1231

###################
## 2)Cutoff based on median intensity for that CV setting proportional to the max and min median intensity (BEST)####
for(i in grep("X[0-9]+", colnames(CleanReport))){
  CleanReport[1:nrow(CleanReport), i] <- as.data.frame(apply(CleanReport[1:nrow(CleanReport),i, drop = F], c(1,2), binarizer, threshold = medianArea_df$Threshold[grep(colnames(CleanReport[i]), medianArea_df$cvsetting)]))
}
##Numbers of each label
## X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
## 5489 13296 24087 29467 26872 19748 25219 28509 28265 24214 21750 16993 12243  8360  5265  3134 
###################
## 3) Cutoff based on the order of the median intensity similar to Spearman####
##Third try with dynamic threshold that is proportional to rank

##for(i in grep("X[0-9]+", colnames(CleanReport))){
##  CleanReport[1:nrow(CleanReport), i] <- as.data.frame(apply(CleanReport[1:nrow(CleanReport),i, drop = F], c(1,2), binarizer, threshold = medianArea_df$BlockThreshold[grep(colnames(CleanReport[i]), medianArea_df$cvsetting)]))
##}
##Numbers of each label
##X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
##2077  7102 15911 21852 26033 20398 23241 24369 21996 12523 17698 12586  8566  5787  3655  2061


###################
###Saving Clean Data#####

#1
#write.csv(CleanReport, "D:/Projects/FAIMS_MachineLearning/2020/January/InitialDataProcessingFromJesse/NewJGMDataWFeatures.csv", row.names = F)
#2
write.csv(CleanReport, "D:/Projects/FAIMS_MachineLearning/2020/March/DynamicThresholdSecondJGMDataset.csv", row.names = F)
#3
#write.csv(CleanReport, "D:/Projects/FAIMS_MachineLearning/2020/March/SpearmanThresholdSecondJGMDataset.csv", row.names = F)

###################







##################          OLD(03/20/19 Pre Multilabel)      STRATEGY    FOR     PREPROCESSING      MAXQUANT  DATA      ###################






##Code used to thin out the evidence file and find max CVs
##Could be made more efficient possibly by using the peptides file
####Preprocessing just to pull out info we want####

##Read in necessary files
myEvidence <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/evidence.txt")
myPeptides <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/peptides.txt")


##Get rid of half of experiments so they can be used for validation
#myEvidenceHalved <- myEvidence[-(grep("_2$", myEvidence$Raw.file)), ]

#Remove modified, bad scoring peptides and those with missing values
workingEvidence <- myEvidence %>% filter(Modifications == "Unmodified", MS.MS.IDs != "-1") %>% filter(is.na(Intensity) == FALSE)

##find max intensity experiment for each Sequence charge combination
maxCVs <- workingEvidence %>% group_by(Sequence, Charge) %>% summarise(myMaxIntensity = max(Intensity))
names(maxCVs) <- c(names(maxCVs[1:2]), "Intensity")


##merge Max CVs with other data to get raw file and amino acid counts
AllFeaturesMaxCV <- merge(maxCVs, workingEvidence, by = c("Sequence", "Charge", "Intensity"))



###############################################
####Filter peptides with limited points####

#Add Total Observation row to Feature list
AllFeaturesMaxCV$Total_Observations = 0

#Iterate through and count total observations

#Add CV_Setting column to both evidence and features dataframe
AllFeaturesMaxCV <- AllFeaturesMaxCV %>% mutate(CV_Setting = gsub("CV_|_2", "", Experiment))
workingEvidence <- workingEvidence %>% mutate(CV_Setting = gsub("CV_|_2", "", Experiment))

#Make data frame with CV experiment observation counts for each peptide/charge combo and give it the same number
#of rows as the feature table

###############################################
####    MAKING OBSERVATION COUNT LAYER  ####
numRows <- nrow(AllFeaturesMaxCV)
Ion_Observations <- data.frame(Sequence=c(rep(NA, numRows)), Charge=c(rep(NA, numRows)), CV_15 = c(rep(NA, numRows)),
                               CV_20=c(rep(NA, numRows)), CV_25=c(rep(NA, numRows)), CV_30=c(rep(NA, numRows)),
                               CV_35=c(rep(NA, numRows)), CV_40=c(rep(NA, numRows)), CV_45=c(rep(NA, numRows)),
                               CV_50=c(rep(NA, numRows)), CV_55=c(rep(NA, numRows)), CV_60=c(rep(NA, numRows)),
                               CV_65=c(rep(NA, numRows)), CV_70=c(rep(NA, numRows)), CV_75=c(rep(NA, numRows)),
                               CV_80=c(rep(NA, numRows)), CV_85=c(rep(NA, numRows)), CV_90=c(rep(NA, numRows)),
                               CV_95=c(rep(NA, numRows)), CV_100=c(rep(NA, numRows)), CV_105=c(rep(NA, numRows)), 
                               CV_110=c(rep(NA, numRows)), CV_115=c(rep(NA, numRows)), CV_120=c(rep(NA, numRows)))

for(i in 1:nrow(AllFeaturesMaxCV)){
  #Extracts peptide charge combo
  MySequence <- AllFeaturesMaxCV$Sequence[i]
  myCharge <- AllFeaturesMaxCV$Charge[i]
  
  #Uses to subset working evidence data frame
  workingFrame <-filter(workingEvidence, Sequence == AllFeaturesMaxCV$Sequence[i], 
                                           Charge == AllFeaturesMaxCV$Charge[i])
  #makes working ion observation row
  Ion_Observations$Sequence[i] <- as.character(MySequence)
  Ion_Observations$Charge[i] <- myCharge
  for(j in 3:ncol(Ion_Observations)){
    Ion_Observations[i, j] <- length(grep(names(Ion_Observations[j]), workingFrame$Experiment))
  }
}


for(i in 1:nrow(AllFeaturesMaxCV)){
  
  AllFeaturesMaxCV$Total_Observations[i] <- nrow(filter(workingEvidence, Sequence == AllFeaturesMaxCV$Sequence[i],
                                                        Charge == AllFeaturesMaxCV$Charge[i]))
}

###############################################
####    MAKING MEAN INTENSITY LAYER     ####
#Make data frame with CV experiment mean intensities for each peptide/charge combo and give it the same number
#of rows as the feature table
Mean_Intensity <- data.frame(Sequence=c(rep(NA, numRows)), Charge=c(rep(NA, numRows)), CV_15 = c(rep(NA, numRows)),
                               CV_20=c(rep(NA, numRows)), CV_25=c(rep(NA, numRows)), CV_30=c(rep(NA, numRows)),
                               CV_35=c(rep(NA, numRows)), CV_40=c(rep(NA, numRows)), CV_45=c(rep(NA, numRows)),
                               CV_50=c(rep(NA, numRows)), CV_55=c(rep(NA, numRows)), CV_60=c(rep(NA, numRows)),
                               CV_65=c(rep(NA, numRows)), CV_70=c(rep(NA, numRows)), CV_75=c(rep(NA, numRows)),
                               CV_80=c(rep(NA, numRows)), CV_85=c(rep(NA, numRows)), CV_90=c(rep(NA, numRows)),
                               CV_95=c(rep(NA, numRows)), CV_100=c(rep(NA, numRows)), CV_105=c(rep(NA, numRows)), 
                               CV_110=c(rep(NA, numRows)), CV_115=c(rep(NA, numRows)), CV_120=c(rep(NA, numRows)))
#Iterating through the features table and calculating the mean at each 
for(i in 1:nrow(AllFeaturesMaxCV)){
  #Extracts peptide charge combo
  MySequence <- AllFeaturesMaxCV$Sequence[i]
  myCharge <- AllFeaturesMaxCV$Charge[i]
  
  #Uses to subset working evidence data frame
  workingFrame <-workingEvidence%>% filter(Sequence == AllFeaturesMaxCV$Sequence[i], 
                        Charge == AllFeaturesMaxCV$Charge[i]) %>% select(Intensity, Experiment)
  #makes working ion observation row
  Mean_Intensity$Sequence[i] <- as.character(MySequence)
  Mean_Intensity$Charge[i] <- myCharge
  for(j in 3:ncol(Mean_Intensity)){
    Mean_Intensity[i, j] <- mean(workingFrame$Intensity[grep(names(Ion_Observations[j]), workingFrame$Experiment)])
  }
}

###############################################
###     CLEANING UP THE FEATURES TABLE     ####

#Select the Sequence, Intensity, Charge, Length and Raw file(as a stand in for max CV)
AllFeaturesMaxCV <- select(AllFeaturesMaxCV, Sequence:Length, Raw.file)

#Get rid of all infor besides CV for the raw file
AllFeaturesMaxCV$Raw.file <- gsub(".*_CV_", "", AllFeaturesMaxCV$Raw.file)

#Make raw file (now basically CV) a factor
AllFeaturesMaxCV$Raw.file <- factor(AllFeaturesMaxCV$Raw.file)

#Merge our peptides with the amino acid count information from the peptide.txt file
AllFeaturesMaxCV <- merge(AllFeaturesMaxCV, myPeptides[, 1:36], by = "Sequence")

#Select columns pulled from evidence(charge, sequence, etc) along with amino acid counts from peptide file
AllFeaturesMaxCV <- select(AllFeaturesMaxCV, Sequence:Raw.file, A.Count:O.Count)

##Change names and write file containing max intensity CV, Charge, Sequence, AA counts and Length to a csv
AllFeaturesMaxCV$Raw.file <- gsub("_2$", "", AllFeaturesMaxCV$Raw.file)
names(AllFeaturesMaxCV) <- gsub("Raw\\.file", "Max Intensity CV", names(AllFeaturesMaxCV))
write.csv(AllFeaturesMaxCV, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv", row.names = F)

###############################################
####   TO NORMALIZE AMINO ACID COUNTS    ####
#One long pipeline operator for making new normalized counts#
normCountsAllFeatures <- normCountsAllFeatures %>% mutate(A.Norm = A.Count / Length.x) %>%
  mutate(R.Norm = R.Count / Length.x) %>%
  mutate(N.Norm = N.Count / Length.x) %>%
  mutate(D.Norm = D.Count / Length.x) %>%
  mutate(C.Norm = C.Count / Length.x) %>%
  mutate(Q.Norm = Q.Count / Length.x) %>%
  mutate(E.Norm = E.Count / Length.x) %>%
  mutate(G.Norm = G.Count / Length.x) %>%
  mutate(H.Norm = H.Count / Length.x) %>%
  mutate(I.Norm = I.Count / Length.x) %>%
  mutate(L.Norm = L.Count / Length.x) %>%
  mutate(K.Norm = K.Count / Length.x) %>%
  mutate(M.Norm = M.Count / Length.x) %>%
  mutate(F.Norm = F.Count / Length.x) %>%
  mutate(S.Norm = S.Count / Length.x) %>%
  mutate(T.Norm = T.Count / Length.x) %>%
  mutate(W.Norm = W.Count / Length.x) %>%
  mutate(Y.Norm = Y.Count / Length.x) %>%
  mutate(V.Norm = V.Count / Length.x) %>%
  mutate(U.Norm = U.Count / Length.x) %>%
  mutate(O.Norm = O.Count / Length.x)

write.csv(normCountsAllFeatures, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/NormCountsFeatures_MaxCV.csv", row.names = F)




###############################################
####    LOOKING AT DISTRIBUTIONS ACROSS RT CHARGE AND PEAK CV    ####
##Shows distribution of the peptides by charge and Max Intensity CV

library(RColorBrewer)

AllFeaturesMaxCV <- read.csv("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Features_MaxCVs.csv")
names(AllFeaturesMaxCV) <- gsub("Max\\.Intensity\\.CV", "Peak_CV", names(AllFeaturesMaxCV))

##Code to make barplot of charge and Peak CV
barplot(table(AllFeaturesMaxCV$Charge, AllFeaturesMaxCV$Peak_CV), col = brewer.pal(7, "Spectral"))
legend("topleft", pch =15, col = brewer.pal(7, "Spectral"), legend = levels(factor(AllFeaturesMaxCV$Charge)))


#Make table showing numeric Charge and Peak CV distribution
CrossComparison <- as.data.frame(table(AllFeaturesMaxCV$Charge, AllFeaturesMaxCV$Peak_CV))
names(CrossComparison) <- c("Charge", "Peak_CV", "Freq")



#Function that takes in charge, CV and dataframe(currently only works with AllFeaturesMaxCV)
#It finds peptide with the most scan observations across all experiments and
#pulls out the score metrics and retention time for those scans
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


#Makes a data frame with most point-rich peptides that are charge 2+ and CV from 30-110 by 10
Charge2by9 <- rbind(FAIMSPeakBuilder(2, 30),FAIMSPeakBuilder(2, 40),FAIMSPeakBuilder(2, 50),
                    FAIMSPeakBuilder(2, 60),FAIMSPeakBuilder(2, 70) ,FAIMSPeakBuilder(2, 80),
                    FAIMSPeakBuilder(2, 90),FAIMSPeakBuilder(2, 100),FAIMSPeakBuilder(2, 110))
#Plot of those peptide/charge combos' intensities across retention time colored by their CV Setting
ggplot(Charge2by9, aes(Retention.time, LogIntensity, color = Raw.file)) + geom_point() + geom_line() +
  facet_wrap(~Peak_CV, ncol = 3, nrow = 3, scales = "free") + 
  theme(axis.text.x = element_text(angle = 90))


#Selects different set of peptide/charge combo
Charge3by4 <- rbind(FAIMSPeakBuilder(3, 40),FAIMSPeakBuilder(3, 60),FAIMSPeakBuilder(3, 80),FAIMSPeakBuilder(3, 100))
#Graps same as above but split into facets by assigned peak CV setting
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

###############################################
###   ANALYSIS OF ERROR IN INITIAL MODELS    ####

#Reads in, isolates and makes density plots showing the error between the model and 
FeatureswError <- read.csv("FeaturesMaxCV_withErr.csv")
AllChargeby9<- FeatureswError[grep("30|40|50|60|70|80|90|110|100", FeatureswError$y2), ]

FeatureswError <- FeatureswError %>% mutate(Signed_Model_Error = y2 - y2_model)

#Generates density plot showing the error distribution colored by charge and split into facets by peak CV
ggplot(AllChargeby9, aes(x=Signd_Model_Error, group = Charge, fill = factor(Charge))) + geom_density(alpha = 0.4) + 
  facet_wrap(~y2, ncol = 3, nrow = 3) + theme_gray()


#Lets make function that takes in the table with error and attaches
#a vector to it describing number of points across chrom peak

#Iterates through Error dataframe and attaches column for number of points at assigned Peak CV setting
for(i in 1:nrow(FeatureswError)){
  
  FeatureswError$No_Points[i] = nrow(filter(workingEvidence, Sequence == FeatureswError$Sequence[i], 
                                            Charge == FeatureswError$Charge[i], 
                                            myCV == FeatureswError$Experiment[i]))
}

#Iterates through Error dataframe and attaches column for number of points across all experiments
for(i in 1:nrow(FeatureswError)){
  
  FeatureswError$TotalObs[i] = nrow(filter(workingEvidence, Sequence == FeatureswError$Sequence[i], 
                                            Charge == FeatureswError$Charge[i]))
}

#Split output with error into quartiles based both on signed and absolute error
FeatureswError <- FeatureswError[order(FeatureswError$model_error), ]
FeatureswError$Abs_Error_Quartile <- c(rep("Q1", 8243), rep("Q2",8243), rep("Q3",8243), rep("Q4",8242))

FeatureswError <- FeatureswError[order(FeatureswError$Signed_Model_Error), ]
FeatureswError$Signed_Error_Quartile <- c(rep("Q1", 8243), rep("Q2",8243), rep("Q3",8243), rep("Q4",8242))


write.csv(FeatureswError, "AllFeatureswErrorObsNumbers.csv", row.names = F)
