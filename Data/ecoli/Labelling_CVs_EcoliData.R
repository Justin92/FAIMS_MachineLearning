#Script for labelling Ecoli data to test final model

#Load necessary packages
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)


#Define binarizing function, applies 0.5 default
binarizer <- function(x, threshold = 0.5){
  
  if(x > threshold){
    x <- 1
  }
  else{
    x <- 0
  }
}


#Read in Ecoli data
Ecoli_Data <- read.delim("C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Data/ecoli/ecoli_maxCVvalues.txt", stringsAsFactors = F)

#Replace some column names to allow for use of some old code
colnames(Ecoli_Data) <- gsub("z_modseq", "SeqCharge", colnames(Ecoli_Data))

#Parse out features such as Charge, Modified Sequence, and True Sequence, Length
Ecoli_Data <- Ecoli_Data %>% mutate(Charge = str_extract(SeqCharge, "[0-9]")) %>% 
  mutate(ModSequence = str_extract(SeqCharge, "[aA-zZ]+")) %>% mutate(Length = nchar(ModSequence))

#Make substitute canonical amino acids with modified amino acids
Ecoli_Data$Sequence <- gsub("a|m", "M", Ecoli_Data$ModSequence)


IntensityDist <- c()
#For each CV setting find the average proportional intensity of second best 
for(cv in unique(Ecoli_Data$maxCV)){
  #Select the peptides with max at that setting
  #Average the proportional intensity of each cv setting
  proportionialIntensities <- Ecoli_Data %>% filter(maxCV == cv) %>%
    select(X20:X95) %>% gather(CV_Setting, NormIntensity) %>% group_by(CV_Setting) %>%
    summarise(Mean_Norm_Int = mean(NormIntensity)) %>% 
    mutate(CV_Setting = gsub("^X", "", CV_Setting), Max = cv) 
  
  IntensityDist <- rbind(IntensityDist, proportionialIntensities)
}

#Plot of the average proportional intensity for each CV setting grouped by their Max CV
ggplot(IntensityDist, aes(CV_Setting, Mean_Norm_Int)) + geom_bar(stat = 'identity') + 
  facet_wrap(~Max, nrow=4, ncol=4) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Start with all Thresholds at 0.5
MaxCV <- unique(Ecoli_Data$maxCV)
Threshold <- rep(0.50, 16)
Threshold50table <- as.data.frame(cbind(MaxCV, Threshold))

#Check the average proportion for second best setting for each max CV group
#If it is >0.5 then that proportion becomes cutoff rather than 0.5
for(cv in Threshold50table$MaxCV){
  subtable <- IntensityDist %>% filter(Max == cv) %>% filter(Mean_Norm_Int < 1)
  newthreshold <- max(subtable$Mean_Norm_Int)
  if(newthreshold > 0.50){
    Threshold50table$Threshold[grep(cv, Threshold50table$MaxCV)] <- newthreshold
  }
}


#Turn CV labels into binary  by row (and in turn by max cv)
for(cv in unique(Ecoli_Data$maxCV)){
  cv_rows <- grep(cv, Ecoli_Data$maxCV)
  Ecoli_Data[cv_rows, grep("X[0-9]+", colnames(Ecoli_Data))] <- apply(Ecoli_Data[cv_rows, grep("X[0-9]+", colnames(Ecoli_Data)), drop=F], c(1,2), 
                                                                        binarizer, threshold = Threshold50table$Threshold[grep(cv, Threshold50table$MaxCV)])
}

#Store overlapping peptides (Need Human data in workspace)
Overlap_peptides <- Ecoli_Data$SeqCharge[Ecoli_Data$SeqCharge %in% CleanReport$SeqCharge]
write.csv(Overlap_peptides, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Data/Overlapping_peptideIons.csv")

#Number of labels(42719 peptide ions overall) (68 overlap with human)
#X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
#304  1632  4120  8092 11152 12831 12225 11268  9198  7118  4969  3053  1861  1169   589   260 

#Remove non-continuous labels
Ecoli_Data <- Ecoli_Data %>% mutate(LabelSequence = paste(X20, X25, X30, X35, X40, X45, X50, X55, X60, X65, X70, X75, X80, X85, X90, X95, sep = ""))
Ecoli_Data2 <- Ecoli_Data[-grep("1[0]+1", Ecoli_Data$LabelSequence), ]


#Number of labels(41063 peptide ions overall) (63 overlap with human)
#X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
#258  1541  4031  7900 10881 12450 11762 10744  8727  6634  4586  2762  1656  1025   503   211 

write.csv(Ecoli_Data2, "C:/Users/jmcketney.AD/Desktop/FAIMS_MachineLearning/Data/ecoli/ProcessedLabelled_Ecoli.csv", row.names = F)

