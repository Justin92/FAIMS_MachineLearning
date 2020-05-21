
###Import necessary libraries####
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)
library(ggplot2)

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

#Read scaled intensity from CSV
CleanReport <- read.delim("P:/JGM_FAIMS_CVprediction/JMM_new_data/NEW_JMMdata_maxCVvalues.txt")

#Replace some column names to allow for use of some old code
colnames(CleanReport) <- gsub("z_modseq", "SeqCharge", colnames(CleanReport))

#removing acetylation notation from the modified sequences and replace with 'a' in modified sequence
CleanReport$SeqCharge <- gsub("\\[\\+42\\]", "a", CleanReport$SeqCharge)

#Parse out features such as Charge, Modified Sequence, and Length
CleanReport <- CleanReport %>% mutate(Charge = str_extract(SeqCharge, "[0-9]")) %>% 
  mutate(ModSequence = str_extract(SeqCharge, "[aA-zZ]+")) %>% mutate(Length = nchar(ModSequence))

#Make substitute canonical amino acids with modified amino acids
CleanReport$Sequence <- gsub("a|m", "M", CleanReport$ModSequence)

CleanReportPerm <- CleanReport

###Exploratory Analysis to decide what we want the algorithm to predict

IntensityDist <- c()
#For each CV setting
for(cv in unique(CleanReportPerm$maxcv_naomit)){
  #Select the peptides with max at that setting
  #Average the proportional intensity of each cv setting
  proportionialIntensities <- CleanReportPerm %>% filter(maxcv_naomit == cv) %>%
    select(X20:X95) %>% gather(CV_Setting, NormIntensity) %>% group_by(CV_Setting) %>%
    summarise(Mean_Norm_Int = mean(NormIntensity)) %>% 
    mutate(CV_Setting = gsub("^X", "", CV_Setting), Max = cv) 
  
  IntensityDist <- rbind(IntensityDist, proportionialIntensities)
}

#Plot of the average proportional intensity for each CV setting grouped by their Max CV
ggplot(IntensityDist, aes(CV_Setting, Mean_Norm_Int)) + geom_bar(stat = 'identity') + 
  facet_wrap(~Max, nrow=4, ncol=4) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#New Threshold Table where 50% of max is default unless the second best settings average is >50% in which case that proportion is the threshold

#Start with all Thresholds at 0.5
MaxCV <- unique(CleanReportPerm$maxcv_naomit)
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


#So now need to do this by row instead of column

CleanReport <- CleanReportPerm

for(cv in unique(CleanReport$maxcv_naomit)){
  
  cv_rows <- grep(cv, CleanReport$maxcv_naomit)
  CleanReport[cv_rows, grep("X[0-9]+", colnames(CleanReport))] <- apply(CleanReport[cv_rows, grep("X[0-9]+", colnames(CleanReport)), drop=F], c(1,2), 
                                                                        binarizer, threshold = Threshold50table$Threshold[grep(cv, Threshold50table$MaxCV)])
}
#colSums(CleanReport[2:17])

#Number of labels(128402 peptides overall)
#X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
#778  5373 15487 25410 34659 35728 34305 30063 25545 14961 17431 11134  6952  3941  2196  1096 

#Remove non-continuous labels
CleanReport <- CleanReport %>% mutate(LabelSequence = paste(X20, X25, X30, X35, X40, X45, X50, X55, X60, X65, X70, X75, X80, X85, X90, X95, sep = ""))
CleanReport <- CleanReport[-grep("1[0]+1", CleanReport$LabelSequence), ]


#Number of labels(122847 peptides overall)
#X20   X25   X30   X35   X40   X45   X50   X55   X60   X65   X70   X75   X80   X85   X90   X95 
#727  5277 15219 24862 33730 34544 32910 28299 22782 14066 14852 10126  6370  3587  1995   972 

write.csv(CleanReport, "D:/Projects/FAIMS_MachineLearning/2020/March/50percentMaxPlusThreshold.csv", row.names = F)

