

#### explore the CV distributions from the re-analysis of stepped CVs
## report exported from skyline to contain the peak area per CV and some quality of identification metrics e.g. q-value
library(stringr)

s<-read.csv("P:/JGM_FAIMS_CVprediction/JMM_new_data/ecoli_Peptide_area_q_idotp.csv", stringsAsFactors = F)

head(s)


## remove peptides matched to decoy proteins
nrow(s)
s<-s[-grep(s[,"Protein"], pattern="rev"),]
nrow(s)
s<-s[-grep(s[,"Protein"], pattern="Decoys"),]
nrow(s)

head(s)
### from modified sequence, replace 'M[+16]' with 'm' 
peptides_fixed<-str_replace_all(s$Peptide.Modified.Sequence, pattern="M\\[\\+16\\]", "m")
### and replace the [+57] with nothing
peptides_fixed<- str_replace_all(peptides_fixed, pattern="\\[\\+57\\]", "")
peptides_fixed
## concatenate charge to peptide
s<-cbind(paste0(s[,"Precursor.Charge"], peptides_fixed, sep=""),s)

head(s)

colnames(s)[1]<-"z_modpep"

###### make matching subset dataframes for idotp, area, identified #######

### find the columns that have "Isotope.Dot.Product"
idotpcols<-grep(colnames(s), pattern="Isotope.Dot.Product")
### find the columns that have "Identified"
qcols<-grep(colnames(s), pattern="Detection.Q.Value")
### find the columns that have "Isotope.Dot.Product"
areacols<-grep(colnames(s), pattern="Area")

s_idotp<-s[,idotpcols]
s_q<-s[,qcols]
s_area<-s[,areacols]

head(s_idotp)
head(s_q)
head(s_area)

### set areas to zero where idotp < 0.9
s_area[s_idotp<0.85]<-0
s_area[s_q>0.01]<-0
s_area[s_q=="#N/A"]<-0
head(s_area)
head(s_q)
# set rownames to z_modpep
rownames(s_area)<- paste0(s[,"Precursor.Charge"], peptides_fixed, sep="")
head(s_area)
max(as.numeric(s_area[1,]))

####### big loop to change data for plotting ########
### 1. combine columns to get average
### 2. get index of which is the
### 3. normalize to max value

### find the column indexes that match to each CV value
cv_values<-seq(from=20, to=95, by=5)
cvcols<-list()
for(x in cv_values){
  cvcols[[as.character(x)]]<-grep(colnames(s_area), pattern=paste("CV",as.character(x), sep=""))
}

cvcols

itern<-nrow(s_area)
itern
i<-1

maxcolindex<-rep(0, times=itern)  ## empty vector to fill with index of largest value

nrow(s_area)
head(s_area)
s_area[1,]
rowSums(is.na(s_area[1,]))

#### average each of the duplicate measures, and then normalize each column to the maximum value 
for( i in 1:itern){
  tmp_vec <- as.numeric(s_area[i,])
  ## combine columns to get average
  j=1
  new_ave <- c(rep(0, times=length(cvcols)))
  for(x in cvcols){
    new_ave[j]<-ave(tmp_vec[c(x)])[1]
    j=j+1
  }
  ## get index of the largest value
  maxcolindex[i]<-which.max(new_ave)
  
  ### normalize to the max value
  newline<-new_ave/new_ave[maxcolindex[i]]
  
  if(i==1){
    s_area_norm<- rbind(newline)
    colnames(s_area_norm)<- names(cvcols)
  }
  if(i>1){
    s_area_norm<-rbind(s_area_norm, newline)
  }
  print(i)  
}

rownames(s_area_norm) <- rownames(s_area)
## get groups of rows where each value is the max
cvmaxrows<-list()
#x<- names(cvcols)[1]
#which(names(cvcols)==x)
maxCVvec<-c(rep(0, times=itern))
#as.numeric(x)

for( x in names(cvcols)){
  cvmaxrows[[x]] <- which(maxcolindex== which(names(cvcols)==x) )
  maxCVvec[ cvmaxrows[[x]] ] <- as.numeric(x)
}

s_area_norm

### 
s_area_norm_naomit<- s_area_norm[rowSums(is.na(s_area_norm)) != ncol(s_area_norm),]
maxcv_naomit<-maxCVvec[rowSums(is.na(s_area_norm)) != ncol(s_area_norm)]
output4machinelearning<-cbind(s_area_norm_naomit, maxcv_naomit)



colnames(output4machinelearning)[ncol(output4machinelearning)]<-"maxCV"
head(output4machinelearning, 10)

write.table(output4machinelearning, file="P:/JGM_FAIMS_CVprediction/JMM_new_data/ecoli_maxCVvalues.txt",
            sep="\t",quote=F, row.names = T, col.names = T)
nrow(output4machinelearning)
hist(maxcv_naomit)

head(output4machinelearning)
cvmaxrows[[1]]

colnames(s_area_norm)
head(s_area_norm)

### count the number that max at each value for histogram making
num_per_cv<-c()
for(x in names(cvcols)){
  num_per_cv<-c(num_per_cv,length(cvmaxrows[[x]]))
}
?barplot
barplot(num_per_cv, names.arg =names(cvcols))

### which are max @ 15?
s_area_norm_max15<-s_area_norm[cvmaxrows[[1]],]
head(s_area_norm_max15)


### which are max @ 60?
s_area_norm_max60<-s_area_norm[cvmaxrows[[10]],]

### loop through each column and get the mean and standard deviation
meanlist<-list()
sdlist<-list()
minlist<-list()
maxlist<-list()
topquantile<-list()
#cvcols
#x<-"70"
for(x in names(cvcols)){
  meanlist[[x]]<-mean(s_area_norm_max60[,x])
  sdlist[[x]]<-sd(s_area_norm_max60[,x])
  minlist[[x]]<-min(s_area_norm_max60[,x])
  maxlist[[x]]<-max(s_area_norm_max60[,x])
  topquantile[[x]]<-quantile(s_area_norm_max60[,x],probs=c(.00,.95))[2]
}


bot1sd<-unlist(meanlist)- unlist(sdlist)
top1sd<-unlist(meanlist)+ unlist(sdlist)
topquantile<-unlist(topquantile)
par(cex=1.5)
plot(names(cvcols), unlist(meanlist), type="l", col="red", lwd=4, ylim=c(-0.2,1.1), 
     main="maxCV=60")
lines(names(cvcols),  bot1sd, col="grey" , lwd=2)
lines(names(cvcols),  top1sd, col="grey" , lwd=2)
lines(names(cvcols),  topquantile, col="blue", lwd=2 )



##### same but put the whole thing into a loop and go through the CV values 
### loop through each column and get the mean and standard deviation
getwd()  
setwd("P:/JGM_FAIMS_CVprediction/")
dev.off()
#par(cex=0.5, mfcol=c(4,6))

for(y in names(cvcols)){
  tmp<-s_area_norm[cvmaxrows[[y]],] ## get the rows that are max @ value 'y'
  
  meanlist<-list()
  sdlist<-list()
  minlist<-list()
  maxlist<-list()
  topquantile<-list()
  #cvcols
  #x<-"70"
  for(x in names(cvcols)){
    meanlist[[x]]<-mean(tmp[,x])
    sdlist[[x]]<-sd(tmp[,x])
    minlist[[x]]<-min(tmp[,x])
    maxlist[[x]]<-max(tmp[,x])
    topquantile[[x]]<-quantile(tmp[,x],probs=c(.00,.95), na.rm=T)[2]
  }
  
  bot1sd<-unlist(meanlist)- unlist(sdlist)
  top1sd<-unlist(meanlist)+ unlist(sdlist)
  topquantile<-unlist(topquantile)
  tiff(file=paste(y, ".tiff", sep=""))
  #par(cex=1.5)
  plot(names(cvcols), unlist(meanlist), type="l", col="red", lwd=4, ylim=c(-0.2,1.1), 
       main=paste("MaxCV=", y, sep=" "), sub=paste("n=", nrow(tmp),sep=" "))
  lines(names(cvcols),  bot1sd, col="grey" , lwd=2)
  lines(names(cvcols),  top1sd, col="grey" , lwd=2)
  lines(names(cvcols),  topquantile, col="blue", lwd=2 )
  dev.off()
}

nrow(tmp)
??plot
