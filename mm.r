# fetch data
tempDNA <- read.table('data.txt', header=TRUE);
DNA <- mutate(tempDNA, DNA = substr(DNA, 4, 1000));
rm(tempDNA);


# the list for modeling
originalDNA <- DNA$DNA[1:20];
# the list wait for clustering
futureDNA <- DNA$DNA[21:40];


library(stringr);
library(dplyr);
#take the frequency matrix
takefreqMatrix <- function(DNA){
  length = nchar(DNA[1]);
  a <- str_count(originalDNA, pattern = "a")/length;
  c <- str_count(originalDNA, pattern = "c")/length;
  t <- str_count(originalDNA, pattern = "t")/length;
  g <- str_count(originalDNA, pattern = "g")/length;
  
  return(cbind(a,c,t,g));
}

# STD transfrom
#afterSTD <- scale(takefreqMatrix(originalDNA));


#take the max-min standard matrix
takeMMS <- function (col) {
  return ((col-min(col)) / (max(col)-min(col)))
}

#max-min transform
#fuzzyMatrix <- apply(afterSTD, takeMMS, MARGIN = 2);

#get fuzzy similar martix
#FSM <- cor(t(fuzzyMatrix));
takeFSMfinal <- function(m) {
  for(i in 1:20) {
    for(j in 1:20) {
      if(m[i,j] < 0) {
        m[i,j] = (m[i,j] + 1) /2;
      }
    }
  }
  return (m)
}
# now I have got a fuzzy similar matrix

#FSMfinal <- takeFSMfinal(FSM)
# to generate a fuzzy equivalent matrix



multiSelf <- function(m, n=m) {
  new <- matrix(0, 20, 20)
  for(i in 1:20) {
    for(j in 1:20){
      new[i,j] = max(pmin(m[i, ], n[,j]));
    }
  }
  return(new)
}

takeFEM <- function (m) {
  mm = multiSelf(m);
  if(all(mm ==m)) {
    return (mm)
  } else {
    takeFEM(mm)
  }
}

takeLamdaMatrix <- function(m) {
  for(i in 1:20) {
    for(j in 1:20) {
      if(m[i,j] >= 0.8) {
        m[i,j] =1;
      } else {
        m[i,j] =0;
      }
    }
  }
  return(m)
}

takeData <- function(m){
  s <- matrix(NA, 20, 1);
  count = 0;
  for(i in 1:20) {
    for(j in 1:20) {
      if ((m[i,j] == 1) && (i != j)) {
        if((is.na(s[j, 1])) && (is.na(s[i, 1]))){
          s[i, 1] = count;
          count = count+1
        } else if(is.na(s[i, 1])) {
          s[i, 1] = s[j, 1];
        }
      }
    }
  }
  return(s)
}

analysis <- function(m) {
  afterSTD <- scale(takefreqMatrix(m));
  fuzzyMatrix <- apply(afterSTD, takeMMS, MARGIN = 2);
  FSM <- cor(t(fuzzyMatrix));
  FSMfinal <- takeFSMfinal(FSM);
  FEM <- takeFEM(FSMfinal);
  lamdaM <- takeLamdaMatrix(FEM);
  return(takeData(lamdaM))
}
takefreqTestMatrix <- function(DNA){
  length = nchar(DNA[1]);
  acg <- str_count(originalDNA, pattern = "acg")/length;
  cag <- str_count(originalDNA, pattern = "cag")/length;
  tag <- str_count(originalDNA, pattern = "tag")/length;
  atg <- str_count(originalDNA, pattern = "atg")/length;
  atc <- str_count(originalDNA, pattern = "atc")/length;
  acc <- str_count(originalDNA, pattern = "acc")/length;
  cac <- str_count(originalDNA, pattern = "cac")/length;

  #seqA <- (a-mean(a))/sd(a);
  #seqC <- (c-mean(c))/sd(c);
  #seqT <- (t-mean(t))/sd(t);
  #seqG <- (g-mean(g))/sd(g);
  
  #afterScale <- scale(cbind(a,c,t,g));
  #return(cbind(seqA,seqC, seqT, seqG));
  return(cbind(acg,cag,tag,atg, atc, acc, cac));
}