# install.packages(seqinr)
library(seqinr)

traindata = read.csv(file="training_data.csv",header=TRUE,sep=",")
# traindata and testdata are "data frames" that contain variables of different types (numeric, character, etc)
# dim(traindata) = 1000 by 6
prseq = traindata[,3]
# length(prseq) = 920 
rtseq = traindata[,4]
# 1000 
train0 = traindata[traindata[,2]==0,]
# 794 rows (subjects)
train1 = traindata[traindata[,2]==1,]
# 206 rows (subjects)
levels(traindata$PR.Seq)
# 920
traindata[,3]=as.character(traindata[,3])
traindata[,4]=as.character(traindata[,4])
# convert factors in column 3 and 4 into character strings

testdata = read.csv(file="test_data.csv",header=TRUE,sep=",")
# dim(testdata) = 692 by 6

# seqinr package has write.fasta() function
# write.fasta(sequences, names, nbchar = 60, file.out, open = "w")

# writing fasta file for all subjects
prlist=traindata[1:920,3]
rtlist=traindata[,4]
for (i in 1:920){
  if (i==1){
    write.fasta(prlist[1],1,nbchar=60,file.out="pr",open="w")
  } else{
    write.fasta(prlist[i],i,nbchar=60,file.out="pr",open="a")
  }
}
for (i in 1:1000){
  if (i==1){
    write.fasta(rtlist[1],1,nbchar=60,file.out="rt",open="w")
  } else{
    write.fasta(rtlist[i],i,nbchar=60,file.out="rt",open="a")
  }
}

# Get row index of 0 patients and 1 patients
pr0index=c()
pr1index=c()
for (i in 1:920){
  if (traindata[i,2]==0){
    pr0index=c(pr0index,i)
  } else {
    pr1index=c(pr1index,i)
  }
}

rt0index=c()
rt1index=c()
for (i in 1:1000){
  if (traindata[i,2]==0){
    rt0index=c(rt0index,i)
  } else {
    rt1index=c(rt0index,i)
  }
}

# write fasta file for 0 and 1 patients separately
pr0=traindata[pr0index,3]
pr1=traindata[pr1index,3]
rt0=traindata[rt0index,4]
rt1=traindata[rt1index,4]
# pr0 and pr1
for (i in 1:length(pr0index)){
  if (i==1){
    write.fasta(pr0[1],pr0index[1],nbchar=60,file.out="pr0",open="w")
  } else{
    write.fasta(pr0[i],pr0index[i],nbchar=60,file.out="pr0",open="a")
  }
}
for (i in 1:length(pr1index)){
  if (i==1){
    write.fasta(pr1[1],pr1index[1],nbchar=60,file.out="pr1",open="w")
  } else{
    write.fasta(pr1[i],pr1index[i],nbchar=60,file.out="pr1",open="a")
  }
}

# rt0 and rt1
for (i in 1:length(rt0index)){
  if (i==1){
    write.fasta(rt0[1],rt0index[1],nbchar=60,file.out="rt0",open="w")
  } else{
    write.fasta(rt0[i],rt0index[i],nbchar=60,file.out="rt0",open="a")
  }
}
for (i in 1:length(pr1index)){
  if (i==1){
    write.fasta(rt1[1],rt1index[1],nbchar=60,file.out="rt1",open="w")
  } else{
    write.fasta(rt1[i],rt1index[i],nbchar=60,file.out="rt1",open="a")
  }
}



