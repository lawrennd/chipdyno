#Sample#1
rm(list=ls())

load("/home/muhammad/H-drive/CElegans/Results_cElegans_Optim_Sample1.RData")

tfId=8 # To examin the activity of TransNames[tfId,]

source("chipDynoTransFactNoise.R")
source("chipDynoMaxDiff.R")

nTrans=nrow(TransNames);
lst=list();
newX=array(0, dim <-c(dim(X)));
newXVals=array(0, dim <-c(dim(X)));

i=tfId
expectations =chipDynoTransFactNoise(data,X,Sigma,beta, precs, gamma,mu, TransNames, annotation, TransNames[i,]);
S1_TF = expectations[[1]]
S1_TFError = expectations [[2]]
S1_TFErrorDiff = expectations [[3]]

###################

i=tfId
name = TransNames[i,]
index=which(name == TransNames);
genesIn=which(X[,index]!=0);
actGeneNames = annotation[genesIn]


# The graphical representation of "chipDynoActTransFact's" TF activity! 

plotErrorBar <- function(y,SE, actGeneName){

add.error.bars <- function(x,y,SE,w,col=1){
x0 = x; y0 = (y-SE); x1 =x; y1 = (y+SE);
arrows(x0, y0, x1, y1, code=3,angle= 90,length=w,col=col);
}

### Test
x <- c(1:length(y));
#Y <- TF[1,];
#SE <- TFError[1,];
#plot(x,y, type = 'l', col='green4', las=1);
#plot(x,y, type = 'l', col='green4', las=1, main/sub="gene specific TFA", xlab="t", ylab="TFA");

plot(x,y, type = 'l', col='green4', ylim=c(0,13.5), las=1, xlab="t", ylab="TFA", main=bquote("Gene : " ~ .(actGeneName)));
#plot(x,y, type = 'l', col='green4' las=1, xlab="t", ylab="TFA");
add.error.bars(x,y,SE,0.05,col='red');

###
}

##load("test130313_2.RData")
# source("plotErrorBar.R")
#plotErrorBar(TF[1,],TFError[1,]);

M <- matrix(c(rep(1:6)), byrow=TRUE, nrow=2) # Choose the position by matrix setting!
layout(M) 

#jpeg('rplot.jpg')
#bmp("testplot.bmp")
#png(filename="TFA_data_sample1.png")
for (i in 1:6){
plotErrorBar(S1_TF[i,],S1_TFError[i,],actGeneNames[i])
}
#dev.off()