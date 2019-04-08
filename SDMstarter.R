##Species Distribution Model
library(dismo)
data<-read.csv("distribution",header=TRUE)

library(maptools)
data(wrld_simpl,axes=T,col="light yellow")
points(data$longititude,data$latitude,col='red',cex=.75)

#environments
path<-list.files("Bioclimate")
files<-list.files(path, pattern='grd$',full.names=TRUE)

predictors<-stack(files)
predictors
names(predictore)
plot(predictors)

library(maptools)
data(wrld_simpl)
gene<-read.csv("/genotypes.csv",sep="")
gene<-gene[,-1]

plot(predictors,1)
plot(wrld_simpl,add=TRUE)
points(gene,col="red")

presvals <- extract(predictors,gene)
# setting random seed to always create the same
# random set of points for this example
set.seed(2019)
background <- randomPoints(predictors, 500)
absvals <- extract(predictors, background)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata[,'biome'] = as.factor(sdmdata[,'biome'])
head(sdmdata)
summary(sdmdata)

pairs(sdmdata[,-1],cex=0.1)
saveRDS(sdmdata, "sdm.Rds")
saveRDS(presvals, "pvals.Rds")
m1<-glm(pb~bio1+bio5+bio20+bio23,data=sdmdata)
m2<-glm(pb~.,data=sdmdata)
m2

bc<-bioclim(presvals[,c('bio1','bio5','bio20','bio23')])
class(bc)
bc
pairs(bc)

bio1=c(40,150,200)
bio5=c(60,115,290)
bio20=c(80,250,300)
bio23=c(50,200,300)
pd=data.frame(cbind(bio1,bio5,bio20,bio23))
pd
predict(m1,pd)
predict(bc,pd)

response(bc)


names(predictors)
## [1] "bio1"  "bio5" "bio20" "bio23" "biome"
p <- predict(predictors, m1)
plot(p)

p <- rnorm(50, mean=0.7, sd=0.3)
a <- rnorm(50, mean=0.4, sd=0.4)
par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
       pch=c(21,24), col=c('red', 'blue'))
comb <- c(p,a)
group <- c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~phylogroup, col=c('blue', 'red'))
group = c(rep(1, length(p)), rep(0, length(a))) 
cor.test(comb, phylogroup)$estimate
mv <- wilcox.test(p,a)
auc <- as.numeric(mv$statistic) / (length(p) * length(a))
auc
e <- evaluate(p=p, a=a)
class(e)
e
par(mfrow=c(1, 2))
density(e)
boxplot(e, col=c('blue', 'red'))
samp <- sample(nrow(sdmdata), round(0.75 * nrow(sdmdata)))
traindata <- sdmdata[samp,]
traindata <- traindata[traindata[,1] == 1, 2:9]
testdata <- sdmdata[-samp,]
bc <- bioclim(traindata)
e <- evaluate(testdata[testdata==1,], testdata[testdata==0,], bc)
e
plot(e, 'ROC')

pres <- sdmdata[sdmdata[,1] == 1, 2:9]
back <- sdmdata[sdmdata[,1] == 0, 2:9]
k <- 5
group <- kfold(pres, k)
group[1:10]
unique(group)
e <- list()
for (i in 1:k) {
  train <- pres[group != i,]
  test <- pres[group == i,]
  bc <- bioclim(train)
  e[[i]] <- evaluate(p=test, a=back, bc)
}   
auc <- sapply(e, function(x){x@auc}) 
auc
mean(auc)
sapply( e, function(x){ threshold(x)['spec_sens'] } )

set.seed(2019)
backgr <- randomPoints(predictors, 500)
nr <- nrow(gene)
s <- sample(nr, 0.25 * nr)
pres_train <- gene[-s, ]
pres_test <- gene[s, ]
nr <- nrow(backgr)
set.seed(2019)
s <- sample(nr, 0.25 * nr)
back_train <- backgr[-s, ]
back_test <- backgr[s, ]
sb <- ssb(pres_test, back_test, pres_train)
sb[,1] / sb[,2]
i <- pwdSample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2]
bc <- bioclim(predictors, pres_train)
evaluate(bc, p=pres_test, a=back_test, x=predictors)
evaluate(bc, p=pres_test_pwd, a=back_test_pwd, x=predictors)