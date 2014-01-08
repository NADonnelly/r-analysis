# Intro ------------------------------------------------------------------

#This script, or perhaps later a function, will take a set of data from the LFP recordings,
#including LFP power, movement and task data, and then apply lme or lmer and attempt to 
#determine which factors are most important in determining LFP power. I will start by 
#looking at the instantaneous power, but I will think about whether using a short moving
#window might also be helpful.


# load required packages --------------------------------------------------

library(car)
library(gmodels)
library(MASS)
library(reshape2)
library(ggplot2)
library(miscTools)


# load up the data --------------------------------------------------------

event  <- "leav"
freq   <- "Gamma60"
mType  <- "lmer"
struct <- c("Core","Shell","Prl","Il")
flav   <- c("Low","High")

setwd("C://Users/Nick/Documents/Dropbox/PhD/Work/data/LFP_data")

switch(mType,         
       nlme = {
                library(nlme)          
       },       
       lmer = {
                library(lme4)
       }
)


#preallocate a list
ml <- vector(mode = "list", length = length(flav)*length(struct))
nl <- vector(mode = "list", length = length(flav)*length(struct))
vl <- vector(mode = "list", length = length(flav)*length(struct))


cnt <- 0

for (i in 1:length(struct)){  
  for (j in 1:length(flav)){

      
    dataNL <- read.delim(paste(event, struct[i], freq, flav[j], "Neuro.txt",sep = ""),header = F)
    dataNL <- data.matrix(dataNL)
    
    dataVL <- read.delim(paste(event, struct[i], freq, flav[j],  "Velo.txt",sep = ""),header = F)
    dataVL <- data.matrix(dataVL)
    
    dataML      <- read.csv(paste(event, struct[i], freq, flav[j], "Meta.dat",sep = ""),header = T)    
    dataML$sOrd <- factor(dataML$sOrd)
    dataML$tOrd <- factor(dataML$tOrd)
    
    #here I actually need to build in a check to cut out rats affected by the poke and tray artefacts...
    #I have removed this given the new artefact removal system
    
#     if (flav[j] == "Low" & event == "poke" & (freq == "Theta" | freq == "Delta" | freq == "Beta")){
#       
#       dataNL <- dataNL[dataML$ratID != "94_03" & dataML$ratID != "97_14" & dataML$ratID != "97_16",]
#       dataVL <- dataVL[dataML$ratID != "94_03" & dataML$ratID != "97_14" & dataML$ratID != "97_16",]
#       dataML <- subset(dataML,dataML$ratID != "94_03" & dataML$ratID != "97_14" & dataML$ratID != "97_16")            
#     }      
        
    cnt <- cnt + 1
      
    ml[[cnt]] <- dataML
    nl[[cnt]] <- dataNL
    vl[[cnt]] <- dataVL
      
    rm(list = c("dataML","dataNL","dataVL"))          
  }    
}

dataM <- do.call("rbind",ml)
dataN <- do.call("rbind",nl)
dataV <- do.call("rbind",vl)

rm(list = c("ml","nl","vl"))

dataM$ps <- with(dataM, interaction(phenotype, structure))


# basic plotting ----------------------------------------------------------

nResult <- data.frame(x = seq(-5, 5, b = 1/25), 
                      correct_c   = colMedians(dataN[dataM$outcome == " correct"   & dataM$rHX == " c",]),
                      correct_x   = colMedians(dataN[dataM$outcome == " correct"   & dataM$rHX == " x",]),
                      incorrect_c = colMedians(dataN[dataM$outcome == " incorrect" & dataM$rHX == " c",]),
                      incorrect_x = colMedians(dataN[dataM$outcome == " incorrect" & dataM$rHX == " x",]),
                      premature_c = colMedians(dataN[dataM$outcome == " premature" & dataM$rHX == " c",]),
                      premature_x = colMedians(dataN[dataM$outcome == " premature" & dataM$rHX == " x",]))

                      
nResultMelt <- melt(nResult, id = "x")

plotCols2 <- c("dodgerblue1","dodgerblue3",
               "springgreen1","springgreen4",
               "firebrick1","firebrick3")

np.1 <- ggplot(nResultMelt, aes(x = x, y = value, colour = variable)) 
np.1 + geom_line(size = 1) + 
  coord_cartesian(xlim = c(-4.5, 4.5)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "median LFP z score, all data")+   
  scale_colour_manual(values = plotCols2)


# fancy LFP plot ----------------------------------------------------------


#preallocate a list
l <- vector(mode = "list", length = (nlevels(dataM$ps)*nlevels(dataM$outcome)*nlevels(dataM$rHX)))

cnt <- 0

pL <- levels(dataM$ps)
oL <- levels(dataM$outcome)
rL <- levels(dataM$rHX)

for (i in 1:nlevels(dataM$ps)){  
  for (j in 1:nlevels(dataM$outcome)){
    for (k in 1:nlevels(dataM$rHX)){                    
      
      temp <- data.frame(x = seq(-5, 5, b = 1/25), estimate = colMeans(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$ps == pL[i],]))
      temp$lower <- temp$estimate - apply(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$ps == pL[i],],2,sd)/sqrt(dim(dataN)[2])
      temp$upper <- temp$estimate + apply(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$ps == pL[i],],2,sd)/sqrt(dim(dataN)[2])
      temp$type  <- rep(paste(oL[j],rL[k],sep = "_"),dim(temp)[1])
      temp$ps    <- rep(pL[i],dim(temp)[1])
      
      cnt      <- cnt + 1      
      l[[cnt]] <- temp      
      rm(temp)
            
    }
  }    
}


nRibbon <- do.call("rbind",l)


plotCols2 <- c("dodgerblue1","dodgerblue3",
               "springgreen1","springgreen4",
               "firebrick1","firebrick3")

r1 <- ggplot(nRibbon,aes(x = x,y = estimate,colour = type))
r1 <- r1 + geom_ribbon(aes(ymin = lower, ymax = upper,fill = type),alpha = 0.25) + 
           geom_line(size = 1) + 
           geom_hline(yintercept = 0)  + 
           geom_vline(xintercept = 0)  + 
           facet_wrap(~ ps, ncol = 4) +  
           coord_cartesian(xlim = c(-1, 2)) +
           theme(panel.grid.minor = element_line(colour="white", size=0.5)) + 
           scale_x_continuous(minor_breaks = seq(-4.5, 4.5, 0.5)) +
           labs(x = "time", y = paste("z score ", freq, " power",sep = ""), title = paste("z score ", freq, " power over time, all rats",event, paste(struct,collapse = " ") , sep = " ")) + 
           scale_fill_manual(values = plotCols2) + 
           scale_colour_manual(values = plotCols2) + 
           theme(panel.background = element_rect(fill = "white")) + 
           theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) + 
           theme(strip.background = element_rect(fill = "white")) + 
           theme(plot.title = element_text(size = 20)) + 
           theme(axis.title.y = element_text(size = 20)) + 
           theme(axis.title.x = element_text(size = 20))
r1

ggsave(r1, file=paste("C://Users/Nick/Documents/Dropbox/PhD/Work/data/LFP_data/figures/", freq,"Multi", event,"HandL.png",sep = ""), width=16, height=10)


# fancy tracking plot -----------------------------------------------------

#fancier plot
#preallocate a list
l <- vector(mode = "list", length = (nlevels(dataM$phenotype)*nlevels(dataM$outcome)*nlevels(dataM$rHX)))

cnt <- 0

pL <- levels(dataM$phenotype)
oL <- levels(dataM$outcome)
rL <- levels(dataM$rHX)

for (i in 1:nlevels(dataM$phenotype)){  
  for (j in 1:nlevels(dataM$outcome)){
    for (k in 1:nlevels(dataM$rHX)){                    
      
      temp <- data.frame(x = seq(-5, 5, b = 1/25), estimate = colMeans(dataV[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],]))
      temp$lower <- temp$estimate - apply(dataV[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],],2,sd)/sqrt(dim(dataV)[2])
      temp$upper <- temp$estimate + apply(dataV[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],],2,sd)/sqrt(dim(dataV)[2])
      temp$type  <- rep(paste(oL[j],rL[k],sep = "_"),dim(temp)[1])
      temp$flav  <- rep(pL[i],dim(temp)[1])
      
      cnt <- cnt + 1
      
      l[[cnt]] <- temp
      
      rm(temp)
      
    }
  }    
}

vRibbon <- do.call("rbind",l)

plotCols2 <- c("dodgerblue1","dodgerblue3",
               "springgreen1","springgreen4",
               "firebrick1","firebrick3")

r2 <- ggplot(vRibbon,aes(x = x,y = estimate,colour = type))
r2 <- r2 + geom_ribbon(aes(ymin = lower, ymax = upper,fill = type),alpha = 0.25) + 
           geom_line(size = 1) + 
           geom_hline(yintercept = 0)  + 
           geom_vline(xintercept = 0)  + 
           facet_grid(flav ~ .) +  
           coord_cartesian(xlim = c(-2, 2)) +
           theme(panel.grid.minor = element_line(colour="white", size=0.5)) + 
           scale_x_continuous(minor_breaks = seq(-4.5, 4.5, 0.5)) +
           labs(x = "time", y = "velocity", title = paste("velocity over time, all rat/trials",event, sep = " ")) + 
           scale_fill_manual(values = plotCols2) + 
           scale_colour_manual(values = plotCols2) +   
           theme(panel.background = element_rect(fill = "white")) + 
           theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0)) + 
           theme(strip.background = element_rect(fill = "white")) + 
           theme(plot.title = element_text(size = 10)) + 
           theme(axis.title.y = element_text(size = 10)) + 
           theme(axis.title.x = element_text(size = 10))
r2

ggsave(r2, file=paste("C://Users/Nick/Documents/Dropbox/PhD/Work/data/LFP_data/figures/veloMulti", event,"HandL.png",sep = ""), width=16, height=10)



# main loop ---------------------------------------------------------------

#allocate the model formula
# a <- c("outcome","rHX", "velo", "phenotype","structure","outcome:structure","rHX:structure",
#        "phenotype:rHX","phenotype:outcome","phenotype:structure","phenotype:rHX:structure")

a <- c("outcome","rHX", "velo", "phenotype","structure","outcome:structure","rHX:structure","rHX:outcome",
       "phenotype:rHX","phenotype:outcome","phenotype:structure","phenotype:rHX:outcome","phenotype:outcome:structure")

#allocate a list to put the models into
mO <- vector(mode = "list", length = dim(dataN)[2])

#Now loop through every time point
for (i in 1:dim(dataN)[2]){
  
  dataT <- dataM
  
  dataT$LFP  <- dataN[,i]
  dataT$velo <- dataV[,i]
  
  #Get rid of times when the rat has a velocity of zero
  dataT <- subset(dataT, velo > 0)
  
  contrasts(dataT$outcome)   <- contr.Helmert(nlevels(dataT$outcome))
  contrasts(dataT$phenotype) <- contr.Helmert(nlevels(dataT$phenotype))
  contrasts(dataT$rHX)       <- contr.Helmert(nlevels(dataT$rHX))
  
  if (nlevels(dataT$structure) > 1){
    contrasts(dataT$structure) <- contr.Helmert(nlevels(dataT$structure)) 
  }
  
  switch(mType,         
         nlme = {
           #fit a simple model
           m1 <- lme(as.formula(paste("LFP ~ ", paste(a,collapse = " + "),sep = "")), 
                     random =  ~velo|ratID/sOrd/channel, data = dataT, 
                     na.action = na.exclude,control = list(opt="optim"))  
           #,control = list(opt="optim")           
           mO[[i]] <- m1
           
         },
         
         lmer = {
           #Try something more complicated with lmer             
           m1 <- lmer(as.formula(paste("LFP ~ ", paste(a,collapse = " + "), " + (velo|ratID/sOrd/channel)",sep = "")),
                      data = dataT, 
                      na.action = na.exclude,
                      REML = 1)             
           mO[[i]] <- m1
           
         }
  )
  
  print(i)
  
}


#save the list of models
save(mO,file=paste("C://Users/Nick/Documents/Dropbox/PhD/Work/data/LFP_data/modelSet", freq, event,"HandL.RData",sep = ""))


#Now we have all these models we can use lapply to get out the bits we want
fMat <- do.call("rbind", lapply(mO, function(x){
  y <- Anova(x,type = 3)
  y[2:nrow(y),1]
}))

pMat <- do.call("rbind", lapply(mO, function(x){
  y <- Anova(x,type = 3)
  y[2:nrow(y),3]
}))


# Plotting and summarising ------------------------------------------------

#get the model factors out form the output matrixes
fstatOutput <- data.frame(x = seq(-5, 5, by = 1/25)) 
for (i in 1:length(a)){ 
  fstatOutput[,paste(a[i])] <- fMat[,i] 
}

pstatOutput <- data.frame(x = seq(-5, 5, by = 1/25)) 
for (i in 1:length(a)){ 
  pstatOutput[,paste(a[i])] <- pMat[,i] 
}


#melt the fstats for plotting
fstatMelt <- melt(fstatOutput,id = "x")


#Quick temporary plot
qplot(x,value,group = variable,geom = "line", data = fstatMelt, facets = variable ~ .)


#fiddle with the p stats
pstatMelt <- melt(pstatOutput,id = "x")

fstatMelt$id <- "fstat"
pstatMelt$id <- "pstat"

# I need to find a way of sorting out these p values so that we only get a "significant"
# result if p is less than 0.01 and also sustained forl onger than 100 ms

#Get the significance data
pW <- pstatMelt$value < 0.01

#Now we want to loop through each effect we've asked for, get the p values for that effect, and 
#check if sequences of significant effects
aL <- levels(pstatMelt$variable)

sMat <- matrix(data=NA,nrow=dim(dataN)[2],ncol=length(aL))

for (i in 1:length(aL)){
  
  pWT <- pW[pstatMelt$variable == aL[i]]*1
  pWD <- diff(pWT)
  
  onTimes  <- seq(along=pWD)[pWD ==  1]
  offTimes <- seq(along=pWD)[pWD == -1]
  
  
  #if there are no noted on times and the first index is a sig, then
  #all indices are sigs, so we quickly note this and skip to the next
  #variable
  if (length(onTimes) == 0 & pWT[1]){
    sMat[,i] <- 1
    next    
  } else if(length(onTimes) == 0 ){
    next    
  }
  
  #Attempt to address epochs that start/end outside of our window  
  if (length(onTimes) != length(offTimes)){
    if (offTimes[1] < onTimes[1]){
      onTimes <- c(1,onTimes)      
    } else if(tail(onTimes,1) > tail(offTimes,1)){
      offTimes <- c(offTimes,dim(dataN)[2])
    }    
  }
  
  if (offTimes[1] < onTimes[1]){
    onTimes <- c(1,onTimes)      
  }
  
  if(tail(onTimes,1) > tail(offTimes,1)){
    offTimes <- c(offTimes,dim(dataN)[2])
  }  
  
  #Ensure epochs are continuous and free of silly errors
  for (h in 1:length(onTimes)){
    if ((onTimes[h] < offTimes[h]) & all(as.logical(pWT[(onTimes[h]+1):offTimes[h]]))){
      next
    } else {
      onTimes[h]  <- NA
      offTimes[h] <- NA
    }          
  }
  
  onTimes  <- onTimes[!is.na(onTimes)]
  offTimes <- offTimes[!is.na(offTimes)]
  
  if (length(onTimes) == 0){
    next    
  }
  
  #Now we've worked out our epochs are all valid, lets censor those that 
  #are shorter than 200 ms
  for (h in 1:length(onTimes)){
    
    tD <- pstatOutput$x[offTimes[h]]-pstatOutput$x[onTimes[h]+1]
    
    if (tD > 0.2){
      next
    } else {
      onTimes[h]  <- NA
      offTimes[h] <- NA
    }          
  }
  
  onTimes  <- onTimes[!is.na(onTimes)]
  offTimes <- offTimes[!is.na(offTimes)]
  
  if (length(onTimes) == 0){
    next    
  }
    
  #Finally convert the selected indices to 1s
  for (h in 1:length(onTimes)){    
    sMat[(onTimes[h]+1):offTimes[h],i] <- 1   
  }  
}

sMat <- c(sMat)

pstatMelt$value <- sMat*500

sMelt <- rbind(fstatMelt, pstatMelt)

#Finally do a plot
mp.1 <- ggplot(sMelt, aes(x = x, y = value, group = id, colour = variable)) 
mp.1 <- mp.1 + geom_line(size = 1) + 
  facet_grid(variable ~ .) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  coord_cartesian(xlim = c(-2, 2)) +
  theme(panel.grid.minor = element_line(colour="white", size=0.5)) + 
  scale_x_continuous(minor_breaks = seq(-4.5, 4.5, 0.5)) +
  labs(x = "time", y = "F", title = "model statistics over time for different explanatory variables") + 
  theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0)) + 
  theme(strip.background = element_rect(fill = "white")) + 
  theme(panel.background = element_rect(fill = "white")) + 
  theme(plot.title = element_text(size = 10)) + 
  theme(legend.position="none") +  
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10))
mp.1


ggsave(mp.1, file=paste("C://Users/Nick/Documents/Dropbox/PhD/Work/data/LFP_data/figures/pstatMulti", freq, event,"HandL+velo.png",sep = ""), width=16, height=12)





# parameters plot ---------------------------------------------------------


#Plot the parameter estiamtes if this is interesting
# paramOutput <- data.frame(x = seq(-5, 5, by = 1/25)) 
# for (i in 1:length(b)){ 
#   paramOutput[,paste(b[i])] <- paramMat[,i] 
# }
# 
# paramMelt <- melt(paramOutput, id = "x")
# 
# qplot(x,value,group = variable,geom = "line", data = paramMelt, facets = variable ~ . )
# 
# 
# # #Do some faffing to make the ribbon plots over facets work
# 
# l <- vector(mode = "list", length = length(b))
# 
# cnt <- 0
# 
# 
# for (i in 1:length(b)){  
#               
#       
#       temp <- data.frame(x = paramOutput$x)
#       temp[,paste(b[i])] <- paramMat[,i]
#                           estimate = paramOutput$outcomeIe, 
#                                lower    = paramOutput$outcomeIl,
#                                upper    = paramOutput$outcomeIu)
#         
#         x = seq(-5, 5, b = 1/25), estimate = colMeans(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],]))
#       temp$lower <- temp$estimate - apply(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],],2,sd)/sqrt(dim(dataN)[2])
#       temp$upper <- temp$estimate + apply(dataN[dataM$outcome == oL[j] & dataM$rHX == rL[k] & dataM$phenotype == pL[i],],2,sd)/sqrt(dim(dataN)[2])
#       temp$type  <- rep(paste(oL[j],rL[k],sep = "_"),dim(temp)[1])
#       temp$flav  <- rep(pL[i],dim(temp)[1])
#       
#       cnt <- cnt + 1
#       
#       l[[cnt]] <- temp
#       
#       rm(temp)
#  
# }
# 
# nRibbon <- do.call("rbind",l)
# 

# outcomeI.F <- data.frame(x = paramOutput$x, 
#                          estimate = paramOutput$outcomeIe, 
#                          lower    = paramOutput$outcomeIl,
#                          upper    = paramOutput$outcomeIu)
# outcomeI.F$type <-rep("I v C",dim(outcomeI.F)[1])
# 
# outcomeP.F <- data.frame(x = paramOutput$x, 
#                          estimate = paramOutput$outcomePe, 
#                          lower    = paramOutput$outcomePl,
#                          upper    = paramOutput$outcomePu)
# outcomeP.F$type <-rep("P v C",dim(outcomeP.F)[1])
# 
# outcomeIP.F <- data.frame(x = paramOutput$x, 
#                          estimate = paramOutput$outcomeIPe, 
#                          lower    = paramOutput$outcomeIPl,
#                          upper    = paramOutput$outcomeIPu)
# outcomeIP.F$type <-rep("P v I",dim(outcomeIP.F)[1])
# 
# rHXx.F <- data.frame(x = paramOutput$x, 
#                      estimate = paramOutput$rHXxe, 
#                      lower    = paramOutput$rHXxl,
#                      upper    = paramOutput$rHXxu)
# rHXx.F$type <-rep("RH/x v RH/c",dim(rHXx.F)[1])
# 
# 
# paramRibbon <-rbind(outcomeI.F,outcomeP.F,outcomeIP.F,rHXx.F)
# 
# r1 <- ggplot(paramRibbon,aes(x = x,y = estimate,colour = type))
# r1 + geom_ribbon(aes(ymin = lower, ymax = upper,fill = type),alpha = 0.25) + 
#   geom_line(size = 1) + 
#   geom_hline(yintercept = 0)  + 
#   geom_vline(xintercept = 0)  + 
#   facet_grid(type ~ .) +  
#   coord_cartesian(xlim = c(-4.5, 4.5)) +
#   theme(panel.grid.minor = element_line(colour="white", size=0.5)) + 
#   scale_x_continuous(minor_breaks = seq(-4.5, 4.5, 0.5)) +
#   labs(x = "time", y = "parameter estimate", title = "parameter estimates over time for different explanatory variables")

 

# Using path analysis... --------------------------------------------------

#What I'm goign to try and do is use path analysis to, out of interest, see what the direction of causation
#is forthe velocity and gamma- thisi s going to be very very simplistic and at first wrong
#because I'm goign to fit SEM to the whole dataset, disregarding the hieracharchys of rat ID
#and session and channel. Later on I might model these separately


#First take some data from just after the nose poke...

i <- 135 #~400ms post-poke

dataT <- dataM

dataT$LFP  <- dataN[,i]
dataT$velo <- dataV[,i]

#Get rid of times when the rat has a velocity of zero
dataT <- subset(dataT, velo > 0)

#Define a simple path model where outcome causes both velocity and LFP, and velocity also causes LFP

dataT$outcomeN <- rep(0,dim(dataT)[1])
dataT$outcomeN[dataT$outcome == " correct" ] = 1
dataT$outcomeN <- factor(dataT$outcomeN,ordered = T)

library(lavaan)

#Specify the model

mymodel<-'LFP  ~ outcomeN + velo
          velo ~ outcomeN'

mymodel<-'LFP  ~ outcomeN
          velo ~ outcomeN'

mymodel<-'LFP ~ velo'
mymodel<-'velo ~ LFP'

#Fit the model

fit.path<-sem(mymodel,data=dataT)

#Summary of the model fit

summary(fit.path,fit.measures = 1)

#See the coefficients

coef(fit.path)

#In a data frame

parameterEstimates(fit.path)

# autoregressive model... -------------------------------------------------

#this will be the same as the main loop but I'm going to try and use an autoregressive 
#covariance structure and incorporate a sort of moving-window time thing


#preallocate the output frames/matrix for f and p values

#a is the factors used in the model
a  <- c("outcome","rHX", "velo", "phenotype", "phenotype:rHX","phenotype:outcome","phenotype:velo")
a1 <- c("outcome","rHX", "velo", "phenotype")

fMat <- matrix(data=NA,nrow=dim(dataN)[2],ncol=length(a))
pMat <- matrix(data=NA,nrow=dim(dataN)[2],ncol=length(a))


#Now loop through every time point
for (i in 5:dim(dataN)[2]){
  
  dataT <- dataM
  
  dataT$L1  <- dataN[,i]
  dataT$L2  <- dataN[,i-1]
  dataT$L3  <- dataN[,i-2]
  dataT$L4  <- dataN[,i-3]
  dataT$L5  <- dataN[,i-4]
  
  dataT$velo <- dataV[,i]
  
  #Get rid of times when the rat has a velocity of zero
  dataT <- subset(dataT, velo > 0)
  
  dataT <- melt(dataT, id.vars = c(a1,"ratID","channel","sOrd"), measure.vars = c("L1","L2","L3","L4","L5"), variable.name = "tWin", value.name = "LFP")
  
  #fit a simple model
  m1 <- lme(LFP ~ tWin,
            random =  ~tWin|ratID/channel/sOrd, data = dataT, 
            na.action = na.exclude,control = list(opt="optim"), 
            correlation = corAR1())  
  
  
  
  
  ARModel<-lme(Life_Satisfaction~Time, data = restructuredData, 
               random = ~Time|Person, 
               correlation = corAR1(0, form = ~Time|Person), 
               method = "ML", na.action = na.exclude, 
               control = list(opt="optim"))
  summary(ARModel)
  
  
  #harvest model parameters- f stats and p vals
  aM1 <- anova(m1) 
  for (h in 1:length(a)){
    fMat[i,h] <- aM1[h+1,3]
    pMat[i,h] <- aM1[h+1,4]*dim(dataN)[2]
  }
  
  rm(list = c("aM1","dataT"))
  
}


# Look at the relationship between latency and gamma60 leav ---------------

xAx  <- seq(-5, 5, by = 1/25)
bins <- seq(0,5,by = 1/4)

dataM2 <- dataM

dataM2$leavLFP <- rowMeans(dataN[,xAx>0 & xAx <= 0.1])
dataM2$OR       <- with(dataM2, interaction(outcome, rHX))

sp.1 <- ggplot(dataM2, aes(x = T2T, y = log10(leavLFP), colour = OR)) 
sp.1 <- sp.1 + geom_point(size = 1) + 
  facet_grid(ps ~ .) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  coord_cartesian(xlim = c(-2, 2)) +
  theme(panel.grid.minor = element_line(colour="white", size=0.5)) + 
  scale_x_continuous(minor_breaks = seq(-4.5, 4.5, 0.5)) +
  labs(x = "time", y = "F", title = "model statistics over time for different explanatory variables") + 
  theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0)) + 
  theme(strip.background = element_rect(fill = "white")) + 
  theme(panel.background = element_rect(fill = "white")) + 
  theme(plot.title = element_text(size = 10)) + 
  theme(legend.position="none") +  
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10))
mp.1

#Lets look at just the core/correct/low/c trials
dataCLC <- dataM2[dataM2$outcome == " correct" & dataM2$phenotype == " Low" & dataM2$structure == " core" & dataM2$rHX == " c",]

sp.1 <- ggplot(dataCLC, aes(x = T2T, y = log10(leavLFP))) 
sp.1 + geom_point(size = 1) + theme(panel.background = element_rect(fill = "white")) 
  