##----------------------------------------------------------##
## Data analysis of the OECD unemployment rates (see paper)
##----------------------------------------------------------##

rm(list=ls())

(oecd = read.table("oecd_unemployment_rates.txt",header=T,sep="\t",row.names=1)) 
dim(oecd)

## Plot the series (delete NAs)
##-----------------------------##
par(ask=T,las=1)
for(ii in 1:nrow(oecd)){
  norig = na.omit(t(as.vector(oecd[ii,]))); print(1960 + ncol(oecd) - length(norig))
  (country = ts(norig,start=c(1960 + ncol(oecd) - length(norig),1),end=c(2010,1),fre=1))
  #country = log(country)
  plot(country,type="l",main=paste("OECD unemployment rate of ",rownames(oecd)[ii],sep=""),font.main=12,lwd=2,ylab="in %")
}


## Test series for unit root using the RJMCMC scheme for model selection:
##--------------------------------------------------------------------------##

(bayes.urtest = dget("bayes.urtest.R"))
start_bpoints = dget("start_bpoints.R")


## Anpassen vor Beginn der Session!!
##------------------------------------##
folder = "OECD_analysis_transf2"						#Target folder
label = NULL							#Label for simulation design
(path = paste("..\\Outputs\\",folder,"\\",sep=""))
dir.create(path)						#create folder in specified directory with given name
file = folder						#Anfang aller Filenamen, ob Plots oder txt-Files
##------------------------------------##

ii=1
(ii=ii+1)
#(ii=ii-1)

dir.create(paste(path,rownames(oecd)[ii],"\\",sep=""))		#für jedes Land eigenes Verzeichnis anlegen
(path = paste(path,rownames(oecd)[ii],"\\",sep=""))

(norig = na.omit(t(as.vector(oecd[ii,])))); periods = ncol(oecd)			#available observations and all number of years
(series = ts(norig,start=c(1960 + periods - length(norig),1),end=c(2010,1),fre=1))
#(series = series/100);  
#(lseries = log(series))
#(lseries = log(series/(1-series))) 
(lseries = log(series/(100-series)))

cat("\n",rownames(oecd)[ii],"// Start:",1960 + periods - length(norig),"\n\n")

#pdf(paste(path,file,"_",rownames(oecd)[ii],"_orig.pdf",sep="")) 

layout(matrix(c(1,1,2,3),ncol=2,byrow=T)); par(las=1,font.main=11,cex.main=1.2)
#par(mfrow=c(2,1),lwd=1,font.main=4,las=1,lwd=1,ask=F)
plot(series,type="l",main=paste("OECD unemployment rate of ",rownames(oecd)[ii],sep=""),font.main=12,lwd=2,ylab="in %")
#plot(lseries,type="l",main=paste("OECD (log) unemployment rates of",rownames(oecd)[ii]),font.main=12,lwd=2,ylab=expression(logs))

#series=series/100 ;par(mfrow=c(3,1));plot(density(series));plot(density(log(series)));plot(density(log(series/(1-series))))

MHsim = 10000; burn = 2000 ; determc="both"
#add some uniformly distrib. noise via jitter over the starting seq. of break dates (for some variation)
#if(m_max>0){(appetizer = start_bpoints(m_max=m_max,nobs=length(series),enter.dates=floor(jitter(seq(20,length(series)-10,length.out=m_max)))))}
grid = round(seq(from=.01,to=1.2,by=.001),digit=5)

set.seed(12)
(results = bayes.urtest(data=series, PLife=0.2,sigma=sqrt(length(series)),
                        transfdata=lseries, MMm=10^(.5), MMp=10^(0), m_fix=NULL, p_fix=NULL, m_max=5, p_max=5, determc=determc, burn=burn, MHsim=MHsim, sprob = 0.975, grid=grid , v0 = 2.001, lambda0 = 0.001, chm=1.1, chp=1.1, Jm_scale = 5, Jp_scale = 5, thinM = 2, thinP = 2))


## Individual plotting areas (ADF case) for nice plotting:
Lims = switch(ii,
              Lim1 = c(.1,.9),
              Lim2 = c(.1,.8),	
              Lim3 = c(.1,1),
              Lim4 = c(.2,.9),
              Lim5 = c(.5,1),
              Lim6 = c(.1,.9),
              Lim7 = c(.6,1),
              Lim8 = c(.1,1.1),	
              Lim9 = c(.4,1),
              Lim10 = c(.1,.8),
              Lim11 = c(.9,1.05),
              Lim12 = c(.5,1),
              Lim13 = c(.3,.9),
              Lim14 = c(.75,1.05),
              Lim15 = c(.3,1),
              Lim16 = c(.1,1.05),
              Lim17 = c(0.4,1))

ups = switch(ii,
             up1 = c(.15,1.2),
             up2 = c(.15,1.2),	
             up3 = c(.1,1.1),
             up4 = c(.15,1.15),
             up5 = c(.15,1.2),
             up6 = c(.1,1.1),
             up7 = c(.2,1.2),
             up8 = c(.1,1.2),	
             up9 = c(.2,1.2),
             up10 = c(.1,1.2),
             up11 = c(.7,1.2),
             up12 = c(.2,1.2),
             up13 = c(.15,1.2),
             up14 = c(.25,1.2),
             up15 = c(.15,1.1),
             up16 = c(.15,1.2),
             up17 = c(.15,1.2))

plot(results$grid,results$post_Norm,type="l",xlim=Lims,axes=T,ylim=c(min(results$post_Norm)+ups[1],max(results$post_Norm)*ups[2]),lty=1,lwd=1,col="black",font.main=4,xlab="",ylab="density",sub=paste("(p = ",results$p,")",sep=""))
title(paste("Marginal posterior densities\n '",rownames(oecd)[ii],"'",sep=""))
polygon(results$grid,results$post_Norm,col="lightgrey",border="black",lty=1)
lines(results$grid,results$post_Jeff,col="blue4",ylab="",lty=2,lwd=2,xlab="",main="Posterior with Jeffreys Prior")
#points(results$grid,results$post_Jeff,col="blue4",pch=1,cex=1.1)
lines(results$grid,results$post_BY,col="sienna4",ylab="",lty=3,lwd=2,xlab="",main="Posterior with Berger & Yang's Prior")
#points(results$grid,results$post_BY,col="sienna4",pch=4,cex=1.1)
lines(results$grid,results$post_LU,col="red4",ylab="",lty=4,lwd=2,xlab="",main="Posterior with Berger & Yang's Prior")
#points(results$grid,results$post_LU,col="red4",pch=0,cex=1.1) 
legend("topleft",legend=c("Nor.","Jef.","BY","LU"),col=c("black","blue4","sienna4","red4"),lwd=2,bty="n",lty=c(1,2,3,4), pch=c(NA,NA,NA,NA))

#axis(1);axis(2); box(lty='1373')

#graphics.off()

#(results = dget(paste(path,file,"_",rownames(oecd)[ii],".R",sep="")))

(gg = results$m)
(pp = results$p)

#dput(results,paste(path,file,"_",rownames(oecd)[ii],".R",sep=""))

(lags = results$lags) 
(sigma2 = results$s2eps)
(nofb = results$nofb)
(accrateA = results$accrateA)
(accrateB = results$accrateB)
(theta = results$BETAS);		#Draws of all model parameters after burnin (only given that p_fix and m_fix were specified above)
(betas = results$savebetas);
ynew = results$y
(indi = results$indi)
Xnew = results$X
(postdraws_mp = results$postdraws_mp)
(dates = results$bpoints[[gg]][[pp]])
(k_ma = results$kmap1)
(k_map = results$final_kmap)

(modelfreq = round(table(nofb)/length(nofb),4))
(lagsfreq = round(table(lags)/length(lags),4))

lagsfreqA <- lagsfreq
lagsfreqA
modelfreqA <- modelfreq
modelfreqA

#pdf(paste(path,file,"_",rownames(oecd)[ii],"_acf_periodo.pdf",sep=""))
layout(matrix(c(1,2,3,3), ncol=2, nrow=2, byrow = T));acf(series);pacf(series,main="");spectrum(series)
#layout.show(n=3)
graphics.off()

## Falls 'm' und 'p' fix, können hier Trajectorien aller beta komponenten gespeichert werden, z.B. für spätere Plots: 
#dput(theta,paste(path,file,"_",rownames(oecd)[ii],"thetas_traces.R",sep=""))


## Plot marginal! posterior distributions of number of breaks 'm'
## and number of AR lags 'p'
##---------------------------------------------------------------##
#win.graph(width=6.5,height=3.5); 

pdf(paste(path,file,"_",rownames(oecd)[ii],"_post_nofb_lags.pdf",sep=""))

par(ask=F,font.main=11,mfrow=c(1,2));

if(max(modelfreqA)<=.8){ybound = c(0,max(modelfreqA)+.05)} else{ybound=NULL}
barplot(modelfreqA,
        axes = F,axisnames = TRUE,space=1,cex.main=1.4,xlab=expression(paste("Number of breaks 'm'")),
        col="sienna3",density=100,border="black",ylim=ybound,ylab="Probabilities")
title("Posterior probability function")
#title(expression(paste("Posterior probability function of the number of breaks")))
axis(2,las=1) #;axis(1,labels = F)
abline(h=0)

#pdf(paste(path,file,"_",rownames(oecd)[ii],"_hist_lags.pdf",sep=""))
if(max(lagsfreqA)<=.8){ybound = c(0,max(lagsfreqA) + .05)} else{ybound=NULL}
barplot(lagsfreqA,
        axes = F,axisnames = TRUE,space=1,cex.main=1.4,xlab=expression(paste("Number of AR lags 'p'")),
        col="sienna3",density=100,border="black",ylim=ybound,ylab="Probabilities")
title("Posterior probability function")
#title(expression(paste("Posterior probability function of the number of AR lags")))
axis(2,las=1);abline(h=0)

graphics.off()


#(gg=gg+1)
#(pp=pp+1)

## Plot posterior distribution of break points together with time series:
## for all considered number of breaks 'm'
##--------------------------------------------------------------------------##

#(results = dget(paste(path,file,"_",rownames(oecd)[ii],".R",sep="")))			#get the saved results
#for(gg in 1:length(results$bpoints))
#{

#(path_latex =c("C:\\Users\\Alexander\\Documents\\Dissertation\\LaTex\\Paper1\\"))		#Pfad direkt zur Präsentation
#pdf(paste(path_latex,folder,"_",rownames(oecd)[ii],"_",gg,"_",pp,"new.pdf",sep="")) 

#pdf(paste(path,file,"_",rownames(oecd)[ii],"_",gg,"_",pp,"new.pdf",sep="")) 

#graphics.off()
par(ask=F,mfrow=c(2,1),font.main=11,cex.main=1,las=0)

plot(lseries,type="l",axes=F,xlab="",ylab="",lwd=2,col="grey9",ylim=c(min(lseries)*1.1,max(lseries)*1));
#plot(series,type="l",axes=F,xlab="",ylab="",lwd=2,col="grey9",ylim=c(min(series)*1.3,max(series)*1.1));
axis(4,las=0);axis(1); 

hh=2 ; dates = results$bpoints[[gg]][[pp]]		#[[gg]] auskommentieren falls nofb gefixed
while(hh<=(ncol(dates)-1))
{
  par(new=T)
  barplot(height = dates[,hh], las=1, axes=F, axisnames=F, xpd=F,  
          ,xlab="Time", ylab = "probability",beside = TRUE, ylim = c(0,1),col=heat.colors(10),density = 40, main="")
  hh=hh+1
}
axis(2,las=1);
#title(expression(paste("Posterior probability distributions of ",k[i])), cex.main = 1.2,font.main=12)
title(paste("Posterior probability distribution(s) of",gg,"break date(s)"))
#title(paste("Series of",rownames(oecd)[ii],"with the posterior(s) of",gg,"break date(s)"))
#graphics.off()

plot(results$grid,results$post_Norm,xlim=Lims,ylim=c(ups[1],max(results$post_Norm))*ups[2],type="l",las=1,lty=1,lwd=1,col="black",font.main=4,xlab="",ylab="density",sub=paste("(p = ",pp,")",sep=""))
title(paste("Marginal posterior densities '",rownames(oecd)[ii],"'",sep=""))
polygon(results$grid,results$post_Norm,col="lightgrey",border="black",lty=1)
lines(results$grid,results$post_Jeff,type="l",col="blue4",lty=2,ylab="",lwd=2,xlab="",main="Posterior with Jeffreys Prior",xlim=c(min(results$grid)*1.5,max(results$grid)))
legend("topleft",legend=c("Normal","Jeffreys"),col=c("black","blue4"),lwd=2,bty="n",lty=c(1,2), pch=c(NA,NA))

#}
graphics.off()



## Plot joint posterior probability function via stereogram of (m,p):
##--------------------------------------------------------------------##
(zz = postdraws_mp/sum(postdraws_mp))
dim(zz)
persp(x=0:5,y=1:5,z=zz,xlab="breaks 'm'", ylab="lags 'p'", zlab="probability",
      col="lightgrey",theta = 35, phi = 15,shade=NA,ticktype="detailed")
contour(x=0:5,y=1:5,z=zz,xlab="breaks 'm'", ylab="lags 'p'")
graphics.off()
#help(trans3d)










