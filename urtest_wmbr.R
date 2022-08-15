
###-----------------------------------------------------------------------------------###
### MCMC Sampler for an ADF-Test with variable lag order 'p'
### using the first p values of y_{t}; t = y_{-p+1},...,y_{0} as initial values:
###-----------------------------------------------------------------------------------### 

rm(list=ls())

## Anpassen vor Beginn der Session!!
##------------------------------------##
folder = "Example_nobreaks1"						#Target folder
label = NULL							#Label for simulation design
#(path = paste("Z:\\VosselerA\\Eigene Dateien\\Dissertation\\R\\Zivot_Wang_Model\\Outputs\\",folder,"\\",sep=""))
(path = paste("C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Zivot_Wang_Model\\Outputs\\",folder,"\\",sep=""))
dir.create(path)						#create folder in specified directory with given name
file = folder						#Anfang aller Filenamen, ob Plots oder txt-Files
##------------------------------------##


# Simulate data with structural breaks:
#---------------------------------------#
setwd("C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Functions\\")

(ar1break.sim = dget("ar1break.sim.R"))
(arma.break = dget("arma.break.R"))

n = 200                      #sample size; simulations

graphics.off()
#pdf(paste(path,file,"_data.pdf",sep="")) ; file.exists(paste(path,file,"_data.pdf",sep=""))

# Design 1 (Breaks in Level & Trend):
#--------------------------------------#
label = "Design 1: Stationary AR(1) with 3 breaks in level & trend"
k.true = c(33,151)
(y = data = ar1break.sim(n = n, k = k.true, coef = list(ar1 = 1, inter=c(.15,.75,.1), slope = c(0,0,0),errvar = c(.55,.55,.55)), init = 2, graph = T)$y)
title(label)
dev.off()

# Design 1b (Breaks in Level & Trend):
#--------------------------------------#
label = "Design 1: Stationary AR(1) with 5 breaks in level & trend"
k.true = c(33,55,82,121,151)
(y = data = ar1break.sim(n = n, k = k.true, coef = list(ar1 = .55, inter=c(1,1.1,-0.2,.5,1.6,.4), slope = c(.01,.025,.01,.015,.02,.02),errvar = c(.35,.35,.35,.35,.35,.35)), init = 2, graph = T)$y)
title(label)
dev.off()

#Design 2 (Breaks in level & variance):
#-----------------------------------------#
label = "Design 2: Stationary AR(1) with 2 breaks in level & variance"
k.true = c(51,101)
(y = data = ar1break.sim(n = n, coef = list(ar1 = .25, inter = c(0,.1,0),slope = c(0,0,0),errvar = c(.05,0.15,.05)), k = k.true, init = 0,graph = T, col = "sienna4")$y)
title(label)
dev.off()

#Design 3 (Random Walk with breaks in level & variance):
#---------------------------------------------------------#
label = "Random Walk with 2 breaks in level & variance"
k.true = c(33,151)
(y = ar1break.sim(n = n, k = k.true, coef = list(ar1 = 1, inter=c(.15,.75,.1), slope = c(0,0,0),errvar = c(.55,.65,.55)), init = 1, graph = T)$y)
title(label)
graphics.off()

#Design 4 (Stationary ARMAX with level breaks):
#------------------------------------------------------#
label = "Stationary ARMAX with level break"
k.true = c(103)
(y = arma.break(n = n, init=5, k = k.true, coef = list(ar = c(.7,.4,-.5,.4,-.35,-.11), ma=0, inter = c(0.85,.45),slope = c(0,0),errvar = c(.15,.15)), graph=T)$y)
title(label)
graphics.off()


## No-break design 1:
##--------------------##
n=200
set.seed(123)
eps = rnorm(n+1000,0,sd=0.25); 
y = data = arima.sim(list(order = c(2,0,0), ar = c(.8,-.35)), n = n, innov=eps[-c(1:1000)]);plot(y,type="l")  
#(y = arima.sim(list(order = c(2,0,1), ar = c(.4,.24), ma = c(.35)), n = n));plot(y,type="l")  

###########################################################################
###########################################################################
###########################################################################
###########################################################################

setwd("C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Zivot_Wang_Model")         

#rjbreak_ADF6 = dget("rjbreak_ADF6.R")

y=data;
burn=50;
m = 3
MHsim=500;
m_max=7;
p_max=10;
MMm=10^(-2);
MMp=10^(-1);
p_fix = 1;
m_fix = 1;
chm = 1;
chp = 1;
Jm_scale=5;
determc = "drift"
p_init=NULL; m_init=NULL; initial.breaks = NULL
Jp_scale = 5;
v0 = 2.001 ; lambda0 = 2.001; CC=100


rjbreak_ADF6 = function(data, burn = 300, MHsim = 1000, m_max = 5, p_max = 5, p_fix = NULL, m_fix = NULL, p_init=NULL, m_init=NULL, MMp = 10^(-2), MMm = 10^(-2),  chm = 1.1, chp = 1.1, Jm_scale = 5, Jp_scale = 5, v0 = 2.001 , lambda0 = 0.001, thinPAR=2, thinS2=1, thinM=1, thinP=1, determc = c("both","drift"), initial.breaks = NULL, Mlife=0.1, Plife=0.1, CC=100, ppL=10^(-2), ppU=10^(-10), tune=50, tria.m = c(A=0,C=1), tria.p = c(A=1,C=2))
{
  require(mvtnorm); require(coda); require(VGAM); require(MASS)  #VGAM for laplace distr.; MASS for General. Moore/Penrose inverse
  
  indBrkL=dget("indBrkL.R"); gibbsdraw_k=dget("gibbsdraw_k.R");wetbreak=dget("wetbreak.R"); 
  lagw=dget("lagw.R");inv=dget("inv.R"); dual2dec = dget("dual2dec.R") ; 
  mlogpost = dget("mlogpost2.R") 
  
  # Some useful functions:
  #--------------------------------------------------#
  Xupdate = function(y,p,k=NULL,type,y0=F, sbreaks = m>0)   
  {     
    lagw=dget("lagw.R");indBrkL=dget("indBrkL.R")
    n = length(y) ; 
    if(y0){ 
      y0=rep(0,p); (ylag1 = matrix(c(y0[1],y[1:(n-1)]),ncol=1))		#y0=TRUE means set zeros for initial values (T obs.) vs. use approx. Likelihood  T-p obs.
      y01 = y0[1] ; y0p = y0[-p]
    } else{y0=y01=y0p=NULL; ylag1 = matrix(y[p:(n-1)],ncol=1)}       
    if(p > 1){
      (d1.y = matrix(c(y01,diff(y)),ncol=1))				
      (d1.ylagk = lagw(d1.y,k=p-1,initial=y0p)[,-1])		
      (RiteX = cbind(ylag1,d1.ylagk))
      colnames(RiteX) = c("y1",paste("dy",1:(p-1),sep=""))
    } else if(p==1){(RiteX = cbind(ylag1)) ; colnames(RiteX)="y1" 
    } else{stop("AR order 'p' must be at least 1!\n\n")}
    #Break Case:
    if(sbreaks && !is.null(k)){
      (Dt = switch(type[1],drift=indBrkL(k)$ei1, both=indBrkL(k)$ei))		##Determ. Comp.; ei: Drift + Trend ##ei1: Drift
      X = cbind(Dt,RiteX)
    } else{X = switch(type[1],
                      drift=matrix(cbind(1,RiteX),ncol=ncol(RiteX)+1,dimnames=list(NULL,c("drift",colnames(RiteX)))),
                      both=matrix(cbind(1,(p+1):n,RiteX),ncol=ncol(RiteX)+2,dimnames=list(NULL,c("drift","trend",colnames(RiteX)))))
           cat("No structural breaks model!\n")
    } 
    return(list(ynew=y[-c(1:p)],X=X,RiteX=RiteX))
  }
  #dput(Xupdate,"Xupdate.R")
  #xy = Xupdate(y=y,p=5,k=c(1,10,length(y)+1),type="both") ;  head(xy$X) ;tail(xy$X)
  
  modelID = function(p,p_max,type,m_max,m)
  {
    if(p>p_max || m>m_max){stop("Wrong number of breaks or number of lags specified!\n")}
    if(type[1]=="drift"){
      dt = rep(0,m_max+1) ; dt[1:(m+1)] = rep(1,m+1)
    }
    if(type[1]=="both"){
      dt = rep(0,2*(m_max+1)) ; dt[1:(2*(m+1))] = rep(1,2*(m+1))
    }
    st = numeric(p_max) ; st[1:p]= rep(1,p)		#stochastic comp.
    return(c(dt,st))
  }
  #(mbits = modelID(p=p,p_max=p_max,type=type,m_max=m_max,m=m))
  #(mdec = dual2dec(mbits))					#One-to-one Mapping: Binary -> Decimal 
  #(mbits = dec2dual(dec=mdec))				#Decimal -> Binary system
  
  #Computes first two moments of theta_p in (19) and (20), p.5:
  moments12 = function(y,X,s2eps){
    Sigma2 = inv(as.double(s2eps)^(-1) * t(X)%*%X + as.double(s2eps)^(-1) * diag(ncol(X)))
    mu = (as.double(s2eps)^(-1)) * Sigma2%*%t(X)%*%y ; rownames(mu)=rownames(Sigma2)
    return(list(Sigma2=Sigma2,mu=mu))
  }
  
  ##Initializes a list or matrix object to save posterior draws of break dates 'k' depending on 'm' and 'p':
  ##---------------------------------------------------------------------------------------------------------##
  initial_pm = function(p_fix,m_fix,m_max,p_max)
  {
    save.bpoints=NULL
    #Case 1: (all fix)
    if(!is.null(p_fix) && !is.null(m_fix) && m_fix>0)
    {
      p_max=p_fix ; m_max=m_fix; save.bpoints = matrix(,ncol=m_max+2)
      #cat("\nCase1\n")
      nNew= n-p_fix ; k = matrix(c(1,(1:m_fix)*floor(nNew/(m_fix +1)),nNew+1),nrow=1)	
      save.bpoints[1,] = k
    }
    #Case 2: (only p fix)
    if(!is.null(p_fix) && is.null(m_fix) && m_max>0)
    { 
      save.bpoints = vector("list",m_max) 
      for(m in 1:m_max)
      {
        datearr = array(,dim=c(1,m+2,p_fix))
        nNew= n - p_fix ; k = matrix(c(1,(1:m)*floor(nNew/(m+1)),nNew+1),nrow=1)	 	#vector of initial break dates; Dimension: m+2!!
        #cat("\nBreak dates:",k,"\n\n") ; 
        save.bpoints[[m]] = k
      }  #end for
      #cat("\nCase2\n")
    }	#end if
    #Case 3: (only m fix)
    if(is.null(p_fix) && !is.null(m_fix) && m_fix>0)
    {  save.bpoints = vector("list",p_max)
       for(p in 1:p_max){
         nNew= n-p ; k = matrix(c(1,(1:m_fix)*floor(nNew/(m_fix + 1)),nNew+1),nrow=1)	 	#vector of initial break dates; Dimension: m+2!!
         save.bpoints[[p]] = k
       }  #end for
       #cat("\nCase3\n")
    }
    #Case 4: both flexible 
    if(is.null(p_fix) && is.null(m_fix) && m_max>0)
    {
      save.bpoints = vector("list",m_max) 
      #Initial values of break dates for each m=1,...,m_max:
      if(m_max>0)
      {
        for(m in 1:m_max){
          datearr = array(,dim=c(1,m+2,p_max))
          for(p in 1:p_max){
            nNew= n-p ; k = matrix(c(1,(1:m)*floor(nNew/(m+1)),nNew+1),nrow=1)	 	#vector of initial break dates; Dimension: m+2!!
            save.bpoints[[m]][[p]] = k
          } #end for
        } #end for
      } #end if
      #cat("\nCase4\n")
    } #end if
    return(save.bpoints)
  }
  #initial_pm(p_fix=5,m_fix=NULL, m_max=4, p_max=4) 
  #dput(initial_pm,"initial_pm.R")
  
  ## Function to correct the break dates finally (adding p-1) and to compute the break frequencies of each (m,p) combination:
  ##---------------------------------------------------------------------------------------------------------------------------##
  correct_datesfreq = function(p_fix,m_fix,m_max,p_max,save.bpoints)
  {		
    brfreq=NULL
    #Case 1: (all fix)
    if(!is.null(p_fix) && !is.null(m_fix) && m_fix>0){
      (save.bpoints[,-1] = save.bpoints[,-1]+(p_fix -1))
      brfreq = matrix(,ncol=ncol(save.bpoints),nrow=nrow(save.bpoints))
      brfreq = round(wetbreak(kMat = save.bpoints, n=n),6)
    } #end if 
    #Case 2: (only p fix)
    if(!is.null(p_fix) && is.null(m_fix) && m_max>0){ 
      brfreq = vector(class(save.bpoints),length=length(save.bpoints))
      for(m in 1:m_max){
        (save.bpoints[[m]][,-1] = save.bpoints[[m]][,-1]+(p_fix -1))
        brfreq[[m]] = round(wetbreak(kMat = save.bpoints[[m]], n=n),6)
      }
    } #end if
    #Case 3: (only m fix)
    if(is.null(p_fix) && !is.null(m_fix) && m_fix>0){
      brfreq = vector(class(save.bpoints),length=length(save.bpoints))
      for(p in 1:p_max){
        (save.bpoints[[p]][,-1] = save.bpoints[[p]][,-1]+(p-1))
        brfreq[[p]] = round(wetbreak(kMat = save.bpoints[[p]], n=n),6)
      }
    }
    #Case 4: both flexible 
    if(is.null(p_fix) && is.null(m_fix) && m_max>0){
      brfreq = vector(class(save.bpoints),length=length(save.bpoints))
      for(m in 1:m_max){
        for(p in 1:p_max){
          (save.bpoints[[m]][[p]][,-1] = save.bpoints[[m]][[p]][,-1]+(p-1))
          brfreq[[m]][[p]] = round(wetbreak(kMat = save.bpoints[[m]][[p]], n=n),6)
        }
      }
    }
    return(brfreq)
  }
  #correct_dates(p_fix=NULL,m_fix=2, m_max=3, p_max=4,save.bpoints) 
  
  ##Function to append ('save') new break dates or to get saved break dates as initial values:  
  ##-------------------------------------------------------------------------------------------##
  inout_bdates = function(p_fix,m_fix,m,p,save.bpoints,k_in=NULL)
  {
    #Put out some initial break dates from the 'save.bpoints' object:
    if(is.null(k_in))
    {
      if(!is.null(p_fix) && !is.null(m_fix) && m_fix>0)         #Case 1: (all fix)
      {(k = as.numeric(save.bpoints[nrow(save.bpoints),]))}
      if(is.null(p_fix) && !is.null(m_fix) && m_fix>0)		#Case 2: (only m fix)
      {(k = as.numeric(save.bpoints[[p]][nrow(save.bpoints[[p]]),]))} 
      if(!is.null(p_fix) && is.null(m_fix) && m>0)			#Case 3: (only p fix)
      {(k = as.numeric(save.bpoints[[m]][nrow(save.bpoints[[m]]),]))} 
      if(is.null(p_fix) && is.null(m_fix) && m>0)			#Case 4: both flexible 
      {(k = as.numeric(save.bpoints[[m]][[p]][nrow(save.bpoints[[m]][[p]]),]))
      }   
      return(k)
      #Or append some new drawn break dates 'kstar' on the 'save.bpoints' object
    } else{
      if(!is.null(p_fix) && !is.null(m_fix) && m_fix>0)		#Case 1: (all fix)
      {(save.bpoints = rbind(save.bpoints,k_in))} 
      if(!is.null(p_fix) && is.null(m_fix) && m>0) 		#Case 2: (only p fix)
      {(save.bpoints[[m]] = rbind(save.bpoints[[m]],k_in))}
      if(is.null(p_fix) && !is.null(m_fix) && m_fix>0) 		#Case 3: (only m fix)
      {(save.bpoints[[p]] = rbind(save.bpoints[[p]],k_in))}
      if(is.null(p_fix) && is.null(m_fix) && m>0)		#Case 4: both flexible 
      {(save.bpoints[[m]][[p]] = rbind(save.bpoints[[m]][[p]],k_in))}
      return(save.bpoints) 
    }
  }
  #inout_bdates(p_fix=NULL,m_fix=NULL,m=2,p=3,save.bpoints=save.bpoints)
  #dput(inout_bdates,"inout_bdates.R")
  
  ## Function to initialise the object to save the posterior draws of beta for each 'p \times m' combination:
  ##----------------------------------------------------------------------------------------------------------##
  initialise_savebetas = function(p_max, m_max, p_fix, m_fix)
  {
    if(is.null(m_fix) && is.null(p_fix))
    {
      (savebetas = vector(length=p_max*(m_max+1),"list")); (indi_ii = matrix(0,ncol=p_max*(m_max+1),nrow=1))
      hh=1; nam=NULL; 
      for(mm in 0:m_max){ for(pp in 1:p_max){
        (mbits = modelID(p=pp,p_max=p_max,type=type,m_max=m_max,m=mm))
        (nam=c(nam,paste(mbits,collapse=""))); savebetas[[hh]]=0; hh=hh+1
      }}   #loop
    }
    if(!is.null(m_fix) && is.null(p_fix))
    {
      (savebetas = vector(length=p_max,"list")); (indi_ii = matrix(0,ncol=p_max,nrow=1))
      hh=1; nam=NULL; 
      for(pp in 1:p_max){   
        (mbits = modelID(p=pp,p_max=p_max,type=type,m_max=m_max,m=m_fix))
        (nam=c(nam,paste(mbits,collapse="")))
        #(nam=c(nam,paste("m=",m_fix,"//p=",pp,sep="")))
        savebetas[[hh]]=0; hh=hh+1
      }
    }
    if(is.null(m_fix) && !is.null(p_fix))
    {
      (savebetas = vector(length=m_max+1,"list")); (indi_ii = matrix(0,ncol=m_max,nrow=1))  #m_max + 1 wegen 0,1,..
      hh=1; nam=NULL; 
      for(mm in 0:m_max){
        (mbits = modelID(p=p_fix,p_max=p_max,type=type,m_max=m_max,m=mm))
        (nam=c(nam,paste(mbits,collapse="")))
        savebetas[[hh]]=0; hh=hh+1
      }
    }
    if(!is.null(m_fix) && !is.null(p_fix))
    {
      (savebetas = vector(length=1,"list")) ;  indi_ii=0
      (mbits = modelID(p=p_fix,p_max=p_max,type=type,m_max=m_max,m=m_fix))
      (nam=paste(mbits,collapse=""))
      savebetas[[1]]=0
    }
    names(savebetas)=nam
    return(list(savebetas=savebetas,indi_ii=indi_ii))
  }
  #(sb = initialise_savebetas(p_max=5, m_max=5, p_fix=2, m_fix=2))
  
  # Triangular distribution as a (penalty/dilution) prior for the model indicators:
  #--------------------------------------------------------------------------------#
  triangular = function(x,A,B,C,g=F,logS=F,tune=50){
    ff=numeric(length(x));
    for(ii in 1:length(x)){
      if(x[ii] < A){ff[ii]=0} 
      if(A<=x[ii] && x[ii]<=C){ff[ii]=2*(x[ii]-A)/((B-A)*(C-A))}    #Note that C equals the mode!     
      if(C<x[ii] && x[ii]<=B){ff[ii]=2*(B-x[ii])/((B-A)*(B-C))} 
      if(B < x[ii]){ff[ii]=0}  
    }
    #ff = ff/sum(ff);
    if(logS){ ff[which(ff==0)] <- 10^(-tune) ; ff = log(ff)}
    if(g){plot(x,ff,type="b",xlab="",ylab="",main="Triangular prior",font.main=11);abline(v=c(A,C,B),lty=3,col="lightblue")}
    return(ff)
  }
  #------------------------------------------------------------------------------# 
  
  if(any(is.na(data))==TRUE){data = ts(na.omit(data))}  
  (n = length(data)) ; (type = match.arg(arg = determc, choices = c("drift","both"), several.ok = TRUE))
  
  # Initialisations:
  #------------------# 
  save.lags = nofb = mdec = save_theta = NULL
  if(is.null(m_init)){m = sample(1:m_max,1)} else{m = ifelse(m_max>0,yes=m_init, no=0)}	#start with 'm' breaks
  if(is.null(p_init)){p = sample(1:p_max,1)} else{p = ifelse(p_max>1,yes=p_init, no=1)}
  if(!is.null(p_fix)){p = p_max = p_fix} ; if(!is.null(m_fix)){m = m_max = m_fix}
  save.lags = c(save.lags,p) ; nofb = c(m,nofb)
  
  postdraws_mp = bics = matrix(0,ncol=p_max,nrow=m_max+1,dimnames=list(paste("m",0:m_max,sep=""),paste("p",1:p_max,sep="")));
  postdraws_mp[m+1,p]=1			#count absolute occurences of (m,p) dupels
  
  if(m>0){
    (save.bpoints = initial_pm(p_fix=p_fix, m_fix=m_fix, m_max=m_max, p_max=p_max)) 
    (k = inout_bdates(p_fix=p_fix,m_fix=m_fix,m=m,p=p,save.bpoints=save.bpoints))    #initial break dates vector 'k'
  }	
  (nNew = n-p)
  xy = Xupdate(y=data,p=p,k=k,type=type) ;  
  X = xy$X ; RiteX = xy$RiteX ; y=xy$y           	#update Designmatrix X according to new posterior draw of break date vector k
  
  npar = switch(type[1], drift = (2*m_max + p_max + 1), both = (3*m_max + p_max + 2)) 	     #number of unknown parameters estimated (break dates k_i included!) 
  if(npar>=n){stop("\nThere are more parameters to estimate than observations available!\n")}
  
  accrateA = accrateB = numeric(1); s2eps = numeric(MHsim+burn); 
  
  (sb = initialise_savebetas(p_max=p_max, m_max=m_max, p_fix=p_fix, m_fix=m_fix))
  indi_ii = sb$indi ; (savebetas = sb$save)
  
  #(savebetas = vector(length=p_max*(m_max+1),"list")) ; (indi_ii = matrix(1,ncol=p_max*(m_max+1),nrow=1))
  #hh=1;nam=NULL; for(mm in 0:m_max){ for(pp in 1:p_max){(nam=c(nam,paste("m=",mm,"//p=",pp,sep="")))
  #savebetas[[hh]]=0; hh=hh+1}} ;  names(savebetas)=nam  					#Initialise savebetas with zeros
  
  #------- AR parameters and resid.variance ------------------------------#
  
  (theta = solve(t(X)%*%X)%*%t(X)%*%y)						#Ols estimates of all coeff. (beta_hats)!!!!
  id = which(colnames(X)=="y1") ; pAR = theta[id] 
  err = y-X%*%theta ; d = length(theta)
  s2eps[1] = t(err)%*%err/(nNew-d)    				# Residual variance (no break) 
  # Note: setting v0=lambda0=0 yields the usual (uninformative) Jeffreys prior for \sigma^2 
  #Compute moments of theta (Step1):(Current values):
  mom.theta = moments12(y=y,X=X,s2eps=s2eps[1])				#compute first two Moments for model p
  (mu_theta = mom.theta$mu) ; var_theta = mom.theta$Sigma2
  
  #npar = m+p 		
  npar = switch(type[1], drift = (2*m + p + 1), both = (3*m + p + 2))   #all parameters including the m break dates
  #npar = switch(type[1], drift = (m + p + 1), both = (2*m + p + 2))	 #only the estimated number of regression parameters excluding the number of break dates	
  
  #logpost0p = logpost0m = mlogpost(yy=y, Xgam=X, a=v0, b=lambda0, cc=CC) -(npar/2)*log(nNew,base=exp(1))  		#T-p, statt T
  (logpost0p = logpost0m = mlogpost(yy=y, Xgam=X, a=v0, b=lambda0, cc=CC) -npar*log(n,base=exp(1)))
  #logpost0p = logpost0m = mlogpost(yy=y, Xgam=X, a=v0, b=lambda0, cc=CC) -(npar/2)*log(n,base=exp(1))			
  
  ii = jj = hh = 2; 				#probabilty of a life move for m,p
  repeat{
    
    ## Step 1: Model move p -> p*; with p=1,...,p_max (AR lag order):
    ##----------------------------------------------------------------##
    cat("\nDraw ",ii," (",length(save.lags),"/",length(nofb),"/",MHsim+burn,"):\n",sep="")  	#Hyperparameters
    
    #Propose new state for lag order 'p':  
    
    #Version1:
    #----------# 
    (pst = ceiling(rlaplace(1, location=p, scale=Jp_scale)))				#Proposal J(p,p*): draw candidate for p from J() double exp.
    #(pstar = ifelse(pst<1 || pst>p_max, yes=sample(1:p_max,1), no=pst))   	#left/right censoring p \in [1;p_max]   
    pstar=pst ; if(pst<1){(pstar=1)} ; if(pst>p_max){(pstar=p_max)}
    
    #Proposal at the margins, but only allowed if previous draw wasn't either:
    #if(pstar==1 && p==1){pstar=pstar+1}   
    #if(pstar==p_max && p==p_max){pstar=pstar-1}   
    
    #Version 2:
    #-----------#
    #if(p==1){pstar = p+1} 			#birth
    #if(p==p_max){pstar = p-1}			#death
    #if(p>1 && p<p_max){pstar = ifelse(runif(1)<=0.5, yes=p+1, no=p-1)}   #death/birth
    #if(runif(1)<Plife){pstar=p}			#life move
    
    #For penalty term:
    #-------------------#
    (nNew.star = n-pstar) ; 
    #npar.star = pstar   
    #npar.star = switch(type[1], drift = (2*m + pstar + 1), both = (3*m + pstar + 2))    
    npar.star = switch(type[1], drift = (m + pstar + 1), both = (2*m + pstar + 2))
    
    if(runif(1)<chp && is.null(p_fix))
    {
      if(m>0){
        (kstar = inout_bdates(p_fix=p_fix,m_fix=m_fix,m=m, p=pstar, save.bpoints=save.bpoints))   #fetch k_{m,pstar}  
      } else{kstar=NULL}
      
      xy.new = Xupdate(y=data,p=pstar,k=kstar,type=type)				#Update X matrix and y vector (Rows: n-pstar!)
      Xstar = xy.new$X ; RiteX.star = xy.new$RiteX.star ; ystar = xy.new$y    		# X_{p,m,k}  ; RiteX conatins only stochastic terms (vs. determ. trend and so on)
      
      ##Update moments under proposed values of 'p'(candidate!):
      mom.theta_star = moments12(y=ystar, X=Xstar, s2eps=s2eps[ii-1])
      mu_theta_star = mom.theta_star$mu ; var_theta_star = mom.theta_star$Sigma2
      
      #Compute acceptance probability:
      #----------------------------------#
      #Numerator and denominator in Logs: 
      
      JA10 = dlaplace(p,location=pstar,scale=Jp_scale)  		#Version 1
      JA01 = dlaplace(pstar,location=p,scale=Jp_scale)
      # JA10=JA01=0.5								#Version 2
      
      if(p!=pstar){
        #logpost1p = mlogpost(yy=ystar, Xgam=Xstar, a=v0, b=lambda0, cc=CC) -log(nNew.star,base=exp(1))*(npar.star/2)
        #logpost1p = mlogpost(yy=ystar, Xgam=Xstar, a=v0, b=lambda0, cc=CC) -log(n,base=exp(1))*(npar.star/2)	
        
        if(is.null(ppL) && is.null(ppU)){(prob = rep(1/p_max,times=p_max))}
        if(!is.null(ppL) && is.null(ppU)){(prob = c(ppL,rep((1-(ppL))/(p_max-1),p_max-1)))}
        if(is.null(ppL) && !is.null(ppU)){(prob = c(rep((1-(ppU))/(p_max-1),p_max-1),ppU))}
        if(!is.null(ppL) && !is.null(ppU)){(prob = c(ppL,rep((1-(ppL+ppU))/(p_max),p_max),ppU))}
        
        logpost1p = mlogpost(yy=ystar, Xgam=Xstar, a=v0, b=lambda0, cc=CC) #+ log(prob)[pstar]
        
        #logpost1p = mlogpost(yy=ystar, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + triangular(x=pstar,A=tria.p[1],B=p_max, C=tria.p[2], logS=T,tune=tune)      
  
        #logpost1p = mlogpost(yy=ystar, Xgam=Xstar, a=v0, b=lambda0, cc=CC) -log(n,base=exp(1))*npar.star
      } else{logpost1p=logpost0p}
      
      cat("proposal 'p*'(lpost):",logpost1p,"\n") ; cat("current 'p'(lpost):",logpost0p,"\n")
      (log_u = log(runif(1)*MMp)) ; accA=0
      
      #Some MH monitoring:
      cat("proposal 'p*':",pstar,"\n") ; cat("current 'p':",p,"\n")  
      
      ## Take equation (23),acceptance probability in logs:
      ##-----------------------------------------------------##
      if(logpost1p-logpost0p + log(JA10)-log(JA01) > log_u)
      {
        X = Xstar ; p = pstar; var_theta = var_theta_star ; mu_theta = mu_theta_star	  #accept proposal!
        RiteX = RiteX.star ; y = ystar; k = kstar ; nNew = nNew.star
        logpost0p = logpost1p ; save.lags = c(save.lags,p) ; 
        accA=1 ; cat("'p*' accepted!\n") ; postdraws_mp[m+1,pstar] = postdraws_mp[m+1,pstar] + 1		#given the current value of m
        #bics[m+1,pstar] = bics[m+1,pstar] + BIC(lm(y~X-1))     #sum and finally average (see below)
      } else{if(runif(1)<.05){logpost0p = -10^10}}
      
      accrateA[jj] = (accrateA[jj-1]*(jj-1) + accA)/jj
      if(jj>1){cat("accept.rate 'p':", accrateA[jj],"\n\n")}
      jj=jj+1
    }  #end if  
    
    ## Step 2: Propose new number of structural breaks 'mstar':
    #--------------------------------------------------------------# 
    #Version 1:
    mst = ceiling(rlaplace(1, location=m, scale=Jm_scale))	   		#Proposal J(p,p*): draw candidate for p from J() double exp.
    mstar=mst ; if(mst<0){(mstar=0)} ; if(mst>m_max){(mstar=m_max)}
    
    #Proposal at the margins, but only allowed if previous draw wasn't either:
    #if(mstar==1 && m==1){mstar=mstar+1}   
    #if(mstar==m_max && m==m_max){mstar=mstar-1} 
    
    #Version 2:
    #   if(m==0){mstar = m+1} 
    #   if(m == m_max){mstar = m-1}
    #   if(m>0 && m<m_max){mstar = ifelse(runif(1)<=0.5, yes=m+1, no=m-1)}
    #   if(runif(1)<Mlife){mstar=m}			#life move
    
    #npar.star = mstar
    #npar.star = switch(type[1], drift = (2*mstar + p + 1), both = (3*mstar + p + 2))  
    npar.star = switch(type[1], drift = (mstar + p + 1), both = (2*mstar + p + 2))    
    
    if(runif(1)<chm && m_max>0 && is.null(m_fix))
    {
      if(mstar>0){
        (start.breaks = inout_bdates(p_fix=p_fix, m_fix=m_fix, m=mstar, p=p,save.bpoints=save.bpoints))
        kstar = gibbsdraw_k(y=y, RiteX=RiteX, m=mstar,type=type,s2eps=s2eps[ii-1], start.breaks=start.breaks)  			
        (save.bpoints = inout_bdates(p_fix=p_fix, m_fix=m_fix, m=mstar, p=p,save.bpoints=save.bpoints, k_in = kstar))
      } else{kstar=NULL}   									#update X matrix according to proposal!
      
      xy.new = Xupdate(y=data,p=p,k=kstar,type=type)	
      Xstar = xy.new$X ; RiteX.star = xy.new$RiteX  							# update X_{p,m,k}
      #theta.star = try(ginv(t(Xstar)%*%Xstar)%*%t(Xstar)%*%y,silent=T)
      
      ##Update moments of theta under proposed values of 'm'(Candidate!):
      #--------------------------------------------------------------------#
      mom.theta_star = moments12(y=y,X=Xstar,s2eps=s2eps[ii-1])					#for gibbs step -> draw from full cond. of theta
      (mu_theta_star = mom.theta_star$mu) ; (var_theta_star = mom.theta_star$Sigma2)
      
      #Compute acceptance probability of candidate 'mstar' (new number of breaks):
      #------------------------------------------------------------------------------#
      JB10 = dlaplace(m,location=mstar,scale=Jm_scale)		#transition probabilitiy
      JB01 = dlaplace(mstar,location=m,scale=Jm_scale)
      # JB10=JB01=0.5							#Version 2
      
      if(m!=mstar){
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) -log(nNew.star,base=exp(1))*(npar.star/2)
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) -log(n,base=exp(1))*(npar.star/2)              #with penalty term (from prior)
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + 1/(1+ exp(tune*(npar.star)))
        
        #(const = pnorm(m_max-1,mean=trunc(.5*(m_max)),sd=tune,log=F)-pnorm(1,mean=trunc(.5*(m_max)),sd=tune,log=F))
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + log(dnorm(mstar,mean=trunc(.5*(m_max)),sd=tune,log=F)/const)
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + dnorm(mstar,mean=trunc(.5*(m_max+1)),sd=tune,log=T)
        
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + triangular(x=mstar,A=tria.m[1],B=m_max,C=tria.m[2], logS=T,tune=tune)
        
        if(is.null(ppL) && is.null(ppU)){(prob = rep(1/(m_max+1),times=m_max + 1))}
        if(!is.null(ppL) && is.null(ppU)){(prob = c(ppL,rep((1-(ppL))/(m_max),m_max)))}
        if(is.null(ppL) && !is.null(ppU)){(prob = c(rep((1-(ppU))/(m_max),m_max),ppU))}
        if(!is.null(ppL) && !is.null(ppU)){(prob = c(ppL,rep((1-(ppL+ppU))/(m_max-1),m_max-1),ppU))}
        
        logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + log(prob)[mstar+1]
        
        #logpost1m = mlogpost(yy=y, Xgam=Xstar, a=v0, b=lambda0, cc=CC) + dpois(x=mstar,lambda=trunc(.5*(m_max))-1,log=T) 
        
      } else{logpost1m=logpost0m}
      
      cat("proposal 'm*'(lpost):",logpost1m,"\n") ; cat("current 'm'(lpost):",logpost0m,"\n")
      (log_v = log(runif(1)*MMm)) ; accB=0		#reject
      
      #Some MH monitoring:
      cat("proposal 'm*':",mstar,"\n") ; cat("current 'm':",m,"\n") ; 
      
      ## Take equation (23) with acceptance probability in logs:
      ##----------------------------------------------------------##
      if(logpost1m-logpost0m+log(JB10)-log(JB01) > log_v)		#accept!
      {
        m = mstar; X = Xstar ; k = kstar; var_theta = var_theta_star ; mu_theta = mu_theta_star		#accept proposal!
        RiteX = RiteX.star; logpost0m = logpost1m; 
        nofb = c(m,nofb) ; accB=1 ; postdraws_mp[mstar+1,p] = postdraws_mp[mstar+1,p] + 1
        #bics[mstar+1,p] = bics[mstar+1,p] + BIC(lm(y~X-1))     #sum and finally average (see below)
        
        cat("'m*' accepted!\n")
      } else{ if(runif(1)<.05){logpost0m = -10^10}}		#ensure acceptance of candidate from time to time (for better mixing)
      
      if(ii>burn){
        mbits = modelID(p=p, p_max=p_max, type=type, m_max=m_max, m=m)   #take accepted model combination (p,m)
        mdec = c(mdec, dual2dec(mbits))}
      
      accrateB[hh] = (accrateB[hh-1]*(hh-1) + accB)/hh
      if(hh>1){cat("accept.rate 'm':",accrateB[hh],"\n")}
      
      hh=hh+1
    } #end if
    
    ## Step 3: Draw new break dates 'k' from full conditional distr. (given 'm' and 'p') 
    ##-------------------------------------------------------------------------------------##
    if(m>0){
      (start.breaks = inout_bdates(p_fix=p_fix, m_fix=m_fix, m=m, p=p,save.bpoints=save.bpoints))
      k = gibbsdraw_k(y=y, RiteX=RiteX, m=m,type=type,s2eps=s2eps[ii-1], start.breaks=start.breaks) 
      (save.bpoints = inout_bdates(p_fix=p_fix, m_fix=m_fix, m=m, p=p,save.bpoints=save.bpoints, k_in = k))
      xy.new = Xupdate(y=data,p=p,k=k,type=type) ; 		
      X = xy.new$X ; RiteX = xy.new$RiteX   
      
      mom.theta = moments12(y=y,X=X,s2eps=s2eps[ii-1])			##Update moments under new value of 'm'(Candidate!)
      (mu_theta = mom.theta$mu) ; (var_theta = mom.theta$Sigma2)
    }  	#end if
    
    # Step 4: Draw AR parameters :
    ##-------------------------------## 
    (theta = t(rmvnorm(1,mean = mu_theta, sigma = var_theta)))    #Beta parameters 	
    rownames(theta)=colnames(X)
    id = which(colnames(X)=="y1") 
    err = y-X%*%theta 			#posterior draws of theta; update prediction errors!!
    
    if(!is.null(p_fix) && !is.null(m_fix))
    {(save_theta = mcmc(cbind(save_theta, theta)))}    #in work
    
    # Compute sample means of betas after some burnin:
    if(ii>burn){
      if(is.null(m_fix) && is.null(p_fix)){combi = m*p_max + p}
      if(!is.null(m_fix) && is.null(p_fix)){combi = p}
      if(is.null(m_fix) && !is.null(p_fix)){combi = m + 1}
      if(!is.null(m_fix) && !is.null(p_fix)){combi = 1}
      zz = indi_ii[combi] + 1							
      savebetas[[combi]] = (savebetas[[combi]]*(zz-1) + theta)/zz    #compute the posterior means -> Bayes estimator under quadratic loss function
      indi_ii[combi] = indi_ii[combi] + 1
    }
    
    # For plotting the marginal posterior of the long run comp.:
    if(!is.null(p_fix) && !is.null(m_fix)){
      pAR = cbind(pAR,theta) ; rownames(pAR)=colnames(X)		#..or saving the posterior draws of the whole theta vector
    } else{pAR = cbind(pAR,theta[id])}
    
    # Step 5: Draw from cond. Posterior of s2e:
    #--------------------------------------------#  
    v.e = v0 + nNew/2
    lambda.e = lambda0 + t(err)%*%err/2                          #with B0==0
    (s2eps[ii] = 1/rgamma(1,shape = v.e, rate = lambda.e))	#Draw from Inverse Gamma Distr. (s2e)
    
    #Stopping Rules:
    if(is.null(m_fix) || is.null(p_fix)){
      if(length(nofb)<(MHsim+burn) && length(save.lags)<(MHsim+burn)){ii=ii+1} else{break}
    } else{
      if(ii<(MHsim+burn)){ii=ii+1} else{break}
    }
    cat("\nCurrent (p/m): (",p,",",m,")\n",sep="")
  }  #MC end loop
  #-----------------------------------------------------------------------------------#
  
  brfreq = correct_datesfreq(p_fix=p_fix,m_fix=m_fix, m_max=m_max, p_max=p_max,save.bpoints=save.bpoints) 
  
  if(m_max>0){ if(is.null(m_fix)){ nofbs = window(mcmc(nofb[-c(1:burn)]),thin=thinM)} else{nofbs=NULL}
  } else{nofbs=brfreq=NULL}
  
  if(is.null(p_fix)){lags=window(mcmc(save.lags[-c(1:burn)]),thin=thinP)} else{lags=NULL}
  
  bics = (1/postdraws_mp) * bics		#finally compute arithmeic means for each matrix entry, i.e. elementwise
  
  return(invisible(list(bp=brfreq,lags=lags, 
                        postdraws_mp=postdraws_mp,bics=bics, thetas=save_theta,
                        nofb=nofbs, accrateA=accrateA, mdec = mdec,
                        accrateB=accrateB, s2eps=window(mcmc(s2eps[-c(1:burn)]),thin=thinS2), 
                        indi_ii=indi_ii, savebetas=savebetas, pAR=window(mcmc(t(pAR[,-c(1:burn)])),thin=thinPAR)
  )))
}
##############################################################################################################################
##############################################################################################################################
dput(rjbreak_ADF6,"rjbreak_ADF6.R")

## Call:
##-------##
out = rjbreak_ADF6(data=y, determc="drift", burn=1000, MHsim=10000, MMm=10^(0.1), MMp=10^(0.3), Jm_scale = 2, Jp_scale = 5, m_max=5, p_max=10, p_fix=NULL, m_fix=NULL, chm=1.1, chp=1.1, v0 = 2.001 , lambda0 = 0.001, CC=100, ppL=NULL, ppU=10^(-10), tune=20, tria.m = c(A=-1,C=0), tria.p = c(A=-1,C=0))


(mdec = out$mdec) 
(lags = out$lags) 
(sigma2 = out$s2eps)
(bpoints = out$bp)
(nofb = out$nofb)
(accrateA = out$accrateA)
(accrateB = out$accrateB)
(theta = out$pAR); dim(theta)
betatraces = out$theta; dim(betatraces)
(betas = out$savebetas); 
(indi = out$indi)
(postdraws_mp = out$postdraws_mp)
(bics = out$bics)

(modelfreq = round(table(nofb)/length(nofb),4))
(lagsfreq = round(table(lags)/length(lags),4))

#mom.beta = apply(mcmc(betatraces),1,summary); colnames(mom.beta)=rownames(betatraces);mom.beta


## Save summary statistics of beta coefficients:
##----------------------------------------------------##
#sink(paste(path,"summariesBETAS_Design1.out",sep=""))
cat("\n Variables: ",rownames(betatraces),'sigma2',"\n",rep("=",20),cat("\n"));
for(ii in 1:nrow(betatraces)){ 
  print(summary(mcmc(betatraces[ii,]))) ; print(HPDinterval(mcmc(betatraces[ii,])))
}
print(summary(mcmc(sigma2))) ; print(HPDinterval(mcmc(sigma2)))

#sink()

## Traceplots of hyperparameters 'lags' and 'nofb':
##--------------------------------------------------------##
#pdf(paste(path,file,"_traceplots.pdf",sep=""))

par(ask=F,mfrow=c(2,1),font.main=11,cex.main=1.4)

IND = 100*c(1:floor(min(c(length(lags),length(nofb)))/100));
(ii=length(IND)-1)
#for(ii in 1:(length(IND)-1)){
plot(IND[ii]:IND[ii+1],mcmc(lags[IND[ii]:IND[ii+1]]),type="l",xlab="",ylab="",main="AR lag order")
plot(IND[ii]:IND[ii+1],mcmc(nofb[IND[ii]:IND[ii+1]]),type="l",xlab="",ylab="",main="Number of structural breaks")
#}
#graphics.off()


## Plot posterior probabilities of number of breaks 'm'
## and number of AR lags 'p'
##------------------------------------------------------------##
#win.graph(width=6.5,height=3.5); 
par(ask=F,font.main=11,mfrow=c(1,2));

#pdf(paste(path,file,"_hist_nofb.pdf",sep=""))
if(max(modelfreq)<=.8){ybound = c(0,max(modelfreq) + .2)} else{ybound=NULL}
barplot(modelfreq,main=expression(paste("Posterior probability function of the number of breaks")),
        axes = F,axisnames = TRUE,space=1,cex.main=1.4,xlab=expression(paste("Number of breaks ",m)),
        col="sienna3",density=100,border="black",ylim=ybound,ylab="Probabilities")
axis(2,las=1) #;axis(1,labels = F)
abline(h=0)

#pdf(paste(path,file,"_hist_lags.pdf",sep=""))
if(max(lagsfreq)<=.8){ybound = c(0,max(lagsfreq) + .2)} else{ybound=NULL}
barplot(lagsfreq,main=expression(paste("Posterior probability function of the number of AR lags")),
        axes = F,axisnames = TRUE,space=1,cex.main=1.4,xlab=expression(paste("Number of AR lags ",p)),
        col="sienna3",density=100,border="black",ylim=ybound,ylab="Probabilities")
axis(2,las=1);abline(h=0)

graphics.off()


dec2dual = dget("C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Zivot_Wang_Model\\dec2dual.R")
dual2dec = dget("C:\\Users\\Alexander\\Documents\\Dissertation\\R\\Zivot_Wang_Model\\dual2dec.R") 

## Compute model frequencies and best model specifications :
##-------------------------------------------------------------##
(hh = sort(round(table(mdec)/length(mdec),5)))
best = length(hh)
Mdec = as.double(names(hh)) ; Mfreq = as.double(hh) ; mbits = vector("list",length(Mdec))
for(pp in 1:length(Mdec)){mbits[[pp]] = dec2dual(dec=Mdec[pp])} ; #mbits
ww=1; for(ss in length(mbits):(length(mbits)-best+1)){
  cat("Model ",ww," (",dual2dec(mbits[[ss]]),")\n\n",sep=""); cat(mbits[[ss]],"\n\n");ww=ww+1}
##-------------------------------------------------------------##

hh

for(oo in which(indi!=0)){print(names(betas)[oo]);print(betas[[oo]])}
dual2dec(100)
as.numeric(names(betas)[oo])


## Plot MH - acceptance rates:
##-----------------------------##
#pdf(paste(path,file,"_MHaccrates.pdf",sep=""))
par(ask=F,mfrow=c(2,1),las=1,font.main=12)
plot(accrateA,xlab="Iterations",ylab=expression(paste("in %")),main="Acceptance rates 'p'",type="l",xlim=c(0,1000))
plot(accrateB,xlab="Iterations",ylab=expression(paste("in %")),main="Acceptance rates 'm'",type="l",xlim=c(0,1000))

graphics.off()


## s2eps: Residual variance
##-----------------------------##
#win.graph(width=7,height=4); 
#pdf(paste(path,file,"_post_s2eps.pdf",sep=""))
par(ask=F,font.main=13,cex.main=1.4,las=1);
plot(sigma2,main="")
hist(sigma2,prob=T,nclass=120,col="wheat2",xlab=expression(sigma[epsilon]^2),xlim=c(.5*min(sigma2),1.1*max(sigma2)),main="")
lines(density(sigma2),lwd=2)
title(expression(paste("Posterior density of ",sigma[epsilon]^2)))
#acf(sigma2)
#cumuplot(sigma2,main="Chain of sigma2",probs=c(0.025,0.5,0.975))

graphics.off()

effectiveSize(sigma2)
(effs = effectiveSize(sigma2)/length(sigma2))  


## Theta: Long term coefficient 
##-----------------------------------##
#win.graph(width=7,height=4); 
#pdf(paste(path,file,"_post_AR.pdf",sep=""))

par(ask=F,font.main=13,cex.main=1.4,las=1);
hist(theta,prob=T,nclass=120,col="wheat2",xlab=expression(theta),xlim=c(.5*min(theta),1.2*max(theta)),main="")
lines(density(theta),lwd=2)
title(expression(paste("Posterior density of ",theta)))
#acf(theta)
#cumuplot(theta,main="Chain of theta",probs=c(0.025,0.5,0.975))
summary(as.vector(theta))

dev.off()

effectiveSize(theta)						#Effizienz: effectiveSize/GibbsRuns * 100  (RC S.256)
(effb = effectiveSize(theta)/length(theta)) 			#falls 1, MC=iid


## All beta coefficients
##-----------------------------------##
#pdf(paste(path,file,"_post_betas.pdf",sep=""))

par(ask=T,font.main=13,cex.main=1.4,las=1);
win.graph(width=7,height=4); par(ask=T,font.main=13);

for(ii in 1:nrow(betatraces))
{
  plot(mcmc(betatraces[ii,]),main=paste("Chain of '",rownames(betatraces)[ii],"'",sep=""))
  #cumuplot(betatraces[,ii],main=paste("Chain of '",colnames(betatraces)[ii],"'",sep=""))     #Berechnet "Running CDF" für jede Kette; Oben 97.5%-Quantil unten 2.5%-Quantil (Robert/Casella, S.243)
  acf(betatraces[ii,],main=paste("Chain of '",rownames(betatraces)[ii],"'",sep="")) 
  #print(autocorr.diag(betatraces[,ii]))									
}
graphics.off()


gg=1
(gg=gg+1)
(gg = as.numeric(names(which.max(modelfreq))))   #MAP number of breaks
dd=1
(dd=dd+1)
(dd = as.numeric(names(which.max(lagsfreq))))	#MAP number of lags


## Plot posterior distribution of break points together with time series:
## for all considered number of breaks 'm'
##--------------------------------------------------------------------------##

par(ask=T)
for(gg in 1:length(bpoints))			#overall nr. of breaks given MAP of lags ('dd')
{
  #pdf(paste(path,file,"_pdistr_breakp",gg,".pdf",sep="")) 
  #graphics.off()
  plot(y,type="l",axes=F,xlab="",ylab="",lwd=3,col="grey9",ylim=c(min(0,min(y)),max(y)));
  axis(4,las=1);axis(1); abline(v = k.true,lty=2)
  ii=2
  while(ii<=(ncol(bpoints[[gg]][[dd]])-1))
  {
    par(new=TRUE)
    barplot(height = bpoints[[gg]][[dd]][,ii], las=1, axes=F, axisnames=F, xpd=F,  
            ,xlab="Time", ylab = "Probability",beside = TRUE, ylim = c(0,1),col=heat.colors(10),density = 50, main="")
    ii=ii+1
  }
  axis(2,las=1);
  title(expression(paste("Posterior probability distributions of ",k[i])), cex.main = 1.2,font.main=12)
  print(gg)
  #graphics.off()
}


?spectrum0
spectrum0(out$betas[,1])$spec/length(out$betas[,1])
spectrum(out$betas[,1])$freq
spectrum(out$betas[,1])$spec

sqrt(sdf(betas))
sqrt(spectrum0(out$betas)$spec)
summary(out$betas)

#sink("Outputs\\gibbs_log2.txt") 
#sink()
#dput(mbreak_ADF,"mbreak_ADF.R")



# (log)POSTERIOR PROBABILITY OF MODEL 'n_i' UP TO A NORMALIZING CONSTANT:
#-------------------------------------------------------------------------#
mlogpost = function(yy, Xgam, betatil.gam=NULL,cc=1,base=exp(1),a=NULL,b=NULL)
{  
  require(MASS)		#Calculates the Moore-Penrose generalized inverse of a matrix X
  if(is.null(betatil.gam)){betatil.gam = matrix(rep(0,ncol(Xgam)),ncol=1)}
  (n = nrow(Xgam)) ; Xgam = as.matrix(Xgam) ; betatil.gam = as.matrix(betatil.gam)
  #Log Posterior of gamma (conjugate prior):
  Mgam = diag(cc,nrow=ncol(Xgam),ncol=ncol(Xgam))
  (S = diag(n) + Xgam%*%ginv(Mgam)%*%t(Xgam)) ;		#mit betatilde.gam Hyperparameter Prior!!!!
  (logpost_ni = -.5*log(det(S),base=base) -(n/2)*log(t(yy-Xgam%*%betatil.gam)%*%ginv(S)%*%(yy-Xgam%*%betatil.gam),base=base))   
  
  #Log Posterior of gamma (Zellner's G-prior):
  #S = diag(n) + n*Xgam%*%ginv(t(Xgam)%*%Xgam)%*%t(Xgam)		#under G Prior
  #äquivalente Darstellung im G-Prior case via quadr.Form, vgl. auch conjugate case!! 
  #(logpost_ni = -.5*log(det(S)) -(n/2)*log(t(yy - Xgam%*%betatil.gam)%*%ginv(S)%*%(yy - Xgam%*%betatil.gam)))
  #(logpost_ni = -((ncol(Xgam)+1)/2)*log(n + 1) -(n/2)*log(t(yy-Xgam%*%betatil.gam)%*%ginv(S)%*%(yy-Xgam%*%betatil.gam)))
  return(as.numeric(logpost_ni))
}
dput(mlogpost,"mlogpost.R")
#mlogpost(yy=y, Xgam=Xgam, betatil.gam=betatil.gam)



# Alternativ: (log)POSTERIOR wie in Paper beschrieben:
# Unterschied zu oben: 
# Oben wird \sigma^2 ~ \sigma^{-2} 
# (Jeffreys prior) verwendet statt Inverse Gammavert.
# Beachte Hyperparameter a,b>0 aus IG(a,b)
#----------------------------------------------------------#

mlogpost2 = function(yy, Xgam, betatil.gam=NULL,cc=100,base=exp(1),a,b,kernel.only=F)
{  
  require(MASS);	
  if(is.null(betatil.gam)){betatil.gam = matrix(rep(0,ncol(Xgam)),ncol=1)}
  (n = nrow(Xgam)) ; Xgam = as.matrix(Xgam) ; betatil.gam = as.matrix(betatil.gam)
  #Log Posterior of gamma (conjugate prior):
  
  Mgam = diag(1/cc,nrow=ncol(Xgam),ncol=ncol(Xgam))
  #(S = ginv(diag(n) + Xgam%*%ginv(Mgam)%*%t(Xgam))) ;	
  (S = diag(n) - Xgam%*%ginv(Mgam + t(Xgam)%*%Xgam)%*%t(Xgam)) ;
  (QF = t(yy-Xgam%*%betatil.gam)%*%S%*%(yy-Xgam%*%betatil.gam))	
  if(kernel.only){ l_nconst = 0;
  } else{ (Det = det(S)^{1/2});
          (l_nconst = log(Det))
          #(l_nconst = as.numeric( -(n/2)*log(pi) + log(Det) + log(gamma((n+a)/2)) - log(gamma(a/2)) + (a/2)*log(b)))   #normalizing constant
  }
  (logpost_ni = l_nconst -.5*(n+a)*log(b + QF))   
  
  #Log Posterior of gamma (Zellner's G-prior):
  #S = diag(n) + n*Xgam%*%ginv(t(Xgam)%*%Xgam)%*%t(Xgam)		#under G Prior
  #äquivalente Darstellung im G-Prior case via quadr.Form, vgl. auch conjugate case!! 
  #(logpost_ni = -.5*log(det(S)) -(n/2)*log(t(yy - Xgam%*%betatil.gam)%*%ginv(S)%*%(yy - Xgam%*%betatil.gam)))
  #(logpost_ni = -((ncol(Xgam)+1)/2)*log(n + 1) -(n/2)*log(t(yy-Xgam%*%betatil.gam)%*%ginv(S)%*%(yy-Xgam%*%betatil.gam)))
  
  return(as.numeric(logpost_ni))
}
dput(mlogpost2,"mlogpost2.R")
#mlogpost2(yy=yy, Xgam=X, betatil.gam=NULL,a=4,b=4)


#--------Only Zellner's G prior------------------------------------------------------------------#
#a=0;b=1; cc=100
#betatil.gam=NULL;cc=1; yy=y; Xgam=X

mlogpost3 = function(yy, Xgam, betatil.gam=NULL,cc=1,a=0,b=0)     #note a,b -> 0, yields Jeffreys ('diffuse') prior
{  
  require(MASS);		
  (n = nrow(Xgam)) ; d = ncol(Xgam);
  if(is.null(betatil.gam)){betatil.gam = matrix(rep(0,d),ncol=1)}
  Xgam = as.matrix(Xgam) ; betatil.gam = as.matrix(betatil.gam)  
  
  #Log Posterior of gamma under Zellner's G-prior:
  #-------------------------------------------------#
  scale = t(Xgam)%*%Xgam; Det = (cc+1)^(-d/2); gam=1;
  Pw = Xgam%*%ginv(scale)%*%t(Xgam);
  (QF = t(yy)%*%yy -(cc/(cc+1))*t(yy)%*%Pw%*%yy -(1/(cc+1))*t(betatil.gam)%*%scale%*%betatil.gam) 
  #S = diag(length(yy))-(cc/(cc+1))*Pw                
  #(QF = t(yy-Xgam%*%betatil.gam)%*%S%*%(yy-Xgam%*%betatil.gam))    
  (l_nconst = as.numeric( -(length(yy)/2)*log(pi) + log(Det) + log(gamma((length(yy)+a)/2)) - log(gam) + (a/2)*log(b)))
  l_nconst = log(Det)
  (logpost_ni = l_nconst -.5*(length(yy)+a)*log(b + QF))
  return(as.numeric(logpost_ni))
}
dput(mlogpost3,"mlogpost3.R")

mlogpost3(yy=y, Xgam=Xgam, betatil.gam=betatil.gam,a=0,b=1,cc=100)














