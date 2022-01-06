


#######################################################################################################################################
########################################################### Load  #####################################################################
#######################################################################################################################################


library(bayesm)
library(evd)
library(copula)
library(stats)
library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(mvtnorm)
options(scipen=500)
options(stringsAsFactors = FALSE)


root_dir = ""

#######################################################################################################################################
########################################################### Settings ##################################################################
#######################################################################################################################################

############################################################## Model

#### 1: Independent
#### 2: Power
#### 3: Sine
#### 4: Guassian

draw_model = 4
est_model = 4
est_M = FALSE
M_min =  20
M_max =  30
z_var = "z"
marg = "L"

Settings = list(draw_model = draw_model,est_model = est_model,z_var = z_var,est_M = est_M,M_min = M_min,M_max = M_max,marg=marg)

setwd(root_dir)
source("directcopula_functions_r.r")
sourceCpp("directcopula_functions_cpp.cpp")

#######################################################################################################################################
########################################################### Simuluation ###############################################################
#######################################################################################################################################


################################################### for simulation
nlgt = 100
nSKU = 2
nobs = 100

nbeta = (nSKU) * 2

ntau = nSKU * (nSKU - 1)/2

psi_out = 1
ntheta = 1
nzbeta=1
nztau=1
nztheta = 1
ZB=matrix(rep(1,nlgt*nzbeta),nrow=nlgt,ncol=nzbeta)
ZT=matrix(rep(1,nlgt*nztau),nrow=nlgt,ncol=nztau)
ZQ=matrix(rep(1,nlgt*nztheta),nrow=nlgt,ncol=nztheta)


Delta=matrix(c(-1,-0.5,-0.5,0.5),nrow=nzbeta,ncol=nbeta)
Gamma=matrix(c(log(20*nSKU/4)),nrow=nztheta,ncol=ntheta)

if(Settings$draw_model==4){

  Omega=matrix(c(-0.6),nrow=nztau,ncol=ntau)

}else{

  Omega=matrix(c(1.6),nrow=nztau,ncol=ntau)

}



iota0=matrix(1,nrow=nbeta,ncol=1)
Vbeta=(diag(nbeta))

iota1=matrix(1,nrow=ntau,ncol=1)
Vtau=(diag(ntau))*0.25

iota2=matrix(1,nrow=ntheta,ncol=1)
Vtheta=(diag(ntheta))*0.25


set.seed(proc.time()[3])
price=cbind(matrix(runif(nobs*(nSKU),min=2,max=4),ncol=nSKU))

powern = 100

lgtdata=NULL
itime = proc.time()[3]

for(i in 1:nlgt){


     tempattri=cbind(rep(1,nobs)%x% diag(nSKU))

     beta=as.vector(t(Delta)%*%ZB[i,]+as.vector(t(chol(Vbeta))%*%rnorm(nbeta)))
     theta=as.vector(t(Gamma)%*%ZT[i,]+as.vector(t(chol(Vtheta))%*%rnorm(ntheta)))


     if(Settings$draw_model==4){
           doloop2=TRUE
           while(doloop2){

              nrep = 0
              repeat{
                    nrep = nrep + 1
                    tau=as.vector(t(Omega)%*%ZT[i,]+as.vector(t(chol(Vtau))%*%rnorm(ntau)))
                    temp_cov=try (gen_cov(tau,nSKU),silent=TRUE)
                    if(!is(temp_cov,"try-error")){
                      break
                    }
                    if(nrep>100000){
                      break
                    }
              }

              if(det(temp_cov)>0){
                   doloop2 = FALSE

                   doloop3=TRUE
                   decimalt = 10
                   while(doloop3){

                      if(isTRUE(all.equal(round(temp_cov,decimalt), t(round(temp_cov,decimalt))))){
                           doloop3 = FALSE
                      }else{
                           decimalt = decimalt - 1
                      }

                   }

              }
           }

     }else{
           tau=as.vector(t(Omega)%*%ZT[i,]+as.vector(t(chol(Vtau))%*%rnorm(ntau)))
     }

     if(Settings$est_M){
        M = rep(exp(theta[ntheta]),nobs)
     }else{
        M = runif(nobs,min=Settings$M_min,max=Settings$M_max)
     }


      beta_attri = beta[1:nSKU]
      temp_gamma = exp(beta[(nSKU+1):nbeta])


     tempX=NULL
     tempPsi=NULL
     tempKKT = NULL
     tempKKT2 = NULL

     if(Settings$draw_model==1){

         errordraw_ = matrix(runif(nSKU*nobs),ncol=nSKU)

     }else if(Settings$draw_model==2){

         errordraw_ = rmycopula_power(nSKU,mylogit2(tau),nobs,powern)

     }else if(Settings$draw_model==3){

         errordraw_ = rmycopula_sine(nSKU,mylogit2(tau),nobs)

     }else if(Settings$draw_model==4){

         errordraw_ = pnorm(rmvnorm(nobs,mean=rep(0,nSKU),sigma=round(temp_cov,decimalt)))

     }



     if(Settings$marg == "N"){
           errordraw = qnorm(errordraw_)
           errordraw[errordraw == Inf]=qnorm(0.99)
     }else if(Settings$marg == "L"){
           errordraw = mylogitinv(errordraw_)
           errordraw[errordraw == Inf]=mylogitinv(0.99)
     }




     for(j in 1:nobs){

        tempattri_t=tempattri[((j-1)*nSKU+1):(j*nSKU),]
        psi=exp(tempattri_t%*%beta_attri + errordraw[j,] )

        ui=rbind(-c(price[j,],1),diag(nSKU+1))
        ci=c(-M[j],rep(0,(nSKU+1)))

        nrep = 0
        repeat{
          nrep = nrep + 1
          if(nrep%%2 == 0){
              iniv = runif(nSKU+1,0.1,0.5)
          }else{
              weight =runif(nSKU+1,1,4)
              iniv = ((exp(weight)/sum(exp(weight))) * (M[j]-0.1))/c(price[j,],1)
          }

          outopt=try (

                         constrOptim(iniv, obj_gamma , gradi_gamma, ui=ui, ci=ci, method = "BFGS",
                                     control = list(fnscale = -1,trace = 0, reltol = 1e-40),
                                     psi=psi,psi_out = psi_out,nSKU = nSKU,gamma = temp_gamma)




                      ,silent=TRUE)
          if(!is(outopt,"try-error")){
            break
          }
        }

        x_star=c(outopt$par)
        opt_value=outopt$value
        KKT = gradi(x_star,psi,psi_out,nSKU)/c(price[j,],1)
        lambda = KKT[which(x_star>0.1)[1]]

        KKT2 = x_star * (gradi_gamma(x_star,psi,psi_out,nSKU,temp_gamma)-lambda * c(price[j,],1))


        tempX=rbind(tempX,c(outopt$par))
        tempPsi=rbind(tempPsi,as.vector(psi))
        tempKKT = rbind(tempKKT,KKT)
        tempKKT2 = rbind(tempKKT2,KKT2)

     }

     tempX2 = tempX
     tempX2[tempX2<0.001]=0

     lgtdata[[i]]=list(X=tempX2[,1:nSKU],oriX=tempX,attri=tempattri,beta=beta,tau=tau,theta=theta,Psi=tempPsi,KKT=tempKKT,KKT2=tempKKT2,errordraw=errordraw,z=tempX2[,nSKU+1],price=price)

     if ((i%%5) == 0) {
        ctime = proc.time()[3]
        timetoend = ((ctime - itime)/i) * (nlgt - i)
        cat(" ", i, " (", round(timetoend/60, 1), ")", fill = TRUE)


     }

}

ctime = proc.time()[3]
cat(" Total Time Elapsed: ", round((ctime - itime)/60, 2),fill = TRUE)


total_npur_obs = NULL
for(i in 1:length(lgtdata)){
     if(nobs==1){
         lgtdata[[i]]$X = matrix(lgtdata[[i]]$X,nrow=1)
         lgtdata[[i]]$oriX = matrix(lgtdata[[i]]$oriX,nrow=1)

     }
     total_npur_obs = c(total_npur_obs,apply(lgtdata[[i]]$X>0,1,sum))
}
mean(total_npur_obs)
table(total_npur_obs)

for(i in 1:length(lgtdata)){
	hh_spending=rowSums(lgtdata[[i]]$price*lgtdata[[i]]$X)
	hh_spending_max=max(hh_spending)
  if(hh_spending_max==0){
     hh_spending_max = 0.1
  }
	lgtdata[[i]]$max_spending = hh_spending_max
}



for(i in 1:length(lgtdata)){

  lgtdata[[i]]$pur_l = apply(lgtdata[[i]]$X>0,1,function(x) which(x)-1)
  lgtdata[[i]]$npur_l = apply(!(lgtdata[[i]]$X>0),1,function(x) which(x)-1)

}

Data=list(lgtdata=lgtdata,psi_out=psi_out,ZB=ZB,ZT=ZT,ZQ=ZQ,nlgt=nlgt,nSKU=nSKU,nbeta = nbeta, ntau = ntau ,ntheta=ntheta,nzbeta=nzbeta,nztau=nztau,nztheta=nztheta)


## Prior

nu0=nbeta+3
V0 = nu0 * diag(rep(1, nbeta))
ADelta = 0.01 * diag(nzbeta)
Deltabar = matrix(rep(0, nzbeta * nbeta), ncol = nbeta)

## hierarchi2

nu1=ntau+3
V1 = nu1 * diag(rep(1, ntau))
AOmega = 0.01 * diag(nztau)
Omegabar = matrix(rep(0, nztau * ntau), ncol = ntau)

## hierarchi2

nu2=ntau+3
V2= nu2 * diag(rep(1, ntheta))
AGamma = 0.01 * diag(nztheta)
Gammabar = matrix(rep(0, nztheta * ntheta), ncol = ntheta)



Prior=list(nu0=nu0,V0=V0,ADelta=ADelta,Deltabar=Deltabar,nu1=nu1,V1=V1,AOmega=AOmega,Omegabar=Omegabar,nu2=nu2,V2=V2,AGamma=AGamma,Gammabar=Gammabar)


## Mcmc

sbeta = 0.05
stheta = 0.005
stau = 0.3

keep = 10
R=10000
Mcmc=list(sbeta=sbeta,stau = stau,stheta=stheta,keep=keep,R=R)


Ini=NULL
set.seed(66)
out=DirectCopula(Data,Prior,Mcmc,Ini,Settings)



summary(out$Deltadraw,tvalues=as.vector(Delta))
plot(out$Deltadraw,tvalues=as.vector(Delta))

summary(out$Gammadraw,tvalues=as.vector(Gamma))
plot(out$Gammadraw,tvalues=as.vector(Gamma))

summary(out$Omegadraw,tvalues=as.vector(Omega))
plot(out$Omegadraw,tvalues=as.vector(Omega))












