


conv_pram_cor = function(temp_Omega,nSKU,temp_model){

    temp_correl = NA
    if(temp_model==1){
          temp_correl = as.vector(temp_Omega)
    }else if(temp_model==2){
          temp_cor_s = integrate(gen_power, lower = 0, upper = 1,n=powern)[1]$value
          temp_correl = mylogit2(t(mat_conv(as.vector(temp_Omega),nSKU))[lower.tri(diag(nSKU),diag=FALSE)]) * 12 * (temp_cor_s)^2

    }else if(temp_model==3){
          temp_correl = mylogit2(t(mat_conv(as.vector(temp_Omega),nSKU))[lower.tri(diag(nSKU),diag=FALSE)]) * 48/pi^4
    }else if(temp_model==4){
          temp_correl = 6/pi*asin(gen_cov(as.vector(temp_Omega),nSKU)[lower.tri(diag(nSKU),diag=FALSE)]/2)
    }

    return(temp_correl)

}

myrEV = function (k) {
    rvec = as.vector(-log(-log(runif(k, 0, 1))))
    return(rvec)
}


mat_conv = function(tau,nvar){
    tau_mat = matrix(rep(0,nvar*nvar),nvar,nvar)
    tau_mat[upper.tri(tau_mat,diag=FALSE)] = tau
    return(tau_mat)
}


gen_cov = function(tau,nSKU){
    temp_r = mat_conv(tau,nSKU)
    diag(temp_r) = 1

    temp_sig_inv = t(temp_r) %*% temp_r
    temp_sig = solve(temp_sig_inv,diag(nSKU))

    temp_sig_diag = diag(nSKU)
    diag(temp_sig_diag) = diag(temp_sig)^(-1/2)
    temp_cov = (temp_sig_diag) %*% temp_sig %*% (temp_sig_diag)
    return(temp_cov)

}


P_make_sine = function(tau,u,m,nvar){

    tau_tilde = mat_conv(tau,nvar)
    tau_tilde[c(m:nvar),]=0
    tau_tilde[,c(m:nvar)]=0

    ret = (dgen_sine(u) %*% tau_tilde %*% t(dgen_sine(u))) + 1
    return(ret)
}



A_make_sine = function(tau,u,m,nvar){
    tau_mat = mat_conv(tau,nvar)
    ret = dgen_sine(u) %*% tau_mat[,m]
    return(ret)

}



P_make_power = function(tau,u,m,nvar,n){

    tau_tilde = mat_conv(tau,nvar)
    tau_tilde[c(m:nvar),]=0
    tau_tilde[,c(m:nvar)]=0

    ret = (dgen_power(u,n) %*% tau_tilde %*% t(dgen_power(u,n))) + 1
    return(ret)
}

A_make_power = function(tau,u,m,nvar,n){
    tau_mat = mat_conv(tau,nvar)
    ret = dgen_power(u,n) %*% tau_mat[,m]
    return(ret)

}


rmycopula_sine = function(nvar,true_tau,nobs){


  errordraw = NULL
  for(i in 1:nobs){
    W = matrix(runif(nvar),nrow=1)
    u_old = matrix(rep(0,nvar),nrow=1)
    u_old[1,1] = W[1,1]

    for(j in 2:nvar){
      AA = A_make_sine(true_tau,u_old,j,nvar)
      BB = P_make_sine(true_tau,u_old,j,nvar)

      fun = function(x){
        ret = x * BB + gen_sine(x) * AA - BB*W[1,j]
        return(ret)
      }
      u_new = uniroot(fun,lower=0,upper=1)$root
      u_old[1,j] = u_new
    }

    errordraw = rbind(errordraw,u_old)

  }
  return(errordraw)

}



rmycopula_power = function(nvar,true_tau,nobs,n){


  errordraw = NULL
  for(i in 1:nobs){
    W = matrix(runif(nvar),nrow=1)
    u_old = matrix(rep(0,nvar),nrow=1)
    u_old[1,1] = W[1,1]

    for(j in 2:nvar){
      AA = A_make_power(true_tau,u_old,j,nvar,n)
      BB = P_make_power(true_tau,u_old,j,nvar,n)

      fun = function(x){
        ret = x * BB + gen_power(x,n) * AA - BB*W[1,j]
        return(ret)
      }
      u_new = uniroot(fun,lower=0,upper=1)$root
      u_old[1,j] = u_new
    }

    errordraw = rbind(errordraw,u_old)

  }
  return(errordraw)

}



mycopula_sine = function(tau,errordraw,nvar){

  tau_mat = mat_conv(tau,nvar)
  temp_ret = diag(dgen_sine(errordraw) %*% tau_mat %*% t(dgen_sine(errordraw)))+1
  ret = temp_ret
  return(sum(log(ret)))
}



mycopula_power = function(tau,errordraw,nvar,n){

  tau_mat = mat_conv(tau,nvar)
  temp_ret = diag(dgen_power(errordraw,n) %*% tau_mat %*% t(dgen_power(errordraw,n)))+1
  ret = temp_ret
  return(sum(log(ret)))
}

mylogit = function(x){

      x[x==Inf]=600
      y=exp(x)/(1+exp(x))
      return(y)
}


mylogitinv = function(y){
  x = log(y) - log(1-y)
  return(x)

}

mylogitpdf = function(x){
  y =  exp(x) / (( 1 + exp(x))^2)
  return(y)

}

mylogit2 = function(x){
      y=(-1+exp(x))/(1+exp(x))
      return(y)
}


gen_power = function(u,n){
  return(1-(u^n+(1-u)^n)^(1/n))
}

dgen_power = function(u,n){
  return(((1-u)^(n-1)-u^(n-1))*(u^n+(1-u)^n)^(1/n-1))
}



gen_sine = function(u){
  return( sin(pi*u)/pi )
}

dgen_sine = function(u){
  return( cos(pi*u) )
}



obj_gamma = function(X,psi,psi_out,nSKU,gamma){

  ret = sum( psi/gamma * log(gamma* X[1:nSKU] + 1) , psi_out * log(X[nSKU+1]) )
  return( ret )
}

gradi_gamma=function(X,psi,psi_out,nSKU,gamma){

  ret = c (psi/(gamma*X[1:nSKU] + 1), psi_out/(X[nSKU+1]) )

  return( ret )
}



tau_cano_det = function(tau,nSKU,decimalt){
            #decimalt = 6
            temp_cov = gen_cov(tau,nSKU)
            yscena = as.matrix(expand.grid(as.data.frame(matrix(rep(c(FALSE,TRUE),nSKU),ncol=nSKU))))
            ret_sigma = list()
            ret_TF = NULL
            for(i in 1:nrow(yscena)){
                    if(sum(yscena[i,])==nSKU){
                        #ret_TF_temp = (det(temp_cov)>0)
                        ret_TF_temp = isTRUE(all.equal(round(temp_cov,decimalt), t(round(temp_cov,decimalt))))

                        ret_TF = c(ret_TF,ret_TF_temp)
                        ret_sigma[[i]] = temp_cov
                    }else if(sum(yscena[i,])==0){
                        #ret_TF_temp = (det(temp_cov)>0)
                        ret_TF_temp = isTRUE(all.equal(round(temp_cov,decimalt), t(round(temp_cov,decimalt))))
                        ret_TF = c(ret_TF,ret_TF_temp)
                        ret_sigma[[i]] = temp_cov
                    }else{

                        sigma11 = temp_cov[!yscena[i,],!yscena[i,]]
                        sigma12 = temp_cov[!yscena[i,],yscena[i,]]
                        sigma22 = temp_cov[yscena[i,],yscena[i,]]
                        sigma21 = temp_cov[yscena[i,],!yscena[i,]]
                        isigma22 = solve(sigma22,diag(sum(yscena[i,])))

                        temp_cov_sub = sigma11 - sigma12 %*% isigma22 %*% sigma21
                        #ret_TF_temp = (det(temp_cov_sub)>0)
                        ret_TF_temp = isTRUE(all.equal(round(temp_cov_sub,decimalt), t(round(temp_cov_sub,decimalt))))
                        ret_TF = c(ret_TF,ret_TF_temp)

                        ret_sigma[[i]] = temp_cov_sub
                    }

             }

             return(list(yscena=yscena,ret_TF=ret_TF,ret_sigma=ret_sigma))

}

DirectCopula=function(Data,Prior,Mcmc,Ini,Settings){

    ###########
    lgtdata=Data$lgtdata
    ZB=Data$ZB
    ZT=Data$ZT
    ZQ=Data$ZQ

    gammas = Data$gammas
    psi_out = Data$psi_out
    nlgt = Data$nlgt
    nSKU = Data$nSKU

    nbeta = Data$nbeta
    ntau = Data$ntau
    ntheta = Data$ntheta
    nzbeta = Data$nzbeta
    nztau = Data$nztau
    nztheta = Data$nztheta

    nu0 = Prior$nu0
    V0 = Prior$V0
    ADelta = Prior$ADelta
    Deltabar = Prior$Deltabar

    nu1=Prior$nu1
    V1=Prior$V1
    AOmega=Prior$AOmega
    Omegabar=Prior$Omegabar

    nu2=Prior$nu2
    V2=Prior$V2
    AGamma=Prior$AGamma
    Gammabar=Prior$Gammabar

    sbeta = Mcmc$sbeta
    stau = Mcmc$stau
    stheta = Mcmc$stheta
    keep = Mcmc$keep
    R = Mcmc$R

    beta_ini = Ini$beta_ini
    tau_ini = Ini$tau_ini
    theta_ini = Ini$theta_ini

    draw_model = Settings$draw_model
    est_model = Settings$est_model
    z_var = Settings$z_var
    est_M = Settings$est_M
    M_min = Settings$M_min
    M_max = Settings$M_max
    marg = Settings$marg

    Vbetadraw = matrix(double(floor(R/keep) * nbeta * nbeta), ncol = nbeta * nbeta)
    betadraw = array(double(floor(R/keep) * nlgt * nbeta), dim = c(nlgt,nbeta, floor(R/keep)))

    Vtaudraw = matrix(double(floor(R/keep) * ntau * ntau), ncol = ntau * ntau)
    taudraw = array(double(floor(R/keep) * nlgt * ntau), dim = c(nlgt,ntau, floor(R/keep)))

    Vthetadraw = matrix(double(floor(R/keep) * ntheta * ntheta), ncol = ntheta * ntheta)
    thetadraw = array(double(floor(R/keep) * nlgt * ntheta), dim = c(nlgt,ntheta, floor(R/keep)))

    Deltadraw = matrix(double(floor(R/keep) * nbeta * nzbeta), ncol = nbeta * nzbeta)
    Omegadraw = matrix(double(floor(R/keep) * ntau * nztau), ncol = ntau * nztau)
    Gammadraw = matrix(double(floor(R/keep) * ntheta * nztheta), ncol = ntheta * nztheta)

    oldbetas = matrix(rep(0,nlgt * nbeta), ncol = nbeta)
    oldtaus = matrix(rep(0,nlgt * ntau), ncol = ntau)
    oldthetas = matrix(rep(0,nlgt * ntheta), ncol = ntheta)
    oldlikes = double(nlgt)


    cat("initializing calgd, caljacod \n")
    for(i in 1:nlgt){

          doloop=TRUE
          temp_oldbetas = rep(0,nbeta)
          temp_oldtaus = rep(0,ntau)
          temp_oldthetas = rep(0,ntheta)
          while(doloop){

            if(is.null(beta_ini)){
              temp_oldbetas[1:nSKU] = runif(nSKU,-1,0)
            }else{
              temp_oldbetas = beta_ini[i,]
            }

            if(Settings$est_model==1){

            }else{
                if(is.null(tau_ini)){

                     if(Settings$est_model==4){
                           doloop2=TRUE
                           decimalt = 10
                           while(doloop2){

                              temp_oldtaus = runif(ntau,-0.3,0.3)

                              temp_cov = gen_cov(temp_oldtaus,nSKU)


                              temp_cano_cov = tau_cano_det(temp_oldtaus,nSKU,decimalt)
                              if(prod(temp_cano_cov$ret_TF)==1){
                                   doloop2 = FALSE
                              }else{
                                   decimalt = decimalt - 1
                              }

                           }

                     }else if(Settings$est_model==2){
                           temp_oldtaus = runif(ntau,-0.3,0.3)
                           temp_tau_mat = mat_conv(mylogit2(temp_oldtaus),nSKU)
                     }else if(Settings$est_model==3){
                           temp_oldtaus = runif(ntau,-0.3,0.3)
                           temp_tau_mat = mat_conv(mylogit2(temp_oldtaus),nSKU)
                     }

                }else{
                  temp_oldtaus = tau_ini[i,]
                }
            }


            if(Settings$est_M){
                if(is.null(theta_ini)){
                  temp_oldthetas[ntheta] = log(lgtdata[[i]]$max_spending) + 0.01
                }else{
                  temp_oldthetas = theta_ini[i,]
                }
            }else{

            }

            if(Settings$est_M){
                 if(Settings$est_model==1){

                      if(marg == "N"){
                              temp_like=likelihood_estM_indi_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,nSKU,nbeta,ntheta)
                      }else if(marg == "L"){
                              temp_like=likelihood_estM_indi_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,nSKU,nbeta,ntheta)
                      }


                 }else if(Settings$est_model==2){

                      if(marg == "N"){
                            temp_like=likelihood_estM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,temp_tau_mat,nSKU,nbeta,ntheta,powern)
                      }else if(marg == "L"){
                            temp_like=likelihood_estM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,temp_tau_mat,nSKU,nbeta,ntheta,powern)
                      }

                 }else if(Settings$est_model==3){

                      if(marg == "N"){
                              temp_like=likelihood_estM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,temp_tau_mat,nSKU,nbeta,ntheta)
                      }else if(marg == "L"){
                              temp_like=likelihood_estM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,temp_tau_mat,nSKU,nbeta,ntheta)
                      }

                 }else if(Settings$est_model==4){

                      if(marg == "N"){
                            invisible(capture.output(temp_like <-likelihood_estM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                      }else if(marg == "L"){
                            invisible(capture.output(temp_like <-likelihood_estM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_oldthetas,round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                      }


                 }
            }else{
                 if(Settings$est_model==1){
                      if(marg == "N"){
                            temp_like=likelihood_setM_indi_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,nSKU,nbeta)
                      }else if(marg == "L"){
                            temp_like=likelihood_setM_indi_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,nSKU,nbeta)
                      }

                 }else if(Settings$est_model==2){

                      if(marg == "N"){
                           temp_like=likelihood_setM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_tau_mat,nSKU,nbeta,powern)
                      }else if(marg == "L"){
                           temp_like=likelihood_setM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_tau_mat,nSKU,nbeta,powern)
                      }

                 }else if(Settings$est_model==3){

                      if(marg == "N"){
                           temp_like=likelihood_setM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_tau_mat,nSKU,nbeta)

                      }else if(marg == "L"){
                           temp_like=likelihood_setM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,temp_tau_mat,nSKU,nbeta)
                      }


                 }else if(Settings$est_model==4){

                      if(marg == "N"){

                            invisible(capture.output(temp_like <-likelihood_setM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))

                      }else if(marg == "L"){

                            invisible(capture.output(temp_like <-likelihood_setM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,temp_oldbetas,round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                      }

                 }
            }


            if ( sum (temp_like < 0, na.rm=TRUE) > 0 ){

            }
            else{
              doloop = FALSE
              oldbetas[i,]=temp_oldbetas
              oldtaus[i,]=temp_oldtaus
              oldthetas[i,]=temp_oldthetas
              oldlikes[i] = sum(log(temp_like),na.rm=T)


            }



          }
    }

    oldVbeta = diag(nbeta)
    oldVbetai = diag(nbeta)

    oldVtau = diag(ntau)
    oldVtaui = diag(ntau)

    oldVtheta = diag(ntheta)
    oldVthetai = diag(ntheta)

    oldDelta = matrix(double(nbeta * nzbeta), ncol = nbeta)
    oldOmega = matrix(double(ntau * nztau), ncol = ntau)
    oldGamma = matrix(double(ntheta * nztheta), ncol = ntheta)

    betad = array(0, dim = c(nbeta))
    betan = array(0, dim = c(nbeta))

    taud = array(0, dim = c(ntau))
    taun = array(0, dim = c(ntau))

    thetad = array(0, dim = c(ntheta))
    thetan = array(0, dim = c(ntheta))

    reject0 = array(0, dim = c(R/keep))
    reject1 = array(0, dim = c(R/keep))
    reject2 = array(0, dim = c(R/keep))

    llike = array(0, dim = c(R/keep))

    itime = proc.time()[3]

    cat("MCMC Iteration (est time to end - min)", fill = TRUE)


    for (j in 1:R) {

        rej0 = 0
        rej1 = 0
        rej2 = 0
        logl = 0

        sVbeta = sbeta * diag(nbeta)
        rootsVbeta = t(chol(sVbeta))

        sVtau = stau * diag(ntau)
        rootsVtau = t(chol(sVtau))

        sVtheta = stheta * diag(ntheta)
        rootsVtheta = t(chol(sVtheta))

        for (i in 1:nlgt) {

            nobs=nrow(lgtdata[[i]]$X)

            betad = oldbetas[i, ]
            liked = oldlikes[i]

            if(Settings$est_model==1){

            }else if(Settings$est_model==4){

                 doloop2=TRUE
                 decimalt = 10
                 while(doloop2){

                    temp_cov = gen_cov(oldtaus[i,],nSKU)
                    temp_cano_cov = tau_cano_det(oldtaus[i,],nSKU,decimalt)

                    if(prod(temp_cano_cov$ret_TF)==1){
                         doloop2 = FALSE
                    }else{
                         decimalt = decimalt - 1
                    }

                 }


            }else if(Settings$est_model==2){
                temp_tau_mat = mat_conv(mylogit2(oldtaus[i,]),nSKU)
            }else if(Settings$est_model==3){
                temp_tau_mat = mat_conv(mylogit2(oldtaus[i,]),nSKU)
            }

            doloop = TRUE
            while(doloop){

               betan = as.vector(betad + rootsVbeta %*% rnorm(nbeta))

               if(Settings$est_M){

                   if(Settings$est_model==1){

                        if(marg == "N"){
                                 liken=likelihood_estM_indi_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],nSKU,nbeta,ntheta)
                        }else if(marg == "L"){
                                 liken=likelihood_estM_indi_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],nSKU,nbeta,ntheta)
                        }

                   }else if(Settings$est_model==2){

                        if(marg == "N"){
                                 liken = likelihood_estM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }else if(marg == "L"){
                                 liken = likelihood_estM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }

                   }else if(Settings$est_model==3){

                        if(marg == "N"){
                                 liken=likelihood_estM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta)
                        }else if(marg == "L"){
                                 liken=likelihood_estM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta)
                        }

                   }else if(Settings$est_model==4){

                        if(marg == "N"){
                            invisible(capture.output(liken <-likelihood_estM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))

                        }else if(marg == "L"){
                            invisible(capture.output(liken <-likelihood_estM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,oldthetas[i,],round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }

                   }

               }else{

                   if(Settings$est_model==1){

                        if(marg == "N"){
                                  liken=likelihood_setM_indi_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,nSKU,nbeta)
                        }else if(marg == "L"){
                                  liken=likelihood_setM_indi_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,nSKU,nbeta)
                        }

                   }else if(Settings$est_model==2){

                          if(marg == "N"){
                               liken = likelihood_setM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,temp_tau_mat,nSKU,nbeta,powern)
                          }else if(marg == "L"){
                               liken = likelihood_setM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,temp_tau_mat,nSKU,nbeta,powern)
                          }

                   }else if(Settings$est_model==3){

                          if(marg == "N"){
                               liken = likelihood_setM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,temp_tau_mat,nSKU,nbeta)
                          }else if(marg == "L"){
                               liken = likelihood_setM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,temp_tau_mat,nSKU,nbeta)
                          }

                   }else if(Settings$est_model==4){

                          if(marg == "N"){
                                invisible(capture.output(liken <-likelihood_setM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                          }else if(marg == "L"){
                                invisible(capture.output(liken <-likelihood_setM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,betan,round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                          }

                   }
               }



               if ( sum (liken < 0, na.rm=TRUE) > 0 ){

               }
               else{
                 doloop = FALSE
               }



            }

            priord = -0.5 * (t(betad) - ZB[i, ] %*% oldDelta) %*% oldVbetai %*% (betad - t(ZB[i, ] %*% oldDelta))
            priorn = -0.5 * (t(betan) - ZB[i, ] %*% oldDelta) %*% oldVbetai %*% (betan - t(ZB[i, ] %*% oldDelta))

              if(Settings$est_model==5){
                    c_liken = sum((liken),na.rm=T)
              }else{
                    c_liken = sum(log(liken),na.rm=T)
              }

            if(c_liken==Inf){
               break
            }
            ralpha = exp( c_liken  + priorn - liked - priord)

            if (ralpha == "NaN"){
                ralpha = -1
            }

            u = runif(n = 1, min = 0, max = 1)

            if (u < ralpha) {
                oldbetas[i,] = betan
                oldlikes[i] = c_liken
                logl = logl + c_liken
            }else {

                rej0 = rej0 + 1
                logl = logl + liked
            }

            if(Settings$est_M){

                thetad = oldthetas[i, ]
                liked = oldlikes[i]

                if(Settings$est_model==1){

                }else if(Settings$est_model==4){

                     doloop2=TRUE
                     decimalt = 10
                     while(doloop2){

                        temp_cov = gen_cov(oldtaus[i,],nSKU)
                        temp_cano_cov = tau_cano_det(oldtaus[i,],nSKU,decimalt)

                        if(prod(temp_cano_cov$ret_TF)==1){
                             doloop2 = FALSE
                        }else{
                             decimalt = decimalt - 1
                        }

                     }

                }else if(Settings$est_model==2){
                    temp_tau_mat = mat_conv(mylogit2(oldtaus[i,]),nSKU)
                }else if(Settings$est_model==3){
                    temp_tau_mat = mat_conv(mylogit2(oldtaus[i,]),nSKU)
                }

                doloop = TRUE
                while(doloop){

                   doloop_sub = TRUE
                   while(doloop_sub){

                      thetan = as.vector(thetad + rootsVtheta %*% rnorm(ntheta))

                      if(exp(thetan[ntheta]) > lgtdata[[i]]$max_spending){
                         doloop_sub = FALSE
                      }
                   }

                   if(Settings$est_model==1){

                        if(marg == "N"){
                                 liken=likelihood_estM_indi_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,nSKU,nbeta,ntheta)
                        }else if(marg == "L"){
                                 liken=likelihood_estM_indi_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,nSKU,nbeta,ntheta)
                        }


                   }else if(Settings$est_model==2){

                        if(marg == "N"){
                                 liken = likelihood_estM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }else if(marg == "L"){
                                 liken = likelihood_estM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }

                   }else if(Settings$est_model==3){

                        if(marg == "N"){
                                liken=likelihood_estM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,temp_tau_mat,nSKU,nbeta,ntheta)
                        }else if(marg == "L"){
                                liken=likelihood_estM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,temp_tau_mat,nSKU,nbeta,ntheta)
                        }

                   }else if(Settings$est_model==4){

                        if(marg == "N"){
                               invisible(capture.output(liken <-likelihood_estM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }else if(marg == "L"){
                               invisible(capture.output(liken <-likelihood_estM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],thetan,round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }


                   }


                   if ( sum (liken < 0, na.rm=TRUE) > 0 ){

                   }
                   else{
                     doloop = FALSE
                   }



                }

                priord = -0.5 * (t(thetad) - ZQ[i, ] %*% oldGamma) %*% oldVthetai %*% (thetad - t(ZQ[i, ] %*% oldGamma))
                priorn = -0.5 * (t(thetan) - ZQ[i, ] %*% oldGamma) %*% oldVthetai %*% (thetan - t(ZQ[i, ] %*% oldGamma))



                c_liken = sum(log(liken),na.rm=T)

               if(c_liken==Inf){
                   break
                }

                ralpha = exp( c_liken  + priorn - liked - priord)

                if (ralpha == "NaN"){
                    ralpha = -1
                }

                u = runif(n = 1, min = 0, max = 1)

                if (u < ralpha) {
                    oldthetas[i, ] = thetan
                    oldlikes[i] = c_liken

                }else {

                    rej2 = rej2 + 1

                }

            }else{


            }

            ######### taudraw

            if(Settings$est_model==1){

            }else{

              taud = oldtaus[i, ]
              liked = oldlikes[i]

              doloop=TRUE

              while(doloop){

                 if(Settings$est_model==4){
                     doloop2=TRUE
                     decimalt = 10
                     taun = as.vector(taud + rootsVtau %*% rnorm(ntau))
                     while(doloop2){

                        temp_cov = gen_cov(taun,nSKU)
                        #if(det(temp_cov)>0){
                        #     doloop2 = FALSE
                        #}
                        temp_cano_cov = tau_cano_det(taun,nSKU,decimalt)
                        if(prod(temp_cano_cov$ret_TF)==1){
                             doloop2 = FALSE
                        }else{
                             decimalt = decimalt - 1
                        }
                        if(decimalt < (-1)){
                             break
                        }
                     }
                 }else{
                      taun = as.vector(taud + rootsVtau %*% rnorm(ntau))
                      if(Settings$est_model==3){
                          temp_tau_mat = mat_conv(mylogit2(taun),nSKU)
                      }else if(Settings$est_model==2){
                          temp_tau_mat = mat_conv(mylogit2(taun),nSKU)
                      }

                 }

                 if(Settings$est_M){

                     if(Settings$est_model==2){

                        if(marg == "N"){
                               liken = likelihood_estM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }else if(marg == "L"){
                               liken = likelihood_estM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta,powern)
                        }

                     }else if(Settings$est_model==3){

                        if(marg == "N"){
                               liken=likelihood_estM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta)
                        }else if(marg == "L"){
                               liken=likelihood_estM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],temp_tau_mat,nSKU,nbeta,ntheta)
                        }

                     }else if(Settings$est_model==4){

                        if(marg == "N"){
                                invisible(capture.output(liken <-likelihood_estM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }else if(marg == "L"){
                                invisible(capture.output(liken <-likelihood_estM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],oldthetas[i,],round(temp_cov,decimalt),nSKU,nbeta,ntheta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }

                     }

                 }else{

                     if(Settings$est_model==2){

                          if(marg == "N"){
                              liken = likelihood_setM_power_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],temp_tau_mat,nSKU,nbeta,powern)
                          }else if(marg == "L"){
                              liken = likelihood_setM_power_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],temp_tau_mat,nSKU,nbeta,powern)
                          }

                     }else if(Settings$est_model==3){

                          if(marg == "N"){
                              liken = likelihood_setM_sine_cpp_norm(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],temp_tau_mat,nSKU,nbeta)
                          }else if(marg == "L"){
                              liken = likelihood_setM_sine_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],temp_tau_mat,nSKU,nbeta)
                          }

                     }else if(Settings$est_model==4){

                        if(marg == "N"){
                              invisible(capture.output(liken <-likelihood_setM_gau_cpp(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }else if(marg == "L"){
                              invisible(capture.output(liken <-likelihood_setM_gau_cpp_logit(lgtdata[[i]]$X,(lgtdata[[i]]$X>0)*1,lgtdata[[i]][[Settings$z_var]],lgtdata[[i]]$attri,lgtdata[[i]]$price,oldbetas[i,],round(temp_cov,decimalt),nSKU,nbeta,lgtdata[[i]]$pur_l,lgtdata[[i]]$npur_l)))
                        }

                     }
                 }


                 if ( sum (liken < 0, na.rm=TRUE) > 0 ){

                 }
                 else{
                   doloop = FALSE
                 }



              }

              priord = -0.5 * (t(taud) - ZT[i, ] %*% oldOmega) %*% oldVtaui %*% (taud - t(ZT[i, ] %*% oldOmega))
              priorn = -0.5 * (t(taun) - ZT[i, ] %*% oldOmega) %*% oldVtaui %*% (taun - t(ZT[i, ] %*% oldOmega))

              c_liken = sum(log(liken),na.rm=T)

              if(c_liken==Inf){
                 break
              }
              ralpha = exp( c_liken + priorn - liked - priord)

              if (ralpha == "NaN"){
                  ralpha = -1
              }

              u = runif(n = 1, min = 0, max = 1)

              if (u < ralpha) {
                  oldtaus[i, ] = taun
                  oldlikes[i] = c_liken

              }else {

                  rej1 = rej1 + 1

              }
           }

        }

        out0 = rmultireg(oldbetas, ZB, Deltabar, ADelta, nu0, V0)
        oldDelta = out0$B
        oldVbeta = out0$Sigma
        oldVbetai = chol2inv(chol(oldVbeta))

        if(Settings$est_M){

            out2 = rmultireg(oldthetas, ZQ, Gammabar, AGamma, nu2, V2)
            oldGamma = out2$B
            oldVtheta = out2$Sigma
            oldVthetai = chol2inv(chol(oldVtheta))

        }else{

        }

        if(Settings$est_model==1){

        }else{
            out1 = rmultireg(oldtaus, ZT, Omegabar, AOmega, nu1, V1)
            oldOmega = out1$B
            oldVtau = out1$Sigma
            oldVtaui = chol2inv(chol(oldVtau))

        }
        if ((j%%100) == 0) {
            ctime = proc.time()[3]
            timetoend = ((ctime - itime)/j) * (R - j)
            cat(" ", j, " (", round(timetoend/60, 1), ")", fill = TRUE)
            cat("iteration: ",j, "Delta :",as.vector(oldDelta),"Gamma :",as.vector(oldGamma),"Omega :",as.vector(oldOmega)," rej0 = ", round(rej0/nlgt,2)," rej1 = ", round(rej1/nlgt,2)," rej2 = ", round(rej2/nlgt,2)," logl = ",round(logl,2),fill=TRUE)

        }

        mkeep = j/keep

        if (mkeep * keep == (floor(mkeep) * keep)) {
            Deltadraw[mkeep, ] = as.vector(oldDelta)
            Vbetadraw[mkeep, ] = as.vector(oldVbeta)
            betadraw[, , mkeep] = oldbetas
            Omegadraw[mkeep, ] = as.vector(oldOmega)
            Vtaudraw[mkeep, ] = as.vector(oldVtau)
            taudraw[, , mkeep] = oldtaus
            Gammadraw[mkeep, ] = as.vector(oldGamma)
            Vthetadraw[mkeep, ] = as.vector(oldVtheta)
            thetadraw[, , mkeep] = oldthetas

            llike[mkeep] = logl
            reject0[mkeep] = rej0/nlgt
            reject1[mkeep] = rej1/nlgt
            reject2[mkeep] = rej2/nlgt
        }



    }

    ctime = proc.time()[3]

    cat(" Total Time Elapsed: ", round((ctime - itime)/60, 2),fill = TRUE)
    attributes(betadraw)$class = c("bayesm.hcoef")
    attributes(Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(Deltadraw)$mcpar = c(1, R, keep)
    attributes(Vbetadraw)$class = c("bayesm.var", "bayesm.mat","mcmc")
    attributes(Vbetadraw)$mcpar = c(1, R, keep)
    attributes(taudraw)$class = c("bayesm.hcoef")
    attributes(Omegadraw)$class = c("bayesm.mat", "mcmc")
    attributes(Omegadraw)$mcpar = c(1, R, keep)
    attributes(Vtaudraw)$class = c("bayesm.var", "bayesm.mat","mcmc")
    attributes(Vtaudraw)$mcpar = c(1, R, keep)
    attributes(thetadraw)$class = c("bayesm.hcoef")
    attributes(Gammadraw)$class = c("bayesm.mat", "mcmc")
    attributes(Gammadraw)$mcpar = c(1, R, keep)
    attributes(Vthetadraw)$class = c("bayesm.var", "bayesm.mat","mcmc")
    attributes(Vthetadraw)$mcpar = c(1, R, keep)

    return(list(betadraw = betadraw,Vbetadraw = Vbetadraw,Deltadraw = Deltadraw,
                 taudraw = taudraw,Vtaudraw = Vtaudraw,Omegadraw = Omegadraw,
                 thetadraw = thetadraw,Vthetadraw = Vthetadraw,Gammadraw = Gammadraw,
                llike = llike, reject0 = reject0,reject1 = reject1,reject2 = reject2,
                esttime=round((ctime - itime)/60, 2)))
}


