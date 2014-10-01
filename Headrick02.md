

f_skew <- function(x){
  
  sd.x <- sd(x)
  mu3.x <- mean((x-mean(x))^3)
  mu3.x/sd.x^3
  
}

f_kurt <- function(x){
  
  sd.x <- sd(x)
  mu4 <- mean((x-mean(x))^4)
  mu4/sd.x^4 - 3
  
}


f_gamma3 <- function(x){
  sd.x <- sd(x)
  mu3 <- mean((x-mean(x))^3)
  mu5 <- mean((x-mean(x))^5)
  mu5/sd.x^5-10*mu3/sd.x^3
}

f_gamma4 <- function(x){
  sd.x <- sd(x)
  mu3 <- mean((x-mean(x))^3)
  mu4 <- mean((x-mean(x))^4)
  mu6 <- mean((x-mean(x))^6)
  mu6/sd.x^6-15*mu4/sd.x^4 + 45-10*mu3^2/sd.x^3-15
}

headrick02.poly.coeff <- function(skewness, kurtosis, gam3, gam4, control = list(trace = T, max.ntry = 10, obj.tol = 1e-10, n.valid.sol = 2)){
  
  gam1 <- skewness
  gam2 <- kurtosis
  
  gam <- c(gam1, gam2, gam3, gam4)
  
  obj.fun <- function(x, gam){
    
    if(length(x) != 6){
      stop("coefficients of fifth-order polynomial should be length-six")
    }
    
    c0 <- x[1]
    c1 <- x[2]
    c2 <- x[3]
    c3 <- x[4]
    c4 <- x[5]
    c5 <- x[6]
    
    gam1 <- gam[1]
    gam2 <- gam[2]
    gam3 <- gam[3]
    gam4 <- gam[4]
    
    eq.18 <- 0 + c0 + c2 + 3 * c4
    eq.22 <- -1 + c1^2 + 2 * c2^2 + 24 * c2 * c4 + 
      6 * c1 * (c3 + 5 * c5) + 
      3 * (5 * c3^2 + 32 * c4^2 + 70 * c3 * c5 + 315 * c5^2)
    
    eq.B1 <- -gam1 + 2 * (
      4 * c2^3 + 108 * c2^2 * c4 + 3 * c1^2 * (c2 + 6 * c4) + 
        18 * c1 * (2 * c2 * c3 + 16 * c3 * c4 + 15 * c2 * c5 + 150 * c4 * c5) + 
        9 * c2 * (15 * c3^2 + 128 * c4^2 + 280 * c3 * c5 + 1575 * c5^2) +
        54 * c4 * (25 * c3^2 + 88 * c4^2 + 560 * c3 * c5 + 3675 * c5^2)
      )
    eq.B2 <- -gam2 + 24 * (
      2 * c2^4 + 96 * c2^3 * c4 + c1^3 * (c3 + 10 * c5) + 
        30 * c2^2 * (6 * c3^2 + 64 * c4^2 + 140 * c3 * c5 + 945 * c5^2) +
        c1^2 * (2 * c2^2 + 18 * c3^2 + 36 * c2 * c4 + 192 * c4^2 + 375 * c3 * c5 + 2250 * c5^2) + 
        36 * c2 * c4 * (125 * c3^2 + 528 * c4^2 + 3360 * c3 * c5 + 25725 * c5^2) + 
        3 * c1 * (45 * c3^3 + 1584 * c3 * c4^2 + 1590 * c3^2 * c5 + 21360 * c4^2 * c5 + 21525 * c3 * c5^2 + 
                    110250 * c5^3 + 12 * c2^2 * (c3 + 10 * c5) + 8 * c2 * c4 * (32 * c3 + 375 * c5)) + 
        9 * (45 * c3^4 + 8704 * c4^4 + 2415 * c3^3 * c5 + 932400 * c4^2 * c5^2 + 3018750 * c5^4 + 
               20 * c3^2 * (178 * c4^2 + 2765 * c5^2) + 35 * c3 * (3104 * c4^2 * c5 + 18075 * c5^3))
      )
    eq.B3 <- -gam3 + 24 * (
      16 * c2^5 + 5 * c1^4 * c4 + 1200 * c2^4 * c4 + 10 * c1^3 * (3 * c2 * c3 + 42 * c3 * c4 + 40 * c2 * c5 + 570 * c4 * c5) + 
        300 * c2^3 * (10 * c3^2 + 128 * c4^2 + 280 * c3 * c5 + 2205 * c5^2) + 
        1080 * c2^2 * c4 * (125 * c3^2 + 3920 * c3 * c5 + 28 * (22 * c4^2 + 1225 * c5^2)) + 
        10 * c1^2 * (2 * c2^3 + 72 * c2^2 * c4 + 3 * c2 * (24 * c3^2 + 320 * c4^2 + 625 * c3 * c5 + 4500 * c5^2) + 
                       9 * c4 * (109 * c3^2 + 528 * c4^2 + 3130 * c3 * c5 + 24975 * c5^2)) + 
        30 * c1 * (8 * c2^3 * (2 * c3 + 25 * c5) + 40 * c2^2 * c4 * (16 * c3 + 225 * c5) + 
                     3 * c2 * (75 * c3^3 + 3168 * c3 * c4^2 + 3180 * c3^2 * c5 + 49840 * c4^2 * c5 + 50225 * c3 * c5^2 + 294000 * c5^3) + 
                     6 * c4 * (555 * c3^3 + 8704 * c3 * c4^2 + 26225 * c3^2 * c5 + 152160 * c4^2 * c5 + 459375 * c3 * c5^2 + 2963625 * c5^3)) + 
        90 * c2 * (270 * c3^4 + 16905 * c3^3 * c5 + 280 * c3^2 * (89 * c4^2 + 1580 * c5^2) + 
                     35 * c3 * (24832 * c4^2 * c5 + 162675 * c5^3) + 
                     4 * (17408 * c4^4 + 2097900 * c4^2 * c5^2 + 7546875 * c5^4)) + 
        27 * c4 * (14775 * c3^4 + 1028300 * c3^3 * c5 + 50 * c3^2 * (10144 * c4^2 + 594055 * c5^2) + 
                     700 * c3 * (27904 * c4^2 * c5 + 598575 * c5^3) + 
                     3 * (316928 * c4^4 + 68908000 * c4^2 * c5^2 + 806378125 * c5^4))
      )
    eq.B4 <- -gam4 + 120 * (
      32 * c2^6 + 3456 * c2^5 * c4 + 6 * c1^5 * c5 + 
        3 * c1^4 * (9 * c3^2 + 16 * c2 * c4 + 168 * c4^2 + 330 * c3 * c5 + 2850 * c5^2) + 
        720 * c2^4 * (15 * c3^2 + 224 * c4^2 + 490 * c3 * c5 + 4410 * c5^2) + 
        6048 * c2^3 * c4 * (125 * c3^2 + 704 * c4^2 + 4480 * c3 * c5 + 44100 * c5^2) + 
        12 * c1^3 * (4 * c2^2 * (3 * c3 + 50 * c5) + 60 * c2 * c4 * (7 * c3 + 114 * c5) + 
                       3 * (24 * c3^3 + 1192 * c3 * c4^2 + 1170 * c3^2 * c5 + 20440 * c4^2 * c5 + 
                              20150 * c3 * c5^2 + 124875 * c5^3)) + 
        216 * c2^2 * (945 * c3^4 + 67620 * c3^3 * c5 + 
                        560 * c3^2 * (178 * c4^2 + 3555 * c5^2) + 
                        315 * c3 * (12416 * c4^2 * c5 + 90375 * c5^3) + 
                        6 * (52224 * c4^4 + 6993000 * c4^2 * c5^2 + 27671875 * c5^4)) + 
        6 * c1^2 * (8 * c2^4 + 480 * c2^3 * c4 + 
                    180 * c2^2 * (4 * c3^2 + 64 * c4^2 + 125 * c3 * c5 + 1050 * c5^2) + 
                    72 * c2 * c4 * (327 * c3^2 + 1848 * c4^2 + 10955 * c3 * c5 + 99900 * c5^2) + 
                    9 * (225 * c3^4 + 22824 * c3^2 * c4^2 + 69632 * c4^4 + 15090 * c3^3 * c5 + 
                           830240 * c3 * c4^2 * c5 + 412925 * c3^2 * c5^2 + 
                           8239800 * c4^2 * c5^2 + 5475750 * c3 * c5^3 + 29636250 * c5^4)) + 
        1296 * c2 * c4 * (5910 * c3^4 + 462735 * c3^3 * c5 + 
                            c3^2 * (228240 * c4^2 + 14851375 * c5^2) + 
                            175 * c3 * (55808 * c4^2 * c5 + 1316865 * c5^3) + 
                            3 * (158464 * c4^4 + 37899400 * c4^2 * c5^2 + 483826875 * c5^4)) + 
        27 * (9945 * c3^6 + 92930048 * c4^6 + 1166130 * c3^5 * c5 + 
                35724729600 * c4^4 * c5^2 + 977816385000 * c4^2 * c5^4 + 
                1907724656250 * c5^6 + 180 * c3^4 * (16082 * c4^2 + 345905 * c5^2) + 
                140 * c3^3 * (1765608 * c4^2 * c5 + 13775375 * c5^3) + 
                15 * c3^2 * (4076032 * c4^4 + 574146160 * c4^2 * c5^2 + 
                               2424667875 * c5^4) + 
                210 * c3 * (13526272 * c4^4 * c5 + 687499200 * c4^2 * c5^3 + 
                              1876468125 * c5^5)) + 
        18 * c1 * (80 * c2^4 * (c3 + 15 * c5) + 160 * c2^3 * c4 * (32 * c3 + 525 * c5) + 
                     12 * c2^2 * (225 * c3^3 + 11088 * c3 * c4^2 + 11130 * c3^2 * c5 + 
                                    199360 * c4^2 * c5 + 200900 * c3 * c5^2 + 1323000 * c5^3) + 
                     24 * c2 * c4 * (3885 * c3^3 + 69632 * c3 * c4^2 + 209800 * c3^2 * c5 + 
                                       1369440 * c4^2 * c5 + 4134375 * c3 * c5^2 + 29636250 * c5^3) + 
                     9 * (540 * c3^5 + 48585 * c3^4 * c5 + 
                            20 * c3^3 * (4856 * c4^2 + 95655 * c5^2) + 
                            80 * c3^2 * (71597 * c4^2 * c5 + 513625 * c5^3) + 
                            4 * c3 * (237696 * c4^4 + 30726500 * c4^2 * c5^2 + 
                                        119844375 * c5^4) + 
                            5 * c5 * (4076032 * c4^4 + 191074800 * c4^2 * c5^2 + 
                                        483826875 * c5^4)))
    )
    
    eqs <- c(eq.18, eq.22, eq.B1, eq.B2, eq.B3, eq.B4)
    obj <- sum(eqs^2)
    
    obj
    
  }
  
  OPT <- list()
  ntry <- 0
  cnt <- 0
  while(ntry < control[["max.ntry"]]){
    
    ntry <- ntry + 1
    start <- rnorm(6, sd = .5)
    opt <- nlminb(start = start, objective = obj.fun, scale = 10, 
                  lower = -2, upper = 2, 
                  control = list(trace = F, abs.tol = 1e-20, rel.tol = 1e-15, eval.max = 1e6, iter.max = 1e6), gam = gam)
    #print(opt$objective)
    if(opt$convergence == 0 && opt$objective <= control[["obj.tol"]]){
      cnt <- cnt + 1
      OPT[[cnt]] <- opt
      if(control[["trace"]]){
        #cat(cnt, "/", ntry, "\n", sep="")
      }
    }

    
    if(length(OPT) >= control[["n.valid.sol"]] || (opt$objective <= min(1e-15, control[["obj.tol"]]) && opt$convergence == 0)){
      break
    }
  }
  
  if(length(OPT) == 0){
    return(NULL)
    stop(paste0("cannot find the coefficients of polynomial after ", control[["max.ntry"]], " attempts"))
  }
  
     }
  }
  if(control[["trace"]]){
    #cat("minimum objective: ", min.obj, "\n", sep="")
  }
  
  coeff <- OPT[[idx]]$par
  list(coeff = coeff, min.obj = min.obj)
  
  
}

headrick02.corr.match <- function(poly.coeff, corr){
  
  obj.fun2 <- function(x, c0, c1, c2, c3, c4, c5, i, j, rho.Y){
    
    eq <- -rho.Y + 
      3*c0[j]*c4[i] + 3*c2[j]*c4[i] + 9*c4[i]*c4[j] + c0[i]*(c0[j] + c2[j] + 3*c4[j]) + 
      c1[i]*c1[j]*x + 3*c1[j]*c3[i]*x + 3*c1[i]*c3[j]*x + 9*c3[i]*c3[j]*x + 
      15*c1[j]*c5[i]*x + 45*c3[j]*c5[i]*x + 15*c1[i]*c5[j]*x + 
      45*c3[i]*c5[j]*x + 225*c5[i]*c5[j]*x + 12*c2[j]*c4[i]*x^2 + 
      72*c4[i]*c4[j]*x^2 + 6*c3[i]*c3[j]*x^3 + 60*c3[j]*c5[i]*x^3 + 
      60*c3[i]*c5[j]*x^3 + 600*c5[i]*c5[j]*x^3 + 24*c4[i]*c4[j]*x^4 + 
      120*c5[i]*c5[j]*x^5 + 
      c2[i]*(c0[j] + c2[j] + 3*c4[j] + 2*c2[j]*x^2 + 12*c4[j]*x^2)
    obj <- eq^2
    obj
    
  }
  
  k <- ncol(poly.coeff)
  
  c0 <- as.vector(poly.coeff[1,], mode = "numeric")
  c1 <- as.vector(poly.coeff[2,], mode = "numeric")
  c2 <- as.vector(poly.coeff[3,], mode = "numeric")
  c3 <- as.vector(poly.coeff[4,], mode = "numeric")
  c4 <- as.vector(poly.coeff[5,], mode = "numeric")
  c5 <- as.vector(poly.coeff[6,], mode = "numeric")  
  
   
  l <- 0
  inter.corr <- diag(1, k)
  obj <- matrix(NA, k , k)
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      l <- l + 1
      rho.Y <- corr[i, j]
      opt <- nlminb(start = rho.Y, objective = obj.fun2, scale = 10, lower = -1, upper = 1, 
                    control = list(trace = F, abs.tol = 1e-20, eval.max = 1e5, iter.max = 1e3), 
                    c0 = c0, c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, i = i, j = j, rho.Y = rho.Y)
      if(opt$convergence == 0){
        inter.corr[i, j] <- opt$par
        inter.corr[j, i] <- opt$par
        obj[i, j] <- opt$objective
      }else{
        stop("error in solving intermediate correlation")
      }
      
    }
  }
  
  list(inter.corr = inter.corr, obj = obj)
  
}

headrick02 <- function(n, mean, sd, corr, skewness, kurtosis, gam3=NaN, gam4=NaN, replication = 1, control = list(seed = NULL, trace = T, max.ntry = 5, obj.tol = 1e-10, n.valid.sol = 1)){
  
  ##setting up
  
  if (! file.exists("compiled.txt")){
	print("Error: compiled.txt not found. Please change the working directory and try again.")
	return
  } 
  
   
  if(!is.null(control[["seed"]])){
    set.seed(control[["seed"]])
  }
  

  if(is.null(control[["trace"]])){
    control[["trace"]] <- T
  }
  
  if(is.null(control[["max.ntry"]])){
    control[["max.ntry"]] <- 5
  }
  
  if(is.null(control[["obj.tol"]])){
    control[["obj.tol"]] <- 1e-10
  }
  
  if(is.null(control[["n.valid.sol"]])){
    control[["n.valid.sol"]] <- 1
  }
  
  start = Sys.time()  
  
  k <- nrow(corr)

  
  if (is.nan(gam3[1]) && !is.nan(gam4[1])){
	print("Error: Please provide both gam3 and gam4, or neither.")
	return
  }
  
  if (is.nan(gam4[1]) && !is.nan(gam3[1])){
	print("Error: Please provide both gam3 and gam4, or neither.")
	return
  }
  
  default_gam3 = F
  default_gam4 = F
  
  if(is.nan(gam3[1])){
    gam3 = pmax(skewness, kurtosis)
	default_gam3 = T
  }
  
  if(is.nan(gam4[1])){
    gam4 = pmax(skewness,kurtosis)^2
	default_gam4 = T
  }
  
  len <- c(length(mean), length(sd), length(skewness), length(kurtosis), length(gam3), length(gam4))
  if(var(len) != 0){
    stop("Lengths of mean, std, skewness, kurtosis, gam3 and gam4 must be equal")
  }
  
  if(len[1] != 1 && len[1] != k){
    stop("Inconsistent length/dim of moments and correlation")
  }
  
  if(len[1] == 1){
    mean <- rep(mean, k)
    sd <- rep(sd, k)
    skewness <- rep(skewness, k)
    kurtosis <- rep(kurtosis, k)
    gam3 <- rep(gam3, k)
    gam4 <- rep(gam4, k)
  }
  
    for (i in 1:k){
	if (kurtosis[i] <= skewness[i]^2 - 2){
		cat("Error: the ", i," th component of kurtosis is not bigger than skewness squared minus 2.\n")
		return
	}
  }
  
  
  ##Solve for coefficients c0-c5 using equation 18, 22, B1-B4
  
  coeff <- NULL
  obj.poly.coeff <- NULL
  poly.coeff <- NULL
  gam4_fit = rep(0,k)
  gam3_fit = rep(0,k)
  
  for(i in 1:k){   
		if (control[["trace"]]){
			cat("Time elapsed ", as.numeric(Sys.time()-start, units="secs") ," seconds. Start fitting c0 - c5 for distribution ", i, ".\n", sep="")
		}	
	
		compiled = read.table("compiled.txt", header = T)
	
		if (default_gam3 && default_gam4){
			matched = compiled[compiled["g1"]==skewness[i] & compiled["g2"]==kurtosis[i] & compiled["tol"]<=control[["obj.tol"]],]
			if (nrow(matched)>0){
				if (control[["trace"]]){
					cat("Configuration found in compiled list. Compiled coefficients will be used. \n")
				}
				matched = matched[order(-matched$tol),]
				curr.coeff = c(as.vector(matched[1,c("c0", "c1", "c2", "c3", "c4", "c5")]))
				curr.obj = matched[1, "tol"]
				gam3_fit[i] = matched[1, "g3"]
				gam4_fit[i] = matched[1, "g4"]
			}
			else{
				j <- 1
				j3 <- 1
				poly.coeff <- NULL
				upper = 4
				step_size = 4
				iterations = 0
				while(j<=15 && j3<=15 && is.null(poly.coeff)){
					iterations = iterations + 1
					tic = Sys.time()
					gam3_fit[i] = gam3[i]/2*2^j3
					gam4_fit[i] = gam4[i]/2*2^j
					poly.coeff <- headrick02.poly.coeff(skewness[i], kurtosis[i], gam3_fit[i], gam4_fit[i], control = control)
					if(is.null(poly.coeff)){
						if (control[["trace"]]){
							cat("Trial ",iterations," unsuccessful. Time spent: ", as.numeric(Sys.time()-tic, units="secs") , " seconds.\n", sep = "")
						}
						j3<-j3+1						
						if(j3==upper+1 && j<upper){
							j3 <- 1
							j <- j+1
						}
						else if(j3==upper+1 && j==upper){
							input = "y"
							cat("No solutions found after ", iterations," iterations. Do you want to continue searching? (y/n)\n", sep = "")
							input = readline()
							while(input!="y" && input !="n" && input != "Y" && input != "N"){
								cat("Invalid input. Please try again.\n")
								input = readline()
							}
							if(input=="y" || input == "Y"){
								j3 = 1
								j = j + 1
								upper = upper + step_size								
							}
							if(input == "n" || input == "N"){
								break
							}
						}
					}
				}		
				if(!is.null(poly.coeff)){
					if (control[["trace"]]){
						cat("Trial ",iterations," successful. Time spent: ", as.numeric(Sys.time()-tic, units="secs") , " seconds.\n", sep = "")
					}
					curr.coeff = poly.coeff$coeff
					curr.obj = poly.coeff$min.obj
					to_append = matrix(c(skewness[i],kurtosis[i],curr.obj,gam3_fit[i], gam4_fit[i], curr.coeff), nrow=1)
					write.table(to_append, "compiled.txt",append = T, row.names = F, col.names = F)
				}
				
				if (is.null(poly.coeff)){
					cat("Error: no solution found for the combination of skewness: ", skewness[i], "; kurtosis: ",
						kurtosis[i], ".\n", sep = "")
					return(NULL)
				}
			}
		}
		
		else{
			matched = compiled[compiled["g1"]==skewness[i] & compiled["g2"]==kurtosis[i] 
				& compiled["g3"] == gam3[i] & compiled["g4"] == gam4[i]
				& compiled["tol"]<=control[["obj.tol"]],]
			if (nrow(matched)>0){
				if (control[["trace"]]){
					cat("Configuration found in compiled list. Compiled coefficients will be used. \n")
				}
				matched = matched[order(-matched$tol)]
				curr.coeff = c(as.vector(matched[1,c("c0", "c1", "c2", "c3", "c4", "c5")]))
				curr.obj = matched[1, "tol"]
				gam3_fit[i] = matched[1, "g3"]
				gam4_fit[i] = matched[1, "g4"]
			}
			else {
				gam3_fit[i]  = gam3[i]
				gam4_fit[i] = gam4[i]
				poly.coeff <- headrick02.poly.coeff(skewness[i], kurtosis[i], gam3_fit[i], gam4_fit[i], control = control)
				if(is.null(poly.coeff)){
					cat("Error: no solution found for the combination of skewness: ", skewness[i], "; kurtosis: ",
						kurtosis[i], " gam3: ", gam3_fit[i], "; gam4: ", gam4_fit[i], ".\n", sep = "")
						return
				}
				else{
					curr.coeff = poly.coeff$coeff
					curr.obj = poly.coeff$min.obj
					to_append = matrix(c(skewness[i],kurtosis[i],curr.obj,gam3_fit[i], gam4_fit[i], curr.coeff), nrow=1)
					write.table(to_append, "compiled.txt",append = T, row.names = F, col.names = F)
				}
			}
		}
		
		coeff <- c(coeff, curr.coeff)
		obj.poly.coeff <- c(obj.poly.coeff, curr.obj)
	}
	
	coeff = matrix(coeff, nrow = 6)
	
	desired.moments <-data.frame(mean = mean, sd=sd, skewness = skewness, kurtosis = kurtosis, gam3 = gam3_fit, gam4=gam4_fit)
	rownames(desired.moments) <- paste0("Y", 1:nrow(desired.moments))

	if (control[["trace"]]){
			cat("Finished fitting c0 - c5. Time elapsed ", as.numeric(Sys.time()-start, units="secs") ," seconds. \n", sep="")
		}
  
	  summary.poly.coeff <- rbind(obj.poly.coeff, coeff)
	  colnames(summary.poly.coeff) <- paste0("Distribution ", 1:ncol(summary.poly.coeff))
	  rownames(summary.poly.coeff) <- c("obj value @ convergence", paste0("c", 0:5))
	  colnames(coeff) <- paste0("Distribution ", 1:ncol(summary.poly.coeff))
	  rownames(coeff) <- paste0("c", 0:5)
	  
	  ##Solve for intermediate correlation using equation 26
	  if(k>1){
	  if (control[["trace"]]){
			cat("\nStart solving for intermediate correlation matrix...\n")
	    }
	  
	  corr.match <- headrick02.corr.match(coeff, corr)
	  inter.corr <- corr.match$inter.corr
	  obj.corr.match <- corr.match$obj
	  colnames(inter.corr) <- paste0("Z", 1:ncol(inter.corr))
	  rownames(inter.corr) <- paste0("Z", 1:nrow(inter.corr))
	  
	  colnames(obj.corr.match) <- paste0("Z", 1:ncol(obj.corr.match))
	  rownames(obj.corr.match) <- paste0("Z", 1:nrow(obj.corr.match))
	  
	  if (control[["trace"]]){
			cat("Finished solving for intermediate correlation matrix. Time elapsed ", as.numeric(Sys.time()-start, units="secs") ," seconds.\n", sep="")
	    }
	  }
	  else{
	  inter.corr <- corr
	  colnames(inter.corr) <- paste0("Z", 1:ncol(inter.corr))
	  rownames(inter.corr) <- paste0("Z", 1:nrow(inter.corr))
	 
	  }
	  
		  c0 <- as.vector(coeff[1,], mode = "numeric")
		  c1 <- as.vector(coeff[2,], mode = "numeric")
		  c2 <- as.vector(coeff[3,], mode = "numeric")
		  c3 <- as.vector(coeff[4,], mode = "numeric")
		  c4 <- as.vector(coeff[5,], mode = "numeric")
		  c5 <- as.vector(coeff[6,], mode = "numeric")
		  library("MASS")
		  
		  obs.mean = NULL
		  obs.sd = NULL
		  obs.skew = NULL
		  obs.kurt = NULL
		  obs.gam3 = NULL
		  obs.gam4 = NULL
		  
	  for (replica in 1:replication){
		  ## Generate intermediate normal distribution with desired intermediate correlation	

		  Z <- mvrnorm(n, mu = rep(0, k), Sigma = inter.corr)
		  Z2 <- Z^2
		  Z3 <- Z^3
		  Z4 <- Z^4
		  Z5 <- Z^5	  
		  
		  ## Generate multivariate distribution with desired property
		  
		  Y <- matrix(0, nrow = n, ncol = k)
		  for(i in 1:k){
				Y[, i] <- mean[i] + sd[i]*(c0[i] + c1[i] * Z[, i] + c2[i] * Z2[, i] + c3[i] * Z3[, i] + c4[i] * Z4[, i] + c5[i] * Z5[, i])
			}
		  
		  obs.mean = rbind(obs.mean, apply(Y, 2, mean))
		  obs.sd = rbind(obs.sd, apply(Y, 2, sd))
		  obs.skew = rbind(obs.skew, apply(Y, 2, f_skew))
		  obs.kurt = rbind(obs.kurt, apply(Y, 2, f_kurt))
		  obs.gam3 = rbind(obs.gam3, apply(Y, 2, f_gamma3))
		  obs.gam4 = rbind(obs.gam4, apply(Y, 2, f_gamma4))
		}  
		  
		  
		obs.moments <- data.frame(mean = apply(obs.mean, 2, mean), sd = apply(obs.sd, 2, mean), skewness = apply(obs.skew, 2, mean), 
						kurtosis = apply(obs.kurt, 2, mean), gam3 = apply(obs.gam3,2,mean), gam4 = apply(obs.gam4,2,mean))
		obs.moments.sd <- data.frame(mean = apply(obs.mean, 2, sd), sd = apply(obs.sd, 2, sd), skewness = apply(obs.skew, 2, sd), 
						kurtosis = apply(obs.kurt, 2, sd), gam3 = apply(obs.gam3,2,sd), gam4 = apply(obs.gam4,2,sd))
		rownames(obs.moments) <- paste0("Y", 1:nrow(obs.moments))
		rownames(obs.moments.sd) <- paste0("Y", 1:nrow(obs.moments.sd))
		
		obs.corr <- cor(Y)
		rownames(obs.corr) <- paste0("Y", 1:nrow(obs.corr))
		colnames(obs.corr) <- paste0("Y ", 1:ncol(obs.corr))
		
		if (replication>1){
			obs.corr = NULL
		}
		
	  if (control[["trace"]]){
			cat("\nDesired moments:\n")
			print(desired.moments)
	  
			cat("\nSampling moments:\n")
			print(obs.moments)
			
			if (replication >1){
				cat("\nSampling moment standard deviations:\n")
				print(obs.moments.sd)
			}
	  
			if(replication == 1){
				cat("\nDesired correlation matrix:\n")
				print(corr)
	  
				cat("\nSampling correlation matrix:\n")
				print(obs.corr)
			}
			cat("\nTotal time elapsed ", as.numeric(Sys.time()-start, units="secs"), " seconds.\n", sep="")
		}
		
	  list(obs = Y, obs.corr = obs.corr, obs.moments = obs.moments, obs.moments.sd = obs.moments.sd, desired.corr = corr, desired.moments = desired.moments,
		   summary.poly.coeff = summary.poly.coeff, inter.corr = inter.corr)
	  
}






