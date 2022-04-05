library(foreach)
library(doParallel)

#####################################################################
# weighted sum of normal distributions for uniroot
#####################################################################
wdf <- function(x, w=w, m=m, s=s, a=0 ) {
  return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
}

#####################################################################
# 
#####################################################################

# non-singlet
Z_QQ1  <- c(1.151, 0.004 )
Z_QQ2  <- c(1.160, 0.003 )

# singlet
Z_QQs1 <- c(1.161, 0.024 )
Z_QQs2 <- c(1.163, 0.012 )

# gluon
Z_GG2 <- c(1.08, 0.17 )

# mixing coefficients

Z_GQ1 <- c(  0.232, 0 )
Z_GQ2 <- c(  0.083, 0 )
Z_QG1 <- c( -0.027, 0 )
# Z_QG2 <- c(  0.0  , 0 )
Z_QG2 <- Z_QG1

################################################################
#
################################################################
aic_ren <- function( ns = 600, nprobe=2^16 , singlet=F , nonsinglet=F, stats=F) {

  ensemble <- "cB211.072.64"
  lvl <- 0

  nstout <- 10
  rstout <- 0.1290
  chisq_red_min = 0.9
  chisq_red_max = 1.1


#  sZ_QQ1  <- rnorm( n = ns, mean = Z_QQ1[1],  sd = Z_QQ1[2]  )
#  sZ_QQs1 <- rnorm( n = ns, mean = Z_QQs1[1], sd = Z_QQs1[2] )
#  sZ_GQ1  <- rnorm( n = ns, mean = Z_GQ1[1],  sd = Z_GQ1[2]  )
#  sZ_QG1  <- rnorm( n = ns, mean = Z_QG1[1],  sd = Z_QG1[2]  )

  sZ_QQ2  <- rnorm( n = ns, mean = Z_QQ2[1],  sd = Z_QQ2[2]  )

  sZ_QQs2 <- rnorm( n = ns, mean = Z_QQs2[1], sd = Z_QQs2[2] )

  sZ_GG2  <- rnorm( n = ns, mean = Z_GG2[1],  sd = Z_GG2[2]  )

  sZ_GQ2  <- rnorm( n = ns, mean = Z_GQ2[1],  sd = Z_GQ2[2]  )

  sZ_QG2  <- rnorm( n = ns, mean = Z_QG2[1],  sd = Z_QG2[2]  )


  #####################################################
  # read aic data
  #####################################################
  read_aic <- function(q) {
    if ( is.na ( q$nstout ) ) {
      if ( !is.na(q$op ) ) {
        op_tag <- paste ( "-", q$op, sep="" )
      } else {
        op_tag <- "" 
      }
      f <- paste ( q$name, op_tag, "/lvl", q$lvl, "/fit.", ensemble, ".", q$name, ".aic", sep="" )
    } else {
      f <- paste ( q$name, "/nstout", q$nstout, "/lvl", q$lvl, "/fit.", ensemble, ".", q$name, "-", q$op, ".nstout", q$nstout, "_", formatC(q$rstout, width=6, digits=4, format="f" ) , ".aic", sep="" )
    }

    if ( !file.exists ( f ) ) stop ( "[read_aic] no file ", f )
    else message ( "# [read_aic] reading from file ", f )
    q$f <- read.table ( f ) 
    return ( q )
  }

  #####################################################
  # filter by chisq / dof
  #####################################################
  cr_filter <- function(q, col_cr=6, col_dof=9 ) {
    idx <- which( ( q$f[,col_cr] / q$f[,col_dof] >= q$cr[1] ) & ( q$f[,col_cr] / q$f[,col_dof] <= q$cr[2]  ) )
    q$f <- q$f[idx,]
    q$n <- length ( q$f$V1 )
    message( "# ", q$name, "  n = ", q$n )
    rm(idx)
    return(q)
  }

  #####################################################
  # start filling list
  #####################################################

  q_c_ll <- list()
  q_c_ll$name <- "xq-conn-kaon-ll"
  q_c_ll$op <- "4k"
  q_c_ll$cr <- c( chisq_red_min, chisq_red_max )
  q_c_ll$lvl <- 0
  q_c_ll$nstout <- NA
  q_c_ll$f <- NULL

  q_c_ll <- read_aic (q = q_c_ll )
  q_c_ll <- cr_filter ( q = q_c_ll )

  q_c_ss <- list()
  q_c_ss$name <- "xq-conn-kaon-ss"
  q_c_ss$op <- "4k"
  q_c_ss$cr <- c( chisq_red_min, chisq_red_max )
  q_c_ss$lvl <- 0
  q_c_ss$nstout <- NA
  
  q_c_ss <- read_aic (q = q_c_ss )
  q_c_ss <- cr_filter ( q = q_c_ss )

  q_d_ll <- list()
  q_d_ll$name <- "xq-disc-kaon-ll"
  q_d_ll$op <- NA
  q_d_ll$cr <- c( chisq_red_min, chisq_red_max )
  q_d_ll$lvl <- 0
  q_d_ll$nstout <- NA
 
  q_d_ll <- read_aic (q = q_d_ll )
  q_d_ll <- cr_filter ( q = q_d_ll )

  q_d_ss <- list()
  q_d_ss$name <- "xq-disc-kaon-ss"
  q_d_ss$op <- NA
  q_d_ss$cr <- c( chisq_red_min, chisq_red_max )
  q_d_ss$lvl <- 0
  q_d_ss$nstout <- NA
 
  q_d_ss <- read_aic (q = q_d_ss )
  q_d_ss <- cr_filter ( q = q_d_ss )

  q_d_cc <- list()
  q_d_cc$name <- "xq-disc-kaon-cc"
  q_d_cc$op <- NA
  q_d_cc$cr <- c( chisq_red_min, chisq_red_max )
  q_d_cc$lvl <- 0
  q_d_cc$nstout <- NA
 
  q_d_cc <- read_aic (q = q_d_cc )
  q_d_cc <- cr_filter ( q = q_d_cc )

  g_d    <- list()
  g_d$name <- "xg-disc-kaon"
  g_d$op <- "clover"
  g_d$cr <- c( chisq_red_min, chisq_red_max )
  g_d$lvl <- 0
  g_d$nstout <- nstout
  g_d$rstout <- 0.1290

  g_d <- read_aic (q = g_d )
  g_d <- cr_filter ( q = g_d )

  #####################################################
  # read selected bootstrap data
  #####################################################
  read_bs <- function (q) {
    q$data <- array ( dim=c( q$n, ns) )

    for( i in 1:q$n ) {
      lvl <- q$f$V1[i]
      tf  <- c( q$f$V2[i],  q$f$V3[i] )
      tc  <- q$f$V4[i]
      ts  <- as.character(q$f$V5)[i]

      if ( is.na ( q$nstout ) ) {
        if ( !is.na(q$op ) ) {
          op_tag <- paste ( "-", q$op, sep="" )
        } else {
          op_tag <- "" 
        }

        f <- paste ( q$name, op_tag, "/lvl", q$lvl, "/fit.", ensemble, ".", q$name, ".", ts, ".tc", tc, ".tf", tf[1], "_", tf[2], ".bs", sep="" )

      } else {

        f <- paste ( g_d$name, "/nstout", q$nstout, "/lvl", q$lvl, "/fit.", ensemble, ".", q$name, "-", g_d$op, ".nstout", nstout, "_", formatC(rstout, width=6, digits=4, format="f"), ".", ts, ".tc", tc, ".tf", tf[1], "_", tf[2], ".bs", sep="" )
      }

      if ( !file.exists(f) ) {
        stop ( "[read_bs] cannot find file ", f )
      } else {
        message( "# [read_bs] reading from file ", f)
      }
      q$data[i,] <- read.table ( f )$V3
    }
    return ( q )
  }

  #####################################################
  # apply read_bs
  #####################################################
  q_c_ll <- read_bs (q = q_c_ll )
  q_c_ss <- read_bs (q = q_c_ss )
  q_d_ll <- read_bs (q = q_d_ll )
  q_d_ss <- read_bs (q = q_d_ss )
  q_d_cc <- read_bs (q = q_d_cc )
  g_d    <- read_bs (q = g_d )
  
  #####################################################
  # combined AIC analysis for sum of flavor, singlet
  #   and gluon
  #
  # SINGLET CASE
  #
  #####################################################
  if ( singlet ) {

    w <- list()
    # w$res <- array ( dim=c( nprobe, 5 ) )
  
    w$n <-   q_c_ll$n * q_c_ss$n * q_d_ll$n * q_d_ss$n * q_d_cc$n * g_d$n
    LL <- c(            q_c_ss$n * q_d_ll$n * q_d_ss$n * q_d_cc$n * g_d$n,
                                   q_d_ll$n * q_d_ss$n * q_d_cc$n * g_d$n,
                                              q_d_ss$n * q_d_cc$n * g_d$n,
                                                         q_d_cc$n * g_d$n,
                                                                    g_d$n )
    cat("# [aic_ren] w n = ", w$n, "\n",
        "# [aic_ren] LL = ", formatC(LL, width=6, format="d"), "\n" , sep="")
  
    start_time <- Sys.time()

    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    
    end_time <- Sys.time()
    message ( "# [aic_ren] time for idx = ", end_time - start_time )
  
    # return(idx)

    start_time <- Sys.time()

    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
  
    # for ( k in 1:nprobe ) 
    w$res <- foreach ( k = 1:nprobe , .combine=rbind ) %dopar% 
    {
  
      i1 <- ( (idx[k]-1)          ) %/% LL[1] + 1
      i2 <- ( (idx[k]-1) %% LL[1] ) %/% LL[2] + 1 
      i3 <- ( (idx[k]-1) %% LL[2] ) %/% LL[3] + 1
      i4 <- ( (idx[k]-1) %% LL[3] ) %/% LL[4] + 1
      i5 <- ( (idx[k]-1) %% LL[4] ) %/% LL[5] + 1
      i6 <- ( (idx[k]-1) %% LL[5] ) + 1
  
      # message ( " idx ", idx[k] , " ---> ", formatC( c(i1, i2, i3, i4, i5, i6), width=4, format="d" ) )
  
      # per quark flavor and gluon
      q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
      q_d <-                    q_d_ll$data[i3,]
      
      q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]
      
      q_c <-                    q_d_cc$data[i5,]
  
      g   <-                                      g_d$data[i6,]
  
      # sum of quark flavors
      q_f <- q_u + q_d + q_s + q_c
  
      # renormalized combinations
      q_r <- sZ_QQs2 * q_f + sZ_QG2 * g
  
      g_r <- sZ_GG2  * g   + sZ_GQ2 * q_f
  
      qg_r <- q_r + g_r
  
      sZ_QQd2 <- sZ_QQs2 - sZ_QQ2
  
      qg_aux <- ( sZ_QQd2 * q_f + sZ_QG2 * g ) / 4.
  
      q_u_r <- sZ_QQ2 * q_u + qg_aux
  
      q_d_r <- sZ_QQ2 * q_d + qg_aux
  
      q_s_r <- sZ_QQ2 * q_s + qg_aux
  
      q_c_r <- sZ_QQ2 * q_c + qg_aux
  
  
      a <- q_c_ll$f$V10[i1] * q_c_ss$f$V10[i2] * q_d_ll$f$V10[i3] * q_d_ss$f$V10[i4] * q_d_cc$f$V10[i5] * g_d$f$V10[i6]
  
      # w$res[k,] <- c( mean(q_r), sqrt(var(q_r)), mean(g_r), sqrt(var(g_r)), a )
      c( a, 
         mean(q_u_r), sqrt(var(q_u_r) ),
         mean(q_d_r), sqrt(var(q_d_r) ),
         mean(q_s_r), sqrt(var(q_s_r) ),
         mean(q_c_r), sqrt(var(q_c_r) ),
         mean(q_r),   sqrt(var(q_r)   ), 
         mean(g_r),   sqrt(var(g_r)   ),
         mean(qg_r),  sqrt(var(qg_r)  ) 
      )
    }
  
    #stop cluster
    stopCluster(cl)
 
    end_time <- Sys.time()
    message( "# [aic_red] time for aic bs combination = ", end_time - start_time )

    start_time <- Sys.time()
  
    output_filename <- paste( "x-kaon.singlet.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
    saveRDS( object = w$res, file = output_filename )
              
    end_time <- Sys.time()
    message( "# [aic_red] time for saveRDS = ", end_time - start_time )

    # return(w)
  
    if ( stats ) {
  
      start_time <- Sys.time()
     
      cores=detectCores()
      cl <- makeCluster(cores[1]-1)
      registerDoParallel(cl)
    
      bres <- foreach ( b = 1:7, .combine = rbind ) %dopar% {
    
        wdf <- function(x, w=w, m=m, s=s, a=0 ) {
          # return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
  
          return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
        }
    
        lower <- c( -0.5, -0.5, -0.5 )
        upper <- c( +1.5, +1.5, +1.5 )
          
        ba  <- w$res[,1]
        bm  <- w$res[,(2*b  )]
        bs  <- w$res[,(2*b+1)]
    
        res <- array( dim=c( log2(nprobe) - 3 ,3 ) )
    
        for ( i in 1:(log2(nprobe)-3) ) {
          k <- 2^(i+3)
    
          a   <- ba[1:k] / sum ( ba[1:k] )
    
          m  <- bm[1:k]
          s  <- bs[1:k]
    
          u_m <- uniroot( wdf, lower=lower[1], upper=upper[1] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.50, s=s, w=a, m=m )
          u_l <- uniroot( wdf, lower=lower[2], upper=upper[2] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.16, s=s, w=a, m=m )
          u_u <- uniroot( wdf, lower=lower[3], upper=upper[3] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.84, s=s, w=a, m=m )
    
          res[i, ] <- c( u_m$root, u_l$root , u_u$root )
    
    ##      cat ( "# g convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_filename, append=TRUE )
    
          ### res[i,] <- c( 2*as.numeric(b), 2*as.numeric(i)+1, as.numeric(k) )
        }
    
        res
    
      }  # end of loop on observables
      stopCluster(cl)
    
      end_time <- Sys.time()
      message( "# [aic_red] time for quantiles = ", end_time - start_time )
  
      # return (bres)
    
      bnames <- c( "u", "d", "s", "c", "q", "g", "qg" )
    
      output_filename <- paste( "x-kaon.singlet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
      cat ( "# ", date(), "\n", file = output_filename, append=F )
    
      bi <- 0
      for ( b in 1:7 ) {
        for ( i in 1:(log2(nprobe)-3) ) {
          bi <- bi + 1 
          k <- 2^(i+3)
          cat( formatC( bnames[b] , width=4, format="s" ),
               formatC( k , width=10, format="d"),
               formatC( bres[bi,1], width=8, digits=4, format="f" ),
               formatC( bres[bi,1] - bres[bi,2], width=6, digits=4, format="f" ),
               formatC( bres[bi,3] - bres[bi,1], width=6, digits=4, format="f" ),
               "\n", file=output_filename, append=T )
        }
      }
      
    }  # end of if stats

  }  # end of singlet

  #####################################################
  # combined AIC analysis for sum of flavor non-singlet
  #
  # NON-SINGLET CASE
  #
  #####################################################
  if ( nonsinglet ) {

    #####################################################
    #
    # u, d, s, c
    #
    # 4 flavor
    #
    #####################################################


    w <- list()
  
    w$n <-   q_c_ll$n * q_c_ss$n * q_d_ll$n * q_d_ss$n * q_d_cc$n

    LL <- c(            q_c_ss$n * q_d_ll$n * q_d_ss$n * q_d_cc$n, 
                                   q_d_ll$n * q_d_ss$n * q_d_cc$n, 
                                              q_d_ss$n * q_d_cc$n, 
                                                         q_d_cc$n )
    cat("# [aic_ren] w n = ", w$n, "\n",
        "# [aic_ren] LL = ", formatC(LL, width=6, format="d"), "\n" , sep="")
  
    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    message ( "# [aic_ren] idx done " )
  
    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
  
    # for ( k in 1:nprobe ) 
    w$res <- foreach ( k = 1:nprobe , .combine=rbind ) %dopar% 
    {
  
      i1 <- ( (idx[k]-1)          ) %/% LL[1] + 1
      i2 <- ( (idx[k]-1) %% LL[1] ) %/% LL[2] + 1 
      i3 <- ( (idx[k]-1) %% LL[2] ) %/% LL[3] + 1
      i4 <- ( (idx[k]-1) %% LL[3] ) %/% LL[4] + 1
      i5 <- ( (idx[k]-1) %% LL[4] ) + 1
  
      q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
      q_d <-                    q_d_ll$data[i3,]
      
      q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]
      
      q_c <-                    q_d_cc$data[i5,]
  
      # sum of quark flavors
      q_f <- q_u + q_d + q_s - 3 * q_c
  
      # renormalized combination
      q_r <- sZ_QQ2 * q_f
  
      a <- q_c_ll$f$V10[i1] * q_c_ss$f$V10[i2] * q_d_ll$f$V10[i3] * q_d_ss$f$V10[i4] * q_d_cc$f$V10[i5]
  
      # w$res[k,] <- c( mean(q_r), sqrt(var(q_r)), mean(g_r), sqrt(var(g_r)), a )
      c( a, mean(q_r),   sqrt(var(q_r)) )
    }
  
    #stop cluster
    stopCluster(cl)
  
    output_filename <- paste( "x-kaon.nonsinglet.q4.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
#    saveRDS( object = w$res, file = output_filename )
              
    # return(w)
  
    #cores=detectCores()
    #cl <- makeCluster(cores[1]-1)
    #registerDoParallel(cl)
  
    bres <- foreach ( b = 1, .combine = rbind ) # %dopar% 
    {
  
      wdf <- function(x, w=w, m=m, s=s, a=0 ) {
        return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
      }
  
      lower <- c( -0.5, -0.5, -0.5 )
      upper <- c( +1.5, +1.5, +1.5 )
        
      ba  <- w$res[,1]
      bm  <- w$res[,(2*b  )]
      bs  <- w$res[,(2*b+1)]
  
      res <- array( dim=c( log2(nprobe) - 3 ,3 ) )
  
      for ( i in 1:(log2(nprobe)-3) ) {
        k <- 2^(i+3)
  
        a   <- ba[1:k] / sum ( ba[1:k] )
  
        m  <- bm[1:k]
        s  <- bs[1:k]
  
        u_m <- uniroot( wdf, lower=lower[1], upper=upper[1] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.50, s=s, w=a, m=m )
        u_l <- uniroot( wdf, lower=lower[2], upper=upper[2] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.16, s=s, w=a, m=m )
        u_u <- uniroot( wdf, lower=lower[3], upper=upper[3] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.84, s=s, w=a, m=m )
  
        res[i, ] <- c( u_m$root, u_l$root , u_u$root )
  
  ##      cat ( "# g convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_filename, append=TRUE )
  
      }
  
      res
  
    }  # end of loop on observables
    # stopCluster(cl)
  
    # return (bres)
  
    bnames <- c( "q4" )
  
    output_filename <- paste( "x-kaon.nonsinglet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
    cat ( "# ", date(), "\n", file = output_filename, append=F )
  
    bi <- 0
    for ( b in 1 ) {
      for ( i in 1:(log2(nprobe)-3) ) {
        bi <- bi + 1 
        k <- 2^(i+3)
        cat( formatC( bnames[b] , width=4, format="s" ),
             formatC( k , width=10, format="d"),
             formatC( bres[bi,1], width=8, digits=4, format="f" ),
             formatC( bres[bi,1] - bres[bi,2], width=6, digits=4, format="f" ),
             formatC( bres[bi,3] - bres[bi,1], width=6, digits=4, format="f" ),
             "\n", file=output_filename, append=T )
      }
    }

    rm(w)

    #####################################################
    #
    # u, d, s
    #
    # 3 flavor
    #
    #####################################################

    w <- list()
  
    w$n <-   q_c_ll$n * q_c_ss$n * q_d_ll$n * q_d_ss$n

    LL <- c(            q_c_ss$n * q_d_ll$n * q_d_ss$n,
                                   q_d_ll$n * q_d_ss$n, 
                                              q_d_ss$n )

    cat("# [aic_ren] w n = ", w$n, "\n",
        "# [aic_ren] LL = ", formatC(LL, width=6, format="d"), "\n" , sep="")
  
    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    message ( "# [aic_ren] idx done " )
  
    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
  
    # for ( k in 1:nprobe ) 
    w$res <- foreach ( k = 1:nprobe , .combine=rbind ) %dopar% 
    {
  
      i1 <- ( (idx[k]-1)          ) %/% LL[1] + 1
      i2 <- ( (idx[k]-1) %% LL[1] ) %/% LL[2] + 1 
      i3 <- ( (idx[k]-1) %% LL[2] ) %/% LL[3] + 1
      i4 <- ( (idx[k]-1) %% LL[3] ) %/% + 1
  
      q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
      q_d <-                    q_d_ll$data[i3,]
      
      q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]
      
      # sum of quark flavors
      q_f <- q_u + q_d - 2 * q_s
  
      # renormalized combination
      q_r <- sZ_QQ2 * q_f
  
      a <- q_c_ll$f$V10[i1] * q_c_ss$f$V10[i2] * q_d_ll$f$V10[i3] * q_d_ss$f$V10[i4]
  
      c( a, mean(q_r),   sqrt(var(q_r)) )
    }
  
    stopCluster(cl)
  
    output_filename <- paste( "x-kaon.nonsinglet.q3.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
#    saveRDS( object = w$res, file = output_filename )
              
    # return(w)
  
    #cores=detectCores()
    #cl <- makeCluster(cores[1]-1)
    #registerDoParallel(cl)
  
    bres <- foreach ( b = 1, .combine = rbind ) # %dopar% 
    {
  
      wdf <- function(x, w=w, m=m, s=s, a=0 ) {
        return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
      }
  
      lower <- c( -0.5, -0.5, -0.5 )
      upper <- c( +1.5, +1.5, +1.5 )
        
      ba  <- w$res[,1]
      bm  <- w$res[,(2*b  )]
      bs  <- w$res[,(2*b+1)]
  
      res <- array( dim=c( log2(nprobe) - 3 ,3 ) )
  
      for ( i in 1:(log2(nprobe)-3) ) {
        k <- 2^(i+3)
  
        a   <- ba[1:k] / sum ( ba[1:k] )
  
        m  <- bm[1:k]
        s  <- bs[1:k]
  
        u_m <- uniroot( wdf, lower=lower[1], upper=upper[1] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.50, s=s, w=a, m=m )
        u_l <- uniroot( wdf, lower=lower[2], upper=upper[2] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.16, s=s, w=a, m=m )
        u_u <- uniroot( wdf, lower=lower[3], upper=upper[3] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.84, s=s, w=a, m=m )
  
        res[i, ] <- c( u_m$root, u_l$root , u_u$root )
  
  ##      cat ( "# g convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_filename, append=TRUE )
  
      }
  
      res
  
    }  # end of loop on observables
    # stopCluster(cl)
  
    bnames <- c( "q3" )
  
    output_filename <- paste( "x-kaon.nonsinglet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
  
    bi <- 0
    for ( b in 1 ) {
      for ( i in 1:(log2(nprobe)-3) ) {
        bi <- bi + 1 
        k <- 2^(i+3)
        cat( formatC( bnames[b] , width=4, format="s" ),
             formatC( k , width=10, format="d"),
             formatC( bres[bi,1], width=8, digits=4, format="f" ),
             formatC( bres[bi,1] - bres[bi,2], width=6, digits=4, format="f" ),
             formatC( bres[bi,3] - bres[bi,1], width=6, digits=4, format="f" ),
             "\n", file=output_filename, append=T )
      }
    }
  
  }  # end of nonsinglet

  return ( NULL )

}  # end of aic_ren
