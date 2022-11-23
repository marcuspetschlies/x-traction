#####################################################################
#####################################################################
##
## CLUSTER / SOCKET PARALLELIZED VERSION
##
#####################################################################
#####################################################################

library(foreach)
library(doParallel)

source("/home/marcuspe/software/x-traction/zfactors.R" )

#####################################################################
# weighted sum of normal distributions for uniroot
#####################################################################
wdf <- function(x, w=w, m=m, s=s, a=0 ) {
  return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
}

################################################################
#
################################################################
aic_ren_cl <- function( ns = 600, q_c_pvec=c(0,0,1), q_c_op="4k", nprobe=2^16 , singlet=F , nonsinglet=F, with_quantiles =F, ens, field="kaon", chisq_red_range = c( 0.9, 1.1) ) {

  if ( missing ( ens   ) ) stop ( "need ensemble name" )

  # if ( missing ( field ) ) stop ( "need field name" )

  ensemble_full <- list()
  ensemble_full[["cB64"]] <- "cB211.072.64"
  ensemble_full[["cC80"]] <- "cC211.06.80"
  ensemble_full[["cD96"]] <- "cD211.054.96"


  ensemble <- ensemble_full[[ens]]
  lvl <- 0

  nstout <- 10
  rstout <- 0.1290
  chisq_red_min = chisq_red_range[1]
  chisq_red_max = chisq_red_range[2]

  Z <- init_z ( ens )

#  sZ_QQ1  <- rnorm( n = ns, mean = Z_QQ1[1],  sd = Z_QQ1[2]  )
#  sZ_QQs1 <- rnorm( n = ns, mean = Z_QQs1[1], sd = Z_QQs1[2] )
#  sZ_GQ1  <- rnorm( n = ns, mean = Z_GQ1[1],  sd = Z_GQ1[2]  )
#  sZ_QG1  <- rnorm( n = ns, mean = Z_QG1[1],  sd = Z_QG1[2]  )

  sZ_QQ2  <- rnorm( n = ns, mean = Z$QQ2[1],  sd = Z$QQ2[2]  )

  sZ_QQs2 <- rnorm( n = ns, mean = Z$QQs2[1], sd = Z$QQs2[2] )

  sZ_GG2  <- rnorm( n = ns, mean = Z$GG2[1],  sd = Z$GG2[2]  )

  sZ_GQ2  <- rnorm( n = ns, mean = Z$GQ2[1],  sd = Z$GQ2[2]  )

  sZ_QG2  <- rnorm( n = ns, mean = Z$QG2[1],  sd = Z$QG2[2]  )

  #####################################################
  # read aic data
  #
  read_aic <- function(q)
  #####################################################
  {
    if ( is.na ( q$nstout ) ) {
      if ( !is.na(q$op ) ) {
        op_tag <- paste ( "-", q$op, sep="" )
      } else {
        op_tag <- "" 
      }
      ps_tag <- ""
      if ( !is.na(q$p ) ) {
        ps_tag <- paste ( q$p, "-s", ns, sep="" )
      } else {
        ps_tag <- paste( "s" , ns, sep="" )
      }
      f <- paste ( q$name, op_tag, "/lvl", q$lvl, "/", ps_tag, "/fit.", ensemble, ".", q$name, ".aic", sep="" )
    } else {
      f <- paste ( q$name, "/nstout", q$nstout, "/lvl", q$lvl, "/s", ns, "/fit.", ensemble, ".", q$name, "-", q$op, ".nstout", q$nstout, "_", formatC(q$rstout, width=6, digits=4, format="f" ) , ".aic", sep="" )
    }

    if ( !file.exists ( f ) ) stop ( "[read_aic] no file ", f )
    else message ( "# [read_aic] reading from file ", f )
    q$f <- read.table ( f ) 
    return ( q )
  }

  #####################################################
  #####################################################

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
  # filling list for connected ll and ss
  #####################################################

  q_c_ll <- list()
  q_c_ss <- list()
  for ( f in c( "ll", "ss" ) ) {
    x.name <-  paste( "q_c_", f, sep="" )
    message ( "x = ", x.name )
    x.val <- list()
    x.val$name <- paste( "xq-conn-", field, "-", f , sep="" )
    x.val$op <- q_c_op # "4k"
    x.val$cr <- c( chisq_red_min, chisq_red_max )
    x.val$lvl <- 0
    x.val$nstout <- NA
    x.val$f <- NULL
    x.val$p <- paste( "p", q_c_pvec[1], q_c_pvec[2], q_c_pvec[3], sep="" )
    x.val$norm <- 0.5  # multiply factor 0.5 to get the single flavor contribution

    x.val <- read_aic (q = x.val )
    x.val <- cr_filter ( q = x.val )
    assign( x.name, x.val )

#    eval(as.name(x)) <- read_aic (q = v )
#    eval(as.name(x)) <- cr_filter ( q = v )

#    eval(as.name(x))[["name"]] <- paste( "xq-conn-", field, "-", f , sep="" )
#    eval(as.name(x) )$op <- "4k"
#    eval(as.name(x) )$cr <- c( chisq_red_min, chisq_red_max )
#    eval(as.name(x) )$lvl <- 0
#    eval(as.name(x) )$nstout <- NA
#    eval(as.name(x) )$f <- NULL
#    eval(as.name(x) )$p <- paste( "p", p[1], p[2], p[3], sep="" )

#    v <- eval(as.name(x) )
#    eval(as.name(x)) <- read_aic (q = v )
#    eval(as.name(x)) <- cr_filter ( q = v )
  }  # end of loop on flavor combinations

#  return ( list(q_c_ll=q_c_ll, q_c_ss = q_c_ss ) )


#  q_c_ss$name <- "xq-conn-kaon-ss"
#  q_c_ss$op <- "4k"
#  q_c_ss$cr <- c( chisq_red_min, chisq_red_max )
#  q_c_ss$lvl <- 0
#  q_c_ss$nstout <- NA
#  
#  q_c_ss <- read_aic (q = q_c_ss )
#  q_c_ss <- cr_filter ( q = q_c_ss )

  #####################################################
  # filling list for disconnected ll, ss and cc
  #####################################################
  q_d_ll <- list()
  q_d_ss <- list()
  q_d_cc <- list()

  for ( f in c( "ll", "ss" , "cc" ) ) {
    x.name <-  paste( "q_d_", f, sep="" )

    x.val <- list()
    x.val$name <- paste( "xq-disc-", field, "-", f, sep="" )
    x.val$op <- NA
    x.val$cr <- c( chisq_red_min, chisq_red_max )
    x.val$lvl <- 0
    x.val$nstout <- NA
    x.val$f <- NULL
    x.val$p <- NA
    x.val$norm <- 1.0  # NO factor 0.5 here; factor of 1/2 for single flavor has been already included in first analysis step

    x.val <- read_aic (q = x.val )
    x.val <- cr_filter ( q = x.val )

    assign ( x.name, x.val )
  }

#  return ( list(q_d_ll=q_d_ll, q_d_ss = q_d_ss, q_d_cc = q_d_cc ) )

#  q_d_ss <- list()
#  q_d_ss$name <- "xq-disc-kaon-ss"
#  q_d_ss$op <- NA
#  q_d_ss$cr <- c( chisq_red_min, chisq_red_max )
#  q_d_ss$lvl <- 0
#  q_d_ss$nstout <- NA
 
#  q_d_ss <- read_aic (q = q_d_ss )
#  q_d_ss <- cr_filter ( q = q_d_ss )

#  q_d_cc <- list()
#  q_d_cc$name <- "xq-disc-kaon-cc"
#  q_d_cc$op <- NA
#  q_d_cc$cr <- c( chisq_red_min, chisq_red_max )
#  q_d_cc$lvl <- 0
#  q_d_cc$nstout <- NA
 
#  q_d_cc <- read_aic (q = q_d_cc )
#  q_d_cc <- cr_filter ( q = q_d_cc )

  #####################################################
  # filling list for (disconnected) gluon
  #####################################################
  g_d    <- list()
  g_d$name <- paste( "xg-disc-", field, sep="" )
  g_d$op <- "clover"
  g_d$cr <- c( chisq_red_min, chisq_red_max )
  g_d$lvl <- 0
  g_d$nstout <- nstout
  g_d$rstout <- 0.1290
  g_d$p <- NA
  g_d$f <- NA
  g_d$norm <- 1.0  # no further multiplicative normalization here

  g_d <- read_aic (q = g_d )
  g_d <- cr_filter ( q = g_d )

#  return ( list(q_d=g_d ))

  #####################################################
  # read selected bootstrap data
  #
  read_bs <- function (q) 
  #####################################################
  {
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

        ps_tag <- ""
        if ( !is.na(q$p ) ) {
          ps_tag <- paste ( q$p, "-s", ns, sep="" )
        } else {
          ps_tag <- paste( "s" , ns, sep="" )
        }
        f <- paste ( q$name, op_tag, "/lvl", q$lvl, "/", ps_tag, "/fit.", ensemble, ".", q$name, ".", ts, ".tc", tc, ".tf", tf[1], "_", tf[2], ".bs", sep="" )

      } else {

        f <- paste ( g_d$name, "/nstout", q$nstout, "/lvl", q$lvl, "/s", ns, "/fit.", ensemble, ".", q$name, "-", g_d$op, ".nstout", nstout, "_", formatC(rstout, width=6, digits=4, format="f"), ".", ts, ".tc", tc, ".tf", tf[1], "_", tf[2], ".bs", sep="" )
      }

      if ( !file.exists(f) ) {
        stop ( "[read_bs] cannot find file ", f )
      } else {
        message( "# [read_bs] reading from file ", f)
      }
      q$data[i,] <- read.table ( f )$V3 * q$norm
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
    cat ( "n values = ", formatC( c(q_c_ll$n, q_c_ss$n, q_d_ll$n, q_d_ss$n, q_d_cc$n, g_d$n), width=6, format="d"), "\n" , sep="")

    cat("# [aic_ren_cl] w n = ", w$n, "\n",
        "# [aic_ren_cl] LL = ", formatC(LL, width=10, format="d"), "\n" , sep="")
  
    start_time <- Sys.time()

    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    
    end_time <- Sys.time()
    message ( "# [aic_ren_cl] time for idx = ", end_time - start_time )
  
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

      #####################################################
      # per quark flavor and gluon
      #####################################################
      q_u <- NA
      q_d <- NA
      q_s <- NA
      q_c <- NA
      g   <- NA

      if ( field == "kaon" ) 
      {
        q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
        q_d <-                    q_d_ll$data[i3,]
      
        q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]
      
        q_c <-                    q_d_cc$data[i5,]
  
        g   <-                                      g_d$data[i6,]
      }
  
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
    message( "# [aic_ren_cl] time for aic bs combination = ", end_time - start_time )

    start_time <- Sys.time()
  
    output_filename <- paste( "x-", field, ".", ensemble, ".singlet.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
    saveRDS( object = w$res, file = output_filename )
              
    end_time <- Sys.time()
    message( "# [aic_ren_cl] time for saveRDS = ", end_time - start_time )

    # return(w)

    #####################################################
    # determine quantiles
    #####################################################
    if ( with_quantiles  ) 
    {
  
      start_time <- Sys.time()
     
      cores=detectCores()
      cl <- makeCluster(cores[1]-1)
      registerDoParallel(cl)
    
      bres <- foreach ( b = 1:7, .combine = rbind ) %dopar% {
    
        wdf <- function(x, w=w, m=m, s=s, a=0 ) {
          # return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
  
          return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
        }
     
        lower <- c( -2.5, -2.5, -2.5 )
        upper <- c( +2.5, +2.5, +2.5 )
           
        ba  <- w$res[,1]
        bm  <- w$res[,(2*b  )]
        bs  <- w$res[,(2*b+1)]
     
        # res <- array( dim=c( log2(nprobe) - 3 ,3 ) )
        res <- array( 0, dim=c( log2(nprobe) - 3 ,3 ) )
     
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
      message( "# [aic_ren_cl] time for quantiles = ", end_time - start_time )
  
      # return (bres)
    
      bnames <- c( "u", "d", "s", "c", "q", "g", "qg" )
    
      output_filename <- paste( "x-", field, ".", ensemble, ".singlet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
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
      
    }  # end of if with_quantiles

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
    cat("# [aic_ren_cl] w n = ", w$n, "\n",
        "# [aic_ren_cl] LL = ", formatC(LL, width=6, format="d"), "\n" , sep="")
  
    nprobe_orig <- nprobe
    if ( nprobe > w$n ) nprobe <- 2^( as.integer( log2(w$N) ) )


    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    message ( "# [aic_ren_cl] idx done " )
  
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
      i5 <- ( (idx[k]-1) %% LL[4] )           + 1
  
      q_u <- NA
      q_d <- NA
      q_s <- NA
      q_c <- NA
      if ( field == "kaon" ) 
      {
        q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
        q_d <-                    q_d_ll$data[i3,]
      
        q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]
      
        q_c <-                    q_d_cc$data[i5,]
  
      }

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
  
    output_filename <- paste( "x-", field, ".", ensemble, ".nonsinglet.q4.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
    saveRDS( object = w$res, file = output_filename )
              
    # return(w)
  
    #####################################################
    if ( with_quantiles  )
    {

      #cores=detectCores()
      #cl <- makeCluster(cores[1]-1)
      #registerDoParallel(cl)
    
      b <- 1
      # bres <- foreach ( b = 1:1, .combine = rbind ) # %dopar% 
      #{
    
        wdf <- function(x, w=w, m=m, s=s, a=0 ) {
          return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
        }
    
        lower <- c( -2.5, -2.5, -2.5 )
        upper <- c( +2.5, +2.5, +2.5 )
          
        ba  <- w$res[,1]
        bm  <- w$res[,(2*b  )]
        bs  <- w$res[,(2*b+1)]
    
        bres <- array( dim=c( log2(nprobe) - 3 ,3 ) )
    
        for ( i in 1:(log2(nprobe)-3) ) {
          k <- 2^(i+3)
    
          a   <- ba[1:k] / sum ( ba[1:k] )
    
          m  <- bm[1:k]
          s  <- bs[1:k]
    
          u_m <- uniroot( wdf, lower=lower[1], upper=upper[1] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.50, s=s, w=a, m=m )
          u_l <- uniroot( wdf, lower=lower[2], upper=upper[2] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.16, s=s, w=a, m=m )
          u_u <- uniroot( wdf, lower=lower[3], upper=upper[3] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.84, s=s, w=a, m=m )
    
          bres[i, ] <- c( u_m$root, u_l$root , u_u$root )
    
    ##      cat ( "# g convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_filename, append=TRUE )
    
        }
        #
        #res
    
      # }  # end of loop on observables
      # stopCluster(cl)
      
      
      # return (bres)
    
      bnames <- c( "q4" )
    
      output_filename <- paste( "x-", field, ".", ensemble, ".nonsinglet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
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

    }  # end of if with_quantiles

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

    cat("# [aic_ren_cl] w n = ", w$n, "\n",
        "# [aic_ren_cl] LL = ", formatC(LL, width=6, format="d"), "\n" , sep="")
  
    nprobe <- nprobe_orig
    if ( nprobe > w$n ) nprobe <- 2^( as.integer ( log2 ( w$n) ) )


    idx <- sample.int ( n=w$n, size = nprobe, replace = FALSE )
    message ( "# [aic_ren_cl] idx done " )
  
    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)


    # for ( k in 1:nprobe ) 
    w$res <- foreach ( k = 1:nprobe , .combine=rbind ) %dopar% 
    {
  
      i1 <- ( (idx[k]-1)          ) %/% LL[1] + 1
      i2 <- ( (idx[k]-1) %% LL[1] ) %/% LL[2] + 1 
      i3 <- ( (idx[k]-1) %% LL[2] ) %/% LL[3] + 1
      i4 <- ( (idx[k]-1) %% LL[3] )           + 1
  
      q_u <- q_c_ll$data[i1,] + q_d_ll$data[i3,]
  
      q_d <-                    q_d_ll$data[i3,]
      
      q_s <- q_c_ss$data[i2,] + q_d_ss$data[i4,]

      # sum of quark flavors
      q_f <- q_u + q_d - 2 * q_s

      # renormalized combination
      q_r <- sZ_QQ2 * q_f

      a <- q_c_ll$f$V10[i1] * q_c_ss$f$V10[i2] * q_d_ll$f$V10[i3] * q_d_ss$f$V10[i4]


      str(a)
      # str(q_r)

      c( a, mean(q_r),   sqrt(var(q_r)) )

    }
  
    stopCluster(cl)
  
    output_filename <- paste( "x-", field, ".", ensemble, ".nonsinglet.q3.nprobe", as.integer(log2(nprobe)), ".aic.prb", sep="" )
    # cat ( "# ", date(), "\n", file = output_filename, append=F )
    # write.table ( w$res, file=output_filename, col.names=F, row.names=F, append=T )
    saveRDS( object = w$res, file = output_filename )
              
    # return(w)

    #####################################################
    if ( with_quantiles  )
    {

      #cores=detectCores()
      #cl <- makeCluster(cores[1]-1)
      #registerDoParallel(cl)
 
      b <- 1   
      #bres <- foreach ( b = 1:1, .combine = rbind ) # %dopar% 
      #{
    
        wdf <- function(x, w=w, m=m, s=s, a=0 ) {
          return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
        }
    
        lower <- c( -0.5, -0.5, -0.5 )
        upper <- c( +1.5, +1.5, +1.5 )
          
        ba  <- w$res[,1]
        bm  <- w$res[,(2*b  )]
        bs  <- w$res[,(2*b+1)]
    
        #res <- array( dim=c( log2(nprobe) - 3 ,3 ) )
        bres <- array( dim=c( log2(nprobe) - 3 ,3 ) )
    
        for ( i in 1:(log2(nprobe)-3) ) {
          k <- 2^(i+3)
    
          a   <- ba[1:k] / sum ( ba[1:k] )
    
          m  <- bm[1:k]
          s  <- bs[1:k]
    
          u_m <- uniroot( wdf, lower=lower[1], upper=upper[1] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.50, s=s, w=a, m=m )
          u_l <- uniroot( wdf, lower=lower[2], upper=upper[2] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.16, s=s, w=a, m=m )
          u_u <- uniroot( wdf, lower=lower[3], upper=upper[3] , check.conv = T, tol = 1.e-8, maxiter = 100000, a=0.84, s=s, w=a, m=m )
    
          bres[i, ] <- c( u_m$root, u_l$root , u_u$root )
    
    ##      cat ( "# g convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_filename, append=TRUE )
    
        #}
        #
        #res
    
      }  # end of loop on observables
      # stopCluster(cl)
 
      bnames <- c( "q3" )
    
      output_filename <- paste( "x-", field, ".", ensemble, ".nonsinglet.nprobe", as.integer(log2(nprobe)), ".aic.stats", sep="" )
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
  
    }  # end of if with_quantiles

  }  # end of nonsinglet

  return ( NULL )

}  # end of aic_ren_cl
