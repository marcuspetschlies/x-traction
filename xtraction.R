library(hadron)

#######################################################
# plot functions including error bars from the script
#######################################################
plotwitherror <- function(x, y, dy, col="black", ...) {
      plot(x=x, y=y, col=col, ...)
  arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01,
                           angle=90, code=3, col=col)
}

#######################################################
# clean up open devices
#######################################################
cleandev <- function(d) {
    if(missing(d)) {
          while(dev.cur()>1) dev.off(dev.cur())
  } else for(i in 1:length(d)) dev.off(d[i])
}

#############################################################
# 
# p = Z0, E0,        ... lvl 0
# p = Z0, E0, Z1, E1 ... lvl 1 
#############################################################
ftwop <- function( p, tf, fitinfo ) {
                  
  lvl <- fitinfo$lvl
  TT  <- fitinfo$TT

  # str( tf )

  #message( "# [ftwop] TT = ", TT )
  #message( "# [ftwop] lvl= ", lvl )

  res <- p[1]^2 / (2 * p[2] ) * ( exp ( -p[2] * tf )  + exp ( -p[2] * ( TT - tf ) ) )
  
  if ( lvl > 0 ) {
    res <- res +  p[3]^2 / (2 * p[4] ) * ( exp ( -p[4] * tf )  + exp ( -p[4] * ( TT - tf ) ) )
  }
 
  #message( "# [ftwop] done ")
  return ( res )  
}

#############################################################
# Operator 4k 3-point function time-dependence
# 
# p = Z0, E0, x,                  ... lvl 0
# p = Z0, E0, Z1, E1, x, M01, M11 ... lvl 1
#############################################################
fthreep4k <- function ( p, tc, tf, fitinfo) {

  lvl <- fitinfo$lvl
  TT  <- fitinfo$TT
  res <- numeric()

  if ( lvl == 0 ) {
    res <- rep( -p[1]^2 / ( 2 * p[2] ) * p[3] * exp ( -p[2] * tf ), times=length(tc) )
  }

  if ( lvl == 1 ) {
    res <- ( -p[1]^2 / ( 2 * p[2] ) * p[5] * exp ( -p[2] * tf )
         + p[1]*p[3] / ( 4 * p[2] * p[4] ) * p[6] * ( exp ( -p[2] * (tf -tc) -p[4] * tc ) + exp ( -p[4] * (tf -tc) -p[2] * tc ) )
         + p[3]^2 / ( 4 * p[4]^2 ) * p[7] * exp( -p[4] * tf ) )
  }

  return ( res )
}


#############################################################
# Operator 44 3-point function time-dependence
# 
# p = Z0, E0, x,                  ... lvl 0
# p = Z0, E0, Z1, E1, x, M01, M11 ... lvl 1
#
# !!! NOTE: only for zero 3-momentum here !!!
#############################################################
fthreep44 <- function ( p, tc, tf, fitinfo ) {

  lvl <- fitinfo$lvl
  TT  <- fitinfo$TT

  c <- 1 + sum( fitinfo$p^2 ) * (2*pi/fitinfo$LL)^2 / ( 3 * p[2]^2 )

  res <- 0

  # str( tf )
  # str( tc )
  # message( "# [threep44] tc = ", tc, ", tf = ", tf )

  if ( lvl == 0 ) {
    res <- rep( -3./8. * c * p[1]^2 * p[3] * exp ( -p[2] * tf ), times=length(tc) )
  }

  if ( lvl == 1 ) {
    res <- ( -3./8. * c * p[1]^2 * p[5] * exp ( -p[2] * tf ) 
           + p[1]*p[3] / ( 4 * p[2] * p[4] ) * p[6] * ( exp ( -p[2] * (tf -tc) -p[4] * tc ) + exp ( -p[4] * (tf -tc) -p[2] * tc ) ) 
           + p[3]^2 / ( 4 * p[4]^2 ) * p[7] * exp( -p[4] * tf ) ) # + p[8] * exp(-p[2] * (TT - tf)) )
  }

  # message( "# [fthreep44] done ")
  return ( res )
}

#############################################################
# chisq function
#
# p = parameters
# d = data
# ...$twop   = 2pt data
# ...$threep = 3pt data
#############################################################

fchisq <- function (p, d, fitinfo ) {

  ntf <- length( fitinfo$tf_range_threep )
  ntc <- length( fitinfo$tc_range_threep )

  # message( "# [fchisq] ntf = ", ntf )
  # message( "# [fchisq] ntc = ", ntc )

  # enforce ordering of masses
  if ( fitinfo$lvl == 1 ) {
    if ( p[4] < p[2] ) {
      aux <- p[1:2]
      p[1:2] <- p[3:4]
      p[3:4] <- aux
      rm(aux)
    }
  }


  # 2-point function parametrization
  a <- fitinfo$ftwop ( p, fitinfo$tf_range_twop, fitinfo )

  # return ( a )

  b <- array ( dim=c(length( fitinfo$tf_range_threep ), length( fitinfo$tc_range_threep ) ) )

  for ( itf in 1:ntf ) {
    tf <- fitinfo$tf_range_threep[itf]

    tc <- fitinfo$tc_range_threep + tf %/% 2

    # message( "# [fchisq] tc = ", tc )

    b[itf,] <- fitinfo$fthreep ( p, tc, tf, fitinfo )
  }

  fvec <- c(a, t(b) )

  # return ( fvec )

  yvec <- c( d$twop, t( d$threep ) )

  # return ( yvec )

  covi <- d$covi
  if ( is.null ( covi ) ) {
    if ( is.null ( d$cov ) ) {
      # diagonal 1 
      covi <- diag ( rep( 1, times=length( fvec ) ) )
    } else {
      # svd-based inverse of cov
      a <- svd ( d$cov )
      covi <- a$u %*% diag( 1./a$d ) %*% t( a$u )
    }
  # } else {
  #   message( "using existing covi" )
  }

  dfy <- yvec - fvec

  # return( dfy )

  res <- t(dfy) %*% covi %*% dfy 

  return ( res )

}  # end of fchisq

#############################################################
# minimization
# 
# IN:
#  kind of fit
#  fit time ranges tf for 2pt, tf, tc for 3pt
#  
#############################################################
fminimize <- function (par0, fitinfo, fitdata, bs = NULL ) {
  
  par <- par0

  d <- list()

  d[["twop"]] <- apply ( fitdata$twop, c(1), mean )
  
  d[["threep"]] <- apply ( fitdata$threep, c(1,2), mean )

  #############################################################
  # fill covariance matrix
  #############################################################

  twop_dim   <-  dim( fitdata$twop )
  threep_dim <-  dim( fitdata$threep )

  ddim <- c( twop_dim[1] +  prod( threep_dim[1:2] ), twop_dim[2] )

  # message( "# [fminimize] ddim = ", ddim[1], ", ", ddim[2] )

  # full covariance matrix
  b <- NULL
  if ( threep_dim[1] * threep_dim[2] == 1 ) {
    b <- cbind( t( fitdata$twop ),   matrix( fitdata$threep, threep_dim[3], 1) )
  } else {
    b <- cbind( t( fitdata$twop ),   t( apply( fitdata$threep, c(3) , t ) ) )
  }

  d[["cov"]] <- cov( b ) / ddim[2]

  # block diagonal covariance matrix
  #d[["cov"]] <- array ( 0, c(ddim[1], ddim[1]) )
  #d$cov[1:twop_dim[1], 1:twop_dim[1]] <- cov ( t( fitdata$twop ) ) / ddim[2]
  #b <- cbind( t( apply( fitdata$threep, c(3) , t ) ) )
  #d$cov[(twop_dim[1]+1):ddim[1], (twop_dim[1]+1):ddim[1]] <- cov ( b ) / ddim[2]

  rm(b)

  b <- svd ( d$cov )
  d[["covi"]] <- b$u %*% diag( 1./b$d ) %*% t( b$u )

  rm(b)

  # return ( d )
  # return( fchisq(p=par0, d=d, fitinfo ) )
  # return ( list(d=d, i=fitinfo) )

  #############################################################
  # fit on ORIGINAL data
  #############################################################
  uorig <- optim ( par, fchisq, gr = NULL, method="BFGS", control = list(maxit = 100000, abstol=1.e-14, reltol=1.e-14 ), d=d, fitinfo =fitinfo)

  # uorig <- optim ( par, fchisq, gr = NULL, 
  #                 control = list(maxit = 100000, abstol=1.e-14, reltol=1.e-14 ),
  #                 method="BFGS", d=d, fitinfo =fitinfo)

  # return ( uorig )

  if ( uorig$convergence != 0 ) {
    stop( "[fminimize] optim orig exit status ", uorig$convergence )
  }

  #############################################################
  # fit on SAMPLED data
  #############################################################
  if ( !is.null ( bs ) && bs$nsample > 0 ) {
    bs[["fitres"]]      <- array ( dim=c(bs$nsample, length(par) ) )
    bs[["chisq"]]       <- numeric()
    bs[["data"]]        <- array( dim=c( bs$nsample, ddim[1] ) )
    bs[["cov"]]         <- array( dim=c( bs$nsample, ddim[1], ddim[1] ) )
    bs[["counts"]]      <- integer()
    bs[["convergence"]] <- integer()

    for ( s in 1:bs$nsample ) {
  
      ds <- list()

      idx <- sample.int ( n=ddim[2], size=ddim[2], replace = TRUE )
      # cat( "[", s, "]", formatC( idx, width=5, format="d" ) , "\n")
  
      ds[["twop"]] <- apply ( fitdata$twop[,idx, drop=F], c(1), mean )
    
      ds[["threep"]] <- apply ( fitdata$threep[,,idx, drop=F], c(1,2), mean )
  
      # add current data to bs
      # bs$data[s,] <- apply( cbind( t( fitdata$twop[,idx, drop=F] ),   t( apply( fitdata$threep[,,idx, drop=F], c(3) , t ) ) ) , c(2) , mean )
  
  
      # cov constant
      # ds[["cov"]] <- d$cov
  
      # full cov per sample
      b <- NULL
      if ( threep_dim[1] * threep_dim[2] == 1 ) {
        b <- cbind( t( fitdata$twop[,idx, drop=F] ),   matrix( fitdata$threep[,,idx, drop=F] ) )
      } else {
        b <- cbind( t( fitdata$twop[,idx, drop=F] ),   t( apply( fitdata$threep[,,idx, drop=F], c(3) , t ) ) )
      }

      ds[["cov"]] <- cov( b ) / ddim[2]
  
      # block diagonal cov per sample
      #ds[["cov"]] <- array ( 0, c(ddim[1], ddim[1]) )
      #ds$cov[1:twop_dim[1], 1:twop_dim[1]] <- cov ( t( fitdata$twop[,idx, drop=F] ) ) / ddim[2]
      #b <- cbind( t( apply( fitdata$threep[,,idx,drop=F], c(3) , t ) ) )
      #ds$cov[(twop_dim[1]+1):ddim[1], (twop_dim[1]+1):ddim[1]] <- cov ( b ) / ddim[2]
  
      bs$cov[s,,] <- ds$cov
  
      rm(b)

      # add covariance matrix to bs
      # b <- svd ( ds$cov )
      # ds[["covi"]] <- b$u %*% diag( 1./b$d ) %*% t( b$u )
      # rm(b)
  
      par <- par0
  
      u <- optim ( par, fchisq, gr = NULL, method="BFGS", control = list(maxit = 100000, abstol=1.e-14, reltol=1.e-12), d=ds, fitinfo =fitinfo )
  
      if ( u$convergence != 0 ) {
        stop( "[fminimize] optim exit status ", u$convergence )
      }
  
      bs$fitres[s,] <- u$par
      bs$chisq[s]   <- u$value
  
      bs$counts[s]  <- u$counts[1]
      bs$convergence[s] <- u$convergence
  
    }  # end of bootstrap sampling

  }  # end of if bs not null and samples > 0

  #############################################################
  # analyse bootstrap
  #############################################################
  res <- list()
  res$par_value <- uorig$par
  res$chisq     <- uorig$value
  res$dof       <- ddim[1] - length(par)
  res$ndata     <- ddim[1]
  res$npar      <- length(par)
  res$aic       <- exp ( -0.5 * ( res$chisq  + 2 * res$npar - res$ndata ) )
  res$par_cov   <- array( 0, dim=c( length(par), length(par) ) )
  if ( !is.null ( bs ) && bs$nsample > 0 ) {
    res$par_cov   <- cov( bs$fitres )
  }
  res$orig      <- uorig
  res$bs        <- bs
  res$info      <- fitinfo

  return ( res )

}  # end of fminimize

#############################################################
#############################################################

#############################################################
# show fit info content
#############################################################
show_fitinfo <- function ( f ) {
  if ( missing ( f ) ) {
    message( "# [show_fitinfo] no fit info")
    return ( NULL )
  }

  cat ( "\n",
       "# [show_fitinfo] TT = ", f$TT, "\n",
       "# [show_fitinfo] tf_range_twop = ", f$tf_range_twop, "\n",
       "# [show_fitinfo] tf_range_threep = ", f$tf_range_threep, "\n",
       "# [show_fitinfo] tc_range_threep = ",  f$tc_range_threep, "\n",
       "# [show_fitinfo] operator        = ",  f$operator, "\n",
       "# [show_fitinfo] p               = ",  formatC(f$p , width=3, format="d"), "\n",
       "# [show_fitinfo] lvl             = ",  f$lv, "\n" )
}  # end of show_fitinfo

#############################################################
#############################################################

#############################################################
# run the minimize function
#############################################################
run_min <- function( ens="cB211.072.64", obs="xq-conn", nconf=790, nsrc=8, TT=128,
                     p=c(0,0,0), operator="g4_D4",
                     path_to_data="..",
                     twop_tf_range, threep_tf_range, threep_tc_range,
                     lvl=0, par0 =c(1,1,1),
                     nsample = 100, seed,
                     output_tag="fit",
                     twop_prefix="twop.pseudoscalar.orbit",
                     threep_prefix="threep.conn" , 
                     twop_col=2,
                     threep_col=2) {

  twop_filename <- paste ( path_to_data, "/", ens, "/", obs, "/", twop_prefix, ".PX", p[1], "_PY", p[2], "_PZ", p[3], ".re.corr", sep="" )
  if ( !file.exists( twop_filename) ) stop( "Could not find ", twop_filename )

  message( "# [run_min] reading 2pt from file ", twop_filename )
  # d <- apply( array( read.table ( twop_filename )$V2, dim=c(TT, nsrc, nconf) ), c(1,3), mean )
  d <- apply( array( read.table ( twop_filename )[,twop_col], dim=c(TT, nsrc, nconf) ), c(1,3), mean )
  twop_data <- ( d[1:TT, ] + d[ ((TT:1)%%TT + 1), ] ) * 0.5

  # return ( twop_data)

  threep_data <- array( dim=c( length( threep_tf_range ), TT, nconf ) )

  for (idt in 1:length(threep_tf_range) )  {
    dt <- threep_tf_range[idt]

    if ( dt <= TT / 2 ) {
  
      threep_filename <- paste ( path_to_data, "/", ens, "/", obs, "/", threep_prefix, ".", operator, ".dtsnk", dt, ".PX", p[1], "_PY", p[2], "_PZ", p[3], ".re.corr", sep="" )
      if ( !file.exists( threep_filename) ) stop( "Could not find ", threep_filename )
      message( "# [run_min] reading 3pt from file ", threep_filename )

      # threep_data[idt,,] <- apply( array( read.table ( threep_filename )$V2, dim=c(TT, nsrc, nconf) ), c(1,3), mean )
      threep_data[idt,,] <- apply( array( read.table ( threep_filename )[,threep_col], dim=c(TT, nsrc, nconf) ), c(1,3), mean )

    } else {
      threep_filename <- paste ( 
                                path_to_data, "/", ens, "/", obs, "/", threep_prefix, ".", operator, ".dtsnk", TT-dt, ".PX", p[1], "_PY", p[2], "_PZ", p[3], ".re.corr", sep=""
      )

      if ( !file.exists( threep_filename) ) stop( "Could not find ", threep_filename )
      message( "# [run_min] reading 3pt from file ", threep_filename )
      
      idx <- ( TT : 1 ) %% TT + 1

      if        ( operator == "g4_D4" ) {
        threep_data[idt,,] <-  apply( array( read.table ( threep_filename )[,threep_col], dim=c(TT, nsrc, nconf) ), c(1,3), mean )[idx,,drop=F]
      } else if ( operator == "g4_Dk" ) {
        threep_data[idt,,] <- -apply( array( read.table ( threep_filename )[,threep_col], dim=c(TT, nsrc, nconf) ), c(1,3), mean )[idx,,drop=F]
      }
    }
  }

  # return ( threep_data )

  #############################################################
  # prepare fit data set
  #############################################################
  fitdata <- list()
  fitdata[["twop"]]   <- twop_data[ ((twop_tf_range[1]:twop_tf_range[2])+1),, drop=F]
  fitdata[["threep"]] <- array( dim=c(length(threep_tf_range), length( (threep_tc_range[1]:threep_tc_range[2]) ),  nconf ) )
  for (idt in 1:length(threep_tf_range) )  {
    dt <- threep_tf_range[idt]
    fitdata[["threep"]][idt,,] <- threep_data[idt,((threep_tc_range[1]:threep_tc_range[2])+(dt%/%2)+1),]
  }
 
  # return ( fitdata )

  #############################################################
  # prepare fit into
  #############################################################
  fitinfo <- list()
  fitinfo[["TT"]]              <- TT
  fitinfo[["LL"]]              <- TT %/% 2
  fitinfo[["tf_range_twop"]]   <- ( twop_tf_range[1] : twop_tf_range[2] )
  fitinfo[["tf_range_threep"]] <- threep_tf_range
  fitinfo[["tc_range_threep"]] <- ( threep_tc_range[1] : threep_tc_range[2] )
  fitinfo[["ftwop"]]           <- ftwop
  if ( operator == "g4_D4" ) {
    fitinfo[["fthreep"]]       <- fthreep44
  } else if ( operator == "g4_Dk" ) {
    fitinfo[["fthreep"]]       <- fthreep4k
  }
  fitinfo[["lvl"]]             <- lvl
  fitinfo[["operator"]]        <- operator
  fitinfo[["p"]]               <- p

  show_fitinfo ( fitinfo )
  # return (fitinfo)

  #############################################################
  # prepare fit parameter list
  #############################################################
  fitparam <- par0

  # return ( invisible(list(data=fitdata, info=fitinfo ) ) )

  #############################################################
  # prepare boostrap configuration
  #############################################################
  bs <- list()
  bs[["nsample"]] <- nsample

  #############################################################
  # initialize rng
  #############################################################
  if ( missing ( seed ) ) {
    stop( "need seed value" ) 
  } else {
    message ( "# setting seed value ", seed )
  }
  bs[["seed"]]    <- seed
  set.seed ( seed )


  #############################################################
  # call minizer
  #############################################################
  fitres <- fminimize ( fitparam, fitinfo, fitdata, bs )

  # return (fitres)

  plt <- list()

  #############################################################
  # evaluate the fit data and function
  #############################################################

  # data
  plt$data <- list()

  # twop
  plt$data$twop <- array ( dim=c( length( fitinfo$tf_range_twop ), 6 ) )
  for ( i in 1:length( fitinfo$tf_range_twop ) ) {
    t <- fitinfo$tf_range_twop[i]
    u <- uwerrprimary ( fitdata$twop[i,] )
    plt$data$twop [i,] <- c ( t, u$value, u$dvalue, u$ddvalue, u$tauint, u$dtauint )
  }

  # threep
  plt$data$threep <- array ( dim=c( length(fitinfo$tf_range_threep), length( fitinfo$tc_range_threep ), 7 ) )
  for ( i in 1:length( fitinfo$tf_range_threep ) ) {
    tf <- fitinfo$tf_range_threep[i]

    for ( k in 1:length( fitinfo$tc_range_threep ) ) {
      
      tc <- fitinfo$tc_range_threep[k]

      u <- uwerrprimary ( fitdata$threep[i,k,] )
      plt$data$threep[i,k,] <- c ( tf, tc, u$value, u$dvalue, u$ddvalue, u$tauint, u$dtauint )
    }
  }

  # fit function

  plt$fit <- list()

  # twop
  plt$fit$twop <- array ( 0, dim=c( length( fitinfo$tf_range_twop ), 4 ) )

  for ( i in 1:length( fitinfo$tf_range_twop ) ) {

    tf <- fitinfo$tf_range_twop[i]
 
    u <- 0
    if ( bs$nsample > 0 ) {
      u <- apply ( fitres$bs$fitres, c(1), fitinfo$ftwop, tf = tf, fitinfo = fitinfo )
    }
    plt$fit$twop [i,] <- c ( tf, fitinfo$ftwop ( fitres$par_value, tf, fitinfo ),  mean( u ), sqrt( var ( u ) ) )
  }

  # threep
  plt$fit$threep <- array ( 0, dim=c( length( fitinfo$tf_range_threep ), length( fitinfo$tc_range_threep ), 5 ) )

  for ( i in 1:length( fitinfo$tf_range_threep ) ) {
    tf <- fitinfo$tf_range_threep[i]

    for ( k in 1:length( fitinfo$tc_range_threep ) ) {
      tc <- fitinfo$tc_range_threep[k] + tf %/% 2

      u <- 0
      if ( bs$nsample > 0 ) {
        u <- apply ( fitres$bs$fitres, c(1), fitinfo$fthreep, tf = tf, tc = tc, fitinfo = fitinfo )
      }

      plt$fit$threep [i,k,] <- c ( tf, fitinfo$tc_range_threep[k], fitinfo$fthreep( fitres$par_value, tc, tf, fitinfo ),  mean( u ), sqrt( var ( u ) ) )
    }
  }

  #############################################################
  # write fit results to file
  #############################################################
  output_filename <- paste( "fit.", output_tag, ".res", sep="" )
  cat ( "# ", date(), "\n", file=output_filename, append=F )

  cat ( "# tf_range_twop ", fitinfo$tf_range_twop, "\n", file=output_filename, append=T )
  cat ( "# tf_range_threep ", fitinfo$tf_range_threep, "\n", file=output_filename, append=T )
  cat ( "# tc_range_threep ", fitinfo$tc_range_threep, "\n", file=output_filename, append=T )
  cat ( "# lvl ", fitinfo$lvl, "\n", file=output_filename, append=T )
  cat ( "# operator ", fitinfo$operator, "\n", file=output_filename, append=T )
  cat ( "# obs ", obs, "\n", file=output_filename, append=T )

  cat ( "# threep_prefix ", threep_prefix, "\n", file=output_filename, append=T )
  cat ( "# twop_prefix ", twop_prefix, "\n", file=output_filename, append=T )


  for ( i in 1:length(fitres$par_value) ) {
    cat ( "par", formatC(i, width=3, format="d" ), 
         formatC(fitres$par_value[i], width=16,  digits=7),
         formatC( sqrt( fitres$par_cov[i,i]), width=16,  digits=7), "\n",
         file = output_filename, append=T )
  }
  cat ("chisq ", fitres$chisq, "\ndof ", fitres$dof, "\nndata ", fitres$ndata, "\nnpar ", fitres$npar, "\naic ", fitres$aic, "\n",
         file = output_filename, append=T, sep="" )

  for ( i in 1:(length(fitres$par_value)-1) ) {
    for ( k in (i+1):(length(fitres$par_value)) ) {
      c <- fitres$par_cov[i,k] / sqrt ( fitres$par_cov[i,i] ) / sqrt( fitres$par_cov[k,k] )
      cat ( "corr", 
           formatC( i, width=3, format="d" ),
           formatC( k, width=3, format="d" ),
           formatC( c, width=16,  digits=7), "\n",
           file = output_filename, append=T )
    }
  }

  output_filename <- paste( "fit.", output_tag, ".bs", sep="" )
  cat ( "# ", date(), "\n", file=output_filename, append=F )
  cat ("# samples ", fitres$bs$nsample, "\n", file=output_filename, append=T )
  cat ("# seed    ", fitres$bs$seed, "\n", file=output_filename, append=T )
  write.table ( fitres$bs$fitres, file=output_filename, col.names=F, row.names=F, append=TRUE )

  #############################################################
  # write plot data to file
  #############################################################

  # twop
  output_filename <- paste( "twop.", output_tag, ".plt", sep="" )
  cat ( "# ", date(), "\n", file=output_filename, append=F )
  for ( i in 1:dim( plt$data$twop)[1] ) {

    cat ( 
         formatC ( plt$data$twop[i,1], width=6, format="d" ),
         formatC ( plt$data$twop[i,2], width=16, digits=7, format="e" ),
         formatC ( plt$data$twop[i,3], width=16, digits=7, format="e" ),
         formatC ( plt$data$twop[i,4], width=16, digits=7, format="e" ),
         formatC ( plt$data$twop[i,5], width=16, digits=7, format="e" ),
         formatC ( plt$data$twop[i,6], width=16, digits=7, format="e" ),
         #
         formatC ( plt$fit$twop[i,2], width=16, digits=7, format="e" ),
         formatC ( plt$fit$twop[i,3], width=16, digits=7, format="e" ),
         formatC ( plt$fit$twop[i,4], width=16, digits=7, format="e" ),
         "\n", sep="", file=output_filename, append=T )
  }

  # threep
  output_filename <- paste( "threep.", output_tag, ".plt", sep="" )
  cat ( "# ", date(), "\n", file=output_filename, append=F )
  for ( i in 1:dim( plt$data$threep)[1] ) {
    for ( k in 1:dim( plt$data$threep)[2] ) {

      cat ( 
         formatC ( plt$data$threep[i,k,1], width=6, format="d" ),
         formatC ( plt$data$threep[i,k,2], width=6, format="d" ),
         formatC ( plt$data$threep[i,k,3], width=16, digits=7, format="e" ),
         formatC ( plt$data$threep[i,k,4], width=16, digits=7, format="e" ),
         formatC ( plt$data$threep[i,k,5], width=16, digits=7, format="e" ),
         formatC ( plt$data$threep[i,k,6], width=16, digits=7, format="e" ),
         formatC ( plt$data$threep[i,k,7], width=16, digits=7, format="e" ),
         #
         formatC ( plt$fit$threep[i,k,3], width=16, digits=7, format="e" ),
         formatC ( plt$fit$threep[i,k,4], width=16, digits=7, format="e" ),
         formatC ( plt$fit$threep[i,k,5], width=16, digits=7, format="e" ),
         "\n", sep="", file=output_filename, append=T )
    }
  }


 


  #############################################################
  # return data, info and res
  #############################################################
  return ( invisible(list(data=fitdata, info=fitinfo, res=fitres, plt=plt) ) )

}  # end of run_min

#############################################################
#############################################################
 
#############################################################
# sequence of fits to check systematics
#############################################################
fit_sequence_conn <- function (seed, aic_file, obs, threep_prefix,   nsample = 600 , operator = "g4_D4", mom = c(0,0,0), ens, nconf, TT ) {

  if ( missing ( obs   ) | missing ( threep_prefix ) ) stop ( "need obs and threep_prefix" )
  if ( missing ( ens   ) ) stop ( "need ens" )
  if ( missing ( nconf ) ) stop ( "need nconf" )
  if ( missing ( TT    ) ) stop ( "need TT" )

  # TT    <- 128
  # ens   <- "cB211.072.64"
  # nconf <- 745

  nsrc  <- 8
  
  path_to_data  <- "../../../"
  # nsample       <- 0

  lvl <- 0

  

  par0 <- c( 5., 0.19, 1. )
  # par0 <- c( 50, 0.06, 0, 1, 0.5, 0, 0 )

  # threep_tf_list <- c( 24, 36, 48, 56, 64, 72, 80, 92, 104 )
  # threep_tf_list <- c( 36, 48, 56, 64, 72 )

  # cC80 
  # threep_tf_list <- c( 30, 48, 64, 72, 80 )

  # cD96 
  threep_tf_list <- c( 36, 58, 76, 86, 96 )


  # twop_tf_range <- c( 48, 64 ) 
  twop_tf_range <- c( 36, 72 ) 


  #obs   <- "xq-conn-kaon-ll"
  #threep_prefix <- "threep.l-gd-ls-gi.conn"

  #obs   <- "xq-conn-kaon-ss"
  #threep_prefix <- "threep.s-gd-sl-gi.conn"

  twop_prefix   <- "twop.s-gf-l-gi.pseudoscalar.orbit"
  twop_col      <- 2
  threep_col    <- 2


  if ( !missing(aic_file ) ) {

    d <- read.table ( aic_file ) 

    wt <- sum ( d$V10 )

    r <- d$V10 / wt

    idx <- which(r>0.00001 )

    message ( "# [fit_sequence_conn] including fraction ",   sum ( r [ idx ] ) , " with ", length(idx), " entries" )

    for ( i in idx ) {

      lvl <- d[i,1]

      twop_tf_range <- c( d[i,2], d[i,3])

      tc <- d[i,4]

      threep_tc_range <- c(-tc, tc )

      s <- unlist( strsplit(x=as.character(d[i,5]), split="[ts_]") )

      threep_tf_list <- as.integer( s[3:length(s)] )

      threep_tf_range <- c( threep_tf_list )

      ts_tag <- "ts"
      for ( j in 1:length(threep_tf_list) ) {
        if ( j == 1 ) {
          ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
        } else {
          ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
        }
      }

      tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )
      # message ( "# [fit_sequence_conn] tag = ", tag )

      # check existence of res file

      file_res <- paste ( "fit.", tag, ".res", sep="" )
      if ( file.exists ( file_res ) ) {
        message ( "# [fit_sequence_conn]  SKIP  ", tag )
        next
      } else {
        message ( "# [fit_sequence_conn]  START ", tag )
      }

      message ("# [fit_sequence_conn] threep_tf_range = ", formatC(threep_tf_range , width=4, format="d" ))
      message ("# [fit_sequence_conn] threep_tc_range = ", formatC(threep_tc_range , width=4, format="d" ))
      message ("# [fit_sequence_conn] twop_tf_range   = ", formatC(twop_tf_range , width=4, format="d" ))

      r <- run_min ( ens          = ens,
                     obs          = obs,
                     nconf        = nconf,
                     nsrc         = nsrc,
                     TT           = TT,
                     p            = mom,
                     operator     = operator,
                     path_to_data = path_to_data,
                     lvl          = lvl,
                     par0         = par0,
                     nsample      = nsample,
                     seed         = seed,
                     output_tag   = tag ,
                     twop_tf_range, threep_tf_range, threep_tc_range,
                     threep_prefix = threep_prefix,
                     twop_prefix   = twop_prefix,
                     twop_col      = twop_col,
                    threep_col    = threep_col
      )

      message ( "# [fit_sequence_conn]  END   ", tag )

    }

  } else {

    for ( itf in 1:length(threep_tf_list) ) 
    # for ( itf in 1:1 )
    {
      for ( ktf in itf:length(threep_tf_list) ) 
      # for ( ktf in itf:itf)
      # for ( ktf in 4:4)
      {
  
        ts_tag <- "ts"
        for ( j in itf:ktf ) {
          if ( j == itf ) {
            ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
          } else {
            ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
          }
        }
  
      threep_tf_range <- c( threep_tf_list[itf:ktf] )
  
      threep_tc_max <- min( ( threep_tf_range * 3 ) %/% 8 )
  
      for ( tc in 0:threep_tc_max ) 
      {
  
        threep_tc_range <- c(-tc, tc )
  
        tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )
        message ( "# [fit_sequence] tag = ", tag )
  
  
        r <- run_min ( ens          = ens,
                   obs          = obs, 
                   nconf        = nconf,
                   nsrc         = nsrc, 
                   TT           = TT,
                   p            = mom, 
                   operator     = operator,
                   path_to_data = path_to_data,
                   lvl          = lvl,
                   par0         = par0,
                   nsample      = nsample, 
                   seed         = seed,
                   output_tag   = tag ,
                   twop_tf_range, threep_tf_range, threep_tc_range,
                   threep_prefix = threep_prefix,
                   twop_prefix   = twop_prefix,
                   twop_col      = twop_col,
                   threep_col    = threep_col
        )
  
      }
    }}
  }
}  # end of fit_sequence_conn


#############################################################
#############################################################
 
#############################################################
# sequence of fits to check systematics
#############################################################
fit_sequence_xg <- function ( seed , nstout = -1 , path_to_data="./", type = "clover" , aic_file, nsample =600 ) {

  # TT    <- 128
  # ens   <- "cB211.072.64"
  # nconf <- 745

  # TT    <- 160
  # ens   <- "cC211.06.80"
  # nconf <- 400

  TT    <- 192
  ens   <- "cD211.054.96"
  nconf <- 495
  

  nsrc  <- 1
  operator <- "g4_Dk"


  lvl <- 0

  mom <- c(0,0,1)

  par0 <- c( 5., 0.19, 0 )
  # par0 <- c( -0.009409329, 0.1134688, 0.5202603 )

  # threep_tf_list <- c( 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 28, 32)
  threep_tf_list <- c( 12, 16, 20, 24, 28, 32, 36, 40)

  # twop_tf_range <- c( 48, 64 ) 
  # twop_tf_range <- c( 32, 64 ) 
  twop_tf_range <- c( 36, 72 ) 


  obs   <- paste( "xg-disc-kaon-", type, ".nstout", nstout, "_0.1290" , sep ="")


  threep_prefix <- paste( "threep.orbit.src.", type, ".nstout", nstout, "_0.1290", sep="" )
  twop_prefix   <- "twop.orbit.gf5.gi5" 
  twop_col      <- 1
  threep_col    <- 1 

  if ( !missing(aic_file ) ) {

    d <- read.table ( aic_file )

    wt <- sum ( d$V10 )

    r <- d$V10 / wt

    idx <- which(r>0.00001 )

    message ( "# [fit_sequence_xg] including fraction ",   sum ( r [ idx ] ) , " with ", length(idx), " entries" )

    for ( i in idx ) {

      lvl <- d[i,1]

      twop_tf_range <- c( d[i,2], d[i,3])

      tc <- d[i,4]

      threep_tc_range <- c(-tc, tc )

      s <- unlist( strsplit(x=as.character(d[i,5]), split="[ts_]") )

      threep_tf_list <- as.integer( s[3:length(s)] )

      threep_tf_range <- c( threep_tf_list )

      ts_tag <- "ts"
      for ( j in 1:length(threep_tf_list) ) {
        if ( j == 1 ) {
          ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
        } else {
          ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
        }
      }

      tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )
      # message ( "# [fit_sequence_conn] tag = ", tag )

      # check existence of res file

      file_res <- paste ( "fit.", tag, ".res", sep="" )
      if ( file.exists ( file_res ) ) {
        message ( "# [fit_sequence_xg]  SKIP  ", tag )
        next
      } else {
        message ( "# [fit_sequence_xg]  START ", tag )
      }

      message ("# [fit_sequence_xg] threep_tf_range = ", formatC(threep_tf_range , width=4, format="d" ))
      message ("# [fit_sequence_xg] threep_tc_range = ", formatC(threep_tc_range , width=4, format="d" ))
      message ("# [fit_sequence_xg] twop_tf_range   = ", formatC(twop_tf_range , width=4, format="d" ))
  
      r <- run_min ( ens      = ens,
                 obs          = obs,
                 nconf        = nconf,
                 nsrc         = nsrc,
                 TT           = TT,
                 p            = mom,
                 operator     = operator,
                 path_to_data = path_to_data,
                 lvl          = lvl,
                 par0         = par0,
                 nsample      = nsample,
                 seed         = seed,
                 output_tag   = tag ,
                 twop_tf_range, threep_tf_range, threep_tc_range,
                 threep_prefix = threep_prefix,
                 twop_prefix   = twop_prefix,
                 twop_col      = twop_col,
                 threep_col    = threep_col
      )


    }

  } else {

    for ( itf in 1:(length(threep_tf_list)) ) 
    {

    for ( ktf in (itf):length(threep_tf_list) ) 
    {


      ts_tag <- "ts"
      for ( j in itf:ktf ) {
        if ( j == itf ) {
          ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
        } else {
          ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
        }
      }

    threep_tf_range <- c( threep_tf_list[itf:ktf] )

    # cat( "# [fit_sequence] threep_tf_range = ", threep_tf_range, "\n" )

    threep_tc_max <- min( ( threep_tf_range * 3 ) %/% 8 )

    for ( tc in 0:threep_tc_max ) 
    {

      threep_tc_range <- c(-tc, tc )

      tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )

      message ( "# [fit_sequence_xg] tag = ", tag )

      r <- run_min ( ens      = ens,
                 obs          = obs, 
                 nconf        = nconf,
                 nsrc         = nsrc, 
                 TT           = TT,
                 p            = mom, 
                 operator     = operator,
                 path_to_data = path_to_data,
                 lvl          = lvl,
                 par0         = par0,
                 nsample      = nsample,
                 seed         = seed,
                 output_tag   = tag ,
                 twop_tf_range, threep_tf_range, threep_tc_range,
                 threep_prefix = threep_prefix,
                 twop_prefix   = twop_prefix,
                 twop_col      = twop_col,
                 threep_col    = threep_col
      )

    }
    }
    }
  }
}  # end of fit_sequence_xg


#############################################################
#############################################################
 
#############################################################
# sequence of fits to check systematics
#############################################################
fit_sequence_xq_disc <- function ( seed , ens, TT, nconf, obs, threep_prefix, operator, nsample, aic_file ) {

  if ( missing ( TT            ) ) stop ("[fit_sequence_xq_disc] need TT")
  if ( missing ( ens           ) ) stop ("[fit_sequence_xq_disc] need ens")
  if ( missing ( seed          ) ) stop ("[fit_sequence_xq_disc] need seed")
  if ( missing ( nconf         ) ) stop ("[fit_sequence_xq_disc] need nconf")
  if ( missing ( obs           ) ) stop ("[fit_sequence_xq_disc] need obs")
  if ( missing ( operator      ) ) stop ("[fit_sequence_xq_disc] need operator")
  if ( missing ( nsample       ) ) stop ("[fit_sequence_xq_disc] need nsample")
  if ( missing ( threep_prefix ) ) stop ("[fit_sequence_xq_disc] need threep_prefix")

#  TT    <- 128

#  obs   <- "xq-disc-kaon-ll"
#  threep_prefix <- "threep.orbit.src.es3.nev200.Nstoch1"

#  obs   <- "xq-disc-kaon-ss"
#  threep_prefix <- "threep.orbit.src.es1.nev0.Nstoch1"

#  obs   <- "xq-disc-kaon-cc"
#  threep_prefix <- "threep.orbit.src.es1.nev0.Nstoch12"

#  ens   <- "cB211.072.64"

#  nconf <- 745

  nsrc  <- 1
  # operator <- "g4_Dk"
  path_to_data  <- "../../../"
  # nsample       <- 600

  mom <- c(0,0,1)
  lvl <- 0
  par0 <- c( 5., 0.19, 0 )

  # threep_tf_list <- c( 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 28, 32)
  threep_tf_list <- c( 12, 16, 20, 24, 28, 32, 36, 40)

  # twop_tf_range <- c( 48, 64 ) 
  # twop_tf_range <- c( 32, 64 ) 
  twop_tf_range <- c( 36, 72 ) 

  twop_prefix   <- "twop.orbit.gf5.gi5" 
  twop_col      <- 1
  threep_col    <- 1 


  if ( !missing(aic_file ) ) {

    d <- read.table ( aic_file )

    wt <- sum ( d$V10 )

    r <- d$V10 / wt

    idx <- which(r>0.00001 )

    message ( "# [fit_sequence_disc] including fraction ",   sum ( r [ idx ] ) , " with ", length(idx), " entries" )

    for ( i in idx ) {

      lvl <- d[i,1]

      twop_tf_range <- c( d[i,2], d[i,3])

      tc <- d[i,4]

      threep_tc_range <- c(-tc, tc )

      s <- unlist( strsplit(x=as.character(d[i,5]), split="[ts_]") )

      threep_tf_list <- as.integer( s[3:length(s)] )

      ts_tag <- "ts"
      for ( j in 1:length(threep_tf_list) ) {
        if ( j == 1 ) {
          ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
        } else {
          ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
        }
      }

      threep_tf_range <- c( threep_tf_list )

      cat( "# [fit_sequence_disc] threep_tf_range = ", threep_tf_range, "\n" )
    
      tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )
      message ( "# [fit_sequence_disc] tag = ", tag )

      # check existence of res file

      file_res <- paste ( "fit.", tag, ".res", sep="" )
      if ( file.exists ( file_res ) ) {
        message ( "# [fit_sequence_disc]  SKIP  ", tag )
        next
      } else {
        message ( "# [fit_sequence_disc]  START ", tag )
      }

      message ("# [fit_sequence_disc] threep_tf_range = ", formatC(threep_tf_range , width=4, format="d" ))
      message ("# [fit_sequence_disc] threep_tc_range = ", formatC(threep_tc_range , width=4, format="d" ))
      message ("# [fit_sequence_disc] twop_tf_range   = ", formatC(twop_tf_range , width=4, format="d" ))

      r <- run_min ( ens      = ens,
                 obs          = obs, 
                 nconf        = nconf,
                 nsrc         = nsrc, 
                 TT           = TT,
                 p            = mom, 
                 operator     = operator,
                 path_to_data = path_to_data,
                 lvl          = lvl,
                 par0         = par0,
                 nsample      = nsample, 
                 seed         = seed,
                 output_tag   = tag ,
                 twop_tf_range, threep_tf_range, threep_tc_range,
                 threep_prefix = threep_prefix,
                 twop_prefix   = twop_prefix,
                 twop_col      = twop_col,
                 threep_col    = threep_col
      )

    }
  }  else {

    for ( itf in 1:(length(threep_tf_list)) ) 
    {
      for ( ktf in itf:length(threep_tf_list) ) 
      {

        ts_tag <- "ts"
        for ( j in itf:ktf ) {
          if ( j == itf ) {
            ts_tag <- paste ( ts_tag, threep_tf_list[j], sep="" )
          } else {
            ts_tag <- paste ( ts_tag, "_", threep_tf_list[j], sep="" )
          }
        }

        threep_tf_range <- c( threep_tf_list[itf:ktf] )
        cat( "# [fit_sequence_disc] itf = ", itf, ", ktf = ", ktf, ", threep_tf_range = ", threep_tf_range, "\n" )
    
        threep_tc_max <- min( ( threep_tf_range * 3 ) %/% 8 )

        for ( tc in 0:threep_tc_max )
        {

          threep_tc_range <- c(-tc, tc )

          tag <- paste( ens, ".", obs, ".", ts_tag, ".tc", tc, ".tf", twop_tf_range[1], "_", twop_tf_range[2], sep="" )
          message ( "# [fit_sequence_disc] tag = ", tag )

          # check existence of res file

          file_res <- paste ( "fit.", tag, ".res", sep="" )
          if ( file.exists ( file_res ) ) {
            message ( "# [fit_sequence_disc]  SKIP  ", tag )
            next
          } else {
            message ( "# [fit_sequence_disc]  START ", tag )
          }

          message ("# [fit_sequence_disc] threep_tf_range = ", formatC(threep_tf_range , width=4, format="d" ))
          message ("# [fit_sequence_disc] threep_tc_range = ", formatC(threep_tc_range , width=4, format="d" ))
          message ("# [fit_sequence_disc] twop_tf_range   = ", formatC(twop_tf_range , width=4, format="d" ))

          r <- run_min ( ens      = ens,
                 obs          = obs, 
                 nconf        = nconf,
                 nsrc         = nsrc, 
                 TT           = TT,
                 p            = mom, 
                 operator     = operator,
                 path_to_data = path_to_data,
                 lvl          = lvl,
                 par0         = par0,
                 nsample      = nsample, 
                 seed         = seed,
                 output_tag   = tag ,
                 twop_tf_range, threep_tf_range, threep_tc_range,
                 threep_prefix = threep_prefix,
                 twop_prefix   = twop_prefix,
                 twop_col      = twop_col,
                 threep_col    = threep_col
          )

        }  # end of loop on tc
      }  # end of loop on ktf
    }  # end of loop on itf
  }

}  # end of fit_sequence_xq_disc
