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
    res <- res +  p[3]^2 / (2 * p[4] ) * ( exp ( -p[3] * tf )  + exp ( -p[4] * ( TT - tf ) ) )
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
    res <- rep( p[1]^2 / ( 2 * p[2] ) * p[3] * exp ( -p[2] * tf ), times=length(tc) )
  }

  if ( lvl == 1 ) {
    res <- p[1]^2 / ( 2 * p[2] ) * p[5] * exp ( -p[2] * tf )
         + p[1]*p[3] / ( 4 * p[2] * p[4] ) * p[6] * ( exp ( -p[2] * (tf -tc) -p[4] * tc ) + exp ( -p[4] * (tf -tc) -p[2] * tc ) )
         + p[3]^2 / ( 4 * p[4]^2 ) * p[7] * exp( -p[4] * tf )
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
  res <- 0

  # str( tf )
  # str( tc )

  if ( lvl == 0 ) {
    res <- rep( -3./8. * p[1]^2 * p[3] * exp ( -p[2] * tf ), times=length(tc) )
  }

  if ( lvl == 1 ) {
    res <- -3./8. * p[1]^2 * p[5] * exp ( -p[2] * tf )
           + p[1]*p[3] / ( 4 * p[2] * p[4] ) * p[6] * ( exp ( -p[2] * (tf -tc) -p[4] * tc ) + exp ( -p[4] * (tf -tc) -p[2] * tc ) )
           + p[3]^2 / ( 4 * p[4]^2 ) * p[7] * exp( -p[4] * tf )
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

  a <- svd ( d$cov )
  covi <- a$u %*% diag( 1./a$d ) %*% t( a$u )

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

  b <- cbind( t( fitdata$twop ),   t( apply( fitdata$threep, c(3) , t ) ) )

  d[["cov"]] <- cov(b, b) / ddim[2]

  rm(b)

  # return ( d )
  # return( fchisq(p=par0, d=d, fitinfo ) )
  # return ( list(d=d, i=fitinfo) )

  #############################################################
  # fit on ORIGINAL data
  #############################################################
  uorig <- optim ( par, fchisq, gr = NULL, d=d, fitinfo =fitinfo )

  if ( uorig$convergence != 0 ) {
    stop( "[fminimize] optim orig exit status ", uorig$convergence )
  }

  #############################################################
  # fit on SAMPLED data
  #############################################################

  bs[["fitres"]] <- array ( bs$nsample, length(par) )
  bs[["chisq"]] <- numeric()

  for ( s in 1:bs$nsample ) {
  
    ds <- list()

    idx <- sample.int ( n=ddim[2], size=ddim[2], replace = TRUE )

    ds[["twop"]] <- apply ( fitdata$twop[,idx], c(1), mean )
  
    ds[["threep"]] <- apply ( fitdata$threep[,idx], c(1,2), mean )

    b <- cbind( t( fitdata$twop[,idx] ),   t( apply( fitdata$threep[,,idx], c(3) , t ) ) )

    ds[["cov"]] <- cov(b, b) / ddim[2]

    rm(b)

    u <- optim ( par, fchisq, gr = NULL, d=ds, fitinfo =fitinfo )

    if ( u$convergence != 0 ) {
      stop( "[fminimize] optim exit status ", u$convergence )
    }

    bs$fitres[s,] <- u$par
    bs$chisq[s]   <- u$value

  }  # end of bootstrap sampling


  #############################################################
  # analyse bootstrap
  #############################################################
  res <- list()
  res$par_value <- uorig$par
  res$chisq     <- uorig$value
  res$dof       <- ddim[1] - length(par)
  res$par_cov   <- cov( bs$fitres )
  res$orig      <- uorig
  res$bs        <- bs

  return ( res )

}  # end of fminimize



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
                     nsample = 1000 ) {

  twop_filename <- paste ( path_to_data, "/", ens, "/", obs, "/", "twop.pseudoscalar.orbit.PX", p[1], "_PY", p[2], "_PZ", p[3], ".re.corr", sep="" )
  if ( !file.exists( twop_filename) ) stop( "Could not find ", twop_filename )

  message( "# [run_min] reading 2pt from file ", twop_filename )
  twop_data <- apply( array( read.table ( twop_filename )$V2, dim=c(TT, nsrc, nconf) ), c(1,3), mean )

  # return ( twop_data)

  threep_data <- array( dim=c( length( threep_tf_range ), TT, nconf ) )

  for (idt in 1:length(threep_tf_range) )  {
    dt <- threep_tf_range[idt]
  
    threep_filename <- paste ( path_to_data, "/", ens, "/", obs, "/", "threep.conn.", operator, ".dtsnk", dt, ".PX", p[1], "_PY", p[2], "_PZ", p[3], ".re.corr", sep="" )
    if ( !file.exists( threep_filename) ) stop( "Could not find ", threep_filename )
    message( "# [run_min] reading 3pt from file ", threep_filename )

    threep_data[idt,,] <- apply( array( read.table ( threep_filename )$V2, dim=c(TT, nsrc, nconf) ), c(1,3), mean )
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
  fitinfo[["tf_range_twop"]]   <- ( twop_tf_range[1] : twop_tf_range[2] )
  fitinfo[["tf_range_threep"]] <- threep_tf_range
  fitinfo[["tc_range_threep"]] <- ( threep_tc_range[1] : threep_tc_range[2] )
  fitinfo[["ftwop"]]          <- ftwop
  if ( operator == "g4_D4" ) {
    fitinfo[["fthreep"]]       <- fthreep44
  } else if ( operator == "g4_Dk" ) {
    fitinfo[["fthreep"]]       <- fthreep4k
  }
  fitinfo[["lvl"]]             <- 0

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
  # call minizer
  #############################################################
  fitres <- fminimize ( fitparam, fitinfo, fitdata, bs )

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
    plt$data$twop [i,] <- c ( t, u$value, u$dvalue, u$dvalue, u$tauint, u$dtauint )
  }

  # threep
  plt$data$threep <- array ( dim=c( length(fitinfo$tf_range_threep), length( fitinfo$tc_range_threep ), 7 ) )
  for ( i in 1:length( fitinfo$tf_range_threep ) ) {
    tf <- fitinfo$tf_range_threep[i]

    for ( k in 1:length( fitinfo$tc_range_threep ) ) {
      
      tc <- fitinfo$tc_range_threep[k]

      u <- uwerrprimary ( fitdata$threep[i,k,] )
      plt$data$threep[i,k,] <- c ( tf, tc, u$value, u$dvalue, u$dvalue, u$tauint, u$dtauint )
    }
  }

  # fit function

  plt$fit <- list()

  # twop
  plt$fit$twop <- array (  dim=c( length( fitinfo$tf_range_twop ), 3 ) )

  for ( i in 1:length( fitinfo$tf_range_twop ) ) {
    tf <- fitinfo$tf_range_twop[i]

    u <- apply ( fitres$bs$fitres, c(1), fitinfo$ftwop, tf = tf, fitinfo = fitinfo )

    plt$fit$twop [i,] <- c ( tf, mean( u ), sqrt( var ( u ) ) )
  }

  # threep
  plt$fit$threep <- array (  dim=c( length( fitinfo$tf_range_threep ), length( fitinfo$tc_range_threep ), 4 ) )

  for ( i in 1:length( fitinfo$tf_range_threep ) ) {
    tf <- fitinfo$tf_range_threep[i]

    for ( k in 1:length( fitinfo$tc_range_threep ) ) {
      tc <- fitinfo$tc_range_threep[k]

      u <- apply ( fitres$bs$fitres, c(1), fitinfo$fthreep, tf = tf, tc = tc, fitinfo = fitinfo )

      plt$fit$threep [i,k,] <- c ( tf, tc, mean( u ), sqrt( var ( u ) ) )
    }
  }


  #############################################################
  # return data, info and res
  #############################################################
  return ( invisible(list(data=fitdata, info=fitinfo, res=fitres, plt=plt) ) )

}  # end of run_min
