p_list <- list()

p_list[[1]] <- array( dim=c(6,3) )
p_list[[1]][1,] <- c(  0, 0, 1 )
p_list[[1]][2,] <- c(  0, 0,-1 )
p_list[[1]][3,] <- c(  0, 1, 0 )
p_list[[1]][4,] <- c(  0,-1, 0 )
p_list[[1]][5,] <- c(  1, 0, 0 )
p_list[[1]][6,] <- c( -1, 0, 0 )

p_list[[2]] <- array( dim=c(12,3) )
p_list[[2]][1,] <- c(  0, 1, 1 )
p_list[[2]][2,] <- c(  0, 1,-1 )
p_list[[2]][3,] <- c(  0,-1, 1 )
p_list[[2]][4,] <- c(  0,-1,-1 )
p_list[[2]][5,] <- c(  1, 0, 1 )
p_list[[2]][6,] <- c(  1, 0,-1 )
p_list[[2]][7,] <- c( -1, 0, 1 )
p_list[[2]][8,] <- c( -1, 0,-1 )
p_list[[2]][9,] <- c(  1, 1, 0 )
p_list[[2]][10,] <- c(  1,-1, 0 )
p_list[[2]][11,] <- c( -1, 1, 0 )
p_list[[2]][12,] <- c( -1,-1, 0 )

p_list[[3]] <- array( dim=c(8,3) )
p_list[[3]][1,] <- c(  1, 1, 1 )
p_list[[3]][2,] <- c(  1, 1,-1 )
p_list[[3]][3,] <- c(  1,-1, 1 )
p_list[[3]][4,] <- c(  1,-1,-1 )
p_list[[3]][5,] <- c( -1, 1, 1 )
p_list[[3]][6,] <- c( -1, 1,-1 )
p_list[[3]][7,] <- c( -1,-1, 1 )
p_list[[3]][8,] <- c( -1,-1,-1 )


######################################################
#
######################################################
xq_conn_ratio_analyse <- function( dt, p_orbit, op="gddd0123", bootR = 2000, LL,
                                  twop_flavor_tag_list=c("l-gf-l-gi"),
                                  threep_flavor_tag_list=c("DDDl-gc-ll-gi"), threep_reim="im" ) 
{
  if ( missing(dt) )      stop( "need dt input")
  if ( missing(p_orbit) ) stop( "need p_orbit input")

  twop   <- NULL
  threep <- NULL

  TT <- NULL
  nmeas <- NULL
    
  p_orbit_tag <- paste( "PX", p_orbit[1], "_PY", p_orbit[2], "_PZ", p_orbit[3], sep="" )

  ######################################################
  # read twop data
  ######################################################
  for ( fl in twop_flavor_tag_list )
  {

    f <- paste( "twop.", fl, ".pseudoscalar.orbit.", p_orbit_tag, ".re.corr", sep="")
 
    if ( !file.exists(f) ) stop( "cannot find file ", f )
    message("# [xq_conn_ratio_analyse] reading from file ", f )
    d <- read.table ( f )
    if ( is.null ( TT ) )
    {
      TT <- max(d$V1)+1
      if ( missing ( LL ) ) LL <- TT %/% 2
      message( "# [xq_conn_ratio_analyse] T     = ", TT )
      message( "# [xq_conn_ratio_analyse] L     = ", LL )
    }
    if ( is.null ( nmeas ) )
    {
      nmeas <- length(d$V1) %/% TT
      message( "# [xq_conn_ratio_analyse] nmeas = ", nmeas )
    }
    if ( length(d$V1) != TT*nmeas ) stop( "inconsistent data size" )

    if ( any( is.null ( twop) ) )
    {
      twop <- array ( 0, dim=c(TT, nmeas) )
    }

    ######################################################
    # twop is T x nmeas
    ######################################################
    twop <- twop + array( d$V2, dim=c(TT, nmeas) )

  }
  twop <- twop / length( twop_flavor_tag_list )

  ######################################################
  # read threep data
  ######################################################
  for ( fl in threep_flavor_tag_list )
  {
    f <- paste( "threep_orbit.", fl, ".conn.", op, ".dtsnk", dt, ".", p_orbit_tag, ".", threep_reim, ".srcavg.corr", sep="")
    if ( !file.exists(f) ) stop( "cannot find file ", f )
    message("# [xq_conn_ratio_analyse] reading from file ", f )
    d <- read.table ( f )
    if ( length(d$V1) != TT*nmeas ) stop( "inconsistent data size" )

    if ( any ( is.null(threep)) )
    {
      threep <- array ( 0, dim=c( TT, nmeas) )
    }

    norm <- 1.

    ######################################################
    # threep is T x nmeas
    ######################################################

    threep  <- threep + array( d$V2, dim=c(TT, nmeas) ) * norm
  }
  threep <- threep / length ( threep_flavor_tag_list )

  # return(list(twop=twop, threep=threep ) )

  ######################################################
  # build ratio of orbit-averaged twop and threep data sets
  ######################################################
  ratavg <- list()

  itsink <- ( dt + TT ) %% TT + 1

  ratavg$ratio  <- list()
  ratavg$threep <- list()
  ratavg$twop   <- list()

  ratavg$ratio$orig <- apply( threep, c(1), mean ) / mean( twop[itsink,] )
  
  ratavg$threep$orig <- apply( threep, c(1), mean )
  ratavg$twop$orig   <- apply( twop, c(1), mean )

  # bootstrap sampling
  sid <- sample.int ( n=nmeas, size=bootR*nmeas, replace = TRUE )

  bthreep <- apply ( array( threep[,sid], dim=c( TT, bootR, nmeas ) ), c(1,3), mean )
  btwop   <- apply ( array( twop[,sid],   dim=c( TT, bootR, nmeas ) ), c(1,3), mean )

  # str(bthreep)
  # str(btwop)

  ratavg$threep$bs <- array ( c( apply( bthreep, c(1), mean ), sqrt( apply ( bthreep, c(1), var ) ) ), dim=c(TT,2) )
  ratavg$twop$bs <- array ( c( apply( btwop, c(1), mean ), sqrt( apply ( btwop, c(1), var ) ) ), dim=c(TT,2) )


  ratavg$ratio$bs <- array( dim=c(TT,2) )
  for ( it in 1:TT )
  {
    brat <- bthreep[it,] / btwop[itsink,]
    ratavg$ratio$bs[it,] <- c( mean (brat), sqrt(var(brat)) )
  }

  # return(ratavg)

  # write to file

  fout <- paste( "ratio_orbit.", op, ".dtsnk", dt, ".", p_orbit_tag, ".", threep_reim, ".stats", sep="")

  cat( "# ", date(), "\n",
      "# nmeas = ", nmeas, "\n",
      "# bootR = ", bootR, "\n",
      file=fout, append=F )

  for ( i in 1:TT )
  {
    cat( formatC(i-1, width=3, format="d"),
        formatC( c( ratavg$ratio$orig[i], ratavg$ratio$bs[i,] ), width=25, digits=16, format="e" ),
        formatC( c( ratavg$threep$orig[i], ratavg$threep$bs[i,] ), width=25, digits=16, format="e" ),
        formatC( c( ratavg$twop$orig[i], ratavg$twop$bs[i,] ), width=25, digits=16, format="e" ),
        "\n", sep="", file=fout, append=T )
  }

  return(invisible(ratavg))
}

######################################################
#
######################################################
run_xq_conn_ratio_analyse <- function( )
{

  dt_list     <- c(16,20,24,36)
  threep_reim <- "im" 
                 
  p_orbit <- c(1,1,1)
  op      <- "gddd0123"
  bootR   <- 2000
              
  twop_flavor_tag_list   <- c("l-gf-s-gi", "s-gf-l-gi" )
# threep_flavor_tag_list <- c("DDDl-gc-llgi")
# threep_flavor_tag_list <- c("DDDl-gc-ls-gi")
  threep_flavor_tag_list <- c("DDDs-gc-sl-gi")

  for ( dt in dt_list )
  {
    r <- xq_conn_ratio_analyse ( dt, p_orbit=p_orbit, op=op, bootR = bootR,
                                  twop_flavor_tag_list=twop_flavor_tag_list,
                                  threep_flavor_tag_list=threep_flavor_tag_list )
  }
}
