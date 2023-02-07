p_list <- array( dim=c(6,3) )
p_list[1,] <- c(  0, 0, 1 )
p_list[2,] <- c(  0, 0,-1 )
p_list[3,] <- c(  0, 1, 0 )
p_list[4,] <- c(  0,-1, 0 )
p_list[5,] <- c(  1, 0, 0 )
p_list[6,] <- c( -1, 0, 0 )

######################################################
#
######################################################
xg_analyse <- function( dt, nstout=20, p_list, op="g4_Dk", fst="clover" , rstout=0.1290, bootR = 2000, LL, fbwd ) 
{
  if ( missing(dt) )     stop( "need dt input")
  if ( missing(p_list) ) stop( "need p_list input")
  if ( missing(fbwd) )   stop( "need fbwd input")

  p_num <- dim(p_list)[1]

  twop   <- NULL
  threep <- NULL

  for ( i in 1:p_num )
  {
    p_tag <- paste( "PX", p_list[i,1], "_PY", p_list[i,2], "_PZ", p_list[i,3], sep="" )

    f <- paste( "twop.gf5.gi5.", p_tag, ".corr", sep="")
 
    if ( !file.exists(f) ) stop( "cannot find file ", f )
    message("# [] reading from file ", f )
    d <- read.table ( f )
    TT <- max(d$V1)+1
    nmeas <- length(d$V1) %/% TT
    if ( missing ( LL ) ) LL <- TT %/% 2
    message( "T     = ", TT )
    message( "L     = ", LL )
    message( "nmeas = ", nmeas )

    if ( i == 1 )
    {
      twop <- array ( dim=c(p_num, TT, nmeas) )
    }

    ######################################################
    # twop is p x T x nmeas
    ######################################################
    twop[i,,] <- array( d$V2, dim=c(TT, nmeas) )

    f <- paste( "threep.", fst, ".nstout", nstout, "_", formatC(rstout, width=6,digits=4, format="f" ), ".", op, ".dtsnk", dt, ".", p_tag, ".", fbwd, ".corr", sep="")
    if ( !file.exists(f) ) stop( "cannot find file ", f )
    message("# [] reading from file ", f )
    d <- read.table ( f )
    if ( length(d$V1) != TT*nmeas ) stop( "inconsistent data size" )

    if ( i == 1 )
    {
      threep <- array ( dim=c(p_num, (dt+1), nmeas) )
    }

    norm <- 1.
    if ( op == "gi_Dk" ||  op == "g4_Dk" ) norm <- norm / sum( (p_list[i,] * 2*pi / LL )^2 )

    ######################################################
    # threep is p x T x nmeas
    ######################################################

    threep[i,,] <- array( d$V2, dim=c(TT, nmeas) )[1:(dt+1),] * norm

  }

  # return(list(twop=twop, threep=threep ) )

  ######################################################
  ######################################################

  ######################################################
  # average of ratios
  ######################################################
  avgrat <- list()
  avgrat$orig <- numeric()

  itsink <- c( NA, NA )

  if ( fbwd == "fwd" ) 
  {
    itsink <- c( (  dt + TT ) %% TT + 1,
                 (  dt + TT ) %% TT + 1 )
  } else if ( fbwd == "bwd" )
  {
    itsink <- c( ( -dt + TT ) %% TT + 1,
                 ( -dt + TT ) %% TT + 1 )
  } else if ( fbwd == "fbwd" )
  {
    itsink <- c( (  dt + TT ) %% TT + 1,
                 ( -dt + TT ) %% TT + 1 )
  }

  b <- apply ( 0.5 * ( twop[,itsink[1],] + twop[,itsink[2],] ) , c(1), mean )
  for (i in 1:(dt+1) )
  {
    avgrat$orig[i] <- mean ( apply ( threep[,i,], c(1), mean ) / b )
    
    avgrat$othreep[i] <- mean ( threep[,i,] )
  }

  # bootstrap sampling
  sid <- sample.int ( n=nmeas, size=bootR*nmeas, replace = TRUE )

  bthreep <- apply ( array( threep[,,sid], dim=c(p_num, dt+1, bootR, nmeas) ), c(1,2,3), mean )

  btwop   <- apply ( array( 0.5 * ( twop[,itsink[1],sid] + twop[,itsink[2],sid] ) , dim=c(p_num, bootR, nmeas) ), c(1,2), mean )

  #str(bthreep)
  #str(btwop)

  avgrat$bs <- array( dim=c(dt+1, 2 ) )
  avgrat$sthreep <- array( dim=c(dt+1, 2 ) )

  for ( i in 1:(dt+1) )
  {
   brat <- apply ( bthreep[,i,] / btwop, c(2), mean )
  
   avgrat$bs[i,] <- c( mean(brat), sqrt(var(brat) ) )
   
   avgrat$sthreep[i,] <- c( mean ( bthreep[,i,] ),  sqrt( var( apply ( bthreep[,i,], c(2), mean ) ) ) )

  }

  # return( avgrat )

  ######################################################
  # ratio of averages
  ######################################################
  ratavg <- list()
  ratavg$orig <- apply( threep, c(2), mean ) / mean( 0.5 * ( twop[,itsink[1],] + twop[,itsink[2],] ))
  
  ratavg$othreep <- apply( threep, c(2), mean )

  bthreep <- apply ( array( apply ( threep, c(2,3), mean )[,sid], dim=c(dt+1,bootR,nmeas) ), c(1,2), mean )
  d <-  apply ( twop, c(2,3), mean )
  btwop   <- apply ( array( 0.5 * ( d[itsink[1],sid] + d[itsink[2],sid] ), dim=c( bootR, nmeas ) ), c(1), mean )

  # str(bthreep)
  # str(btwop)


  ratavg$bs <- array( dim=c(dt+1,2) )
  ratavg$sthreep <- array( dim=c(dt+1,2) )
  for ( i in 1:(dt+1) )
  {
   brat <- bthreep[i,] / btwop
  
   ratavg$bs[i,] <- c( mean(brat), sqrt(var(brat) ) )
   
   ratavg$sthreep[i,] <- c( mean(bthreep[i,]), sqrt(var(bthreep[i,]) ) )

  }

  # return(ratavg)

  # write to file

  f <- paste( "ratavg.", fst, ".nstout", nstout, "_", rstout, ".", op, ".dtsnk", dt, ".", fbwd, ".stats", sep="")
  cat( "# ", date(), "\n", file=f, append=F )

  for ( i in 1:(dt+1) )
  {
    cat( formatC(i-1, width=3, format="d"),
        formatC( c( ratavg$orig[i], ratavg$bs[i,], ratavg$othreep[i], ratavg$sthreep[i,] ), width=25, digits=16, format="e" ),
        "\n", file=f, append=T )
  }

  f <- paste( "avgrat.", fst, ".nstout", nstout, "_", rstout, ".", op, ".dtsnk", dt, ".", fbwd, ".stats", sep="")
  cat( "# ", date(), "\n", file=f, append=F )

  for ( i in 1:(dt+1) )
  {
    cat( formatC(i-1, width=3, format="d"),
        formatC( c( avgrat$orig[i], avgrat$bs[i,], avgrat$othreep[i], avgrat$sthreep[i,] ), width=25, digits=16, format="e" ),
        "\n", file=f, append=T )
  }

  return(list(ra=ratavg, ar=avgrat) )
}
