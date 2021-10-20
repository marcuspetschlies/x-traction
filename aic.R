#####################################################################
# read data
#####################################################################
read_from_file <- function ( f ) {
  if ( missing ( f ) ) stop ( "[read_from_file] missing parameter values" )
  message ( "# [read_from_file] reading from file ", f )
  return ( read.table ( f ) )
}  # end of read_from_file

#####################################################################
# stats
#####################################################################
stats <- function ( s, t, f, a=T ) {

  if ( missing ( s ) || missing( t ) ) stop ( "[stats] missing argument values" )
  ms <- mean ( s )
  es <- sd   ( s )

  if ( !missing(f) ) {
    cat ( formatC(t, width=-40, format="s"),
         " & ", 
         formatC( ms, width=16, digits=4, format="f" ), 
         " & ", 
         formatC( es, width=16, digits=4, format="f" ),
         " \\\\ ", 
         "\n", file=f, append=a )
  }
  return ( invisible( list ( t=t, m=ms, e=es ) ) )
}  # end of stats

#####################################################################
# weighted sum of normal distributions for uniroot
#####################################################################
wdf <- function(x, w=w, m=m, s=s, a=0 ) {
  return( sum(w * pnorm(x,mean=m,sd=s) ) - a )
}

#####################################################################
# combine RC-factor weighted sample data
#####################################################################
aic_stats <- function(workpath="/data/nf211/x/R/" 
#                        nsample = 1000 , seed 
                       ) {

#  if ( missing ( seed ) ) stop( "[aic_stats] need seed value" )

  prefix_list <- c( 
                   # "xq-conn/lvl0",
                   # "xq-disc/lvl0",
                   # "xq-disc-strange/lvl0",
                   # "xq-disc-charm/lvl0",
                   # "xg-disc/nstout10/lvl0" 
                   "xg-disc/nstout0/lvl0",
                   "xg-disc/nstout1/lvl0",
                   "xg-disc/nstout2/lvl0",
                   "xg-disc/nstout4/lvl0",
                   "xg-disc/nstout8/lvl0",
                   "xg-disc/nstout10/lvl0"
  )

  #####################################################################
  # read data into list v
  #####################################################################
  f_list <- c(
    # "fit.cB211.072.64.xq-conn",
    # "fit.cB211.072.64.xq-disc",
    # "fit.cB211.072.64.xq-disc-strange",
    # "fit.cB211.072.64.xq-disc-charm",
    # "fit.cB211.072.64.xg-disc"
    "fit.cB211.072.64.xg-disc.nstout0.rectangle",
    "fit.cB211.072.64.xg-disc.nstout1.rectangle",
    "fit.cB211.072.64.xg-disc.nstout2.rectangle",
    "fit.cB211.072.64.xg-disc.nstout4.rectangle",
    "fit.cB211.072.64.xg-disc.nstout8.rectangle",
    "fit.cB211.072.64.xg-disc.nstout10.rectangle"
  )

  v <- list()

  for ( i in 1:length(f_list) ) {
    f <- paste( prefix_list[i], "/", f_list[i], ".aic", sep="") 
    v[[i]] <- read_from_file ( f )
  }

  #####################################################################
  # build the sampling
  #####################################################################

  for ( i in 1:length(f_list) ) 
  {

    message( "# [aic_stats] file = ", f_list[i] )
    d <- v[[i]]

    n <- length(d$V1)

    w <- d$V10 / sum( d$V10 )
    m <- numeric()
    s <- numeric()

    #####################################################################
    #
    #####################################################################
#    set.seed ( seed )

    #####################################################################
    #
    #####################################################################
    output_file <- paste( f_list[i], ".aic_stats", sep="" )
    cat( "# ", date(), "\n", file=output_file, append=F)
#    cat( "# nsample  ", nsample , "\n", append=T, file=output_file )

    #####################################################################
    #
    #####################################################################
    for ( k in 1:dim(d)[1] )
    {
      lvl <- d[k,1]
      if ( lvl == 0 ) {
        m[k] <- d[k,15]
        s[k] <- d[k,16]
      }
    }  # end of loop file entries
  
    # return( list(m=m,s=s,w=w))
    lower <- min ( m ) - 5*max ( s )
    upper <- max ( m ) + 5*max ( s )

    u_m <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.50, s=s, w=w, m=m )
    u_l <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.16, s=s, w=w, m=m )
    u_u <- uniroot( wdf, lower=lower, upper=upper , check.conv = T, tol = 1.e-12, maxiter = 100000, a=0.84, s=s, w=w, m=m )


    # return ( list(h=h, l=idxl, u=idxu, m=idxm , a=a) )

    cat ( "# convergence ", u_m$f.root, u_u$f.root, u_l$f.root, "\n", file=output_file, append=TRUE ) 
    cat( 
        formatC( u_m$root, width=6, digits=4, format="f" ),
        formatC( u_m$root - u_l$root, width=6, digits=4, format="f" ),
        formatC( u_u$root - u_m$root, width=6, digits=4, format="f" ),
        "\n", file=output_file, append=T )


  }  # end of loop on files


  return ( NULL )
    
}  # end of aic_stats
