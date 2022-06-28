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
aic_stats <- function( workpath="/data/nf211/x/R/" , lvl=0, 
                       stout_type=c( "clover", "rectangle" ),
                       stout_n=c(0,1,2,4,8,10,20), stout_r = 0.1290,
                       xg_prefix="xg-disc-kaon",
                       op,
                       ens="cB211.072.64",
                       xq_conn_prefix=NULL
                     ) {

#  if ( missing ( seed ) ) stop( "[aic_stats] need seed value" )

  prefix_list <- c( 
                   # "xq-conn/lvl0",
                   # "xq-disc/lvl0",
                   # "xq-disc-strange/lvl0",
                   # "xq-disc-charm/lvl0",
                   #####################################################################
                   # xq disc kaon
                   #####################################################################
                   #"xq-disc-kaon-ll/lvl0",
                   #"xq-disc-kaon-ss/lvl0",
                   #"xq-disc-kaon-cc/lvl0"
                   #####################################################################
                   # xg for nstout pion
                   #####################################################################
                   #"xg-disc/nstout0/lvl0",
                   #"xg-disc/nstout1/lvl0",
                   #"xg-disc/nstout2/lvl0",
                   #"xg-disc/nstout4/lvl0",
                   #"xg-disc/nstout8/lvl0",
                   #"xg-disc/nstout10/lvl0"
                   #####################################################################
                   # xg for nstout kaon
                   #####################################################################
                   #"xg-disc-kaon/nstout0/lvl0",
                   #"xg-disc-kaon/nstout1/lvl0",
                   #"xg-disc-kaon/nstout2/lvl0",
                   #"xg-disc-kaon/nstout4/lvl0",
                   #"xg-disc-kaon/nstout8/lvl0",
                   #"xg-disc-kaon/nstout10/lvl0"
  )


  if ( length(stout_n) != 0 && length(stout_type) != 0 ) {
    for ( t in stout_type ) {
      for ( n in stout_n ) {
        p <- paste( xg_prefix, "/nstout", n, "/lvl", lvl, sep="" )
        prefix_list <- c( prefix_list, p )
      }
    }
  }
  
  if ( length(xq_conn_prefix) != 0 ) {
    for ( s in xq_conn_prefix ) {
      if ( !missing(op) ) {
        p <- paste( s, "-", op,  "/lvl", lvl, sep="" )
      } else {
        p <- paste( s, "/lvl", lvl, sep="" )
      }
      prefix_list <- c( prefix_list, p )
    }
  }

  #####################################################################
  # read data into list v
  #####################################################################
  f_list <- c(
              #####################################################################
              # xq conn pion
              #####################################################################
    # "fit.cB211.072.64.xq-conn",
              #####################################################################
              # xq disc ll, ss, cc pion
              #####################################################################
    # "fit.cB211.072.64.xq-disc",
    # "fit.cB211.072.64.xq-disc-strange",
    # "fit.cB211.072.64.xq-disc-charm",
              #####################################################################
              # xq disc ll, ss, cc kaon
              #####################################################################
     #"fit.cB211.072.64.xq-disc-kaon-ll",
     #"fit.cB211.072.64.xq-disc-kaon-ss",
     #"fit.cB211.072.64.xq-disc-kaon-cc"
              #####################################################################
              # xq disc pion
              #####################################################################
    # "fit.cB211.072.64.xg-disc"
              #####################################################################
              # xg for nstout pion
              #####################################################################
    #"fit.cB211.072.64.xg-disc.nstout0.rectangle",
    #"fit.cB211.072.64.xg-disc.nstout1.rectangle",
    #"fit.cB211.072.64.xg-disc.nstout2.rectangle",
    #"fit.cB211.072.64.xg-disc.nstout4.rectangle",
    #"fit.cB211.072.64.xg-disc.nstout8.rectangle",
    #"fit.cB211.072.64.xg-disc.nstout10.rectangle"
              #####################################################################
              # xg for nstout kaon
              #####################################################################
    #"fit.cB211.072.64.xg-disc-kaon.nstout0.rectangle",
    #"fit.cB211.072.64.xg-disc-kaon.nstout1.rectangle",
    #"fit.cB211.072.64.xg-disc-kaon.nstout2.rectangle",
    #"fit.cB211.072.64.xg-disc-kaon.nstout4.rectangle",
    #"fit.cB211.072.64.xg-disc-kaon.nstout8.rectangle",
    #"fit.cB211.072.64.xg-disc-kaon.nstout10.rectangle"
  )

  if ( length(stout_n) != 0 && length(stout_type) != 0 ) {
    for ( t in stout_type ) {
      for ( n in stout_n ) {
        f <- paste( "fit.", ens, ".", xg_prefix, "-", t, ".nstout", n, "_", formatC(stout_r, width=6, digits=4, format="f"), sep="" )
        f_list <- c( f_list, f )
      }
    }
  }
  
  if ( length(xq_conn_prefix) != 0 ) {
    for ( s in xq_conn_prefix ) {
        f <- paste( "fit.", ens, ".", s, sep="" )
        f_list <- c( f_list, f )
    }
  }

  # return ( list( f=f_list, p=prefix_list ))

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
