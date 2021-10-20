# non-singlet
Z_qq1  <- c(1.151, 0.004 )
Z_qq2  <- c(1.160, 0.003 )

# singlet
Z_qqs1 <- c(1.161, 0.024 )
Z_qqs2 <- c(1.163, 0.012 )

# gluon
Z_gg2 <- c(1.08, 0.17 )

# mixing coefficients
#Z_gq1 <- c(  0.232, 0 )
#Z_gq2 <- c(  0.083, 0 )
#Z_qg1 <- c( -0.027, 0 )

#Z_gq1 <- c( 0, 0 )
#Z_gq2 <- c( 0, 0 )
#Z_qg1 <- c( 0, 0 )

#Z_gq1 <- c(  0.232, 0.232 )
#Z_gq2 <- c(  0.083, 0.083 )
#Z_qg1 <- c( -0.027, 0.027 )

Z_qg1 <- c(  0.232, 0 )
Z_qg2 <- c(  0.083, 0 )
Z_gq1 <- c( -0.027, 0 )
Z_gq2 <- c(  0.0  , 0 )



#####################################################################
# read data
#####################################################################
read_col_from_file <- function ( f, c ) {
  if ( missing ( f ) || missing ( c ) ) stop ( "[read_col_from_file] missing parameter values" )
  message ( "# [read_col_from_file] reading col ", c, " from file ", f )
  return ( read.table ( f )[,c] )
}  # end of read_col_from_file

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
# combine RC-factor weighted sample data
#####################################################################
combine_renormalized <- function(workpath="/data/nf211/x/R/" ) {

  col <- 3
  nf <- 2+1+1

  xL_f <- "xq-conn/lvl0/aic/fit.cB211.072.64.xq-conn.ts56_64_72_80.tc4.tf32_64.bs"

  xl_f <- "xq-disc/lvl0/aic/fit.cB211.072.64.xq-disc.ts12_16_20_24_28_32.tc4.tf32_64.bs"

  xs_f <- "xq-disc-strange/lvl0/aic/fit.cB211.072.64.xq-disc-strange.ts12_16_20_24_28_32.tc4.tf32_64.bs"

  xc_f <- "xq-disc-charm/lvl0/aic/fit.cB211.072.64.xq-disc-charm.ts12_16_20_24_28_32.tc4.tf32_64.bs"

  xg_f <- "xg-disc/nstout10/lvl0/aic/fit.cB211.072.64.xg-disc.ts20_24_28_32.tc7.tf32_64.bs"

  xL_s <- read_col_from_file ( xL_f, col ) / 2 # divide by 2 for single light quark contribution, L <- ( U + D ) / 2
  xl_s <- read_col_from_file ( xl_f, col ) / 2 # divide by 2 for single light quark contribution, l <- ( u + d ) / 2
  xs_s <- read_col_from_file ( xs_f, col )
  xc_s <- read_col_from_file ( xc_f, col )
  xg_s <- read_col_from_file ( xg_f, col )

  ns <- length ( xL_s )
  if ( length ( xl_s ) != ns || length ( xs_s ) != ns || length ( xc_s ) != ns || length ( xg_s ) != ns ) stop ( "different sample sizes" )

  #####################################################################
  # difference singlet - non-singlet
  #####################################################################

  dZ_qq1 <- Z_qqs1 - Z_qq1
  dZ_qq2 <- Z_qqs2 - Z_qq2

  #####################################################################
  # RC factor samples
  #####################################################################
  sZ_qq1  <- rnorm( n = ns, mean = Z_qq1[1], sd = Z_qq1[2] )

  sZ_qq2  <- rnorm( n = ns, mean = Z_qq2[1], sd = Z_qq2[2] )

  sdZ_qq1 <- rnorm( n = ns, mean = dZ_qq1[1], sd = dZ_qq1[2] )
 
  sdZ_qq2 <- rnorm( n = ns, mean = dZ_qq2[1], sd = dZ_qq2[2] )
 
  sZ_gg2  <- rnorm( n = ns, mean = Z_gg2[1], sd = Z_gg2[2] )

  sZ_gq1  <- rnorm( n = ns, mean = Z_gq1[1], sd = Z_gq1[2] )

  sZ_gq2  <- rnorm( n = ns, mean = Z_gq2[1], sd = Z_gq2[2] )

  sZ_qg1  <- rnorm( n = ns, mean = Z_qg1[1], sd = Z_qg1[2] )

  sZ_qg2  <- rnorm( n = ns, mean = Z_qg2[1], sd = Z_qg2[2] )

  #####################################################################
  # combine
  #####################################################################

#  # light
#  rxl_s <- ( sZ_qq1 * xL_s + sZ_qq2 * xl_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg1 * xg_s / nf
#  
#  # strange
#  rxs_s <- (                 sZ_qq2 * xs_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg1 * xg_s / nf
#  
#  # charm
#  rxc_s <- (                 sZ_qq2 * xc_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg1 * xg_s / nf
#  
#  # gluon
#  rxg_s <- sZ_gq1 * 2 * xL_s + sZ_gq2 * ( 2 * xl_s  + xs_s + xc_s )  + sZ_gg2 * xg_s


  # light
  rxl_s <- ( sZ_qq1 * xL_s + sZ_qq2 * xl_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg2 * xg_s / nf
  
  # strange
  rxs_s <- (                 sZ_qq2 * xs_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg2 * xg_s / nf
  
  # charm
  rxc_s <- (                 sZ_qq2 * xc_s ) + ( sdZ_qq1 * 2 * xL_s + sdZ_qq2 * ( 2 * xl_s + xs_s + xc_s ) ) / nf + sZ_qg2 * xg_s / nf
  
  # gluon
  rxg_s <- sZ_gq1 * 2 * xL_s + sZ_gq2 * ( 2 * xl_s  + xs_s + xc_s )  + sZ_gg2 * xg_s

  #####################################################################
  # mean and variance
  #####################################################################
  output_file <- "xq-renormalized.tab"
  cat( "# ", date(), "\n", file=output_file, append=F )

  stats ( s = rxl_s, t = "l", f = output_file )

  stats ( s = rxs_s, t = "s", f = output_file )

  stats ( s = rxc_s, t = "c", f = output_file )

  stats ( s = rxg_s, t = "g", f = output_file )

  stats ( s = 2*rxl_s+rxs_s, t = "u+d+s", f = output_file )

  stats ( s = 2*rxl_s+rxs_s+rxc_s, t = "u+d+s+c", f = output_file )

  stats ( s = 2*rxl_s+rxs_s+rxc_s+rxg_s, t = "u+d+s+c+g", f = output_file )

}  # end of combine_renormalized
