

rm(list=ls())

##Start Script


castPricingDate             <- function( df, dateStr ){
  
  dates                 <- df[ , c( dateStr) ]
  PricingDate           <- chron( as.character( dates ) )
  df[ , c( dateStr) ]   <- PricingDate
  
  df
}

prepareInput                <- function( dfReturns, pReturn, pPerRiskUnit, pSectors, pFlagSimpleName, pOmitNA ){
  
  dfOut               <- data.frame( "PricingDate" = dfReturns$PricingDate )
  
  for( i in 1:length( pSectors ) ){
    
    thisSector <- as.character( pSectors[ i ] )
    print(thisSector)  
    if( pPerRiskUnit ){
      thisString <- paste( c( pReturn, "_", thisSector, "_perRisk" ), collapse = "" )
    } else {
      thisString <- paste( c( pReturn, "_", thisSector), collapse = "" )    
    }
    
    thisVal        <- dfReturns[ , c( thisString ) ]
    
    dfTmp          <- data.frame( "PricingDate" = dfReturns$PricingDate, "thisVal" = thisVal )
    dfOut          <- merge( dfOut, dfTmp, by = c( "PricingDate" ) )
    
    if( pFlagSimpleName ){
      colnames( dfOut )[ grep( "thisVal", colnames( dfOut ), fixed = TRUE ) ]  <- thisSector
    } else {
      colnames( dfOut )[ grep( "thisVal", colnames( dfOut ), fixed = TRUE ) ]  <- thisString      
    }
    
  }
  
  if( pOmitNA ){
    dfOut <- na.omit( dfOut )    
  }
  
  dfOut
  
}


summarizeThresholds <- function( dfReturnData, pPairs, pThresholds ){
  
  count <- 0
  nPairs <- length( pPairs[[ 1 ]])
  firstLeg <- pPairs[[1]]
  secondLeg <- pPairs[[2]]
  legs   <- unique( c( firstLeg, secondLeg))
  nSectors <- length( legs )
  for( i in 1:nSectors ){
    
    thisLeg <- as.character( legs[ i ] )
    
    
    thisVal  <- na.omit( dfReturnData[ , c( thisLeg )] )  
    cdf <- ecdf( thisVal )
    
    nobs <- length( thisVal)
    
    for( j in 1:length( pThresholds ) ){
      
      thisTheta <- pThresholds[ j ]
      if( thisTheta < 0 ){
        idx <- thisVal <= thisTheta 
      } else {
        idx <- thisVal >= thisTheta 
      }
      nobsTheta <- sum(idx)
      pct <- cdf( thisTheta )
      dfOutTmp <- data.frame( "theta" = thisTheta, "sector" = thisLeg, "nobs" = nobsTheta, "pctile" = pct )
      if(count == 0 ){
        dfOut <- dfOutTmp
      } else {
        dfOut <- rbind( dfOut, dfOutTmp )
      }
      count<-count+1
    }
  }
  dfOut
}


replaceWithTheta  <- function( rtrn, theta ){
  
  idx        <- rtrn > theta
  rtrn[idx ] <- theta
  
  rtrn
}


calcTheta         <- function( rtrn, pctile ){
  
  theta <- quantile( rtrn, pctile )[[1]]
  theta
  
}


parameterOverride <- function( R1, R2, params, parameter ){
  
  objValue            <- likelihood(R1, R2, params[1], params[2], params[3], params[4], params[5], params[6], params[7])
  idxMin              <- tail( cumsum( ( objValue == min( objValue ) ) * seq( 1,length(objValue), 1 ) ), 1 )
  r1Min               <- R1[idxMin]
  r2Min               <- R2[idxMin]
  
  objValueMin         <- likelihood(r1Min, r2Min, params[1], params[2], params[3], params[4], params[5], params[6], params[7])
  paramIter           <- params
  
  thisFraction        <- 0.2
  while( objValueMin < 0 ){
    thisFraction <- thisFraction + 0.1   
    paramIter[parameter]<-params[parameter]*thisFraction
    objValueMin         <-likelihood(r1Min, r2Min, paramIter[1], paramIter[2], paramIter[3], paramIter[4], paramIter[5], paramIter[6], paramIter[7])
  }
  paramIter[ parameter ]  
}


errorHandlingMLE  <- function(params){
  tryCatch(
    results            <- optim( params, fn = LogLikelihood, hessian = TRUE ),    
    error=function(e) FALSE
  )
}

returnSeriesAdjustment = function(dfReturnDataPair){
  
  # create monthly dummies
  month = months(dfReturnDataPair$PricingDate)
  jan = 1 * (month == "Jan")
  feb = 1 * (month == "Feb")
  mar = 1 * (month == "Mar")
  apr = 1 * (month == "Apr")
  may = 1 * (month == "May")
  jun = 1 * (month == "Jun")
  jul = 1 * (month == "Jul")
  aug = 1 * (month == "Aug")
  sep = 1 * (month == "Sep")
  oct = 1 * (month == "Oct")
  nov = 1 * (month == "Nov")
  
  # create time trend
  timeTrend = seq(1,nrow(dfReturnDataPair), 1)
  timeTrendSquare = timeTrend^2
  
  dfReturnDataTmp = dfReturnDataPair
  dfReturnDataTmp$PricingDate = NULL
  
  cols = colnames(dfReturnDataTmp)
  for(j in 1:length(cols)){
    
    # regression step 1
    thisYstep1 = dfReturnDataTmp[, c(as.character(cols[j]))]
    regModelStep1 = lm(thisYstep1 ~ jan + feb + mar + apr + may + jun + jul + aug + sep + oct + nov + timeTrend + timeTrendSquare)
    
    # regression step 2
    residModelStep1 = as.array( regModelStep1$residuals)
    thisYstep2 = log(residModelStep1^2)
    regModelStep2 = lm(thisYstep2 ~ jan + feb + mar + apr + may + jun + jul + aug + sep + oct + nov + timeTrend + timeTrendSquare)
    fittedValStep2 = regModelStep2$fitted.values
    adjResid = regModelStep1 / exp(fittedValStep2 / 2)
    bCoeff = sqrt(var(thisYstep1)/var(adjResid))
    aCoeff = mean(thisYstep1) - bCoeff * mean(adjResid)
    
    thisYadj = aCoeff + bCoeff * adjResid
    
    dfReturnDataAdjTmp = data.frame("thisVal" = thisValAdj)
    colnames(dfReturnDataAdjTmp)[grep("thisVal", colnames(dfReturnDataAdjTmp), fixed = TRUE)] = as.character(cols[j])
    
    if(j == 1){
      dfReturnDataAdj = dfReturnDataAdjTmp
    } else {
      dfReturnDataAdj = cbind(dfReturnDataAdj, dfReturnDataAdjTmp)
    }
  }
  
  dfReturnDataAdj$PricingDate = dfReturnDataPair$PricingDate
  dfReturnDataAdj = dfReturnDataAdj[, c("PricingDate", cols)]
  
  dfReturnDataAdj
}


# likelihood sub-functions
FofR1 <- function( r1, p1, epsilon1, sigma1 ){
  
  1 - p1 * ( ( 1 + epsilon1 * ( theta1 - r1 ) / sigma1 ) )^(-1/epsilon1) 
  
}

FofR2 <- function( r2, p2, epsilon2, sigma2 ){
  
  1 - p2 * ( ( 1 + epsilon2 * ( theta2 - r2 ) / sigma2 ) )^(-1/epsilon2) 
  
}

ZofR1 <- function( r1, p1, epsilon1, sigma1 ){
  
  - 1 / ( log( FofR1( r1, p1, epsilon1, sigma1 )))
  
}

ZofR2 <- function( r2, p2, epsilon2, sigma2 ){
  
  - 1 / ( log( FofR2( r2, p2, epsilon2, sigma2 )))
  
}

expDependenceFct <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  exp( -( ZofR1(r1, p1, epsilon1, sigma1)^(-1/alpha) + ZofR2(r2, p2, epsilon2, sigma2)^(-1/alpha))^alpha )
  
}


K1 <- function( r1, p1, epsilon1, sigma1 ){
  
  - p1/sigma1 * ( (1 + epsilon1 * (theta1 - r1)/sigma1)^(-(1+epsilon1)/epsilon1) ) * ZofR1(r1, p1, epsilon1, sigma1)^2 * exp(1 / ZofR1( r1, p1, epsilon1, sigma1))
  
}


K2 <- function( r2, p2, epsilon2, sigma2 ){
  
  - p2/sigma2 * ( (1 + epsilon2 * (theta2 - r2)/sigma2)^(-(1+epsilon2)/epsilon2) ) * ZofR2(r2, p2, epsilon2, sigma2)^2 * exp(1 / ZofR2( r2, p2, epsilon2, sigma2))
  
}

deltaDdeltaR1 = function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha){
  
  ( ZofR1( r1, p1, epsilon1, sigma1)^(-1/alpha) + ZofR2(r2, p2, epsilon2, sigma2)^(-1/alpha) )^(alpha - 1) * ZofR1( r1, p1, epsilon1, sigma1)^(-(1/alpha) - 1)
  
}

deltaDdeltaR2 = function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha){
  
  ( ZofR1( r1, p1, epsilon1, sigma1)^(-1/alpha) + ZofR2(r2, p2, epsilon2, sigma2)^(-1/alpha) )^(alpha - 1) * ZofR2(r2, p2, epsilon2, sigma2)^(-(1/alpha) - 1)
  
}

delta2DdeltaR1deltaR2 = function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha){
  
  ((alpha-1)/alpha) * (( ZofR1( r1, p1, epsilon1, sigma1)^(-1/alpha) + ZofR2(r2, p2, epsilon2, sigma2)^(-1/alpha) )^(alpha - 2)) * ZofR1( r1, p1, epsilon1, sigma1)^(-(1/alpha) - 1) * ZofR2(r2, p2, epsilon2, sigma2)^(-(1/alpha) - 1)
  
}

indicator0_0 <- function( r1, r2 ){
  
  1 * ( ( r1 >= theta1 ) & ( r2 >= theta2 ) )
  
}


indicator0_1 <- function( r1, r2 ){
  
  1 * ( ( r1 >= theta1 ) & ( r2 < theta2 ) )
  
}


indicator1_0 <- function( r1, r2 ){
  
  1 * ( ( r1 < theta1 ) & ( r2 >= theta2 ) )
  
}

indicator1_1 <- function( r1, r2 ){
  
  1 * ( ( r1 < theta1 ) & ( r2 < theta2 ) )
  
}

likelihood0_0 <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  expDependenceFct( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) 
  
}

likelihood0_1 <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  -expDependenceFct( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * deltaDdeltaR2( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * K2( r2, p2, epsilon2, sigma2 )
  
}

likelihood1_0 <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  -expDependenceFct( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * deltaDdeltaR1( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * K1( r1, p1, epsilon1, sigma1 )
  
}


likelihood1_1 <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  expDependenceFct( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * ( deltaDdeltaR1( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) * deltaDdeltaR2( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ) - delta2DdeltaR1deltaR2( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha) ) * K1( r1, p1, epsilon1, sigma1 ) * K2( r2, p2, epsilon2, sigma2 )
  
}

likelihood <- function( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha ){
  
  indicator0_0(r1, r2) * likelihood0_0( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha) + indicator0_1(r1, r2)*likelihood0_1( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha)+indicator1_0(r1, r2)*likelihood1_0( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha)+indicator1_1(r1, r2)*likelihood1_1( r1, r2, p1, epsilon1, sigma1, p2, epsilon2, sigma2, alpha)
  
}

LogLikelihood <- function( param ){
  
  -sum( log( Vectorize( likelihood) ( R1, R2, param[1], param[2], param[3], param[4], param[5], param[6], param[7] )))
  
}



library( chron )
library( evir )



# Set Parameters for input/output data
# file names
# SPECIFY OWN PATH
pPathData                       <- "C:\\Users\\SPECIFY OWN PATH"
pPathOut                    <- pPathData
pFileReturns                <- "R_returnDataPerUnitRisk.csv"        # return data measured as excess returns per unit of risk



# return parameters
pReturn                     <- "rx"   # rx := excess returns, trt := total returns
pPerRiskUnit                <- TRUE
pFlagSimplifyName           <- TRUE
pFlagAdjustReturns          <- FALSE
pSectors                    <- c( "usIGCredit", "mbs", "usHYCredit", "cmbs", "abs", "agency", "em", "ust10" )   
pNSectors                   <- length( pSectors )
pOmitNA                     <- FALSE
pThresholdsUpper            <- c( 0.35, 0.25, 0.15, 0.1 )


# correlation parameters
pPairs <- list(
  "firstLeg" =  c( "usIGCredit", "usIGCredit", "usIGCredit", "mbs" ),
  "secondLeg" = c( "usHYCredit", "cmbs"      , "mbs",        "ust10")
)


# load data
fileName                    <- paste( c( pPathData, "\\", pFileReturns  ), collapse ="" )
dfReturns                   <- read.csv( fileName )
dfReturns                   <- castPricingDate( dfReturns, "PricingDate" )
pPerRiskUnit                <- FALSE
dfReturnData                <- prepareInput( dfReturns, pReturn, pPerRiskUnit, pSectors, pFlagSimplifyName, pOmitNA )


# summary statistics
if(FALSE){
  dfThresholds                <- summarizeThresholds( dfReturnData, pPairs, pThresholdsUpper )
}





# RUN MLE FOR EXCEEDANCES

params <- c()
params[1] <- 0.3                          # p1
params[2] <- 0.6                          # epsilon1
params[3] <- 2.5                          # sigma1
params[4] <- 0.4                          # p2
params[5] <- 0.6                          # epsilon2
params[6] <- 2.5                          # sigma2
params[7] <- 0.9                          # alpha = sqrt(1-rho)
paramOverride <- 7                        # override alpha if log likelihood is negative

# create pairs and prepare input data for MLE
count = 0
for(i in 1:length(pPairs$firstLeg)){
  
  thisFirstLeg = as.character(pPairs$firstLeg[ i ] )
  thisSecondLeg = as.character(pPairs$secondLeg[ i ] )
  
  dfReturnDataTmp <- dfReturnData[ , c( "PricingDate", thisFirstLeg, thisSecondLeg ) ]
  dfReturnDataTmp <- na.omit( dfReturnDataTmp )
  if( pFlagAdjustReturns ){
    dfReturnDataTmp = returnSeriesAdjustment( dfReturnDataTmp )
  }
  rhoUncond          <- cor( dfReturnDataTmp[ , c(thisFirstLeg)], dfReturnDataTmp[ , c(thisSecondLeg)] )
  
  for(j in 1:length(pThresholdsUpper)){
    
    pctile = pThresholdsUpper[j]
    theta1 = calcTheta( dfReturnDataTmp[, c(thisFirstLeg)], pctile)
    theta2 = calcTheta( dfReturnDataTmp[, c(thisSecondLeg)], pctile)
    
    R1 = replaceWithTheta(dfReturnDataTmp[,c(thisFirstLeg)], theta1)
    R2 = replaceWithTheta(dfReturnDataTmp[,c(thisSecondLeg)], theta2)
    
    # get params
    if(TRUE){
      leg1gpd            <- gpd( -R1, -theta1 )
      params[1]          <- leg1gpd$p.less.thresh
      params[2]          <- leg1gpd$par.ests[[1]]
      params[3]          <- leg1gpd$par.ests[[2]]
      
      leg2gpd            <- gpd( -R2, -theta2 )
      params[4]          <- leg2gpd$p.less.thresh
      params[5]          <- leg2gpd$par.ests[[1]]
      params[6]          <- leg2gpd$par.ests[[2]]
      params[7]          <- sqrt(1-rhoUncond)
    }
    
    # run likelihood
    # check positivity of likelihood function
    paramOverrideFlag  <- TRUE
    if(TRUE){
      objValue            <- likelihood(R1, R2, params[1], params[2], params[3], params[4], params[5], params[6], params[7])
      if( sum( objValue < 0 ) > 0 ){
        params[ paramOverride ] <- parameterOverride( R1, R2, params, paramOverride )      
        paramOverrideFlag  <- TRUE
      }
      
    }
    
    print( paste( c( "pctile ", as.character(pctile), ":  ", thisFirstLeg, "  ", thisSecondLeg), collapse = "" ))
    
    out                <- errorHandlingMLE(params)
    if( length(out) > 1  ){
      
      results              <- out
      
      # collect parameter estimates
      dfOutTmp             <- data.frame( "firstLeg" = thisFirstLeg, "secondLeg" = thisSecondLeg, "pctile" = pctile, "theta1" = theta1, "theta2" = theta2, "paramOverride" = paramOverrideFlag, "rhoUncond" = rhoUncond, 
                                          "rho" = 1 - results$par[7]^2, 
                                          "beta" = (1-results$par[7]^2) * results$par[3] / results$par[6], 
                                          "p1"=results$par[1], 
                                          "epsilon1"=results$par[2], 
                                          "sigma1"=results$par[3], 
                                          "p2" = results$par[4],
                                          "epsilon2"=results$par[5], 
                                          "sigma2"=results$par[6], 
                                          "alpha"=results$par[7], 
                                          "iterations" = as.numeric(results$counts["function"]), 
                                          "logLikValue" = results$value, 
                                          "convergence" = results$convergence )
      
      # collect standard errors
      hessian              <- results$hessian
      hessianInv           <- solve( hessian )
      varHat               <- diag( hessianInv )
      
      nobs                 <- length( R1 )
      p1_se                <- sqrt( varHat[ 1 ] ) 
      epsilon1_se          <- sqrt( varHat[ 2 ] ) 
      sigma1_se            <- sqrt( varHat[ 3 ] ) 
      p2_se                <- sqrt( varHat[ 4 ] ) 
      epsilon2_se          <- sqrt( varHat[ 5 ] ) 
      sigma2_se            <- sqrt( varHat[ 6 ] ) 
      alpha_se             <- sqrt( varHat[ 7 ] ) 
      rho_se               <- sqrt( varHat[ 7 ] * 4 * results$par[ 7 ] ^ 2 )
      
      p1_tstat             <- results$par[ 1 ] /  p1_se             
      epsilon1_tstat       <- results$par[ 2 ] /  epsilon1_se          
      sigma1_tstat         <- results$par[ 3 ] /  sigma1_se
      p2_tstat             <- results$par[ 4 ] /  p2_se
      epsilon2_tstat       <- results$par[ 5 ] /  epsilon2_se
      sigma2_tstat         <- results$par[ 6 ] /  sigma2_se
      alpha_tstat          <- results$par[ 7 ] /  alpha_se
      rho_tstat            <- ( 1 - results$par[ 7 ]^ 2 ) / rho_se
      
      tstats               <- c(p1_tstat, epsilon1_tstat, sigma1_tstat, p2_tstat, epsilon2_tstat, sigma2_tstat, alpha_tstat, rho_tstat )
      pValues              <- c()
      for( i in 1:length( tstats )){
        thisTstat <- tstats[ i ]
        if( is.nan( thisTstat )){
          pValues <- c( pValues, 0.9 )
        } else {
          if( thisTstat < 0 ){
            pValues <- c( pValues, pnorm( thisTstat))
          } else {
            pValues <- c( pValues, 1 - pnorm( thisTstat))
          }
        }
      }
      
      dfOutTstatTmp <- data.frame(  "firstLeg" = thisFirstLeg, "secondLeg" = thisSecondLeg, "pctile" = pctile, "nobs" = nobs, "p1_tstat" = p1_tstat, "epsilon1_tstat" = epsilon1_tstat, "sigma1_tstat" = sigma1_tstat, "p2_tstat" = p2_tstat, "epsilon2_tstat" = epsilon2_tstat, "sigma2_tstat" = sigma2_tstat, "alpha_tstat" = alpha_tstat, "rho_tstat" = rho_tstat,
                                    "p1_pvalue" = pValues[1], "epsilon1_pvalue" = pValues[2], "sigma1_pvalue" = pValues[3], "p2_pvalue" = pValues[4], "epsilon2_pvalue" = pValues[5], "sigma2_pvalue" = pValues[6], "alpha_pvalue" = pValues[7], "rho_pvalue" = pValues[8]
      )
      
      if( count == 0 ){
        dfOut              <- dfOutTmp
        dfOutTstat         <- dfOutTstatTmp
      } else {
        dfOut              <- rbind( dfOut, dfOutTmp )
        dfOutTstat         <- rbind( dfOutTstat, dfOutTstatTmp )
      }
      count <- count + 1        
    }

  }
  
}

dfOut$method = "Nelder-Mead"


# write to file
filename             <- paste(c(pPathOut,"\\neg_excedances.csv"), collapse = "")
write.csv(x = dfOut, file = filename, row.names = FALSE, na = "")

filename             <- paste(c(pPathOut,"\\neg_excedances_tstats.csv"), collapse = "")
write.csv(x = dfOutTstat, file = filename, row.names = FALSE, na = "")




