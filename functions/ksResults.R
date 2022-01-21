# This function performs a KS Test and returns a dataframe with the result
# Inputs: 
#     data1 - vector of distrubution data to compare against
#     data2 - vector of distrubution data to be compared 
#     alt - alternative (i.e. test) in function ks.test
#     Dist, ID, Inte are strings containing information about 
#       the Distributions, Animal ID and FTLE integration time 
#       respectively (Note: these are specific to the manuscript 
#       BlueWhales_LagrangianFeatures)
# Returns:
#     A dataframe with the KS Test results and some summary statistics 
#     about each of the distributions

# Ensure the 'stats' package is installed and loaded
if (!require('stats',character.only = TRUE))
{
  install.packages('stats',dep=TRUE)
  if(!require('stats',character.only = TRUE)) stop("Package 'stats' not found")
}


ksResults <- function(data1, data2, alt, dist, ID,Inte) {
  resKS <- ks.test(x=data1,y=data2,alternative=alt)
  res <- data.frame(depid = ID,
                    Integration = Inte,
                    Distribution = dist,
                    test = resKS$alternative,
                    pvalue = as.numeric(resKS$p.value),
                    dvalue = as.numeric(resKS$statistic),
                    sampleSizeX = as.numeric(length(data1)),
                    meanX = mean(data1,na.rm=TRUE),
                    sdX= sd(data1,na.rm = TRUE),
                    sampleSizeY = as.numeric(length(data2)),
                    meanY = mean(data2,na.rm=TRUE),
                    sdY= sd(data2,na.rm = TRUE))
  return(res)
}
