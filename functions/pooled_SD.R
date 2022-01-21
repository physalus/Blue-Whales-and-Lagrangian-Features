## Function to find the Pooled Standard Deviation of a Sample of Samples
# i.e. ILI by dive has a stdev, but a group of dives will have a pooled stdev
# Input: 
#   stdevs = standard deviations of sample
#   numSamples = sample sizes for each observation
# Output:
#   Number (Double)
# 
pooledstd <- function(stdevs,numSamples, avoidNAs=FALSE) {
  if(avoidNAs){
    numSamples <- numSamples[!is.na(stdevs)]
    stdevs <- stdevs[!is.na(stdevs)]
  }
  psd <- sqrt( sum((numSamples-1)*(stdevs^2))/sum(numSamples-1) )
  psd
}