# Importing required files
source("LHS-Support-Function.R")

# Control Parameters
stopTime = 60 #How many cycles to stop each calibration round at
samplesToKeep = 0.5

findMSE = function(tempCohort, LHSGenParameters, targetMatrix){
  timeBefore = proc.time()[[3]] # Used to time calibration
  count = 0
  MSE = 0
  
  while(1){
    count = count + 1
    cat("Cycle Num:", count, "\n")
    
    timeBefore = proc.time()[[3]]
    
    # Applying updates based on LHS generated model parameters
    for(i in 1:nrow(tempCohort)){
      updateState(tempCohort, LHSGenParameters, cycleNum = count)
    }
    
    if(count==60||count>=5){
      #Find MSE Code using percent difference
      var1 = tempCohort$var1
      var2 = tempCohort$var2
      var3 = tempCohort$var3
      
      target1 = targetMatrix[1]
      target2 = targetMatrix[2]
      target3 = targetMatrix[3]
      
      tempMSE = abs(var1 - target1)/((var1+target1)/2)
      tempMSE = tempMSE + abs(var2 - target2)/((var2+target2)/2)
    }
    if(count==stopTime) break
    #Or add any stop parameters
  }
  
  timeAfter = proc.time()[[3]] - timeBefore
  cat("\nMSE Calculation Time:", timeAfter, "\n")
  
  return(MSE)
}


# Function which tests a generated set of LHS parameters and updates the lb + ub
updateBounds = function(numSamples, targetMatrix){
  #Matrix of the best fitting sets of generated data
  topResults = matrix(data = NA, nrow = numSamples*samplesToKeep, ncol = n_param+1)
  
  #temporary matrix to hold results plus MSE
  results = matrix(data = NA, nrow = numSamples, ncol = n_param+1)

  colnames(topResults) = colnames(results) = c(v_param_names, "MSE")
  
  #Getting a random sample
  generatedLHSMatrix = sample_prior(numSamples)
  
  #Either use one cohort or a different cohort each run
  generatedCohort = generateCohort()
  
  # find the MSE associated with each set of samples
  for(i in 1:numSamples){
    #Setting MSE in results matrix
    results[i, 4] = findMSE(generatedCohort, generatedLHSMatrix[i, ], targetMatrix)
  }
  
  #Keeping the top portion of results
  results = results[order(results[,],decreasing=FALSE),] #Sorting by lowest MSE
  
  numToKeep = ncol(results)*samplesToKeep
  numToKeep = round(numToKeep)
  
  # Moving top portion of results into topResults matrix
  for(i in 1:numToKeep){
    topResults[i, ] = results[i, ]
  }
  
  lb = c(HIVToDeath = min(topResults[,"HIVToDeath"]),
         HealthyToDeath = min(topResults[,"HealthyToDeath"]),
         HealthyToHIV = min(topResults[,"HealthyToHIV"]))
  
  ub = c(HIVToDeath = max(topResults[,"HIVToDeath"]),
         HealthyToDeath = max(topResults[,"HealthyToDeath"]),
         HealthyToHIV = max(topResults[,"HealthyToHIV"]))
  
  return(list(lb,ub,topResults))
}
