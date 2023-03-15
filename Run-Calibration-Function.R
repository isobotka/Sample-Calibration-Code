# Run Calibration file
source("LHS-Support-Function.R")
source("Main-Calibration-Function.R")

#Control Parameters
numberOfGenerations = 10

# The Lower and Upper bounds for the random LHS generation
lb = c(1, 2, 3)
ub = c(2, 3, 4)

# Parameters to Calibrate
v_param_names = c("HIVToDeath", "HealthyToDeath", "HealthyToHIV")
n_param = length(v_param_names)

#Target Matrix
targetMatrix = matrix(data = 1, nrow = 1, ncol = 3)

for(i in 1:numberOfGenerations){
  output = updateBounds(numSamples = 10, targetMatrix = targetMatrix)
  lb = output[1]
  ub = output[2]
  topResults = output[3]
  
  write.csv(topResults, paste0("results-gen",i,".csv"))
}
