



add = function(x, y){
  
  return(x+y)
}

add(1,2)


HIVProgression = function(temp_HIV){
  
  # 1. HIV Progression (Comm, RG) =======================================================================
  for (i in 1:n_comm) {
    for (j in 1:n_risk) {
      ## Status quo
      temp_HIV[i, j, ] <- pop_total_t[t, i, j, ] %*% a_HIV_mort_AT[ , , i, j] %*% a_HIV_prog[ , , i, j] %*%
        a_HIV_aware[ , , i, j] %*% a_HIV_trt[ , , i, j] %*% a_HIV_prep[ , , i, j]
    }
  }
  
  if (isFALSE(all.equal(sum(temp_HIV[i, j, ]), sum(pop_total_t[t, i, j, ])))) {   
    print(paste0("X, temp_HIV, cycle t = ", t, " , i = ", i, " , j = ", j))
  } 
  if (temp_HIV[i, j, k] < 0) {   
    print(paste0("X, temp_HIV, Negative Population, cycle t = ", t, " , i = ", i, " , j = ", j, " , k = ", k))
  } 
  
  return(temp_HIV)
}

generateCohort(cohortSize){
  
  # Blank arrays -----------------------------------------------------------------
  ## 1. HIV progression  
  temp_HIV <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))
  
  ## 2. Community transition
  temp_comm <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))
  
  ## 3. Risk group transition
  temp_prison_pop <- temp_prison_RG <- temp_prison_HIV <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))    
  temp_r1_pop <- temp_r1_RG <- temp_r1_HIV <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))    
  temp_RG_final <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))    
  
  ## 4. New HIV infections
  Population_ByCommRisk <- matrix(0, nrow = n_comm, ncol = n_risk, dimnames = list(v_comm_names, v_risk_names))
  
  temp_nHIV <- array(0, dim = c(n_comm, n_risk, n_hiv_neg), dimnames = list(v_comm_names, v_risk_names, v_hiv_neg_names))    
  
  
  ## 5. New population after new HIV infections
  temp_final <- array(0, dim = c(n_comm, n_risk, n_hiv), dimnames = list(v_comm_names, v_risk_names, v_hiv_names))    
  
  
  
  # **** Policy PrEP Implementations -- Initial Population ===========================================================
  ## 1. when t=1, give population in HIV-,PrEP+
  #   1-1. PWID
  pop_total_t[1, 2, 1, 14] <- pop_total_t[1, 2, 1, 1] * iv_pc_comm_nHIV_prep
  pop_total_t[1, 2, 1, 1]  <- pop_total_t[1, 2, 1, 1] * (1 - iv_pc_comm_nHIV_prep)
  
  #   1-2. MSM
  pop_total_t[1, 2, 2, 14] <- pop_total_t[1, 2, 2, 1] * iv_pc_comm_nHIV_prep
  pop_total_t[1, 2, 2, 1]  <- pop_total_t[1, 2, 2, 1] * (1 - iv_pc_comm_nHIV_prep)
  
  #   1-3. MWID
  pop_total_t[1, 2, 3, 14] <- pop_total_t[1, 2, 3, 1] * iv_pc_comm_nHIV_prep
  pop_total_t[1, 2, 3, 1]  <- pop_total_t[1, 2, 3, 1] * (1 - iv_pc_comm_nHIV_prep)
  
  #   1-4. Low
  pop_total_t[1, 2, 4, 14] <- pop_total_t[1, 2, 4, 1] * iv_pc_comm_nHIV_prep
  pop_total_t[1, 2, 4, 1]  <- pop_total_t[1, 2, 4, 1] * (1 - iv_pc_comm_nHIV_prep)
  
  
  # # CHECK) = 'iv_pc_comm_nHIV_prep' = 0.2249
  # for (j in 1:n_risk) {
  #   print(sum(pop_total_t[1, 2, j, 14]) / sum(pop_total_t[1, 2, j, c(1, 14)]))
  # }
  
  return(temp_HIV)
}


simulate = function(temp_HIV, numCycles){
  for (t in 1:numCycles) {
    temp_HIV = HIVProgression(temp_HIV)
    
    temp_HIV = ComTrans(temp_HIV)
  }
  
  return(temp_HIV)
}


updateInputs = function(lhsMatrix){
  a_HIV_Trt[]  = theProperInputs
  
  assign(a_HIV_Trt, "global")
}

findGOF = function(temp_HIV, targets){
  GOF = (temp_HIV[, ,"Death"] - targets["Death"])^2
  
  GOF = 1/GOF
  
  return(GOF)
}

calibrate = function(numSamples, numGenerations = 5){
  for(x in 1:numGenerations){
    
    lhsMatrix = sample_prior(numSamples)
    
    lhsMatrix = cbind(lhsMatrix, vector(mode = "numeric", length = nrow(lhsMatrix)))
    colnames(lhsMatrix) = c("1", "2", "3", "GOF")
    
    for(i in 1:numSamples){
      updateInputs(lhsMatrix[i, ])
      
      temp_HIV = generateCohort()
      
      temp_HIV = simulate(temp_HIV, )
      
      GOF = findGOF(temp_HIV, targets)
    }
    
    #Sort lhsMatrix by GOF
    lhsMatrix = lhsMatrix[order(lhsMatrix[,"GOF"],decreasing=TRUE),]
    
    #Keep top rows of lhsMatrix
    lhsMatrix = lhsMatrix[1:3, ]
    
    write.csv(lb, "lowerBound.csv")
    write.csv(ub, "upperBound.csv")
    
    #update the bounds
    lb = c(HIVToDeath = min(lhsMatrix[,"HIVToDeath"]),
           HealthyToDeath = min(lhsMatrix[,"HealthyToDeath"]),
           HealthyToHIV = min(lhsMatrix[,"HealthyToHIV"]))
    
    ub = c(HIVToDeath = max(lhsMatrix[,"HIVToDeath"]),
           HealthyToDeath = max(lhsMatrix[,"HealthyToDeath"]),
           HealthyToHIV = max(lhsMatrix[,"HealthyToHIV"]))
    
  }
  
  #Percent change between generations is lo
  
}