# LHS Support Code
library(lhs)
#install.packages("lhs")

#Given a number of desired samples (rows) return n_param columns of LHS generated results
sample_prior <- function(n_samp){
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param) # Creating the LHS unit vector
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){ # Randomly generating values which fit with the LHS vector
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = lb[i],
                               max = ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            shape1 = 1,
    #                            shape2 = 1)
  }
  return(m_param_samp)
}

