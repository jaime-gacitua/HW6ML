

PCA_missing <- function(numBasis, data, data.avail, tolerance, startB){

library("corpcor")
  
# On the data matrix, the rows are the features, and columns are the observations.
# Get n = columns of X (number of data points)
n = ncol(data)
# Get d = rows of x (number of features, dimension of x)
d = nrow(data)

# Initialize matrices
Z <- mat.or.vec(numBasis,n)
B <- mat.or.vec(d,numBasis)
B <- startB

cat('Starting Iterations\n')
convergence <- 1000000
iter <- 1
ptm <- 0
ttm <- proc.time()
while (convergence > tolerance) {
  ptm <- proc.time()
  cat("-----------\n")
  cat('Iteration = ', iter, "\n")
  ## First Step, solve i = 1,...,n systems of equations with B fixed. We get Z.
  for (i in 1:n){
  # System of equations M_i * z_i = c_i
    # Define the matrix M_i 
    M <- mat.or.vec(numBasis, numBasis)
    for (k in 1:numBasis){
      for (l in 1:numBasis){
        # Sum over dimension of x
        for (j in 1:d){
          M[k,l] <- M[k,l] + data.avail[j,i]*B[j,l]*B[j,k]
        } 
      } 
    }
    # Define the vector c_i
    c <- mat.or.vec(numBasis,1)
    # Sum over dimension of x
    for (l in 1:numBasis){
        c[l] = 0
      # Sum over dimension of x
      for (j in 1:d){
        if (data.avail[j,i] == 1) {
          c[l] <- c[l] + data.avail[j,i]*data[j,i]*B[j,l]
        } 
      }
    }  
    # Save z.old
    Z.old <- Z

       
    # Solve for z_i
#    print("M = ")
#    print(M)
#    cat("c = ", c, "\n")
    Z[,i] <- pseudoinverse(M) %*% c
#    print("Z = ")
#    print(Z)
    
  }
  
  ## Second Step, solve j = 1,...,d systems of equations with Z fixed. We get B.
  for (j in 1:d){
  
    # System of equations F_j * b_j = m_j
    # Define the matrix F_j 
    Fmat <- mat.or.vec(numBasis, numBasis)
    for (k in 1:numBasis){
      for (l in 1:numBasis){
        # Sum over dimension of x
        for (i in 1:n){
          Fmat[k,l] <- Fmat[k,l] + data.avail[j,i]*Z[l,i]*Z[k,i]
        } 
      } 
    }
    # Define the vector c_i
    m <- mat.or.vec(numBasis,1)
    # Sum over dimension of x
    for (l in 1:numBasis){
      # Sum over dimension of x
      for (i in 1:n){
        if(data.avail[j,i] == 1){
          m[l] <- m[l] + data.avail[j,i]*data[j,i]*Z[l,i]
        }
      } 
        
    }  
    # Save B.old
    B.old <- B
    # Solve for b_j
#    B[,j] <- solve(Fmat, m)
    B[j,] <- pseudoinverse(Fmat) %*% m
#    print("B = ")
#    print(Z)
    
  }
  
  # Calculate convergence criteria
  # Squared difference between iteration results
  convergence <- 0
  for(l in 1:numBasis) {
    for(i in 1:n) {
      convergence <- convergence + (Z[l,n]-Z.old[l,n])^2
    }
    for(j in 1:d) {
      convergence <- convergence + (B[j,l]-B.old[j,l])^2
    }
  }
  iter <- iter+1

  ### Calculate Root-mean-squared-error (RMSE)
  # RMSE Compute objective function
  RMSE <- 0
  temp.RMSE <- 0
  for(i in 1:n){
    for(j in 1:d){
      if(data.avail[j,i] == 1){
        for(k in 1:numBasis){
          temp.RMSE <- temp.RMSE + B[j,k]*Z[k,i]
        } 
        RMSE <- RMSE + (data[j,i]-temp.RMSE)^2
      }
      
    }  
  }
  # RMSE Count number of observations, and divide
  data.avail.count <- 0
  for(i in 1:n){
    for(j in 1:d){
      data.avail.count <- data.avail.count + data.avail[j,i]
    }
  }
  # RMSE Divide and apply square root
  RMSE <- RMSE / data.avail.count
  RMSE <- RMSE^(1/2)
  
  ### Output
  # Convergence Measures
  cat('Tolerance =\t', tolerance, "\n")
  cat('Convergence =\t', convergence, "\n")
  cat('RMSE =\t\t', RMSE, "\n")
  # Timing
  cat("Iteration Time =\t", proc.time() - ptm, "\n")
  cat("Total Time  =\t\t", proc.time() - ttm, "\n\n")
}
returnList <- list("Convergence" = convergence, "RMSE" = RMSE, "Total Time" = proc.time() - ttm, 
                   "Basis" = B)
return(returnList)
}


