
pageRank <- function(connections, epsilon, start.vector, convergence){
  d <- nrow(connections)
  
  ## Create transition matrix as Q_jk = 1/c(j)
  ## Where c(j) is the number of outbout links from page i

  # Initialize transition matrix as zeros.
  trans.mat <- matrix(0, nrow = d, ncol = d)
  
  # Compute c(j)
  c <- rowSums(connections)
  
  # Populate transition matrix  
  for(j in 1:d){
    for(k in 1:d){
      if (connections[j,k] == 1){
        trans.mat[j,k] <- 1 / c[j]
      }      
    }
  }
  
  # Create irreductible transition matrix by adding epsilon to all elements
  ones <- matrix(data = 1, nrow = d, ncol = d)
  trans.mat.irr <- (1-epsilon) * trans.mat + epsilon/d*ones
  
  # Initialize the iterating vector
  mu = start.vector
  mu.old = start.vector
  
  
  check <- 10000
  iter <- 1
  while(check > convergence){
    # Update mu
    mu <- mu.old %*% trans.mat.irr
    # Calculate convergence criteria
    check <- dist(rbind(mu,mu.old))
    # Record old mu
    mu.old <- mu
    # Output stats
    cat('Done with iteration:', iter, '\n')
    cat('Convergence:', check, ' / ', convergence, '\n')
    iter <- iter + 1
  }
return(mu)
}


connections <- matrix(c(0,0,1,1,0,0,
                        1,0,0,0,0,1,
                        1,1,0,1,1,0,
                        0,0,0,0,0,0,
                        0,0,0,0,0,1,
                        0,0,0,0,0,0), 
                      nrow = 6, ncol = 6)
initial <- matrix(1/nrow(connections), nrow = 1, ncol = nrow(connections))
epsilon <- 0.15
convergence <- 0.0001

# Question b
pi <- pageRank(connections,epsilon,initial,convergence) 
pi

# Question c

epsilon.values <- matrix(seq(0.1, 1, 0.05), nrow = 19)
pi.matrix <- matrix(0, nrow = nrow(epsilon.values), ncol = nrow(connections))

for(iter in 1:nrow(epsilon.values)){
  pi.matrix[iter,] <- pageRank(connections,epsilon.values[iter,1],initial,convergence)
}

pi.matrix

a <-image(t(pi.matrix))
title(main = 'Solutions of Figure 14.47 from HTF, for different epsilon', 
      xlab = 'Page #', ylab = 'epsilon')

