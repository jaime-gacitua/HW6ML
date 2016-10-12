source("Q2.R") # Function PCA with missing data.

# Input Parameters
num.users = 200 # Max 943
num.movies = 400 # Max 1682
num.pca.basis = 5
convergence = 1


### Prepare the database (snippet from Q2 HW5)
movies <- read.table("u.item", sep = "|", header = FALSE, stringsAsFactors = FALSE, quote="")
movies <- movies[,c(1,2)]
names(movies) <- c("movieid","movie")
rank   <- read.table("u.data", sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("userid","movieid","rating","ts"))
critics <- merge(movies, rank, by = "movieid")
critics$movie <- NULL
critics$ts <- NULL
names(critics) <- c("movieid","person","rank")

### Prepare matrices
# Data matrix
data = matrix(nrow = num.users, ncol = num.movies)
gamma = matrix(nrow = num.users, ncol = num.movies)
for(i in 1:nrow(critics)) {
  if (critics[i,2] <= num.users & critics[i,1] <= num.movies) {
    data[critics[i,2] , critics[i,1]] <- critics[i,3]
  }
}
# Gamma matrix
for(i in 1:ncol(data)){
  for(d in 1:nrow(data)){
    if(is.na(data[d,i]) ){
      gamma[d,i] = 0
    } else {
      gamma[d,i] = 1
    }
  }
}
# Starting Basis
start.pca.basis = matrix(0, nrow = nrow(data), ncol = num.pca.basis)

# Part A
for(d in 1:nrow(data)){
  start.pca.basis[d,d %% num.pca.basis + 1] = 1
}
finalB <- PCA_missing(num.pca.basis, data, gamma, convergence, start.pca.basis)
finalB

# Part B

K = 15
RMSE <- mat.or.vec(K,1)
iter = 1
for(iter in 1:K){
  start.pca.basis <- matrix( rnorm(nrow(data)*num.pca.basis,mean=0,sd=1), nrow(data), num.pca.basis) 
  output <- PCA_missing(num.pca.basis, data, gamma, convergence, start.pca.basis)
  RMSE[iter] <- output$RMSE
}
RMSE
plot(RMSE)
RMSE.sorted <- sort(RMSE, decreasing = TRUE)
plot(RMSE.sorted, pch = 19, xlab = "Random Starting Points", ylab = "RMSE", 
     main = "PCA with missing values: RMSE different random starting points, decreasing order")


