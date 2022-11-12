#modified script based on Yara Alshwairikh based on source code from https://rdrr.io/bioc/LEA/man/lfmm2.html
#this script needs to be run for each individual environmental variable

library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)
#library(LEA) #library for LEA might work instead of creating function, but previously had issues

print("lfmm2 analysis") 


#Load the allele frequency matrix with 20 simulated individuals per population
rbeta_data <- readRDS("../outputs/sim_data.rds")
rbeta_data$.id <- NULL #remove .id column that gets generated when reading the data

#pull env data in from rda analysis
env_lfmm <- read.csv('../outputs/env_rda.csv')

#duplicate for each of the 20 individuals
# Create an index of the rows you want with duplications
env_lfmm <- env_lfmm[rep(seq_len(nrow(env_lfmm)), each = 20), ]

#remove columns that do not have the env variables
env_lfmm <- env_lfmm[,-c(1:3)]

print("Loaded genomic and environmental data")

# Create the lfmm2 function
setClass("lfmm2Class",
         slots = c(K = "integer", 
                   lambda = "numeric",
                   U = "matrix",
                   V = "matrix"
         )
)

lfmm2 <- function(input,
                  env, 
                  K, 
                  lambda = 1e-5){
  
  ## Check response input matrix 
  ## LEA  
  if (is.character(input)){
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
    }
  }
  
  ## Check independent/covariate env matrix  
  ## LEA 
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  # run SVD of X: X = Q Sigma R
  
  svx <- svd(x = scale(X, scale = FALSE), nu = n)
  Q <- svx$u
  
  d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-d))
  d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-d))
  D_inv <- diag(d_lambda_inv)
  D  <- diag(d_lambda)
  
  # run SVD of modified Y    
  svk <- svd(D %*% t(Q) %*% scale(Y, scale = FALSE), nu = K)
  
  if (K > 1) {
    Sigma_k <- diag(svk$d[1:K])
  } else {
    Sigma_k <- as.matrix(svk$d[1])
  }
  
  # compute the latent matrix W
  W <- Q %*% D_inv %*% tcrossprod(svk$u %*% Sigma_k, svk$v[,1:K])
  
  # compute LFMM factors U and loadings V
  # Non orthogonal factors
  U <- crossprod(t(Q %*% D_inv), svk$u %*% Sigma_k)
  #U <- Q %*% D_inv %*% svk$u %*% Sigma_k
  V <- svk$v[,1:K]
  
  obj <- new("lfmm2Class")
  obj@K <- as.integer(K)
  obj@lambda <- as.numeric(lambda)
  obj@U <- as.matrix(U)
  obj@V <- as.matrix(V)
  
  ## LEA 
  return(obj)
}


setGeneric("lfmm2.test", function(object, input, env, 
                                  genomic.control = TRUE, 
                                  linear = TRUE, 
                                  family  = binomial(link = "logit")) matrix);
setMethod("lfmm2.test", "lfmm2Class",
          function(object, 
                   input,
                   env,
                   genomic.control, 
                   linear,
                   family
          ) {
            
            ## Check input matrix   
            ## LEA  
            if (is.character(input)){
              Y <- read.lfmm(input)
              lst.unique <- unique(as.numeric(Y))
              if (9 %in% lst.unique){
                stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
              }
              if (-9 %in% lst.unique){
                stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
              }
            } else {
              ## Y is an R object       
              if (is.null(input)){
                stop("NULL value for argument 'input'.")
              }
              Y <- as.matrix(input)
              Y[Y == 9] <- NA
              Y[Y == -9] <- NA
              if (anyNA(Y)) {
                stop("The input matrix contains missing values (NA or 9).")
              }
            }
            
            ## Check independent/covariate matrix  
            ## LEA 
            if (is.character(env)){
              X <- read.env(env)
              if (anyNA(X)){
                stop("'env' file contains missing data (NA).")
              }
            } else {
              if (is.null(env)){
                stop("NULL value for argument 'env'.")
              }
              X <- as.matrix(env)
              if (anyNA(X)) {
                stop("The environmental matrix contains NA.")
              }
            }
            
            d <-  ncol(X) #number of environmental variables
            n <-  nrow(X) #number of individuals
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
            }
            
            p <- ncol(Y)
            p_value <- NULL
            z_score <- NULL
            
            if (linear){
              mod_lm <- lm(Y ~ ., data = data.frame(X, object@U)) 
              sm <- summary(mod_lm)
              p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
              z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
            } else {
              for (j in 1:p) {
                mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
                sm <- summary(mod_glm)
                p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
                z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
              }
            }
            if (genomic.control){
              gif <- apply(z_score^2, 2, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
              p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
            } else {
              gif <- NULL
            }
            res <- list(pval = p_value, zscores = z_score, gif = gif)
            return(res)
          }
)

print("finished creating lfmm2 function")

######################################
# Fitting an LFMM with K = 9 factors
######################################

#run the model for each environmental variable
env_lfmm1 <- env_lfmm[,1, drop = FALSE]
env_lfmm2 <- env_lfmm[,2, drop = FALSE]
env_lfmm3 <- env_lfmm[,3, drop = FALSE]
env_lfmm4 <- env_lfmm[,4, drop = FALSE]
env_lfmm5 <- env_lfmm[,5, drop = FALSE]
env_lfmm6 <- env_lfmm[,6, drop = FALSE]
env_lfmm7 <- env_lfmm[,7, drop = FALSE]
env_lfmm8 <- env_lfmm[,8, drop = FALSE]
env_lfmm9 <- env_lfmm[,9, drop = FALSE]


model1 <- lfmm2(input = rbeta_data, env = env_lfmm1, K = 9)
model2 <- lfmm2(input = rbeta_data, env = env_lfmm2, K = 9)
model3 <- lfmm2(input = rbeta_data, env = env_lfmm3, K = 9)
model4 <- lfmm2(input = rbeta_data, env = env_lfmm4, K = 9)
model5 <- lfmm2(input = rbeta_data, env = env_lfmm5, K = 9)
model6 <- lfmm2(input = rbeta_data, env = env_lfmm6, K = 9)
model7 <- lfmm2(input = rbeta_data, env = env_lfmm7, K = 9)
model8 <- lfmm2(input = rbeta_data, env = env_lfmm8, K = 9)
model9 <- lfmm2(input = rbeta_data, env = env_lfmm9, K = 9) ##Set value of K that is appropriate for your data

saveRDS(model1, file = "../outputs/load_LFMM1.RDS")
saveRDS(model2, file = "../outputs/load_LFMM2.RDS")
saveRDS(model3, file = "../outputs/load_LFMM3.RDS")
saveRDS(model4, file = "../outputs/load_LFMM4.RDS")
saveRDS(model5, file = "../outputs/load_LFMM5.RDS")
saveRDS(model6, file = "../outputs/load_LFMM6.RDS")
saveRDS(model7, file = "../outputs/load_LFMM7.RDS")
saveRDS(model8, file = "../outputs/load_LFMM8.RDS")
saveRDS(model9, file = "../outputs/load_LFMM9.RDS")


print("finished running lfmm2 model")



#Code below edited on Nov 11, 2020. Changes: save the RAW output without transposing 
#to check that it's happening ok. Also print the GIF value

# Computing P-values
#save p values, z-scores, and gif into a list of 3 (pv)
df <- lfmm2.test(object = model1, input = rbeta_data, env = env_lfmm1, linear = TRUE)
df_pval <- df$pval #pull out just p value
env <- names(env_lfmm1) #pull in env name
df_pval <- cbind(env, df_pval) #add name to pval
write.csv(df_pval, file = "../outputs/allraw_lfmm1_K9.csv", row.names=FALSE)
 
df <- lfmm2.test(object = model2, input = rbeta_data, env = env_lfmm2, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm2)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm2_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model3, input = rbeta_data, env = env_lfmm3, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm3)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm3_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model4, input = rbeta_data, env = env_lfmm4, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm4)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm4_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model5, input = rbeta_data, env = env_lfmm5, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm5)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm5_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model6, input = rbeta_data, env = env_lfmm6, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm6)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm6_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model7, input = rbeta_data, env = env_lfmm7, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm7)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm7_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model8, input = rbeta_data, env = env_lfmm8, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm8)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm8_K9.csv", row.names=FALSE) 

df <- lfmm2.test(object = model9, input = rbeta_data, env = env_lfmm9, linear = TRUE)
df_pval <- df$pval
env <- names(env_lfmm9)
df_pval <- cbind(df_pval, env)
write.csv(df_pval, file = "../outputs/allraw_lfmm9_K9.csv", row.names=FALSE) 

print("finished running lfmm2 significance test")