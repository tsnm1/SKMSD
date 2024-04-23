# Code to simulate data
#' Title
#'
#' @param n_1
#' @param n_p
#' @param rho_p
#' @param theta
#' @param lambda
#' @param prob_0
#' @param copula
#' @param marginal
#'
#' @return
#' @export
#'
#' @examples
simulate_data_1 <- function(n_1 = 100, n_p = 100, rho_p = 0.6, theta, lambda, prob_0,
                            copula = TRUE, marginal = 'pg'){

  library(mvtnorm)
  library(MASS)

  if(copula){
    Sigma <- matrix( 0, n_p, n_p)
    Sigma = rho_p^abs(row(Sigma) - col(Sigma))
    result1 <- rmvnorm(n = n_1, rep(0.0, n_p), Sigma)
    result1 <- matrix(result1, nrow = n_1)
    result2 <- matrix(apply(result1, 2, pnorm), nrow = n_1)
    # iter <- 2
    # print("result2")
    result3 <- sapply(1:n_p, function(iter){
      p_cdf_0 <- pmax(0.0, result2[, iter] - prob_0[iter]) / (1-prob_0[iter])
      W_simulate <- sapply(1:n_1,function(iter2){
        qnbinom(p_cdf_0[iter2], size = 1/exp(theta[iter2,iter]), mu = exp(lambda[iter2,iter]))
      })
    })
  }
  else{
    if(marginal == 'nb'){
      result3 <- sapply(1:n_p, function(iter){
        W_simulate <- sapply(1:n_1,function(iter2){
          rnbinom(1,size = 1/exp(theta[iter2,iter]), mu = exp(lambda[iter2,iter]))
        })
        zero_sim = rbinom(n_1,size = 1,prob = prob_0[iter])
        W_simulate = W_simulate*(zero_sim==0)
      })
    }
    else if(marginal == 'pg'){
      result3 <- sapply(1:n_p, function(iter){
        W_simulate <- sapply(1:n_1,function(iter2){
          U_simulate = rgamma(1,shape = 1/exp(theta[iter2,iter]),scale = exp(theta[iter2,iter]))
          rpois(1,exp(lambda[iter2,iter])*U_simulate)
        })
        zero_sim = rbinom(n_1,size = 1 ,prob = prob_0[iter])
        W_simulate = W_simulate*(zero_sim==0)
      })
    }
  }
  return(result3)

}


simulate_datasets_3 <- function(times=1, n_data = 2, n_1_all = c(400,400),
                                n_p_all = c(200,400,600,800), n_p = c(20,40,60,80),
                                diff = seq(10,60,10), prob_max = seq(0.1,0.9,0.2),
                                marginal1 = 'pg',copula1 = TRUE){
  data_list <- res_prob_0 <- list()
  for (i in c(1:times)){
    data_list[[i]] <- res_prob_0[[i]] <- list()
    for(j in c(1:length(n_p_all))){
      data_list[[i]][[j]] <- res_prob_0[[i]][[j]] <-  list()
      for(k in c(1:length(diff))){
        data_list[[i]][[j]][[k]] <- res_prob_0[[i]][[j]][[k]] <- list()
        for(p1 in c(1:length(prob_max))){
          data_list[[i]][[j]][[k]][[p1]] <- res_prob_0[[i]][[j]][[k]][[p1]] <- list()
          data_n_data <- c()
          cat(i,j,k,p1,'\n')
          for(i3 in c(1:n_data)){
            sim_data1 <- simulate_PG_6(n_1 = n_1_all[i3], n_p = n_p_all[j], n_p_1 = n_p[j], # n_p[j],
                                       diff = diff[k], k = i3,
                                       prob_up1 = prob_max[p1], prob_up2 = prob_max[p1],
                                       rho_p= rho_p, n_x = n_x1, rho_x = rho_x, M_max = M_max, fc = fc1,
                                       marginal = marginal1, copula = copula1, seed = i*i3)
            data1 <- cbind(i3,sim_data1$W_simulate)

            res_prob_0[[i]][[j]][[k]][[p1]][[i3]] <- sim_data1$prob_0
            data_n_data <- rbind(data_n_data,data1)
          }
          data_list[[i]][[j]][[k]][[p1]] <- data_n_data
        }
      }
    }
  }
  return(list(data_list,res_prob_0))
}
# diff = 1, k = 1;
simulate_PG_6 <- function(n_1 = 20, n_p = 10, n_p_1 = 5,
                          diff = 0, k = 1,
                          prob_up1 = 0.1, prob_up2 = 0.1,
                          rho_p=0.5, n_x = 3, rho_x = 0.1, M_max = 1e2, fc = 2,
                          marginal = c('nb', 'pg'), copula = TRUE, seed = 111){
  # Call libraries
  library(mvtnorm)
  library(MASS)

  if(FALSE){
    n_1 = 20
    n_p = 20
    n_p_1 = 5
    k = 1
    prob_up1 <- 0.5
    prob_up2 <- 0.5
    diff = 0
    rho_p=0.5
    n_x = 3
    rho_x = 0.1
    M_max = 1e2
    seed = 111
    fc <- 2
  }

  M <- round(M_max*n_p*(1 + runif(n_1,-1,1)/3)) #M_max*n_p_1

  prob0_11 <- runif(n_p_1 + k*diff,0,prob_up1)
  prob0_12 <- runif(n_p_1 + k*diff,0,prob_up1)
  prob0_2 <- runif(n_p-n_p_1-k*diff,0,prob_up2)

  prob0_11 <- c(prob0_11,prob0_2)
  prob0_12 <- c(prob0_12,prob0_2)
  prob_0 <- c(prob0_11,prob0_2,prob0_12,prob0_2)


  if(n_x>0){
    # n_x = 3
    # n_1 <- 5
    Sigma_x <- matrix( 0, n_x, n_x )
    Sigma_x = rho_x^abs(row(Sigma_x) - col(Sigma_x))
    data_x <- scale( rmvnorm(n_1,rep(0,n_x),Sigma_x) )
    data_x <- apply(data_x,2,function(x){
      (x - min(x)) / (max(x) - min(x))
    })
    # print(data_x)
    colnames(data_x) <- paste("X_", 1:n_x, sep = "")
    X = as.matrix(cbind(1,data.frame(data_x)))
    # beta <- beta_star <- matrix(1,n_x+1,n_p)
  }else{
    n_x = 3
    data_x <- matrix(0,n_1,n_x)
    colnames(data_x) <- paste("X_", 1:n_x, sep = "")
    X = as.matrix(cbind(1,data.frame(data_x)))
    # beta <- beta_star <- matrix(1,n_x+1,n_p)
  }
  if(TRUE){

    if(diff==0){
      sd <- 0.01
      beta1 <- rbind(1,matrix(1,n_x,n_p_1))
      beta2 <- rbind(1,matrix(0,n_x,n_p-n_p_1))
      beta_1 <- cbind(beta1,beta2)
      theta0 <- X %*% beta_1
      vec_x <- ifelse(rbinom(n_p_1, size=1, prob = 0.5), 2^fc, 2^(-fc))
      theta0[(n_1/2+1):n_1,c(1:n_p_1)] <-
        theta0[(n_1/2+1):n_1,c(1:n_p_1)] * matrix(vec_x, nrow = n_1/2, ncol = length(vec_x), byrow = TRUE)
      theta0 <- theta0 + matrix(rnorm(n_1 * n_p, mean = 0, sd = sd), nrow = n_1)

      beta1 <- rbind(1,matrix(.1,n_x,n_p_1))
      beta2 <- rbind(1,matrix(0,n_x,n_p-n_p_1))
      beta_1 <- cbind(beta1,beta2)
      lambda0 <- X %*% beta_1
      vec_x <- ifelse(rbinom(n_p_1, size=1, prob = 0.5), 2^fc, 2^(-fc))
      lambda0[(n_1/2+1):n_1,c(1:n_p_1)] <-
        lambda0[(n_1/2+1):n_1,c(1:n_p_1)] * matrix(vec_x, nrow = n_1/2, ncol = length(vec_x), byrow = TRUE)
      lambda0 <- lambda0 + matrix(rnorm(n_1 * n_p, mean = 0, sd = sd), nrow = n_1)
    }else{
      sd <- 0.01
      beta1 <- rbind(1,matrix(1,n_x,n_p_1))
      beta3 <- rbind(1,matrix(1,n_x,diff))
      beta2 <- rbind(1,matrix(0,n_x,n_p-n_p_1-diff))
      beta_1 <- cbind(beta1,beta2[,0:((k-1)*diff)],beta3,beta2[,((k-1)*diff+1):(n_p-n_p_1-diff)])
      theta0 <- X %*% beta_1
      vec_x <- ifelse(rbinom(n_p_1+diff, size=1, prob = 0.5), 2^fc, 2^(-fc))
      theta0[(n_1/2+1):n_1,c(1:n_p_1,(n_p_1+(k-1)*diff+1):(n_p_1+k*diff))] <-
        theta0[(n_1/2+1):n_1,c(1:n_p_1,(n_p_1+(k-1)*diff+1):(n_p_1+k*diff))] * matrix(vec_x, nrow = n_1/2, ncol = length(vec_x), byrow = TRUE)
      theta0 <- theta0 + matrix(rnorm(n_1 * n_p, mean = 0, sd = sd), nrow = n_1)

      beta1 <- rbind(1,matrix(.1,n_x,n_p_1))
      beta3 <- rbind(1,matrix(.1,n_x,diff))
      beta2 <- rbind(1,matrix(0,n_x,n_p-n_p_1-diff))
      beta_1 <- cbind(beta1,beta2[,0:((k-1)*diff)],beta3,beta2[,((k-1)*diff+1):(n_p-n_p_1-diff)])
      lambda0 <- X %*% beta_1
      vec_x <- ifelse(rbinom(n_p_1+diff, size=1, prob = 0.5), 2^fc, 2^(-fc))
      lambda0[(n_1/2+1):n_1,c(1:n_p_1,(n_p_1+(k-1)*diff+1):(n_p_1+k*diff))] <-
        lambda0[(n_1/2+1):n_1,c(1:n_p_1,(n_p_1+(k-1)*diff+1):(n_p_1+k*diff))] * matrix(vec_x, nrow = n_1/2, ncol = length(vec_x), byrow = TRUE)
      lambda0 <- lambda0 + matrix(rnorm(n_1 * n_p, mean = 0, sd = sd), nrow = n_1)
    }

  }

  lambda0_1 <- lambda0[1:(n_1/2),]
  lambda0_2 <- lambda0[(n_1/2 +1):(n_1),]

  theta0_1 <- theta0[1:(n_1/2),]
  theta0_2 <- theta0[(n_1/2 +1):(n_1),]
  # copula <- T
  # marginal <- "pg"
  simulate_result1_1 <- simulate_data_2(n_1 = n_1/2, n_p = n_p, rho_p = rho_p, theta = theta0_1, lambda = lambda0_1, prob_0 = prob0_11,
                                        copula = copula, marginal = marginal)
  simulate_result1_2 <- simulate_data_2(n_1 = n_1/2, n_p = n_p, rho_p = rho_p, theta = theta0_2, lambda = lambda0_2, prob_0 = prob0_12,
                                        copula = copula, marginal = marginal)

  # data1 <- rbind(simulate_result1_1,simulate_result1_2)
  W <- rbind(simulate_result1_1,simulate_result1_2)
  W <- round(W/(apply(W,1,sum)+1e-5)*M)
  # table(apply(W,1,sum)-M)
  W_simulate <-  cbind(data_x,M,W)

  return(list(W_simulate = W_simulate, marginal,
              data_x = data_x, M = M, beta = beta, prob_0 = prob_0))

}

# ind <- which(prob_0<0.8)
simulate_data_2 <- function(n_1 = 100, n_p = 100, rho_p = 0.6, theta, lambda, prob_0,
                            copula = TRUE, marginal = 'pg'){

  library(mvtnorm)
  library(MASS)

  if(copula){
    ind <- which(prob_0<0.8)
    ind2 <- (1:n_p)[-ind]
    n_p1 <- length(ind)
    Sigma <- matrix( 0, n_p1, n_p1)
    Sigma = rho_p^abs(row(Sigma) - col(Sigma))
    result1 <- rmvnorm(n = n_1, rep(0.0, n_p1), Sigma)
    result1 <- matrix(result1, nrow = n_1)
    result2 <- matrix(apply(result1, 2, pnorm), nrow = n_1)
    # iter <- 2
    # print("result2")
    result3 <- matrix(0, nrow = n_1, ncol = n_p )
    if(length(ind)>0){
      prob_01 <- prob_0[ind]
      theta1 <- theta[,ind]
      lambda1 <- lambda[,ind]
      result31 <- sapply(1:n_p1, function(iter){
        p_cdf_0 <- pmax(0.0, result2[, iter] - prob_01[iter]) / (1-prob_01[iter])
        W_simulate <- sapply(1:n_1,function(iter2){
          qnbinom(p_cdf_0[iter2], size = 1/exp(theta1[iter2,iter]), mu = exp(lambda1[iter2,iter]))
        })
      })
      result3[, ind] <- result31
    }

    if(length(ind2)>0){
      prob_02 <- prob_0[ind2]
      theta2 <- theta[,ind2]
      lambda2 <- lambda[,ind2]
      result32 <- sapply(1:length(ind2), function(iter){
        W_simulate <- sapply(1:n_1,function(iter2){
          rnbinom(1,size = 1/exp(theta2[iter2,iter]), mu = exp(lambda2[iter2,iter]))
        })
        zero_sim = rbinom(n_1,size = 1,prob = prob_02[iter])
        W_simulate = W_simulate*(zero_sim==0)
      })
      result3[, ind2] <- result32
    }
  }
  else{
    if(marginal == 'nb'){
      result3 <- sapply(1:n_p, function(iter){
        W_simulate <- sapply(1:n_1,function(iter2){
          rnbinom(1,size = 1/exp(theta[iter2,iter]), mu = exp(lambda[iter2,iter]))
        })
        zero_sim = rbinom(n_1,size = 1,prob = prob_0[iter])
        W_simulate = W_simulate*(zero_sim==0)
      })
    }
    else if(marginal == 'pg'){
      result3 <- sapply(1:n_p, function(iter){
        W_simulate <- sapply(1:n_1,function(iter2){
          U_simulate = rgamma(1,shape = 1/exp(theta[iter2,iter]),scale = exp(theta[iter2,iter]))
          rpois(1,exp(lambda[iter2,iter])*U_simulate)
        })
        zero_sim = rbinom(n_1,size = 1 ,prob = prob_0[iter])
        W_simulate = W_simulate*(zero_sim==0)
      })
    }
  }
  return(result3)

}

##################

ZIPG_Estimate_3_1 <- function(data_x,W,M,zp_cutoff=0.8,min_nonzero_num=1,epsilon = 1e-5){

  n1 <- nrow(W)
  p1 <- ncol(W)

  gene_zero_prop <- apply(W, 2, function(y){
    sum(y < 1e-5) / n1
  })
  # plot(sort(gene_zero_prop))
  zp_cutoff <- 0.8
  min_nonzero_num <- 1

  gene_sel1 <- as.integer(which(gene_zero_prop < zp_cutoff))
  gene_sel2 <- as.integer(which(gene_zero_prop <= 1.0 - min_nonzero_num/n1 &
                                  gene_zero_prop >= zp_cutoff))
  gene_sel3 <- (1:p1)[-c(gene_sel1, gene_sel2)]

  # as.numeric(gene_sel1)
  # hist(gene_sel1,breaks = 20)
  # aaa <- W[,gene_sel3]
  # length(aaa[aaa!=0])

  param1 <- c()
  r1 <- c()
  epsilon = 1e-5

  formula_obj <- as.formula( paste("~", paste(colnames(data_x), collapse = " + ")))

  for(i in gene_sel1){
    W1 <- W[,i]              # gene_sel2[2]
    W1[W1==min(W1)] <- 0
    ZIPG_res <- ZIPG_main(data = data_x, X = formula_obj, X_star = formula_obj, W = W1, M = M)
    # ZIPG_res$res$par
    # ZIPG_summary(ZIPG_res)

    # sim_M = M
    X = as.matrix(cbind(1,data.frame(data_x)))
    param_1 = ZIPG_res$res$par  #p = 0.5

    beta = param_1[1:4]
    beta_star = param_1[5:8]
    prob0 <- 1/(1+exp(-param_1[9]))
    lambda =  X %*% beta + log(M)  # exp(X %*% beta )*M
    theta = X %*% beta_star        # exp(X %*% beta_star)
    # hist(lambda/M)
    lambda[is.infinite(lambda)] <- ifelse(lambda[is.infinite(lambda)] == Inf,
                                          max(lambda[!is.infinite(lambda)], na.rm = TRUE), min(lambda[!is.infinite(lambda)], na.rm = TRUE))
    theta[is.infinite(theta)] <- ifelse(theta[is.infinite(theta)] == Inf,
                                        max(theta[!is.infinite(theta)], na.rm = TRUE), min(theta[!is.infinite(theta)], na.rm = TRUE))
    u <- t(sapply(1:n1, function(i_1){
      u1 <- prob0 + (1 - prob0) * pnbinom(W1[i_1], size = 1/exp(theta[i_1]), mu = exp(lambda[i_1]))
      u2 <- (prob0 + (1 - prob0) * pnbinom(W1[i_1] - 1, size = 1/exp(theta[i_1]),mu = exp(lambda[i_1]))) *
        as.integer(W1[i_1] > 0)
      c(u1,u2)
    }))

    # if(TRUE){
    #    v <- runif(n1)
    # }else{
    #   v <- rep(0.5, n1)
    # }
    v <- runif(n1)
    r <- u[,1] * v + u[,2] * (1 - v)
    idx_adjust <- which(1-r < epsilon)
    r[idx_adjust] <- r[idx_adjust] - epsilon
    idx_adjust <- which(r < epsilon)
    r[idx_adjust] <- r[idx_adjust] + epsilon

    param1 <- cbind(param1,param_1)
    r1 <- cbind(r1,r)

    # print(i)
  }

  r1_1 <- apply(r1, 2, function(col) {
    if (any(col > 1)) {
      col[col > 1] <- max(col[col < 1])
    }
    return(col)
  })
  # quantile_normal <- qnorm(r1)
  quantile_normal <- qnorm(r1_1)
  cov_mat <- cor(quantile_normal)

  marginal <- 'zinb'
  marginal_result2 <- fit_marginals(t(W[,gene_sel2]), marginal, DT = FALSE)
  # marginal_result2 <- fit_marginals(t(W[,gene_sel1[c(111,112,141,149)]]), marginal, DT = FALSE)

  # rm(list=c("res","r1","r1_1"))
  return(list(cov_mat = cov_mat,
              marginal_param1 = param1,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_sample = n1, n_read = sum(W)))
}

# scDesign2
fit_marginals <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(gene) > 0)
        {
          mle_NB <- glm.nb(gene ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(gene ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(gene){
      c(0.0, Inf, mean(gene))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- x[iter, ]
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}

simulate_count_copula_3 <- function(copula_result,data_x,M){

  gene_sel1 <- copula_result$gene_sel1
  gene_sel2 <- copula_result$gene_sel2
  gene_sel3 <- copula_result$gene_sel3
  n1 <- copula_result$n_sample
  # param1 <- copula_result$marginal_param1

  X = as.matrix(cbind(1,data.frame(data_x)))

  p1 <- length(gene_sel1)
  result1 <- mvrnorm(n = n1, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
  result1 <- matrix(result1, nrow = n1)
  result2 <- apply(result1, 2, pnorm)
  result2 <- matrix(result2, nrow = n1)
  iter <- 2

  # t6=Sys.time()
  result31 <- sapply(1:p1, function(iter){
    param <- copula_result$marginal_param1[,iter]
    beta = param[1:4]
    beta_star = param[5:8]
    prob0 <- 1/(1+exp(-param[9]))
    lambda =  X %*% beta + log(M)  # exp(X %*% beta )*M
    theta = X %*% beta_star        # exp(X %*% beta_star)
    p_cdf_0 <- pmax(0.0, result2[, iter] - prob0) / (1-prob0)
    # iter2 <- 3
    W_simulate <- sapply(1:n1,function(iter2){
      qnbinom(p_cdf_0[iter2], size = 1/exp(theta[iter2]), mu = exp(lambda[iter2]))
    })
  })
  # t6=Sys.time()
  # print(t6-t5)

  p2 <- length(gene_sel2)
  if(p2 > 0){
    result32 <- sapply(1:p2, function(iter){
      param <- copula_result$marginal_param2[iter, ]
      rbinom(n1, 1, 1-param[1]) * rnbinom(n1, size = param[2], mu = param[3])
    })
  }

  result3 <- matrix(0, nrow = n1, ncol =p1 + p2 + length(gene_sel3) )
  if(p1 > 0){
    result3[, copula_result$gene_sel1] <- result31
  }
  if(p2 > 0){
    result3[, copula_result$gene_sel2] <- result32
  }
  rm(list=c("result2","result31","result32"))
  return(result3)
}

distance_JS <- function(P,Q){
  # P <- odd_vector
  # Q <- even_vector
  # if(length(P)!=length(Q)){
  #   len <- min(length(P),length(Q))
  #   P <- sort(P,decreasing = TRUE)[1:len]
  #   Q <- sort(Q,decreasing = TRUE)[1:len]
  # }
  # P <- sort(P[P!=0])
  # Q <- sort(Q[Q!=0])
  P <- sort(P)
  Q <- sort(Q)
  len <- max(length(P),length(Q))
  P <- c(rep(0,len-length(P)),P)
  Q <- c(rep(0,len-length(Q)),Q)

  avg <- (P+Q)/2
  avg <- avg/(sum(avg)+1e-05)
  P <- P/(sum(P)+1e-05)
  Q <- Q/(sum(Q)+1e-05)
  d1 <- distance(rbind(P,avg), method = "kullback-leibler", unit = "log2",mute.message = TRUE)
  d2 <- distance(rbind(Q,avg), method = "kullback-leibler", unit = "log2",mute.message = TRUE)
  d <- (d1+d2)/2
  d <- ifelse(is.na(d),0,d)
  return(d)
}

clipper_BC <- function(contrastScore, FDR){

  contrastScore[is.na(contrastScore)] = 0 # impute missing contrast scores with 0
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  i = 1
  emp_fdp = rep(NA, length(c_abs))
  emp_fdp[1] = 1
  while(i <= length(c_abs)){
    # print(i)
    t = c_abs[i]
    emp_fdp[i] = min((1 + sum(contrastScore <= -t))/ sum(contrastScore >= t),1)
    if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
    i = i + 1
  }

  c_abs = c_abs[!is.na(emp_fdp)]
  emp_fdp = emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] = 1

  re = lapply(FDR, function(FDR_i){
    thre = c_abs[min(which(emp_fdp <= FDR_i))]
    re_i = list(FDR = FDR_i,
                FDR_control = 'BC',
                thre = thre,
                q = q,
                discovery = which(contrastScore >= thre))
    return(re_i)
  })
  return(re)
}

clipper_BC_2 <- function(contrastScore, FDR){

  # contrastScore <- contrast_score
  contrastScore[is.na(contrastScore)] = 0 # impute missing contrast scores with 0
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  i = 1
  emp_fdp = rep(NA, length(c_abs))
  emp_fdp[1] = 1
  while(i <= length(c_abs)){
    # print(i)
    t = c_abs[i]
    emp_fdp[i] = min((1 + sum(contrastScore <= -t))/ sum(contrastScore >= t),1)
    if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
    i = i + 1
  }

  c_abs = c_abs[!is.na(emp_fdp)]
  emp_fdp = emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] = 1

  for(FDR_i in FDR){
    re_i = list(FDR = FDR_i,
                # thre = thre,
                q = q,
                discovery = NULL)
    if(sum(emp_fdp <= FDR_i)>0){
      thre = c_abs[min(which(emp_fdp <= FDR_i))]
      re_i = list(FDR = FDR_i,
                  thre = thre,
                  q = q,
                  discovery = which(contrastScore >= thre))
      # re_i = c(FDR_i,thre,round(which(contrastScore >= thre)))
      break
    }
  }
  return(re_i)

  # re = lapply(FDR, function(FDR_i){
  #   thre = c_abs[min(which(emp_fdp <= FDR_i))]
  #   re_i = list(FDR = FDR_i,
  #               FDR_control = 'BC',
  #               thre = thre,
  #               q = q,
  #               discovery = which(contrastScore >= thre))
  #   return(re_i)
  # })
  # return(re)
}

# clipper_pooling_intersection
clipper_p_i <- function(contrast_score_k,pooling = FALSE,intersection = FALSE,n_p_j,
                        statistics = 3,FDR = c(0.05,0.1,0.2,0.5,0.9,1)){

  if(pooling){
    result_fdr <- clipper_BC_2(contrast_score_k, FDR)
    a_1 <- result_fdr$discovery
  }else{
    if(intersection){
      result_intersect <- c(1:dim(contrast_score_k)[2])
      for(i in 1:dim(contrast_score_k)[1]){
        result_fdr_i <- clipper_BC_2(contrast_score_k[i,], FDR)
        result_intersect <- intersect(result_intersect, result_fdr_i$discovery)
      }
      a_1 <- result_intersect
    }else{
      contrast_score <- switch(statistics,
                               apply(contrast_score_k, 2, cumprod)[dim(contrast_score_k)[1],],
                               apply(contrast_score_k, 2, max),
                               apply(contrast_score_k, 2, sum))
      result_fdr <- clipper_BC_2(contrast_score, FDR)
      a_1 <- result_fdr$discovery
    }
  }

  b_1 <- 1:n_p_j
  fdr_moni <- length(a_1[!a_1 %in% b_1])/length(a_1)
  power_moni <- length(b_1[b_1 %in% a_1])/length(b_1)
  fdr_moni1 <- c(result_fdr$FDR,fdr_moni,power_moni)
  return(fdr_moni1)

}

clipper_p_i_2 <- function(contrast_score_k,pooling = FALSE,n_p_j,
                          statistics = 1,FDR = c(0.05,0.1,0.2,0.5,0.9,1)){
  b_1 <- 1:n_p_j
  if(pooling){
    result_fdr <- clipper_BC_2(contrast_score_k, FDR)
    a_1 <- result_fdr$discovery
    # fdr_moni <- length(a_1[!a_1 %in% b_1])/length(a_1)
    # power_moni <- length(b_1[b_1 %in% a_1])/length(b_1)
    fdr_moni <- sum(!(a_1 %in% b_1))/max(length(a_1),1)
    power_moni <- sum(b_1 %in% a_1)/length(b_1)
    fdr_moni1 <- list(pooling=c(result_fdr$FDR,fdr_moni,power_moni),pooling_p = a_1)

  }else{
    result_intersect <- c(1:dim(contrast_score_k)[2])
    for(i in 1:dim(contrast_score_k)[1]){
      result_fdr_i <- clipper_BC_2(contrast_score_k[i,], FDR)
      result_intersect <- intersect(result_intersect, result_fdr_i$discovery)
    }
    a_1_i <- result_intersect
    # fdr_moni_i <- length(a_1_i[!a_1_i %in% b_1])/length(a_1_i)
    # power_moni_i <- length(b_1[b_1 %in% a_1_i])/length(b_1)
    fdr_moni_i <- sum(!(a_1_i %in% b_1))/max(length(a_1_i),1)
    power_moni_i <- sum(b_1 %in% a_1_i)/length(b_1)

    contrast_score <- switch(statistics,
                             apply(contrast_score_k, 2, cumprod)[dim(contrast_score_k)[1],],
                             apply(contrast_score_k, 2, max),
                             apply(contrast_score_k, 2, sum))
    result_fdr <- clipper_BC_2(contrast_score, FDR)
    a_1_k <- result_fdr$discovery
    # fdr_moni_k <- length(a_1_k[!a_1_k %in% b_1])/length(a_1_k)
    # power_moni_k <- length(b_1[b_1 %in% a_1_k])/length(b_1)
    fdr_moni_k <- sum(!(a_1_k %in% b_1))/max(length(a_1_k),1)
    power_moni_k <- sum(b_1 %in% a_1_k)/length(b_1)
    # if(length(a_1_k)==0){
    #   fdr_moni_k <- power_moni_k <- 0
    # }

    fdr_moni1 <- list(intersect = c(result_fdr_i$FDR,fdr_moni_i,power_moni_i),intersect_p = a_1_i,
                      simultaneous = c(result_fdr$FDR,fdr_moni_k,power_moni_k),simultaneous_p = a_1_k)
  }
  return(fdr_moni1)
}

clipper_p_i_E <- function(E,pooling = FALSE,n_p_j,
                          statistics = 1,FDR = c(0.05,0.1,0.2,0.5,0.9,1)){
  b_1 <- 1:n_p_j
  if(pooling){
    result_fdr <- ebh_2(E, FDR)
    a_1 <- result_fdr$discovery
    fdr_moni <- sum(!(a_1 %in% b_1))/max(length(a_1),1)
    power_moni <- sum(b_1 %in% a_1)/length(b_1)
    # if(length(a_1)==0){
    #   fdr_moni <- power_moni <- 0
    # }
    fdr_moni1 <- list(pooling=c(result_fdr$FDR,fdr_moni,power_moni),pooling_p = a_1)

  }else{
    result_intersect <- c(1:dim(E)[2])
    for(i in 1:dim(E)[1]){
      result_fdr_i <- ebh_2(E[i,], FDR)
      result_intersect <- intersect(result_intersect, result_fdr_i$discovery)
    }
    a_1_i <- result_intersect
    fdr_moni_i <- sum(!(a_1_i %in% b_1))/max(length(a_1_i),1)
    power_moni_i <- sum(b_1 %in% a_1_i)/length(b_1)
    # if(length(a_1_i)==0){
    #   fdr_moni_i <- power_moni_i <- 0
    # }

    contrast_score <- switch(statistics,
                             apply(E, 2, cumprod)[dim(E)[1],],
                             apply(E, 2, max),
                             apply(E, 2, sum))
    result_fdr <- ebh_2(contrast_score, FDR)
    a_1_k <- result_fdr$discovery
    fdr_moni_k <- sum(!(a_1_k %in% b_1))/max(length(a_1_k),1)
    power_moni_k <- sum(b_1 %in% a_1_k)/length(b_1)
    # if(length(a_1_k)==0){
    #   fdr_moni_k <- power_moni_k <- 0
    # }

    fdr_moni1 <- list(intersect = c(result_fdr_i$FDR,fdr_moni_i,power_moni_i),intersect_p = a_1_i,
                      simultaneous = c(result_fdr$FDR,fdr_moni_k,power_moni_k),simultaneous_p = a_1_k)
  }
  return(fdr_moni1)
}

clipper_p_i_E_1 <- function(E,n_p_j,FDR = c(0.05,0.1,0.2,0.5,0.9,1)){
  b_1 <- 1:n_p_j
  result_fdr <- ebh_2(E, FDR)
  a_1_k <- result_fdr$discovery
  fdr_moni_k <- sum(!(a_1_k %in% b_1))/max(length(a_1_k),1)
  power_moni_k <- sum(b_1 %in% a_1_k)/length(b_1)
  fdr_moni1 <- list(simultaneous = c(result_fdr$FDR,fdr_moni_k,power_moni_k),simultaneous_p = a_1_k)
  return(fdr_moni1)
}


##########

SparseDOSSA_simulation <- function(W){
  # library(SparseDOSSA2)
  Stool_subset <- t(W)   # W[,50:70]
  n_w <- dim(W)[1]
  fitted <- fit_SparseDOSSA2(data = Stool_subset, control = list(verbose = TRUE))
  # fitted$EM_fit$fit$pi0
  # fitted$EM_fit$fit$mu
  # fitted$EM_fit$fit$sigma
  Stool_subset_simulation <- SparseDOSSA2(template = fitted, n_sample = n_w,
                                          new_features = FALSE, verbose = TRUE)
  Stool_simulation_data <- Stool_subset_simulation$simulated_data
  return(Stool_simulation_data)
  # setwd("E:/一些论文/假设检验/Knockoff方式/代码/result")
  # # save.image(file = "231205_1.RData")
  # load("231205_1.RData")
}

scDesign2_simulation <- function(W,class0){
  # library(scDesign2)
  count_mat <- t(W)  # W
  colnames(count_mat) <- class0
  cell_type_sel <- unique(class0)
  # submat <- get_submat(count_mat, cell_type_sel)
  n_cell_new <- ncol(count_mat)
  cell_type_prop <- table(colnames(count_mat))[cell_type_sel]
  copula_result <- fit_model_scDesign2(count_mat, cell_type_sel, sim_method = 'copula', marginal = 'zinb',
                                       ncores =1 ) #length(cell_type_sel)
  sim_count_copula <- simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
                                               # marginal = 'zinb',
                                               cell_type_prop = cell_type_prop)
  return(t(sim_count_copula))
}

konckoffs_simulation <- function(W, class0){
  # library(knockoff)
  # result <- knockoff.filter(W, class0,
  #                           knockoffs = create.second_order, statistic = stat.glmnet_coefdiff,
  #                           fdr = 0.2)
  result = knockoff.filter(W, class0,
                           knockoffs = create.second_order, statistic = stat.random_forest,
                           fdr=0.1)
  # result$selected
  # hist(result$statistic,breaks = 40)
  return(result)
}

contrast_score_computation <- function(W,result3,class0,test1){

  # result3 <- result_konckoffs$Xk

  test2 <- c("wilcox.test",'ks.test','mmd','distance_JS')
  test1 <- match(test1,test2)

  umap_result <- umap(result3, n_components = 3)
  # kmeans_result <- kmeans(umap_result$layout, centers = 2)
  kmeans_result <- kmeans(umap_result, centers = 2)
  class1 <- kmeans_result$cluster

  contrast_score <- sapply(1:dim(result3)[2], function(i_wi){
    # odd_vector <- result3[class1==1,i_Wi]
    # even_vector <- result3[class1==2,i_Wi]
    # KS.test(odd_vector,even_vector)$p.value
    P <- sort(result3[class1==1,i_wi])
    Q <- sort(result3[class1==2,i_wi])

    P1 <- sort(W[class0==1,i_wi])
    Q1 <- sort(W[class0==2,i_wi])

    if(sum(P) + sum(Q)==0 | sum(P1) + sum(Q1)==0){
      result_p <- 0
    }else{
      result_p <- switch(test1,
                         log10(wilcox.test(P,Q)$p.value/wilcox.test(P1,Q1)$p.value),
                         log10(ks.test(P,Q)$p.value/ks.test(P1,Q1)$p.value),
                         mmd(P1,Q1)$stat - mmd(P,Q)$stat,
                         distance_JS(P1,Q1) - distance_JS(P,Q))
    }
    # result_p <- switch(test1,
    #                    format(wilcox.test(P,Q)$p.value, scientific = FALSE),
    #                    format(ks.test(P,Q)$p.value, scientific = FALSE),
    #                    format(mmd(P,Q)$stat, scientific = FALSE),
    #                    format(distance_JS(P,Q), scientific = FALSE))

  })
  # hist(wi_W1,breaks = 40)
  # hist(contrast_score,breaks = 40)
  # result_fdr <- clipper_BC_2(contrast_score,0.5)
  # result_fdr$discovery
  return(contrast_score)
}


FDR_Power <- function(a_1,b_1 = NULL){
  if(is.null(b_1)){
    return(list(FP = c(-1,-1)))
  }else{
    if(length(b_1)==1){
      b_1 <- 1:b_1
    }
    fdr_moni_k <- sum(!(a_1 %in% b_1))/max(length(a_1),1)
    power_moni_k <- sum(b_1 %in% a_1)/length(b_1)
    return(list(FP = c(fdr_moni_k,power_moni_k)))
  }
}

inter_cw <- function(c_w,b_1 = NULL,fdr = 0.2){
  result_intersect <- c(1:dim(c_w)[2])
  for(i in 1:dim(c_w)[1]){
    result_fdr_i <- clipper_BC_2(c_w[i,], fdr)
    result_intersect <- intersect(result_intersect, result_fdr_i$discovery)
  }
  # if(is.null(b_1)){
  #   b_1 = 1:100
  # }
  # a_1_i <- result_intersect
  # fdr_moni_i <- sum(!(a_1_i %in% b_1))/max(length(a_1_i),1)
  # power_moni_i <- sum(b_1 %in% a_1_i)/length(b_1)

  a_1_i <- result_intersect
  FDRPower <- FDR_Power(a_1_i,b_1)
  return(list(S = a_1_i, FDRPower = c(result_fdr_i$FDR,FDRPower$FP)))
}

filter_cw <- function(c_w,filter_statistics = 1,fdr = 0.2,b_1 = NULL){
  c_w_b <- switch(filter_statistics,
                  apply(c_w, 2, cumprod)[dim(c_w)[1],],
                  apply(c_w, 2, max),
                  apply(c_w, 2, sum))
  result_fdr <- clipper_BC_2(c_w_b, fdr)
  S <- result_fdr$discovery
  FDRPower <- FDR_Power(S,b_1)
  FDRPower <- c(result_fdr$FDR,FDRPower$FP)
  return(list(S = S,FDRPower = FDRPower))
}

filter_cw_B <- function(Bstat = 1, combine_1 = NULL, c_w = NULL, e_w = NULL,
                        filter_statistics = 1,B = 30,fdr = 0.2,T_var = NULL,K = 1){
  if(Bstat == 1){
    # B test_statistics, mean; First k and then B; Finally K row E-value;
    if(combine_1 == 'simul'){
      e_value_B <- switch(filter_statistics,
                          apply(e_w, 2, cumprod)[dim(e_w)[1],],
                          apply(e_w, 2, max),
                          apply(e_w, 2, sum))
      result_fdr <- clipper_BC_2(e_value_B, fdr)
      S <- result_fdr$discovery

      # FDRPower = NULL
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S,T_var)
      FDRPower <- c(result_fdr$FDR,FDRPower$FP)
      return(list(fdr = fdr, e_value = e_w, e_value_B = e_value_B,
                  S = S, FDRPower = FDRPower))
    }else if(combine_1 == 'fisher'){
      p_e_w <- e_w
      p_e_w[p_e_w == 0] <- 1
      p_comb <- apply(-2*log(1/p_e_w),2,sum)
      chis <- qchisq(fdr, 2*K, lower.tail = FALSE)
      S <- which(p_comb>chis)

      # FDRPower = NULL
      # FDRPower$FP <- c(-1,-1)
      # if(!is.null(T_var)){
      #   FDRPower <- FDR_Power(S,T_var)
      # }
      FDRPower <- FDR_Power(S,T_var)
      FDRPower <- c(fdr,FDRPower$FP)
      return(list(fdr = fdr, e_value = e_w, p_comb = p_comb,
                  S = S, FDRPower = FDRPower))
    }
  }else if(Bstat == 2){
    # B filter_statistics, mean; First b and then K; Finally B row E-value;
    e_w_B <- c()
    for(b in 1:B){
      c_w_bK <- c_w[seq(b,length(name_data)*B+b-1,B),] # t:11:((n_k*T)+t-1)
      c_w_b <- switch(filter_statistics,
                      apply(c_w_bK, 2, cumprod)[dim(c_w_bK)[1],],
                      apply(c_w_bK, 2, max),
                      apply(c_w_bK, 2, sum))
      gamma <- fdr/2
      offset <- 1
      e_w_b <- ekn(c_w_b, gamma, offset)
      e_w_B <- rbind(e_w_B,e_w_b)
    }
    e_value_B <- apply(e_w_B,2,mean)
    result_ebh <- ebh_2(e_value_B, fdr)
    S <- result_ebh$discovery
    # FDRPower = NULL
    # FDRPower$FP <- c(-1,-1)
    # if(!is.null(T_var)){
    #   FDRPower <- FDR_Power(S,T_var)
    # }
    FDRPower <- FDR_Power(S,T_var)
    FDRPower <- c(result_ebh$FDR,FDRPower$FP)
    return(list(fdr = fdr, c_w_pis = c_w, e_w_B = e_w_B, e_value = e_value_B,
                S = S, FDRPower = FDRPower))
  }
}

# https://github.com/zhimeir/derandomized_knockoffs_fdr/blob/main/utils/utils.R

ebh <- function(E, alpha){
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(rej)
}

ebh_2 <- function(E, alpha){
  for(FDR_i in alpha){
    rej <- ebh(E, FDR_i)
    if(!is.null(rej) | FDR_i ==1 ){
      break
    }
  }
  return(list(discovery = rej,FDR = FDR_i))
}

ekn <- function(C, gamma, offset){

  ts = sort(c(0, abs(C)))
  ratio = sapply(ts, function(t)
    (offset + sum(C <= -t)) / max(1, sum(C >= t)))
  ok = which(ratio <= gamma)
  tau <- ifelse(length(ok) > 0, ts[ok[1]], Inf)
  # tau <- alphakn_threshold(C, fdr =  gamma, offset = offset)

  ord_C <- order(abs(C), decreasing = TRUE)
  sorted_C <- C[ord_C]
  if(sum(C>0) >= 1 / gamma){
    pos_ind <- which(sorted_C > 0)
    tau1 <- sorted_C[pos_ind[ceiling(1/gamma)-1]]
  }else{
    tau1 <- 0
  }
  tau <- min(tau,tau1)

  E  <- length(C)*(C >= tau) / (1 + sum(C <= -tau))


  return(E)
}
