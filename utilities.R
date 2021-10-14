Generate <- function(opt) {
  # point mass
  # bimodal unimod
  # homo vs hetero
  set.seed(opt$seed)
  N <- opt$N
  k <- as.integer(N*(1-opt$w))
  V <- opt$V
  is_homo <- opt$homo
  setting <- opt$setting
  print(opt)
  # unimod
  if (setting == "a") {
    mu <- c(rep(V, k), rep(0, N-k))
    if(is_homo) {
      sds <- rep(1,N)
    } else {
      sds <- invgamma::rinvgamma(N, shape = 2, rate = 1)
    }
    y <- c(rnorm(k, mean = mu[1:k], sd = sds[1:k]), 
           rnorm(N-k, mean = mu[(k+1):N], sd = sds[(k+1):N]))
  }
  
  # multimod
  if (setting == "b") {
    mu <- unlist(mapply(rep, V, k))
    if(is_homo) {
      sds <- rep(1,N)
    } else {
      sds <- invgamma::rinvgamma(N, shape = 2, rate = 1)
    }
    y <- rnorm(N, mean = mu, sd = sds)
  }
  
  if (setting == "c") {
    mu <- c(rnorm(k, rep(V, k), 1), rep(0, N-k))
    if(is_homo) {
      sds <- rep(1,N)
    } else {
      sds <- runif(N, min = 0.5, max = opt$umax)
    }
    y <- rnorm(N, mean = mu, sd = sqrt(sds))
  }
  
  if (setting == "d") {
    
    idx <- sample(1:2, prob=c(0.5, 0.5), size=N, replace=TRUE)
    rate <- c(20,25)
    tau <- rgamma(N, shape = 20, rate = rate[idx])
    sds <- sqrt(1/tau)
    mu <- c(rnorm(k, rep(V, k)), rep(0, N-k))/tau
    
    y <- rnorm(N, mean = mu, sd = sds)
  }
  
  return(list(y = y, mu = mu, sds = sds))
}

data_subset <- function(data, first=1:3, second=4:6, threshold=11, pitcher=NULL) {
  # the function split the data into training and test
  # and exclude players whose Hits less than 11 in either season
  if(!is.data.table(data)) {data <- setDT(data)}
  if(!is.null(pitcher)) {
    if(pitcher) {data <- data[`Pitcher?`==1]}
    else {data <- data[`Pitcher?`==0]}
  }
  first_half  <- data[month%in%first,.(AB=sum(AB), H=sum(H)), by=id]
  second_half <- data[month%in%second,.(AB=sum(AB), H=sum(H)), by=id]
  idx <- intersect(first_half[AB>=threshold, id],  second_half[AB>=threshold, id])
  list(first_half = first_half[AB>=threshold],
       second_half = second_half[id%in%idx])
}

arcsin_transform <- function(H, N, a) {
  asin(sqrt((H + a)/(N + 2*a)))
}

estimate <- function(formula, N, data, estimator, distribution="normal", ...) {
  vars <- as.character(all.vars(formula))
  # ellipsis <- list(...)
  est <- estimate_(data[[vars[1]]], data[[N]], estimator, distribution=distribution, ...)
  if(estimator!="neb") {
    ret <- list(estimate = data.frame(id = data[[vars[2]]], 
                                      estimate = est$estimate,
                                      N_1 = data[[N]],
                                      H_1 = data$H))  
  } else if(estimator=="neb") {
    ret <- list(estimate = data.frame(id = data[[vars[2]]], 
                                      estimate = est$posterior_means,
                                      N_1 = data[[N]],
                                      H_1 = data$H),
                pi_hat = est$pi_hat, 
                tau_grid = est$tau_grid, 
                posterior_hat = est$posterior_hat)
  }
  class(ret) <- "estimate"
  ret
} 

estimate_ <- function(x, N, estimator, distribution, ...) {
  estimate <- switch(estimator,
                     naive      = naive_est(x), 
                     group_mean = group_mean(x), 
                     ebmm       = ebmm(x),
                     ebml       = ebml(x),
                     npeb       = npeb(x, N),
                     npebc      = npebc(x, N),
                     neb        = neb_aux(x, N, 
                                          distribution=distribution, ...),
                     js         = james_stein(x, N))
}

naive_est <- function(data) {
  list(estimate = data)
}

group_mean <- function(data) {
  list(estimate = rep(mean(data), length(data)))
}

# TODO
ebmm <- function(data, N) {
  # parametric empirical Bayes (Method of Moments) estimator
  # init
  mu <- mean(data)
  n <- length(data)
  tmp <- sum((data-mu)^2) * (-(n-1)/n) * sum(1/(4*N)) 
  tau2 <- tmp*(tmp>0)/(n-1)
  tmp <- 
    mu <- sum() / sum ()
  
  mu + tau2/(tau2+sigma2)*(data - mu)
}

# TODO
ebml <- function(data, N) {
  
}

gx <- function(data, i, sigma2, h) {
  si <- sigma2[i]
  idx <- which((1+h)*si>sigma2)
  h <- sqrt((1+h)*pmax(si, sigma2) - sigma2)
  gx <- sum(1/h * dnorm((data[i]- data)/h)) / length(idx)
  gxprime <- sum(1/h * (data-data[i])/h^2 * dnorm((data[i]- data)/h)) / length(idx)
  data[i] + si * gxprime/gx
}

npeb <- function(data, N) {
  n <- length(data)
  # following the 2008 paper
  # h=.25 for n>200 and h=.3 for n=81
  h <- ifelse(n>200, .25, .3)
  sigma2 <- 1/(4*N)
  # sigma2 <- data*(1-data)/N
  list(estimate = sapply(seq_along(data), function(i) gx(data, i, sigma2, h)))
}

npebc <- function(data, N) {
  n <- length(data)
  # following the 2008 paper
  # h=.25 for n>200 and h=.3 for n=81
  h <- ifelse(n>200, .25, .3)
  sigma2 <- data*(1-data)/N
  sigma2[data==0] <- 1/(4*N[data==0])
  list(estimate = sapply(seq_along(data), function(i) gx(data, i, sigma2, h)))
}

func <- function(...) {
  input_list <- list(...)
  sparse <- input_list$sparse
  if(is.null(sparse)) {print("NULL")}
  else if(sparse) {print("yes")}
}

neb_aux <- function(data, N, distribution, ...) {
  if (distribution == "normal") {
    neb(data, N, distribution = "normal", ...)
  } else if (distribution == "binomial") {
    neb(data, N, distribution = "binomial", ...)
  }
}

neb <- function(data, N, posterior="mean",  distribution, ...) {
  input_list <- list(...)
  sparse <- ifelse(is.null(input_list$sparse), F, T)
  full <- ifelse(is.null(input_list$full), F, T)
  esp_gap <- ifelse(is.null(input_list$esp_gap), 0, input_list$esp_gap)
  details <- ifelse(is.null(input_list$details), F, T)
  transformed <- ifelse(is.null(input_list$transformed), F, T)
  const_var <- ifelse(is.null(input_list$const_var), F, T)
  
  ifelse(sparse, print("sparsity"), print("no sparsity"))
  ifelse(full, print("full sample"), print("split sample"))
  ifelse(details, print("return details"), print("return posterior means"))
  ifelse(transformed, print("transformed"), print("original scale"))
  ifelse(const_var, print("const var"), print("plug in var"))
  print(paste0("gap around 0: ", esp_gap, sep=""))
  
  
  n <- length(data)
  
  # Split sample
  set.seed(123)
  itest <- sample(1:n, n/3)
  itrain <- setdiff(1:n, itest)
  N_train <- length(itrain)
  N_test <- length(itest)
  
  # The grid points. The NPMLE is always in the range of the observations.
  if(distribution=="normal") {
    if(transformed) {
      sds <- sqrt(1/(4*N))
    } else {
      if(const_var) {
        sds <- sqrt(1/(4*N))
      } else {
        sds <- sqrt(data*(1-data)/N)
        sds[data==0] <- sqrt(1/(4*N[data==0]))  
      }
    }
    
    a <- min(data)
    b <- max(data)
    R <- b-a
    M <- 10*n
    # lag <- abs(a)/((abs(a)/R*M))
    # tau_grid <- c(seq(a, 0, length.out = ceiling(abs(a)/R*M)),
    #               seq(lag, b, lag))
    # tau_grid <- sort(c(0, seq(a, b, length.out = M-1)))
    if(esp_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(tau_grid[tau_grid>esp_gap],tau_grid[tau_grid< (-1*esp_gap)])
      tau_grid <- sort(c(tau_grid,0))
      M <- length(tau_grid)
    }
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    
    for(i in 1:n) {
      likelihood[, i] <- dnorm(data[i], tau_grid, sds[i])
    }
  } else if (distribution=="binomial") {
    p <- data/N
    
    if(const_var) {
      sds <- sqrt(1/(4*N))
    } else {
      sds <- sqrt(p*(1-p)/N)
      sds[data==0] <- sqrt(1/(4*N[data==0])) 
    }
    
    a <- min(p, na.rm=T)
    b <- max(p, na.rm=T)
    M <- 10*n
    if(esp_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(0,tau_grid[tau_grid>esp_gap])
      M <- length(tau_grid)
    }
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    for(i in 1:n) {
      likelihood[, i] <- dbinom(data[i], size=N[i], tau_grid)
    }
  }
  
  if(!full) { # split sample
    ret <- C_neb(likelihood, itrain, itest, tau_grid, pi_hat, miter = 50, esp=.00001) 
  } else {
    ret <- C_neb_full(likelihood, tau_grid, pi_hat, miter = 100, esp=.00001)
  }
  
  ret$tau_grid <- tau_grid
  ret
  # if(details) {
  #   ret$tau_grid = tau_grid
  #   return(ret)} 
  # else {return(ret$posterior_means)}
}

neb_old <- function(data, N, posterior="mean", sparse=F, distribution) {
  n <- length(data)
  
  # Split sample
  set.seed(123)
  itest <- sample(1:n, n/3)
  itrain <- setdiff(1:n, itest)
  N_train <- length(itrain)
  N_test <- length(itest)
  
  # The grid points. The NPMLE is always in the range of the observations.
  if(distribution=="normal") {
    sds <- sqrt(1/(4*N))
    # sds <- sqrt(data*(1-data)/N)
    
    a <- min(data)
    b <- max(data)
    R <- b-a
    M <- 10*n
    # lag <- abs(a)/((abs(a)/R*M))
    # tau_grid <- c(seq(a, 0, length.out = ceiling(abs(a)/R*M)),
    #               seq(lag, b, lag))
    # tau_grid <- sort(c(0, seq(a, b, length.out = M-1)))
    tau_grid <- seq(a, b, length.out = M)
    tau_grid[which.min(abs(tau_grid))] <- 0
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    
    for(i in 1:n) {
      likelihood[, i] <- dnorm(data[i], tau_grid, sds[i])
    }
  } else if (distribution=="binomial") {
    p <- data/N
    sds <- p*(1-p)/N
    
    a <- min(p, na.rm=T)
    b <- max(p, na.rm=T)
    M <- 10*n
    tau_grid <- seq(a, b, length.out = M)
    tau_grid[which.min(abs(tau_grid))] <- 0
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    for(i in 1:n) {
      likelihood[, i] <- dbinom(data[i], size=N[i], tau_grid)
    }
  }
  
  
  lik_train <- likelihood[, itrain]
  lik_test <- likelihood[, itest]
  loglik_train <- loglik_func(lik_train, pi_hat)/N_train
  loglik_test <- loglik_func(lik_test, pi_hat)/N_test
  loglik_train_hist <- loglik_train
  loglik_test_hist <- loglik_test
  
  loglik(lik_train, pi_hat)/N_train
  loglik(lik_test, pi_hat)/N_train
  # 
  # pdfY <- c(pdfT %*% likelihood * n)
  
  # plot(pdfY, type = 'l')
  
  # results matrix
  ntest <- 100
  MSE <- rep(NA,ntest)
  final_iters <- rep(NA,ntest)
  final_train_loglik <- rep(NA,ntest)
  final_test_loglik <- rep(NA,ntest)
  
  # initialization
  iters <- 1
  total_diff <- 999 
  stop_rule <- 1
  # par(mfrow=c(1,2))
  while((total_diff>.0001) & (stop_rule!=0)){
    # if(iters%%50==0) {print(iters)}
    # unscaled posterior
    posterior_hat_train <- lik_train * many_cols(pi_hat, N_train)
    posterior_hat_test <- lik_test * many_cols(pi_hat, N_test)
    posterior_hat <- likelihood * many_cols(pi_hat, n)
    
    # Bayes estimates and MSE
    # posterior_means_train <- apply(many_cols(tau_grid,N_train)*posterior_hat_train/
    #                                  many_rows(apply(posterior_hat_train,2,sum),M),2,sum)
    # posterior_means_test <- apply(many_cols(tau_grid,N_test)*posterior_hat_test/
    #                                 many_rows(apply(posterior_hat_test,2,sum),M),2,sum)
    posterior_means <- apply(many_cols(tau_grid,n)*posterior_hat/
                               many_rows(apply(posterior_hat,2,sum),M),2,sum)
    # MSE_new <- mean((Th-posterior_means)^2)
    sum(abs(est_posterior_means(posterior_hat, tau_grid) - posterior_means))
    
    
    # update pi
    denom <- colSums(posterior_hat_train)
    denom_mtx <- many_rows(denom, M)
    pi_hat_new <- 1/N_train * rowSums(posterior_hat_train / denom_mtx)
    pi_hat_old <- pi_hat
    pi_hat <- pi_hat_new/sum(pi_hat_new)
    
    loglik_train <- loglik_func(lik_train, pi_hat)/N_train
    loglik_test <- loglik_func(lik_test, pi_hat)/N_test
    loglik_train_hist <- c(loglik_train_hist,loglik_train)
    loglik_test_hist <- c(loglik_test_hist,loglik_test)
    
    # early stopping
    if(((loglik_test <= (loglik_test_hist[length(loglik_test_hist)-1])+.00001))|(iters>50)){
      stop_rule <- 0
      final_iters[iters] <- iters
      final_train_loglik[iters] <- loglik_train_hist[length(loglik_train_hist)-1]
      final_test_loglik[iters] <- loglik_test_hist[length(loglik_test_hist)-1]
      # MSE[k] <- mean((Th-posterior_means)^2)
    }
    
    # MSE_old <- MSE_new
    
    # plot(x=tau_grid,data=pdfT,type="l",main=paste("Estimated versus Actual P(theta), iter", iters), lwd=2)
    # lines(x=tau_grid,data=pi_hat, col=2,lwd=2)
    # legend("bottomright", c("Estimate", "Actual"), lwd=rep(2,2), col=(1:2)) 
    iters <- iters+1
    # print(MSE_new)
  }
  # if(distribution=="binomial") {
  #   posterior_means <- posterior_means/N
  # }
  posterior_means
  
  print(iters)
  
  # sourceCpp("src/neb.cpp")
  # pi_hat <- rep(NA,M)
  # pi_hat[tau_grid==0] <- w0_hat
  # pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
  # ret = C_neb(likelihood, itrain, itest, tau_grid, pi_hat, miter = 50, esp=.00001)
  # sum(abs(ret$posterior_means - posterior_means ))
  ret <- NULL
  ret$pi_hat <- pi_hat
  ret$tau_grid <- tau_grid
  ret$posterior_means <- posterior_means
  ret$posterior_hat <- posterior_hat
  ret$loglik_train <- loglik_train_hist
  ret$loglik_test <- loglik_test_hist
  
  ret
}

james_stein <- function(data, N) {
  n <- length(data)
  sigma2 <- 1/(4*N)
  mu <- sum(data/sigma2) / sum(1/sigma2)
  tmp <- 1 - (n-3)/(sum((data - mu)^2/sigma2))
  list(estimate = mu + tmp*(tmp>0)*(data-mu))  
}

# predict.estimate <- function(est, datanew, transformed=T) {
#   if(class(est)!="estimate") stop("estimator is not of class estimate")
#   pred <- merge(est$estimate, datanew, by="id")
#   tsehat <- tse(pred, transformed)
#   pred_err <- pred_error(pred, transformed)
#   list(prediction = pred,
#        tse = tsehat,
#        pred_error = pred_err)
# }

tse <- function(prediction, transformed) {
  if (transformed) {
    tse <- with(prediction, sum((estimate-x)^2) - sum(1/(4*AB))) 
  } else {
    tse <- with(prediction, sum((estimate-p)^2) - sum(p*(1-p)/AB)) 
  }
  tse
}

pred_error <- function(prediction, transformed) {
  if (transformed) {
    pred_err <- with(prediction, sum((estimate-x)^2)) 
  } else {
    pred_err <- with(prediction, sum((estimate-p)^2)) 
  }
  pred_err
}

predictive_int <- function(est, datanew, nsample=10000, distribution="normal") {
  if(class(est)!="estimate") stop("estimator is not of class estimate")
  pred <- merge(est$estimate, datanew, by="id")
  pred_idx <- which(est$estimate$id %in% pred$id)
  n <- length(pred_idx)
  est_sub <- est$estimate[pred_idx,]
  upper_bound <- rep(0, n)
  lower_bound <- rep(0, n)
  for(i in 1:n) {
    ynew <- sample(est$tau_grid, nsample, prob = est$posterior_hat[,pred_idx[i]], replace = T)
    if (distribution=="normal") {
      # pnew <- rnorm(nsample, ynew, sqrt(1/4/est_sub$N[i]))
      p <- est_sub$H_1[i] / est_sub$N_1[i]
      sd <- ifelse(p, 
                   sqrt((1-p)*p/est_sub$N_1[i]),
                   sqrt(1/(4*est_sub$N_1[i])))
      pnew <- rnorm(nsample, ynew, sd)
      
      # pnew[pnew<0] <- 0
      # if(min(pnew)<0) {
      #   ratio <-(max(pnew)-min(pnew))/max(pnew)
      #   pnew <- max(pnew) - (max(pnew) - pnew)/ratio
      # }
      # if(min(pnew)<0) {
      #   pnew <- sort(pnew)
      #   p0 <- min(which(pnew>=0))
      #   dist <- abs(pnew[1:p0])
      #   pnew[1:p0] <- 0
      #   pnew[nsample:(nsample-p0+1)] <- pnew[nsample:(nsample-p0+1)] - dist
      # }
      if(min(pnew)<0) {
        # mpnew <- median(pnew)
        mpnew <- quantile(pnew, .1)
        ratio <- (mpnew - min(pnew))/mpnew
        pnew[pnew<mpnew] <- mpnew - (mpnew - pnew[pnew<mpnew])/ratio
      }
      # if(min(pnew)<0) {
      #   # mpnew <- median(pnew)
      #   mpnew <- quantile(pnew, .1)
      #   mmpnew <- quantile(pnew, .9)
      #   ratio <- (mpnew - min(pnew))/mpnew
      #   pnew <- sort(pnew)
      #   dist <- (mpnew - pnew[pnew<mpnew])/ratio
      #   pnew[pnew>mmpnew] <- pnew[pnew>mmpnew] + dist
      #   pnew[pnew<mpnew] <- mpnew - (mpnew - pnew[pnew<mpnew])/ratio
      # }
      # hist(pnew, breaks=100)
      
    } else if (distribution=="binomial") {
      # CHEATING! assume we have a perfect estiamte for N_{2i}
      pnew <- rbinom(nsample, datanew$AB[i], ynew)
      # pnew <- rbinom(nsample, est_sub$N[i], ynew)  
    }
    qs <- quantile(pnew, probs = c(0.025, 0.975)) 
    lower_bound[i] <- qs[1]
    upper_bound[i] <- qs[2]
  }
  list(ub = upper_bound,
       lb = lower_bound)
}

posterior_int <- function(est, datanew=NULL, nsample=10000) {
  if(class(est)!="estimate") stop("estimator is not of class estimate")
  if(!is.null(datanew)) {
    pred <- merge(est$estimate, datanew, by="id")
    pred_idx <- which(est$estimate$id %in% pred$id)
    n <- length(pred_idx)
    est_sub <- est$estimate[pred_idx,]
    upper_bound <- rep(0, n)
    lower_bound <- rep(0, n)
    for(i in 1:n) {
      ynew <- sample(est$tau_grid, nsample, prob = est$posterior_hat[,pred_idx[i]], replace = T)
      qs <- quantile(ynew, probs = c(0.025, 0.975)) 
      lower_bound[i] <- qs[1]
      upper_bound[i] <- qs[2]
    }
  } else {
    n <- nrow(est$estimate)
    upper_bound <- rep(0, n)
    lower_bound <- rep(0, n)
    for(i in 1:n) {
      ynew <- sample(est$tau_grid, nsample, prob = est$posterior_hat[,i], replace = T)
      qs <- quantile(ynew, probs = c(0.025, 0.975)) 
      lower_bound[i] <- qs[1]
      upper_bound[i] <- qs[2]
    }
  }
  list(ub = upper_bound,
       lb = lower_bound)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

neb_pi0 <- function(data, sds, opt, sparsity = NULL, 
                    posterior="mean",  distribution="normal", ...) {
  input_list <- list(...)
  eps_gap <- ifelse(is.null(input_list$eps_gap), 0, input_list$eps_gap)
  grid <- ifelse(is.null(input_list$grid), 10*length(data), input_list$grid)
  
  n <- length(data)
  
  # The grid points. The NPMLE is always in the range of the observations.
  if(distribution=="normal") {
    
    a <- quantile(data, 0.01) - 2*sd(data)
    b <- quantile(data, 0.99) + 2*sd(data)
    R <- b-a
    M <- grid
    
    # tau grid
    if(eps_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(tau_grid[tau_grid>eps_gap],tau_grid[tau_grid< (-1*eps_gap)])
      tau_grid <- sort(c(tau_grid,0))
      M <- length(tau_grid)
    }
    
    if(opt$setting == "a") {
      tau_grid[which.min(abs(tau_grid - opt$V))] <- opt$V
    }
    
    # pi0
    pi0 <- rep(0,M)
    pi0[which(tau_grid == 0)] <- opt$w
    if(opt$setting == "a") {
      pi0[which(tau_grid == opt$V)] <- 1-opt$w
    }
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    
    for(i in 1:n) {
      likelihood[, i] <- dnorm(data[i], tau_grid, sds[i])
    }
  } else if (distribution=="binomial") {
    p <- data/N
    
    if(const_var) {
      sds <- sqrt(1/(4*N))
    } else {
      sds <- sqrt(p*(1-p)/N)
      sds[data==0] <- sqrt(1/(4*N[data==0])) 
    }
    
    a <- min(p, na.rm=T)
    b <- max(p, na.rm=T)
    M <- 10*n
    if(eps_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(0,tau_grid[tau_grid>eps_gap])
      M <- length(tau_grid)
    }
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    for(i in 1:n) {
      likelihood[, i] <- dbinom(data[i], size=N[i], tau_grid)
    }
  }
  
  ret <- NULL
  ret$posterior_hat <- est_posterior_hat(likelihood_mat = likelihood, pi_hat = pi0)
  ret$posterior_hat <- sweep(ret$posterior_hat, 2, colSums(ret$posterior_hat), "/")
  ret$posterior_means <- est_posterior_means(ret$posterior_hat, tau_grid)
  ret$tau_grid <- tau_grid
  ret$y <- data
  ret$w <- opt$w
  ret$iter <- 0
  ret$posterior_mode <- ret$tau_grid[sapply(1:n, function(x) which.max(ret$posterior_hat[,x]))]
  
  
  class(ret) <- "neb_pi0" 
  
  ret
  # if(details) {
  #   ret$tau_grid = tau_grid
  #   return(ret)} 
  # else {return(ret$posterior_means)}
}



neb_0 <- function(data, sds, sparsity = NULL, 
                  posterior="mean",  distribution="normal", ...) {
  input_list <- list(...)
  sparse <- ifelse(is.null(input_list$sparse), F, T)
  full <- ifelse(is.null(input_list$full), F, input_list$full)
  split_k <- ifelse(is.null(input_list$split_k), 3, input_list$split_k)
  eps_gap <- ifelse(is.null(input_list$eps_gap), 0, input_list$eps_gap)
  details <- ifelse(is.null(input_list$details), F, T)
  transformed <- ifelse(is.null(input_list$transformed), F, T)
  const_var <- ifelse(is.null(input_list$const_var), F, T)
  miter <- ifelse(is.null(input_list$miter), 50, input_list$miter)
  grid <- ifelse(is.null(input_list$grid), 10*length(data), input_list$grid)
  eps <- ifelse(is.null(input_list$eps), .00001, input_list$eps)
  
  ifelse(sparse, print("sparsity"), print("no sparsity"))
  ifelse(is.null(sparsity), print("unknown sparsity"), print("known sparsity"))
  ifelse(full, print("full sample"), print("split sample"))
  ifelse(details, print("return details"), print("return posterior means"))
  ifelse(transformed, print("transformed"), print("original scale"))
  ifelse(const_var, print("const var"), print("plug in var"))
  print(paste0("Max iterations: ", miter, sep=""))
  print(paste0("Grid points: ", grid, sep=""))
  print(paste0("Split : ", split_k, sep=""))
  print(paste0("gap around 0: ", eps_gap, sep=""))
  
  
  n <- length(data)
  
  # Split sample
  set.seed(123)
  itest <- sample(1:n, n/split_k)
  itrain <- setdiff(1:n, itest)
  N_train <- length(itrain)
  N_test <- length(itest)
  
  # The grid points. The NPMLE is always in the range of the observations.
  if(distribution=="normal") {
    
    a <- quantile(data, 0.01) - 2*sd(data)
    b <- quantile(data, 0.99) + 2*sd(data)
    R <- b-a
    M <- grid
    # lag <- abs(a)/((abs(a)/R*M))
    # tau_grid <- c(seq(a, 0, length.out = ceiling(abs(a)/R*M)),
    #               seq(lag, b, lag))
    # tau_grid <- sort(c(0, seq(a, b, length.out = M-1)))
    if(eps_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(tau_grid[tau_grid>eps_gap],tau_grid[tau_grid< (-1*eps_gap)])
      tau_grid <- sort(c(tau_grid,0))
      M <- length(tau_grid)
    }
    
    # Initialization with fourier method
    if(!is.null(sparsity)) {
      w <- sparsity
    } else {
      # unknown w
      if(length(unique(sds))==1) {
        w <- max(w0_func(data, homosk=TRUE,sigsqy_homo=1),1/M)
      } else {
        w <- max(w0_func(data, homosk=F, sigsq_y = sds^2),1/M) 
      }
    }
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w
    pi_hat[tau_grid!=0] <- rep((1-w)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    
    for(i in 1:n) {
      likelihood[, i] <- dnorm(data[i], tau_grid, sds[i])
    }
  } else if (distribution=="binomial") {
    p <- data/N
    
    if(const_var) {
      sds <- sqrt(1/(4*N))
    } else {
      sds <- sqrt(p*(1-p)/N)
      sds[data==0] <- sqrt(1/(4*N[data==0])) 
    }
    
    a <- min(p, na.rm=T)
    b <- max(p, na.rm=T)
    M <- 10*n
    if(eps_gap==0) {
      tau_grid <- seq(a, b, length.out = M)
      tau_grid[which.min(abs(tau_grid))] <- 0 
    } else {
      tau_grid <- c(0,tau_grid[tau_grid>eps_gap])
      M <- length(tau_grid)
    }
    
    # Initialization with fourier method
    w0_hat <- 0
    if(sparse) {w0_hat <- max(w0_func(data, homosk=F,sigsq_y=sds^2),1/M)}
    # w0_hat <- max(w0_func(data, homosk=T),1/M)
    pi_hat <- rep(NA,M)
    pi_hat[tau_grid==0] <- w0_hat
    pi_hat[tau_grid!=0] <- rep((1-w0_hat)/(M-1),M-1)
    
    # likelihood
    likelihood <- matrix(NA, nrow=M, ncol=n)
    for(i in 1:n) {
      likelihood[, i] <- dbinom(data[i], size=N[i], tau_grid)
    }
  }
  
  if(!full) { # split sample
    ret <- C_neb(likelihood, itrain, itest, tau_grid, pi_hat, miter = miter, esp=eps) 
  } else {
    ret <- C_neb_full(likelihood, tau_grid, pi_hat, miter = miter, esp=eps)
  }
  
  ret$posterior_hat <- sweep(ret$posterior_hat, 2, colSums(ret$posterior_hat), "/")
  ret$tau_grid <- tau_grid
  ret$y <- data
  ret$w <- ret$pi_hat[tau_grid==0]
  ret$posterior_mode <- ret$tau_grid[sapply(1:n, function(x) which.max(ret$posterior_hat[,x]))]
  
  if(!full) {
    class(ret) <- "neb_dnp" 
  } else {
    class(ret) <- "neb_0"
  }
  ret
  # if(details) {
  #   ret$tau_grid = tau_grid
  #   return(ret)} 
  # else {return(ret$posterior_means)}
}


neb_snp <- function(data, sds, miter = 100, 
                    sparsity = NULL, grid = 10*length(data), 
                    eps = 1e-5, ...) {
  
  input_list <- list(...)
  full <- ifelse(is.null(input_list$full), F, input_list$full)
  details <- ifelse(is.null(input_list$details), F, T)
  split_k <- ifelse(is.null(input_list$split_k), 3, input_list$split_k)
  M <- grid
  
  ifelse(is.null(sparsity), print("unknown sparsity"), print("known sparsity"))
  ifelse(full, print("full sample"), print("split sample"))
  ifelse(details, print("return details"), print("return posterior means"))
  print(paste0("Max iterations: ", miter, sep=""))
  print(paste0("Grid points: ", grid, sep=""))
  
  N <- length(data)
  
  # # Split sample
  set.seed(123)
  itest <- sample(1:N, N/split_k)
  itrain <- setdiff(1:N, itest)
  N_train <- length(itrain)
  N_test <- length(itest)
  
  # tau_grid
  a <- quantile(data, 0.01) - 2*sd(data)
  b <- quantile(data, 0.99) + 2*sd(data)
  R <- b-a
  M <- 10*N
  # lag <- abs(a)/((abs(a)/R*M))
  # tau_grid <- c(seq(a, 0, length.out = ceiling(abs(a)/R*M)),
  #               seq(lag, b, lag))
  tau_grid <- seq(a, b, length.out = M)
  tau_grid[which.min(abs(tau_grid))] <- 0 
  
  # known w
  if(!is.null(sparsity)) {
    w <- sparsity
  } else {
    # unknown w
    if(length(unique(sds))==1) {
      w <- max(w0_func(data, homosk=TRUE,sigsqy_homo=1),1/M)
    } else {
      w <- max(w0_func(data, homosk=F, sigsq_y = sds^2),1/M) 
    }
  }
  
  pi_hat <- rep((1-w)/M,M)
  # Start large to make sure sufficient separation
  l0 <- sqrt(N/(w+1))
  
  lik_mat <- dnorm(outer(tau_grid, data, "-"), sd = rep(sds, each = M))
  
  
  if(!full) { # split sample
    ret <- C_snp(lik_mat, itrain, itest, tau_grid, pi_hat, l0, w, miter = miter, eps = eps)
  } else {
    ret <- C_snp_full(lik_mat, tau_grid, pi_hat, l0, w, miter = miter, eps = eps)
  }
  
  ret$y <- data
  ret$posterior_mode <- ret$tau_grid[sapply(1:N, function(x) which.max(ret$posterior_hat[,x]))]
  class(ret) <- "neb_snp"
  ret
}

credint <- function(ret, mu=NULL, alpha=0.05, plot=T, plot_points=100) {
  N <- length(ret$posterior_means)
  
  CI_lo <- rep(NA,N)
  CI_hi <- rep(NA,N)
  # mu_true <- Th
  # x_axis <- 1:N
  prob_zero <- ret$posterior_hat[which(ret$tau_grid==0),]
  for(i in 1:N){
    # Get the CI for each point.  Rescale the posterior distribution point by point.
    # post_dist_i <- ret$posterior_hat[,i]/sum(ret$posterior_hat[,i])
    # Now get the first time the CDF goes above 2.5% and 97.5%.
    post_CDF_i <- cumsum(ret$posterior_hat[,i])
    
    idx_lo <- min(which(post_CDF_i >= alpha/2))
    idx_hi <- min(which(post_CDF_i > 1-alpha/2))
    if(idx_hi < idx_lo) idx_hi <- idx_lo 
    CI_lo[i] <- ret$tau_grid[idx_lo]
    CI_hi[i] <- ret$tau_grid[idx_hi]
    
    # prob_zero[i] <- post_dist_i[which(ret$tau_grid==0)]
  }
  
  spaced_idx <- round(seq(1,N,length.out=plot_points),0)
  
  if(!is.null(mu)) {covered <- (mu <= CI_hi) & (mu >= CI_lo)}
  if(plot) {
    if(is.null(mu)) {
      p <- data.frame(yhat=ret$posterior_means[spaced_idx],
                      y=y[spaced_idx],
                      ui=CI_hi[spaced_idx],
                      li=CI_lo[spaced_idx]) %>%
        arrange(yhat) %>%
        ggplot(aes(x=order(yhat), y=yhat)) + 
        geom_point(aes(shape = "Observed")) + 
        geom_point(aes(x=order(yhat), y=y, shape = "Posterior mean")) +
        geom_errorbar(aes(ymin=li, ymax=ui),
                      position=position_dodge(0.5)) +
        xlab("Observations") + ylab(expression(hat(y))) +
        ggtitle("Posterior means and 95% credible intervals")+
        scale_shape_manual("", values=c("Posterior mean"=4, 
                                        "Observed"=1)) +
        theme(legend.position="bottom")
    } else{
      p <- data.frame(yhat=ret$posterior_means[spaced_idx],
                      y=y[spaced_idx],
                      mu=mu[spaced_idx],
                      ui=CI_hi[spaced_idx],
                      li=CI_lo[spaced_idx],
                      covered = covered[spaced_idx]) %>%
        arrange(yhat) %>%
        ggplot(aes(x=order(yhat), y=yhat, col=covered, group=covered)) + 
        geom_point(aes(shape = "Posterior mean")) + 
        geom_point(aes(x=order(yhat), y=y, shape = "Observed")) +
        geom_point(aes(x=order(yhat), y=mu, shape = "mu"), col = "black") +
        geom_errorbar(aes(ymin=li, ymax=ui),
                      position=position_dodge(0.5))+
        xlab("Observations") + ylab(expression(hat(y))) +
        ggtitle("Posterior means and 95% credible intervals")+
        scale_shape_manual("", values=c("Posterior mean"=4, 
                                        "Observed"=1,
                                        "mu" = 2)) +
        theme(legend.position="bottom")
    }
    print(p)
  }
  
  ret <- list(CI_lwr = CI_lo, 
              CI_upr = CI_hi,
              prob_zero = prob_zero)
  
  if(plot) ret$plot <- p
  if(!is.null(mu)) {
    ret$coverage <- mean(covered)
    ret$coverage_simul <- all(covered)
  }
  
  return(ret)
}

fdr_neb <- function(ret, alpha = 0.05, plot=T) {
  N <- length(ret$y)
  if(class(ret) == "neb_snp") {
    # prob_zero <- ret$p_zero
    l_delta <- ret$delta[1]
    r_delta <- ret$delta[2]
    prob_zero <- colSums(ret$posterior_hat[l_delta:r_delta, ])
  } else {
    prob_zero <- ret$posterior_hat[which(ret$tau_grid==0),]
  }
  
  prob_zero_rank_idx <- order(prob_zero)
  prob_zero_ranked <- sort(prob_zero)
  k <- 1
  for(k in 1:N) {
    if(sum(prob_zero_ranked[1:k]/k) > alpha) {
      break
    } else {
      # print(sum(prob_zero_ranked[1:k]/k))
    }
  }
  maxK <- k-1
  gene_nz_idx <- prob_zero_rank_idx[1:maxK]
  
  if(plot) {
    plot(ret$y, cex = .5, pch = 4)
    points(gene_nz_idx, ret$y[gene_nz_idx], col="red", pch = 16, cex = .5)
  }
  
  return(gene_nz_idx)
}

fnr_neb <- function(ret, alpha = 0.05, plot=T) {
  N <- length(ret$y)
  if(class(ret) == "neb_snp") {
    # prob_zero <- ret$p_zero
    l_delta <- ret$delta[1]
    r_delta <- ret$delta[2]
    prob_zero <- colSums(ret$posterior_hat[l_delta:r_delta, ])
  } else {
    prob_zero <- ret$posterior_hat[which(ret$tau_grid==0),]
  }
  
  prob_zero_rank_idx <- order(prob_zero)
  prob_zero_ranked <- sort(prob_zero)
  k <- 1
  for(k in 1:(N-1)) {
    if( (1-sum(prob_zero_ranked[(k+1):N])/(N-k)) <= alpha ) {
      break
    } else {
      # print(sum(prob_zero_ranked[1:k]/k))
    }
  }
  maxK <- k-1
  gene_nz_idx <- prob_zero_rank_idx[1:maxK]
  
  if(plot) {
    plot(ret$y, cex = .5, pch = 4)
    points(gene_nz_idx, ret$y[gene_nz_idx], col="red", pch = 16, cex = .5)
  }
  
  return(gene_nz_idx)
}


fdr_bh <- function(y, sds, pval = NULL, pi0 = 1, alpha = 0.05, plot = F) {
  N <- length(y)
  if(is.null(pval)) pval <- 2*pnorm(abs(y/sds), lower.tail = F)
  
  pval_ranked <- sort(pval)
  pval_idx <- order(pval)
  
  for(kk in 1:N) {
    if(pval_ranked[kk] > kk/N*alpha/pi0) {break}
  }
  maxK2 <- kk-1
  gene_nz_idx <- pval_idx[1:maxK2]
  
  if(plot) {
    plot(y)
    points(gene_nz_idx, y[gene_nz_idx], col="red")
  }
  
  return(gene_nz_idx)
}

fdr_storey <- function(y, sds, pval = NULL, pi0 = NULL, alpha = 0.05) {
  library(qvalue)
  N <- length(y)
  if(is.null(pval)) pval <- 2*pnorm(abs(y/sds), lower.tail = F)
  qobj <- qvalue::qvalue(p = pval, pi0 = pi0, fdr.level = alpha)
  
  which(qobj$significant)
}

fdr <- function(method = "bh", y, sds, ret = NULL, pi0 = NULL, alpha = 0.05, plot = F) {
  if(method == "bh") idx <- fdr_bh(y = y, sds = sds, alpha = alpha)
  if(method == "storey") idx <- fdr_storey(y = y, sds = sds, alpha = alpha)
  
  if(method == "bh_plugin") idx <- fdr_bh(y = y, sds = sds, pi0 = ret$w, alpha = alpha)
  if(method == "storey_plugin") idx <- fdr_storey(y = y, sds = sds, pi0 = ret$w, alpha = alpha)
  
  if(method == "bh_w0") idx <- fdr_bh(y = y, sds = sds, pi0 = pi0, alpha = alpha)
  if(method == "storey_w0") idx <- fdr_storey(y = y, sds = sds, pi0 = pi0, alpha = alpha)
  
  if(method == "neb") idx <- fdr_neb(ret = ret, alpha = alpha)

  if(method == "hart") idx <- fdr_hart(y = y, sds = sds, fdr_level = alpha)
  
  if(method == "slope") idx <- fdr_slope(y = y, sds = sds, fdr_level = alpha)
    
  idx
} 


plot.neb_snp <- function(ret, mu = NULL, type = "posterior_mean", num.ci = NULL) {
  grid <- c(ret$tau_grid[2:length(ret$tau_grid)]) - c(ret$tau_grid[1:(length(ret$tau_grid)-1)])
  grid <- c(grid[1], grid)
  p1 <- data.frame(tau_grid = rep(ret$tau_grid, 3), 
                   y = 1/grid*c(ret$w* ddlaplace(ret$tau_grid, ret$l0)+(1-ret$w)*ret$pi_hat,
                                ret$w* ddlaplace(ret$tau_grid, ret$l0), 
                                (1-ret$w)*ret$pi_hat),
                   mixture = rep(c("SNP", "Spike", "NP"), each=length(ret$tau_grid))) %>%
    ggplot(aes(x = tau_grid, y = y, col=mixture)) +
    geom_line(aes(linetype=mixture)) +
    ylab(expression(hat(pi)(mu))) + xlab("mu") + ggtitle("Prior estimation")+
    scale_linetype_manual(values=c("solid", "dashed", "dashed"))+
    scale_color_manual(values=c("blue", "cyan3", "red"))+
    theme(legend.position="bottom")
  
  p1 <- data.frame(tau_grid = rep(ret$tau_grid, 2), 
                   y = c(ret$w* ddlaplace(ret$tau_grid, ret$l0), 
                         (1-ret$w)*ret$pi_hat),
                   mixture = rep(c("Spike", "NP"), each=length(ret$tau_grid))) %>%
    ggplot(aes(x = tau_grid, y = y, col=mixture)) +
    geom_line() +
    ylab(expression(hat(pi)(mu))) + xlab("mu") + ggtitle("Prior estimation") +
    # scale_linetype_manual(values=c("solid", "dashed", "dashed"))+
    # scale_color_manual(values=c("blue", "cyan3", "red"))+
    theme(legend.position="bottom")
  
  if(is.null(num.ci)) {
    idx <- 1:length(ret$posterior_means)
  } else {
    idx <- floor(seq(1, length(ret$posterior_means), length.out = num.ci))
  }
  if(type=="MAP") {
    yhat <- ret$posterior_mode
    title_tmp <- "Posterior mode"
  } else {
    yhat <- ret$posterior_means
    title_tmp <- "Posterior mean"
  }
  if(!is.null(mu)) {
    # ci <- credint(ret, mu)
    # covered <- (mu <= ci$CI_upr) & (mu >= ci$CI_lwr)
    # p2 <- data.frame(yhat =yhat[idx], 
    #                  mu=mu[idx],
    #                  lwr = ci$CI_lwr[idx], 
    #                  upper = ci$CI_upr[idx], 
    #                  covered = covered[idx]) %>%
    #   ggplot(aes(x = 1:length(mu), y = yhat, col = covered)) +
    #   geom_point() + 
    #   geom_point(aes(x = 1:length(mu), y = mu), shape = 4) +
    #   # geom_point(aes(x = 1:N, y = posterior_mode), shape = 4, col="green") + 
    #   ylab(expression(hat(y))) + xlab("Observation") +
    #   geom_errorbar(aes(ymin=lwr, ymax=upper), alpha=.2) + ggtitle(title_tmp)+
    #   theme(legend.position="bottom")
    p2 <- credint(ret,mu)$plot
  } else {
    ci <- credint(ret)
    p2 <- data.frame(yhat = yhat[idx], 
                     lwr = ci$CI_lwr[idx], 
                     upper = ci$CI_upr[idx]) %>%
      ggplot(aes(x = 1:length(yhat), y = yhat)) +
      geom_point(col = "#00BFC4") + 
      ylab(expression(hat(y))) + xlab("Observation") +
      geom_errorbar(aes(ymin=lwr, ymax=upper), alpha=.2, col = "#00BFC4") +
      ggtitle(title_tmp) +
      theme(legend.position="bottom")
  }
  
  grid.arrange(p1, p2, ncol=2)
}

plot.neb_dnp <- function(ret, mu = NULL, type = "posterior_mean") {
  p1 <- data.frame(tau_grid = rep(ret$tau_grid, 2), 
                   y = ret$pi_hat) %>%
    ggplot(aes(x = tau_grid, y = y)) +
    geom_line() +
    ylab(expression(hat(pi)(mu))) + xlab("mu") + ggtitle("Prior estimation")+
    theme(legend.position="bottom")
  
  if(type=="MAP") {
    yhat <- ret$posterior_mode
    title_tmp <- "Posterior mode"
  } else {
    yhat <- ret$posterior_means
    title_tmp <- "Posterior mean"
  }
  if(!is.null(mu)) {
    ci <- credint(ret, mu)
    covered <- (mu <= ci$CI_upr) & (mu >= ci$CI_lwr)
    p2 <- data.frame(yhat =yhat, 
                     mu=mu,
                     lwr = ci$CI_lwr, 
                     upper = ci$CI_upr, 
                     covered = covered) %>%
      ggplot(aes(x = 1:length(mu), y = yhat, col = covered)) +
      geom_point() + 
      geom_point(aes(x = 1:length(mu), y = mu), shape = 4) +
      # geom_point(aes(x = 1:N, y = posterior_mode), shape = 4, col="green") + 
      ylab(expression(hat(y))) + xlab("Observation") +
      geom_errorbar(aes(ymin=lwr, ymax=upper), alpha=.2) + ggtitle(title_tmp)+
      theme(legend.position="bottom")
  } else {
    ci <- credint(ret)
    p2 <- 
      data.frame(yhat = yhat, 
                 lwr = ci$CI_lwr, 
                 upper = ci$CI_upr) %>%
      ggplot(aes(x = 1:length(yhat), y = yhat)) +
      geom_point(col = "#00BFC4") + 
      ylab(expression(hat(y))) + xlab("Observation") +
      geom_errorbar(aes(ymin=lwr, ymax=upper), alpha=.2, col = "#00BFC4") +
      ggtitle(title_tmp)+
      theme(legend.position="bottom")
  }
  
  grid.arrange(p1, p2, ncol=2)
}


ddlaplace <- function(y, s) {
  # s/2*exp(-s*abs(y))
  tmp <- exp(log(s/2) - s*abs(y))
  tmp <- tmp/sum(tmp)
  tmp
}

render_report <- function(job_id) {
  rmarkdown::render(
    "simulation.Rmd", params = list(
      job_id = job_id
    ),
    output_file = paste0("../write_up/simulation/simulation-", job_id, ".pdf")
  )
}

# read sim seteup
read_setup <- function(filepath = "../output/simulation/summary/", params) {
  rds_list <- lapply(params$job_id, 
                     function(job) readRDS(paste0(filepath, 
                                                  job, 
                                                  ".RDS")))
  
  # get the union of multiple sim setup
  unionFun <- function(n, obj) {
    unique(unlist(lapply(obj, `[[`, n)))
  }
  
  unions <- lapply(seq_along(rds_list[[1]]), FUN = unionFun, obj = rds_list)
  names(unions) <- names(rds_list[[1]])
  
  unions
  
}

# read multiple sim output
read_output <- function(filepath = "../output/simulation/summary/", params, type = "aggregate") 
{
  filepaths <- unlist(lapply(params$job_id,
                             function(job_id) Sys.glob(paste0(filepath, 
                                                              job_id, 
                                                              "_", type,
                                                              "*.csv"))))
  ret_list <- lapply(filepaths, fread)
  
  rbindlist(ret_list, fill = T, use.names = T)
}


plot_grid <- function(ret, criteria_of_int, methods,
                      x_var, vpanel = NULL, hpanel, 
                      ylabs = criteria_of_int, 
                      ylims = NULL,
                      param_group = c("method", "w0", "V"),
                      title = NULL, FUN = mean,
                      save = F, plotname = paste(criteria_of_int, collapse = "_"),
                      width = 9, height = 12) 
{
  
  id_vars <- param_group[!grepl("method", param_group)]
  dcast_formula <- as.formula(paste0(paste(id_vars, collapse = "+"),
                                     "~ method"))
  
  ret <- ret[method %in% methods]

  labs_cols <- dplyr::inner_join(method_labels_tbl(), 
        data.table(method = unique(ret$method)),
        by = "method")
  
  ret$method <- factor(ret$method, levels = labs_cols$method, ordered = T)
  
  p <- list()
  
  for(i in 1:length(criteria_of_int)) {
    
    crit <- criteria_of_int[i]
    
    if(grepl("ratio|diff", crit)) {
      
      crit_ <- str_replace(crit, "_ratio|_diff", "")
      if(crit_ == "w") {
        ret_tmp <- ret[,.(avg = FUN(get(crit_) - w0)), 
                       by=param_group]
      } else {
        ret_tmp <- ret[,.(avg = FUN(get(crit_))), 
                       by=param_group]
      }
      
      ret_tmp <- dcast(ret_tmp, dcast_formula, value.var = "avg")
      
      if(crit_ == "w") {
        ret_tmp[, (methods) := lapply(.SD, function(x) x), 
                .SDcols = methods]
      } else {
        ret_tmp[, (methods) := lapply(.SD, function(x) (x)/ret_snp), 
                .SDcols = methods]
      }
      
      ret_tmp <- melt(ret_tmp,
                      id.vars = id_vars,
                      variable.name = "method",
                      value.name = "avg")
      
    } else {
      ret_tmp <- ret[,.(avg = FUN(get(crit))), by = param_group]
    }
    
    tmp_p <- ret_tmp %>% 
      # mutate(upper = avg+sd, lower = avg-sd) %>%
      ggplot(aes(x = get(x_var), y = avg, col = method, group = method)) +
      geom_point(aes(shape = method), size = 3) + 
      geom_line(data = ret_tmp[!is.na(avg)]) +
      ylab(TeX(ylabs[i])) +
      xlab(x_var) +
      # ggtitle(paste0("V = ", opt$V, "; N = ", opt$N)) +
      scale_x_continuous(breaks = unique(ret_tmp[[x_var]])) +
      scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
      scale_shape_manual(values = labs_cols$shapes, labels = labs_cols$labels) +
      my_theme + guides(colour = guide_legend(nrow = 1))

    if(!is.null(ylims)) tmp_p <- tmp_p + coord_cartesian(ylim = ylims)
    
    if(crit == "w") tmp_p  <- tmp_p +  geom_abline(slope = 1, intercept = 0)
    if(grepl("diff", crit)) tmp_p  <- tmp_p +  geom_hline(yintercept = 0) 
    if(crit %in% c("fdr", "fnr")) tmp_p  <- tmp_p + geom_hline(yintercept = 0.05)
    if(crit == "mode_zero") tmp_p  <- tmp_p +  geom_abline(slope = ret_sum$N, intercept = 0) 
    if(crit == "mode_nonzero") tmp_p  <- tmp_p +  geom_abline(slope = -ret_sum$N, intercept = ret_sum$N) 
    if(crit == "ci_coverage") tmp_p  <- tmp_p + geom_hline(yintercept = 0.95)
    
    
    if(!is.null(hpanel)) { 
      tmp_p <- tmp_p + 
        facet_wrap(~get(hpanel), ncol = 1)
    }
    if(!is.null(vpanel)) {
      tmp_p <- tmp_p + 
        facet_wrap(get(vpanel)~., nrow = 1)
    }
    if(!is.null(hpanel) & !is.null(vpanel)) { 
      tmp_p <- tmp_p + 
        facet_grid(get(hpanel) ~ get(vpanel))
      tmp_p <- add_general_label(tmp_p, hpanel, vpanel)
    }
    
    # p[[i]] <- tmp_p + theme_bw()
    p[[i]] <- tmp_p
  }
  
  if(i == 1) {
    figure <- p[[1]]
  } else if(!is.null(hpanel) & !is.null(vpanel)) { 
    figure <- ggpubr::ggarrange(plotlist = p, ncol = 1)
  } else {
    figure <- ggpubr::ggarrange(plotlist = p, common.legend = T, 
                                legend = "bottom", ncol = 1)
  }
  
  print(figure)
  
  if(save) {
    ggexport(figure, 
             filename = paste0("../output/simulation/plot/",
                               # paste(ret_sum$job_id, collapse = "-"), 
                               ret_sum$job_id[1],
                               "_", plotname, ".pdf"),
             width = width, height = height)
  }
  
  figure
}


plot_fdr <- function(fdr_mat, criteria_of_int, 
                     methods, method_labels = NULL,
                     x_var = NULL, vpanel = NULL, hpanel, 
                     ylabs = criteria_of_int,
                     ylims = NULL,
                     param_group = c("method", "w0", "V", "fdr_level"),
                     title = NULL, FUN = mean,
                     save = F, plotname = paste(criteria_of_int, collapse = "_"),
                     width = 9, height = 12) 
{
  
  crit <- criteria_of_int
  
  if(is.null(x_var)) {
    x_var <- "fdr_level"
    if(grepl("fnr", crit)) x_var <- "fnr_level"
  }
  
  fdr_mat <- fdr_mat[method %in% methods]
  
  labs_cols <- dplyr::inner_join(method_labels_tbl_fdr(), 
                                 data.table(method = unique(fdr_mat$method)),
                                 by = "method")
  
  fdr_mat$method <- factor(fdr_mat$method, levels = labs_cols$method, ordered = T)
  
  fdr_mat <- fdr_mat[, .(avg = FUN(get(crit))), by=param_group]
  
  tmp_p <- ggplot(fdr_mat[method %in% methods],
                  aes(x = get(x_var), 
                      y = avg, 
                      col = method, 
                      shape = method,
                      group = method)) +
    geom_point(size = 3) +
    geom_line() +
    ylab(ylabs) +
    xlab(x_var) +
    # ggtitle(paste0("V = ", opt$V, "; N = ", opt$N)) +
    scale_x_continuous(breaks = unique(fdr_mat[[x_var]])) +
    scale_color_manual(values = labs_cols$cols, labels = labs_cols$labels) +
    scale_shape_manual(values = labs_cols$shapes, labels = labs_cols$labels) +
    # facet_grid(get(hpanel) ~ get(vpanel), scale = "free") +
    my_theme 
  
  if(grepl("diff", crit)) {
    tmp_p <- tmp_p + geom_hline(yintercept = 0) 
  } else if(grepl("fdr", crit)) {
    tmp_p <- tmp_p + geom_abline(slope = 1, intercept = 0)
  }
  
  if(!is.null(ylims)) tmp_p <- tmp_p + coord_cartesian(ylim = ylims)
  
  if(!is.null(hpanel)) { 
    tmp_p <- tmp_p + 
      facet_wrap(~get(hpanel), ncol = 1, scale = "free")
  }
  if(!is.null(vpanel)) {
    tmp_p <- tmp_p + 
      facet_wrap(get(vpanel)~., nrow = 1, scale = "free")
  }
  if(!is.null(hpanel) & !is.null(vpanel)) { 
    tmp_p <- tmp_p + 
      facet_grid(get(hpanel) ~ get(vpanel), scale = "free")
    tmp_p <- add_general_label(tmp_p, hpanel, vpanel)
  }
  
  if(save) {
    ggsave(plot = tmp_p, 
           filename = paste0("../output/simulation/plot/",
                             ret_sum$job_id[1], "_", plotname, ".pdf"),
           width = width, height = height)
  }
  
  grid.newpage()
  grid.draw(tmp_p)
}


add_general_label <- function(p, labelR, labelT) {
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 15, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 15, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  # Draw it
  # grid.newpage()
  # grid.draw(z)
  
  z
}

method_labels_tbl <- function() {
  cols <- hue_pal()(15)[c(15, 1, 6, 10, 8,
                          4, 9, 13, 2, 3,
                          12, 11, 1, 14 ,5, 5)]
  
  shapes <- c(7, 2, 1, 16, 20,
              17, 25, 0, 
              6, 17,
              3, 3, 0, 0, 4, 4)
  
  methods <- c("ret_pi0", "retf", "ret_dnp", "ret_snp", "ret_snp_split",
               "grouplinear", "XGB.M", "XGB.SG", "sslasso", "sslasso_oracle", 
               "slope", "slope_oracle", "gmleb", "ebthresh", "hart", "ds")
  ret_tbl <- data.table(method = methods,
             labels = c("Oracle", "JZ", "DNP", "SNP", "SNP split",
                        "Grouplinear", "XGB.M", "XGB.SB", 
                        "SSLASSO", "SSLASSO oracle",
                        "SLOPE", "SLOPE oracle",
                        "GMLEB", "EBayesThresh", "HART", "NEST"
             ),
             shapes = shapes,
             cols = cols)
  
  ret_tbl$method <- factor(ret_tbl$method, levels = methods, ordered = T)
  
  ret_tbl
}


method_labels_tbl_fdr <- function() {
  cols <- hue_pal()(15)[c(15, 1, 6, 10, 8,
                          5, 9, 13, 
                          2, 7,
                          15, 11, 3, 14 ,12)]
  shapes <- c(7, 2, 1, 16, 20,
              17, 25, 0, 
              6, 17,
              3, 3, 0, 0, 4)
  
  data.table(method = c("ret_pi0", "bh", "ret_dnp_neb", "ret_snp_neb", "ret_snp_split",
                        "ret_snp_bh_plugin", "ret_snp_storey_plugin", "XGB.SG", 
                        "storey", "bh_w0", 
                        "slope", "slope_oracle", "gmleb", "ebthresh", "hart"),
             labels = c("Oracle", "B&H", "DNP NEB", "SNP NEB", "SNP split",
                        "SNP-BH", "SNP-Storey", "XGB.SB", 
                        "Storey", "SSLASSO oracle",
                        "SLOPE", "SLOPE oracle",
                        "GMLEB", "EBayesThresh", "HART"
             ),
             cols = cols,
             shapes = shapes)
}


w0_func <- function(y, homosk=TRUE,sigsqy_homo=1, sigsq_y="default"){
  # you give me the data and whether the data is homoskedastic or not, and if heterosk., you give me
  #   the vector of variances (otherwise you give me the one variance), and I'll give you back my best
  #   guess of the 0's.
  
  # See Cun Hui Zhang's 'Fourier Methods'
  N = length(y)
  c_n = (log(N)/2)^(1/2)
  
  if(homosk==FALSE){                       # this allows us to default to homoskedastic when it's homosk.
    if(sum(sigsq_y==mean(sigsq_y))==length(sigsq_y)){
      homosk=TRUE
      sigsqy_homo=sigsq_y[1]
    }
  }
  
  if(homosk==TRUE){
    sigsq_y = rep(sigsqy_homo,N)
    erf1 = erfi((c_n*sigsq_y-1i*y)/sqrt(2*sigsq_y))
    erf2 = erfi((c_n*sigsq_y+1i*y)/sqrt(2*sigsq_y))
    allterms = 1/2/c_n*sqrt(pi/2)*exp(y^2/2/sigsq_y)*(erf1+erf2)
    w0 = mean(allterms)
  }
  if(homosk==FALSE){
    c_n = sqrt(1/2*log(N))
    f = function(t) cos(t*y_j)/sum(exp(-rep(t,N)^2*sigsq_y/2))
    K_j = rep(NA,N)
    for(j in 1:N){
      y_j = y[j]
      int_val = integrate(f,lower=-c_n, upper=c_n)$value
      K_j[j] = 1/2*int_val
    }
    multiplicative_factor = 7.7*(9.9/log(N))^(2/3)
    w0 = sum(K_j)*multiplicative_factor
  }
  return(min(max(Re(w0),0),.95))
}
