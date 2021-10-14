# for finalSummary.pdf
# start.time <- Sys.time()
#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman")
if(!require(pacman)){install.packages("pacman")}
pacman::p_load(data.table, Rcpp, ggplot2, optparse, invgamma)
# pacman::p_load(REBayes, Rmosek, rmutil)
# source("src/npeb.r")
# source("src/utilities.R")
# sourceCpp("src/neb.cpp")
source("npeb.r")
source("utilities.R")
sourceCpp("neb.cpp")

source("other_methods.R")

option_list = list(
  make_option(c("-o", "--out"), type="character", default="mboot", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-N", "--N"), type="integer", default="500", 
              help="nrow [default= %default]", metavar="number"),
  make_option(c("--seed"), type="integer", default="1", 
              help="random seed [default= %default]", metavar="number"),
  make_option(c("-x", "--setting"), type="character", default="a", 
              help="y generator options", metavar="character"),
  make_option(c("-w", "--w"), type="double", default="0.95", 
              help="sparsity", metavar="number"),
  make_option(c("-V", "--V"), type="double", default="3", 
              help="sinal strength", metavar="number"),
  make_option(c("--homo"), type="integer", default="0", 
              help="homoscedasticity", metavar="character"),
  make_option(c("--miter"), type="integer", default="500", 
              help="Maximum iteration", metavar="character"),
  make_option(c("--jobid"), type="character", default="1", 
              help="jobid", metavar="number"),
  make_option(c("--umax"), type="double", default="2", 
              help="hetero Uniform max", metavar="number"),
  make_option(c("--methods"), type="character", 
              default="ret_pi0-retf-ret_dnp-ret_snp", 
              help="methods", metavar="number"),
  make_option(c("--fdr"), type="character", 
              default="bh-bh_w0-storey-storey_w0-slope-hart", 
              help="methods", metavar="number")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

# methods <- c("ret_pi0", "retf", "ret_dnp", "ret_snp", "ret_snp_split")
methods <- c("grouplinear", "XGB.M", "XGB.SG", 
            "sslasso", "sslasso_oracle", "slope", "slope_oracle", 
            "gmleb", "ebthresh", "hart")
methods <- strsplit(opt$methods, "-")[[1]]

fdr_methods <- strsplit(opt$fdr, "-")[[1]]
fdr_methods_neb <- c("neb", "bh_plugin", "storey_plugin")

fdr_levels <- c(0.01, seq(0.05, 0.3, 0.05))


data <- Generate(opt)
y <- data$y
sds <- data$sds
mu <- data$mu
plot(y)

system.time({
  if("ret_pi0" %in% methods) {
    ret_pi0 <- neb_pi0(y, sds, opt)
  }
  if("retf" %in% methods) {
    retf <- neb_0(y, sds, sparse=T, full=T, eps = 1e-5, miter = opt$miter)
  }
  if("ret_dnp" %in% methods) {
    ret_dnp <- neb_0(y, sds, sparse=T, eps = 1e-5, split_k = 3, miter = opt$miter)
  }
  if("ret_snp" %in% methods) {
    ret_snp <- neb_snp(y, sds, sparse=T, full=T, eps = 1e-5, miter = opt$miter)
  }
  if("ret_snp_split" %in% methods) {
    ret_snp_split <- neb_snp(y, sds, sparse=T, full=F, split_k = 3, eps = 1e-5, miter = opt$miter)
  }
  if("grouplinear" %in% methods) {
    grouplinear <- est.grouplinear(y, sds)
  }
  if("XGB.M" %in% methods) {
    XGB.M <- thetahat.M(y, sds)
  }
  if("XGB.SG" %in% methods) {
    XGB.SG <- thetahat.SG(y, sds)
  }
  if("sslasso" %in% methods) {
    sslasso <- est.sslasso(y, sds)
  }
  if("sslasso_oracle" %in% methods) {
    sslasso_oracle <- est.sslasso(y, sds, w = opt$w)
  }
  if("slope" %in% methods) {
    slope <- est.slope(y, sds, fdr_level = .1)
  }
  if("slope_oracle" %in% methods) {
    slope_oracle <- est.slope(y, sds, mu = mu, fdr_level = .1)
  }
  if("gmleb" %in% methods) {
    gmleb <- est.gmleb(y, sds)
  }
  if("ebthresh" %in% methods) {
    ebthresh <- est.ebayesthresh(y, sds)
  }
  if("hart" %in% methods) {
    hart <- est.hart(y, sds, fdr_level = .05)
  }
  if("ds" %in% methods) {
    ds <- oracle.ds(y, sds, mu)
  }
})

ret_mat <- NULL
fdr_mat <- NULL
fnr_mat <- NULL

for(fdr_level in fdr_levels) {
  emp_fdr <- sapply(fdr_methods, function(fdr_m) {
    fdr_idx <- fdr(method = fdr_m, y, sds, pi0 = opt$w, alpha = fdr_level, plot = F) 
    c(sum(fdr_idx > (1-opt$w)*opt$N)/length(fdr_idx), 
      length(fdr_idx),
      sum(fdr_idx < (1-opt$w)*opt$N) / ((1-opt$w)*opt$N)
      )
  } )
  
  fdr_mat <- rbind(fdr_mat, 
                   data.table(method = fdr_methods,
                              seed = opt$seed,
                              n = opt$N,
                              V = opt$V,
                              w0 = opt$w,
                              u = opt$umax,
                              homo = opt$homo,
                              fdr_level = fdr_level,
                              emp_fdr = emp_fdr[1,],
                              num_dis = emp_fdr[2,],
                              emp_power = emp_fdr[3,]))
}

for(ret in methods) {
  tmp <- get(ret)
  
  if(grepl("ret", ret)) {
    # NEB methods
    
    # credible interval
    ci <- credint(tmp, mu = mu, plot=F)
    
    # fdr
    fdr_idx <- fdr_neb(tmp, plot=F)
    posterior_means_fdr <- tmp$posterior_means
    posterior_means_fdr[setdiff(1:opt$N, fdr_idx)] <- 0
    
    # fnr
    fnr_idx <- fnr_neb(tmp, plot=F)
    posterior_means_fnr <- tmp$posterior_means
    fnr_nondiscover_idx <- setdiff(1:opt$N, fnr_idx)
    posterior_means_fnr[fnr_nondiscover_idx] <- 0
    
    # fdr levels
    for(fdr_level in fdr_levels) {
      emp_fdr <- sapply(fdr_methods_neb, function(fdr_m) {
        fdr_idx <- fdr(method = fdr_m, y, sds, ret = tmp, pi0 = opt$w, alpha = fdr_level, plot = F) 
        c(sum(fdr_idx > (1-opt$w)*opt$N)/length(fdr_idx), 
          length(fdr_idx),
          sum(fdr_idx < (1-opt$w)*opt$N) / ((1-opt$w)*opt$N))
      } )
      
      fdr_mat <- rbind(fdr_mat, 
                       data.table(method = paste0(ret, "_", fdr_methods_neb),
                                  seed = opt$seed,
                                  n = opt$N,
                                  V = opt$V,
                                  w0 = opt$w,
                                  u = opt$umax,
                                  homo = opt$homo,
                                  fdr_level = fdr_level,
                                  emp_fdr = emp_fdr[1,],
                                  num_dis = emp_fdr[2,],
                                  emp_power = emp_fdr[3,]))
    }
    
    # fnr levels
    for(fnr_level in fdr_levels) {
      fnr_idx <- fnr_neb(tmp, alpha = fnr_level, plot=F)
      fnr_nondiscover_idx <- setdiff(1:opt$N, fnr_idx)
      emp_fnr <- sum(fnr_nondiscover_idx <= (1-opt$w)*opt$N)/length(fnr_nondiscover_idx)
      num_dis <- length(fnr_idx)
      
      fnr_mat <- rbind(fnr_mat, 
                       data.table(method = ret,
                                  seed = opt$seed,
                                  n = opt$N,
                                  V = opt$V,
                                  w0 = opt$w,
                                  u = opt$umax,
                                  homo = opt$homo,
                                  fnr_level = fnr_level,
                                  emp_fnr = emp_fnr,
                                  num_dis = num_dis))
    }
  } else {
    posterior_means_fdr <- posterior_means_fnr <- fdr_idx <- fnr_nondiscover_idx <- NA
    ci <- NULL
  }
  
  
  ret_mat <- rbind(ret_mat,
                   data.table(method = ret,
                              seed = opt$seed,
                              n = opt$N,
                              V = opt$V,
                              w0 = opt$w,
                              u = opt$umax,
                              homo = opt$homo,
                              iter = tmp$iter,
                              MSE = mean((tmp$posterior_means - mu)^2),
                              MSE_mode = mean(abs(tmp$posterior_mode - mu)),
                              MSE_mode_square = mean((tmp$posterior_mode - mu)^2),
                              mode_zero = sum(tmp$posterior_mode == 0),
                              MSE_fdr = mean((posterior_means_fdr - mu)^2),
                              MSE_fnr = mean((posterior_means_fnr - mu)^2),
                              w = tmp$w,
                              ci_len = mean(ci$CI_upr - ci$CI_lwr),
                              ci_coverage = ci$coverage,
                              ci_simul_coverage = ci$coverage_simul,
                              fdr = sum(fdr_idx > (1-opt$w)*opt$N)/length(fdr_idx),
                              fnr = sum(fnr_nondiscover_idx <= (1-opt$w)*opt$N)/length(fnr_nondiscover_idx))
  )
  
}

ret_mat
fdr_mat
fnr_mat


new_dir <- paste0("../hpcc_output/", opt$jobid, "/")
dir.create(new_dir)

optname <- paste(paste0("setting", opt$setting),
                 paste0("n", opt$N),
                 paste0("s", opt$seed),
                 paste0("w", opt$w),
                 paste0("v", opt$V),
                 paste0("u", opt$umax),
                 paste0("h", opt$homo),
                 paste0("m", opt$miter), sep = "_")

filename <- paste0(new_dir, optname, ".csv")
fwrite(ret_mat, filename)

filename <- paste0(new_dir, "fdr_", optname, ".csv")
if(!is.null(fdr_mat)) fwrite(fdr_mat, filename)

filename <- paste0(new_dir, "fnr_", optname, ".csv")
if(!is.null(fnr_mat)) fwrite(fnr_mat, filename)
