#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// unscaled posterior distribution
// [[Rcpp::export]]
arma::mat est_posterior_hat(
	arma::mat likelihood_mat, 
	arma::colvec pi_hat) 
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	// long n = likelihood_mat.ncol(), m = likelihood_mat.nrow();
	arma::mat posterior_hat(m, n);

	for(long i=0; i<n; i++){
		for(long j=0; j<m; j++) {
			posterior_hat(j,i) = likelihood_mat(j,i) * pi_hat(j);
		}
	}
	return posterior_hat;
}

// [[Rcpp::export]]
arma::colvec est_posterior_means(
	arma::mat posterior_hat,
	arma::colvec tau_grid) 
{
	long n = posterior_hat.n_cols, m = posterior_hat.n_rows;
	// long n = posterior_hat.ncol(), m = posterior_hat.nrow();
	arma::colvec posterior_means(n, arma::fill::zeros);
	arma::colvec norm(n, arma::fill::zeros);

	// normalization 
	for(long i=0; i<n; i++){
		for(long j=0; j<m; j++){
			norm(i) += posterior_hat(j,i);
		}
	}

	// posterior means
	for(long i=0; i<n; i++){
		for(long j=0; j<m; j++) {
			// tmp(j,i) = posterior_hat(j,i) * tau_grid(j) / norm(i);
			posterior_means(i) += posterior_hat(j,i)*tau_grid(j)/norm(i);
		}
	}	

	return posterior_means;
}

// [[Rcpp::export]]
arma::colvec update_pi(arma::mat posterior_hat)
{
	long n = posterior_hat.n_cols, m = posterior_hat.n_rows;
	arma::colvec norm(n, arma::fill::zeros), pi_hat(m, arma::fill::zeros);

	for(long i=0; i<n; i++){
		for(long j=0; j<m; j++){
			norm(i) += posterior_hat(j,i);
		}
	}

	for(long j=0; j<m; j++){
		for(long i=0; i<n; i++) {
			pi_hat(j) += posterior_hat(j,i)/norm(i);
		}
		pi_hat(j) /= n;
	}

	double pi_sum = arma::sum(pi_hat);
	for(long j=0; j<m; j++) {
		pi_hat(j) /= pi_sum;
	}

	return pi_hat;
}


// long update_pi(arma::mat posterior_hat, arma::colvec pi_hat)
// {
// 	long n = posterior_hat.n_cols, m = posterior_hat.n_rows;
// 	arma::colvec norm(n, arma::fill::zeros);

// 	for(long i=0; i<n; i++){
// 		for(long j=0; j<m; j++){
// 			norm(i) += posterior_hat(j,i);
// 		}
// 	}

// 	for(long j=0; j<m; j++){
// 		pi_hat(j) = 0;
// 		for(long i=0; i<n; i++) {
// 			pi_hat(j) += posterior_hat(j,i)/norm(i);
// 		}
// 		pi_hat(j) /= n;
// 	}

// 	double pi_sum = arma::sum(pi_hat);
// 	for(long j=0; j<m; j++) {
// 		pi_hat(j) /= pi_sum;
// 	}

// 	return 0;
// }


// [[Rcpp::export]]
double loglik(
	arma::mat likelihood_mat, 
	arma::colvec pi_hat) 
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	double loglik;
	arma::mat posterior_hat(m,n);
	arma::vec norm(n, arma::fill::zeros);

	posterior_hat = est_posterior_hat(likelihood_mat, pi_hat);
	for(long i=0; i<n; i++) {
		for(long j=0; j<m; j++) {
			norm(i) += posterior_hat(j,i);
		}
		norm(i) = log(norm(i));
	}
	loglik = arma::sum(norm);
	
	return loglik;	
}

// [[Rcpp::export]]
List C_neb (
	arma::mat likelihood_mat,
	arma::uvec itrain, arma::uvec itest,
	arma::colvec tau_grid, arma::colvec pi_hat,
	long miter, double esp)
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	long n_train = itrain.n_elem, n_test = itest.n_elem;
	// long n_train = likelihood_mat_train.n_cols, n_test = likelihood_mat_test.n_cols;
	arma::colvec loglik_train(miter+1), loglik_test(miter+1), posterior_means(n), pi_hat_old(m); 
	arma::mat posterior_hat_train(m,n);
	long iter = 0;
	double loglikelihood;
	int stop_rule = 0;

	itrain = itrain - 1;
	itest = itest - 1;

	loglik_train(iter) = loglik(likelihood_mat.cols(itrain), pi_hat)/n_train;
	loglik_test(iter) = loglik(likelihood_mat.cols(itest), pi_hat)/n_test;

	while(iter < miter) {
		posterior_hat_train = est_posterior_hat(likelihood_mat.cols(itrain), pi_hat);
		pi_hat_old = pi_hat;
		pi_hat = update_pi(posterior_hat_train);

		loglik_train(iter+1) = loglik(likelihood_mat.cols(itrain), pi_hat)/n_train;
		loglik_test(iter+1) = loglik(likelihood_mat.cols(itest), pi_hat)/n_test;

		// std::cout << loglik_train(iter) << std::endl;
		// std::cout << loglik_test(iter) << std::endl;

		// stop_rule = (loglik_test(iter+1) <= (loglik_test(iter) + esp));
		stop_rule = (std::abs((loglik_test(iter+1) - loglik_test(iter)) / loglik_test(iter) )) < esp;

		if( ((iter+1)==miter) | stop_rule ) {
			posterior_hat_train = est_posterior_hat(likelihood_mat, pi_hat_old);
			posterior_means = est_posterior_means(posterior_hat_train, tau_grid);
			loglikelihood = loglik(likelihood_mat, pi_hat);
			break;
		}
		iter++;
  	// Rcout << iter << "\t";
	}
	Rcout << iter << "\n";

	return List::create(
		_["pi_hat"] 		 = pi_hat,
		_["posterior_means"] = posterior_means,
		_["posterior_hat"]   = posterior_hat_train,
		_["loglik_train"]    = loglik_train,
		_["loglik_test"]     = loglik_test,
		_["loglik"]          = loglikelihood,
		_["iter"]			= iter);
}

// [[Rcpp::export]]
List C_neb_full (
	arma::mat likelihood_mat,
	arma::colvec tau_grid, arma::colvec pi_hat,
	long miter, double esp)
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	// long n_train = likelihood_mat_train.n_cols, n_test = likelihood_mat_test.n_cols;
	arma::colvec loglik_train(miter+1), posterior_means(n), pi_hat_old(m); 
	arma::mat posterior_hat_train(m,n);
	long iter = 0;
	double loglikelihood, qn = std::min(1.0, exp(1)*sqrt(2*M_PI)/pow(n,2));
	int stop_rule = 0;

	loglik_train(iter) = loglik(likelihood_mat, pi_hat)/n;


	while(iter < miter) {
		posterior_hat_train = est_posterior_hat(likelihood_mat, pi_hat);
		pi_hat_old = pi_hat;
		pi_hat = update_pi(posterior_hat_train);

		loglik_train(iter+1) = loglik(likelihood_mat, pi_hat)/n;

		stop_rule = max(log(pi_hat % (1/pi_hat_old))) <= -log(exp(1)*qn)/n;
		if( ((iter+1)==miter) | stop_rule) {
			posterior_means = est_posterior_means(posterior_hat_train, tau_grid);
		  loglikelihood = loglik(likelihood_mat, pi_hat);
			break;
		}
		iter++;
	  	// Rcout << iter << "\t";
	}
	Rcout << iter << "\n";

	return List::create(
		_["pi_hat"] 		 = pi_hat,
		_["posterior_means"] = posterior_means,
		_["posterior_hat"]   = posterior_hat_train,
		_["loglik_train"]    = loglik_train,
		_["loglik"]          = loglikelihood,
		_["iter"]			= iter);
}

arma::colvec dlaplace(arma::colvec y, double s) {
	long n = y.n_elem;
	arma::colvec tmp(n, arma::fill::zeros);

	for(long i = 0; i < n; i++) {
		tmp(i) = exp(log(s/2) - s*std::abs(y(i)));
	}

	tmp = tmp/sum(tmp);
	return tmp;
}

// [[Rcpp::export]]
arma::colvec sample_C(arma::colvec& y, int n) {
	arma::colvec ret = Rcpp::RcppArmadillo::sample(y,floor(n),false);
	return ret;
}

// [[Rcpp::export]]
arma::colvec find_delta_0(
	arma::colvec tau_grid, 
	arma::colvec pi_hat,
	double l0,
	double w) 
{
	std::vector<int> lengths;
	long m = tau_grid.n_elem;
	arma::colvec laplace_grid = dlaplace(tau_grid, l0), delta(2);
	arma::uvec mixture(m), tau_grid_0 = find(tau_grid == 0);

	if(w > 0.5) {
		mixture  = w*laplace_grid > (1-w)*pi_hat; 
	} else {
		mixture = laplace_grid > pi_hat;
	}

	int i = 0;
	double prev = mixture[0];
	lengths.push_back(0);

	for(long j = 1; j < m; j++) {
		if (prev == mixture[j]) {
	      lengths[i]++;
		} else {
	      lengths.push_back(++lengths[i]);
	      i++;
	      prev = mixture[j];
	    }
	}

	int delta_0 = -1;
	for(int j = 0; j < (i+1); j++) {
		if(lengths[j] > tau_grid_0[0]) {delta_0 = j;break;}
	}

	if(delta_0 == 0) {
		delta(0) = 0; delta(1) = lengths[delta_0];
	} else if (delta_0 == i) {
		delta(0) = lengths[delta_0]; delta(1) = m-1;
	} else {
		delta(0) = lengths[delta_0-1];  delta(1) = lengths[delta_0];
	}

	// return List::create(
	// 	_["lengths"] 		 = lengths,
	// 	_["delta"] = delta,
	// 	_["delta_0"] = delta_0,
	// 	_["tau_grid_0"] = tau_grid_0[0],
	// 	_["i"] = i);
	return delta;
}


// [[Rcpp::export]]
arma::mat pi_mat(
	arma::mat likelihood_mat,
	arma::colvec tau_grid,
	arma::colvec pi_hat,
	double l0,
	double w, 
	long normalized = 1) 
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	arma::mat pi_t(m, n);
	arma::colvec laplace_grid = dlaplace(tau_grid, l0);
	arma::colvec pi_grid(m, arma::fill::zeros); // mixture prior

	pi_grid = w*laplace_grid + (1-w)*pi_hat;
	// for(long j = 0; j < m; j++) {
	// 	pi_grid(j) = w*laplace_grid(j) + (1-w)*pi_hat(j);
	// }

	for(long i = 0; i < n; i++) {
		double pi_j_sum = 0;
		// pi_t.col(i) = likelihood_mat.col(i) * pi_grid;
		for(long j = 0; j < m; j++) {
			pi_t(j,i) = likelihood_mat(j,i)*pi_grid(j);
			pi_j_sum += pi_t(j,i);
		}
		if(normalized) pi_t.col(i) = pi_t.col(i)/pi_j_sum;
	}

	return pi_t;
}

// [[Rcpp::export]]
double loglik_pi_mat(
	arma::mat likelihood_mat,
	arma::colvec tau_grid,
	arma::colvec pi_hat,
	double l0,
	double w) 
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	double loglik;
	arma::mat pi_t(m, n);
	arma::colvec laplace_grid = dlaplace(tau_grid, l0), norm(n, arma::fill::zeros);
	arma::colvec pi_grid(m, arma::fill::zeros); // mixture prior

	pi_grid = w*laplace_grid + (1-w)*pi_hat;
	for(long i = 0; i < n; i++) {
		for(long j = 0; j < m; j++) {
			pi_t(j,i) = likelihood_mat(j,i)*pi_grid(j);
		}
		norm(i) = log(sum(pi_t.col(i)));
	}
	loglik = sum(norm);
	
	return loglik;	
}

// [[Rcpp::export]]
double l_dl0_C (
	arma::mat& pi_t,
	arma::colvec& tau_grid, 
	arma::colvec& pi_hat,
	double l0,
	double w
	)
{
	arma::colvec tmp_grid = (1-arma::abs(tau_grid)*l0) % arma::exp(-arma::abs(tau_grid)*l0) % 
         (1/(w*dlaplace(tau_grid, l0) + (1-w)*pi_hat));
	double ret = sum(sum(pi_t, 1) % tmp_grid);

	return ret;
}

// [[Rcpp::export]]
arma::colvec update_pi_2(
	arma::mat& pi_t, 
	arma::colvec& tau_grid,
	double l0,
	double w
	)
{
	long m = pi_t.n_rows;
	arma::colvec pi_hat(m), rowsum = sum(pi_t, 1);

	pi_hat = rowsum/((1-w)*sum(rowsum)) - w/(1-w)*dlaplace(tau_grid, l0);
	for(auto &it : pi_hat) {if(it < 0) it = 0;}
	if(sum(pi_hat) != 0) pi_hat /= sum(pi_hat);

	return pi_hat;
}


// [[Rcpp::export]]
List C_snp_full (
	arma::mat likelihood_mat,
	arma::colvec tau_grid, 
	arma::colvec pi_hat,
	double l0,
	double w,
	long miter, double eps,
	int split = 1)
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	// long n_train = likelihood_mat_train.n_cols, n_test = likelihood_mat_test.n_cols;
	arma::colvec loglik_train(miter+1), posterior_means(n), pi_hat_old(m), prior_hat(m); 
	arma::colvec delta(2), delta_old(2);
	arma::colvec pi_zero = arma::zeros(m);
	arma::rowvec p_zero;
	arma::mat pi_t(m,n);
	long iter = 0;
	double ll, l0_old, w_old, l0_lower, l0_upper, eps_l = 1e-10;
	int stop_rule = 0;

	delta = find_delta_0(tau_grid, pi_hat, l0, w);

	ll = loglik_pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
	loglik_train(iter) = ll/n;

	while(iter < miter) {
		// Rcout << iter << "\t";
		// Rcout << iter << "\n" << ll << "\n";
		l0_old = l0; w_old = w; pi_hat_old = pi_hat; delta_old = delta;

		// pi_t
		pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);

		// update l0
		l0_lower = eps_l;
		l0_upper = 1;
		while(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * 
			l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) > 0) {
			double tmp = l0_upper;
			l0_upper = 2*l0_lower;
			l0_lower = tmp;
		}
		while(std::abs(l0_upper - l0_lower) > eps_l) {
			double tmp = (l0_upper + l0_lower)/2;
			double tmp_l_dl0 = l_dl0_C(pi_t, tau_grid, pi_hat, tmp, w);
			if(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * tmp_l_dl0 < 0) {
				l0_upper = tmp;
			}
			if (l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) * tmp_l_dl0 < 0) {
				l0_lower = tmp;
			}
		}

		l0 = (l0_upper + l0_lower)/2;

		// update gamma
		pi_hat = update_pi_2(pi_t, tau_grid, l0, w);

		// update delta
		delta = find_delta_0(tau_grid, pi_hat, l0, w);

		// update omega
		pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
		w = accu(pi_t.rows(delta(0), delta(1)))/n;
		arma::mat pi_laplacian = pi_mat(likelihood_mat, tau_grid, pi_zero, l0, w) * w;
		w += accu(pi_laplacian.rows(0,delta(0)))/n;
		w += accu(pi_laplacian.rows(delta(1), m-1))/n;

		pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
		ll = loglik_pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
		loglik_train(iter+1) = ll/n;

		// stop_rule = (ll <= (loglik_train(iter) + eps));
		stop_rule = (std::abs((loglik_train(iter+1) - loglik_train(iter)) / loglik_train(iter) )) < eps;
		
		if((iter+1 == miter) | stop_rule) {
			if(stop_rule) {
				l0 = l0_old; w = w_old; pi_hat = pi_hat_old; delta = delta_old;
			}
			pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
			posterior_means = est_posterior_means(pi_t, tau_grid);
		  	ll = loglik_pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
		  	prior_hat = w*dlaplace(tau_grid, l0) + (1-w)*pi_hat;

		  	p_zero = sum(pi_t.rows(delta(0), delta(1)), 0);
		  	arma::mat pi_laplacian = pi_mat(likelihood_mat, tau_grid, pi_zero, l0, w) * w;
			p_zero += sum(pi_laplacian.rows(0,delta(0)), 0);
			p_zero += sum(pi_laplacian.rows(delta(1), m-1), 0);
			break;
		}
		iter++;
	}

	Rcout << iter << "\n";

	return List::create(
		_["prior_hat"] 		 = prior_hat,
		_["pi_hat"] 		 = pi_hat,
		_["l0"]				 = l0,
		_["w"]				 = w,
		_["delta"]				 = delta,
		_["tau_grid"] 		 = tau_grid,
		_["posterior_means"] = posterior_means,
		_["posterior_hat"]   = pi_t,
		_["p_zero"]   		 = p_zero,
		_["loglik_train"]    = loglik_train.subvec(0, iter),
		_["loglik"]          = ll,
		_["iter"]			= iter);
}

// [[Rcpp::export]]
List C_snp (
	arma::mat likelihood_mat,
	arma::uvec itrain, arma::uvec itest,
	arma::colvec tau_grid, 
	arma::colvec pi_hat,
	double l0,
	double w,
	long miter, double eps,
	int split = 1)
{
	long n = likelihood_mat.n_cols, m = likelihood_mat.n_rows;
	long n_train = itrain.n_elem, n_test = itest.n_elem;
	arma::colvec loglik_test(miter+1),posterior_means(n), pi_hat_old(m), prior_hat(m); 
	arma::colvec delta(2), delta_old(2);
	arma::colvec pi_zero = arma::zeros(m);
	arma::rowvec p_zero;
	arma::mat pi_t(m, n);
	long iter = 0;
	double ll, l0_old, w_old, l0_lower, l0_upper, eps_l = 1e-10;
	int stop_rule = 0;

	itrain = itrain - 1;
	itest = itest - 1;

	delta = find_delta_0(tau_grid, pi_hat, l0, w);

	ll = loglik_pi_mat(likelihood_mat.cols(itest), tau_grid, pi_hat, l0, w);
	loglik_test(iter) = ll/n_test;

	while(iter < miter) {
		// Rcout << iter << "\n" << ll << "\n";
		l0_old = l0; w_old = w; pi_hat_old = pi_hat; delta_old = delta;

		// pi_t
		pi_t = pi_mat(likelihood_mat.cols(itrain), tau_grid, pi_hat, l0, w);

		// update l0
		l0_lower = eps_l;
		l0_upper = 1;
		while(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * 
			l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) > 0) {
			double tmp = l0_upper;
			l0_upper = 2*l0_lower;
			l0_lower = tmp;
		}
		while(std::abs(l0_upper - l0_lower) > eps_l) {
			double tmp = (l0_upper + l0_lower)/2;
			double tmp_l_dl0 = l_dl0_C(pi_t, tau_grid, pi_hat, tmp, w);
			if(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * tmp_l_dl0 < 0) {
				l0_upper = tmp;
			}
			if (l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) * tmp_l_dl0 < 0) {
				l0_lower = tmp;
			}
		}

		l0 = (l0_upper + l0_lower)/2;

		// update gamma
		pi_hat = update_pi_2(pi_t, tau_grid, l0, w);

		// update delta
		delta = find_delta_0(tau_grid, pi_hat, l0, w);

		// update omega
		pi_t = pi_mat(likelihood_mat.cols(itrain), tau_grid, pi_hat, l0, w);
		w = accu(pi_t.rows(delta(0), delta(1)))/n_train;
		arma::mat pi_laplacian = pi_mat(likelihood_mat.cols(itrain), tau_grid, pi_zero, l0, w) * w;
		w += accu(pi_laplacian.rows(0,delta(0)))/n_train;
		w += accu(pi_laplacian.rows(delta(1), m-1))/n_train;

		pi_t = pi_mat(likelihood_mat.cols(itrain), tau_grid, pi_hat, l0, w);
		ll = loglik_pi_mat(likelihood_mat.cols(itest), tau_grid, pi_hat, l0, w);
		loglik_test(iter+1) = ll/n_test;

		// stop_rule = (ll <= (loglik_test(iter) + eps));
		stop_rule = (std::abs((loglik_test(iter+1) - loglik_test(iter)) / loglik_test(iter) )) < eps;
		// stop_rule = (std::abs(ll-loglik_train(iter)/loglik_train(iter)) < eps |
		// 	(ll-loglik_train(iter))/std::abs(loglik_train(iter)) < eps );
		if((iter+1 == miter) | stop_rule) {
			if(stop_rule) {
				l0 = l0_old; w = w_old; pi_hat = pi_hat_old; delta = delta_old;
			}
			pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);

			// update l0
			l0_lower = eps_l;
			l0_upper = 1;
			while(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * 
				l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) > 0) {
				double tmp = l0_upper;
				l0_upper = 2*l0_lower;
				l0_lower = tmp;
			}
			while(std::abs(l0_upper - l0_lower) > eps_l) {
				double tmp = (l0_upper + l0_lower)/2;
				double tmp_l_dl0 = l_dl0_C(pi_t, tau_grid, pi_hat, tmp, w);
				if(l_dl0_C(pi_t, tau_grid, pi_hat, l0_lower, w) * tmp_l_dl0 < 0) {
					l0_upper = tmp;
				}
				if (l_dl0_C(pi_t, tau_grid, pi_hat, l0_upper, w) * tmp_l_dl0 < 0) {
					l0_lower = tmp;
				}
			}

			l0 = (l0_upper + l0_lower)/2;

			// update gamma
			pi_hat = update_pi_2(pi_t, tau_grid, l0, w);

			// update delta
			delta = find_delta_0(tau_grid, pi_hat, l0, w);

			// update omega
			pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
			w = accu(pi_t.rows(delta(0), delta(1)))/n;
			arma::mat pi_laplacian = pi_mat(likelihood_mat, tau_grid, pi_zero, l0, w) * w;
			w += accu(pi_laplacian.rows(0,delta(0)))/n;
			w += accu(pi_laplacian.rows(delta(1), m-1))/n;
			
			pi_t = pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
			posterior_means = est_posterior_means(pi_t, tau_grid);
		  	ll = loglik_pi_mat(likelihood_mat, tau_grid, pi_hat, l0, w);
		  	prior_hat = w*dlaplace(tau_grid, l0) + (1-w)*pi_hat;

		  	p_zero = sum(pi_t.rows(delta(0), delta(1)), 0);
			p_zero += sum(pi_laplacian.rows(0,delta(0)), 0);
			p_zero += sum(pi_laplacian.rows(delta(1), m-1), 0);
			break;
		}
		iter++;
	}

	Rcout << iter << "\n";


	return List::create(
		_["prior_hat"] 		 = prior_hat,
		_["pi_hat"] 		 = pi_hat,
		_["l0"]				 = l0,
		_["w"]				 = w,
		_["delta"]			 = delta,
		_["tau_grid"] 		 = tau_grid,
		_["posterior_means"] = posterior_means,
		_["posterior_hat"]   = pi_t,
		_["p_zero"]   		 = p_zero,
		_["loglik_test"]    = loglik_test.subvec(0, iter),
		_["loglik"]          = ll,
		_["iter"]			= iter);
}