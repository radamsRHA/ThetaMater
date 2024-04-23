//Includes/namespaces
#include <Rcpp.h>
#include <boost/math/special_functions/binomial.hpp>

using namespace Rcpp;

//' @title
//' ThetaDataLike
//' @description
//' Function to compute the probability of a given dataset
//' 
//' @param theta Population parameter to compute the likelihood of the dataset
//' 
//' @param kvec Vector of the mutation counts for each locus
//' 
//' @param lvec Vector of the locus lengths for each locus
//' 
//' @param nvec Vector of the number of samples for each locus
//' 
//' @param cvec Vector of the number of pattern counts
//' 
//' @details
//' \code{ThetaDataLike} Computes the likelihood of a dataset for a given value of theta
//' @export
// [[Rcpp::export]]
double ThetaDataLike(double theta, NumericVector kvec, NumericVector lvec, NumericVector nvec, NumericVector cvec) {
		int numLoci = kvec.size();
		double DataLnL = 0;
		double k;
		double l;
		double n;
		double ThetaL = 0;
		double pSk = 0;
		double a = 0;
		double b = 0;
		double c = 0;
		double d = 0;
		double e = 0;
		double f = 0;
		double cl = 0;

		for (int j = 0; j<numLoci;++j){
				pSk = 0;
				k = kvec[j];
				l = lvec[j];
				n = nvec[j];
				cl = cvec[j];

				ThetaL = theta*l;
				for (int i = 2; i<=n;++i){
					a = pow(-1, i);
					b = boost::math::binomial_coefficient<double>(n-1, i-1);
					c = (i-1)/(ThetaL+i-1);
					d = ThetaL/(ThetaL+i-1);
					e = pow(d,k);
			
					f = a*b*c*e;
					pSk += f;
				}
				DataLnL += cl*log(pSk);
			}
		return DataLnL;
		}
					
