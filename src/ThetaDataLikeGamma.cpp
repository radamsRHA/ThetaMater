//Includes/namespaces
#include <Rcpp.h>
#include <boost/math/special_functions/binomial.hpp>

using namespace Rcpp;

//' @title
//' ThetaDataLikeGamma
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
//' \code{ThetaDataLikeGamma} Computes the likelihood of a dataset for a given value of theta
//' @export
// [[Rcpp::export]]
double ThetaDataLikeGamma(double theta, NumericVector kvec, NumericVector lvec, NumericVector nvec, NumericVector cvec, NumericVector rvec, int Kclasses) {
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
		double avgR = 0;
		double plocus = 0.0;
		double ThetaL_avgR = 0;

		for (int j = 0; j<numLoci;++j){
				k = kvec[j];
				l = lvec[j];
				n = nvec[j];
				cl = cvec[j];

				ThetaL = theta*l;
				plocus = 0.0; 

				for (int r = 0; r<Kclasses; ++r){
					avgR = rvec[r];
					ThetaL_avgR = ThetaL*avgR;
					pSk = 0;
					
					for (int i = 2; i<=n;++i){
						a = pow(-1, i);
						b = boost::math::binomial_coefficient<double>(n-1, i-1);
						c = (i-1)/(ThetaL_avgR+i-1);
						d = ThetaL_avgR/(ThetaL_avgR+i-1);
						e = pow(d,k);
			
						f = a*b*c*e;
						pSk += f;
					}
					plocus += pSk*1/Kclasses;
				}
				DataLnL += cl*log(plocus);
			}
		return DataLnL;
		}
					
