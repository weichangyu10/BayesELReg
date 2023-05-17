#define ARMA_USE_OPENMP

#include<cmath>
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
arma::mat GSLRcpp(const arma::vec& tet, const arma::mat& x){
  
  const int p = tet.n_elem;
  
  const arma::vec y = x.col(p);
  const arma::mat xpred = x.cols(0,p-1);
  const arma::vec u =  y - xpred * tet;
  const arma::mat ans = xpred.t() * diagmat(u);
  
  return ans.t();
}

// [[Rcpp::export]]
arma::mat STARDIntermediateRcpp(const arma::vec& tet, const arma::mat& x){
  
  const int p = x.n_cols - 2;
  const int numCoef = tet.n_elem/2;
  
  const arma::vec x21 = x.col(p);
  const arma::vec x22 = x.col(p+1);
  const arma::mat xpred = x.cols(0,p-1);
  const arma::vec u1 =  x21 - xpred * tet.rows(0,numCoef-1);
  const arma::vec u2 =  x22 - xpred * tet.rows(numCoef,2*numCoef-1);
  const arma::mat ans1 = xpred.t() * diagmat(u1);
  const arma::mat ans2 = xpred.t() * diagmat(u2);
  const arma::mat ans = join_cols(ans1,ans2);
  
  return ans.t();
}




// [[Rcpp::export]]
arma::vec NablaAELCh2(const arma::mat& X, const arma::mat& Hmat, const arma::vec& lambda, const arma::vec& extrah, const double& a, const arma::vec& tet){
  
  const int p = X.n_cols;
  const int n = X.n_rows;
  
  const arma::vec WeightsOEL = ones(n)/(ones(n) + Hmat*lambda);
  const double PerturbWeight = a*(1.0/n)/(1.0 + accu(extrah.t()*lambda));
  const arma::vec FirstPart = X.t() * diagmat(WeightsOEL) * X*lambda;
  const arma::vec SecondPart = PerturbWeight* (X.t() *X*lambda);
  const arma::vec ans = FirstPart - SecondPart - tet/(10000.0*ones(p));
  
  
  return ans;
}

// [[Rcpp::export]]
arma::vec NablaAELChSTARD(const arma::mat& X, const arma::mat& Hmat, const arma::vec& lambda, const arma::vec& extrah, const double& a, const arma::vec& tet){
  
  const int n = X.n_rows;
  
  const arma::vec WeightsOEL = ones(n)/(ones(n) + Hmat*lambda);
  const double PerturbWeight = a*(1.0/n)/(1.0 + accu(extrah.t()*lambda));
  const arma::mat XtransXWeighted = X.t() * diagmat(WeightsOEL) * X;
  const arma::mat XtransX = X.t() * X;
  const arma::mat zeroMAT = zeros(X.n_cols,X.n_cols);
  const arma::vec FirstPart = join_cols(join_rows(XtransXWeighted, zeroMAT), join_rows(zeroMAT,XtransXWeighted))*lambda;
  const arma::vec SecondPart = PerturbWeight*(join_cols(join_rows(XtransX, zeroMAT), join_rows(zeroMAT, XtransX)))*lambda;
  const arma::vec ans = FirstPart - SecondPart - tet/(10000.0*ones(tet.n_elem));

  return ans;
}

arma::mat OptTrtSTARDemplik(const arma::mat& etaHMC, const arma::mat& omega1HMC, const arma::mat& omega2HMC, const arma::vec& x1, const arma::mat& XiMatrix1, const arma::mat& XiMatrix2){
  
  const int B = etaHMC.n_rows;
  const int N = XiMatrix2.n_cols;
  const int A1size = 2;
  const arma::mat eta1HMC = etaHMC.cols(0,  (etaHMC.n_cols/2-1));
  const arma::mat eta0HMC = etaHMC.cols((etaHMC.n_cols/2), (etaHMC.n_cols-1));
  const arma::mat eta13 = eta1HMC.cols((eta1HMC.n_cols-2), (eta1HMC.n_cols-1));
  const arma::mat eta03 = eta0HMC.cols((eta0HMC.n_cols-2), (eta0HMC.n_cols-1));
  arma::mat EYoptMatrix(B,A1size,fill::zeros);
  arma::vec omega11temp, omega21temp, eta11temp, eta01temp;
  arma::mat eta12temp, eta02temp, omega12temp, omega22temp;
  double frontVectemp, backVectemp;
  
  for(int s = 1; s <= A1size; s++){
    
    eta11temp = eta1HMC.col((s-1)*5);
    eta01temp = eta0HMC.col((s-1)*5);
    eta12temp = eta1HMC.cols( ((s-1)*5+1), ((s-1)*5+2) );
    eta02temp = eta0HMC.cols( ((s-1)*5+1), ((s-1)*5+2) );
    omega11temp = omega1HMC.col(((s-1)*3));
    omega12temp = omega1HMC.cols(((s-1)*3+1), ((s-1)*3+2));
    omega21temp = omega2HMC.col(((s-1)*3));
    omega22temp = omega2HMC.cols(((s-1)*3+1), ((s-1)*3+2));
    backVectemp = 0.0;
    for(int b = 0; b < B; b++){
      
      frontVectemp = 0.5*( (eta11temp(b) +  eta01temp(b)) + (eta12temp(b,0)+ eta02temp(b,0))*x1(0) +  (eta12temp(b,1)+ eta02temp(b,1))*x1(1) + (eta13(b,0) + eta03(b,0))*( omega11temp(b) + omega12temp(b,0)*x1(0) + omega12temp(b,1)*x1(1) ) + (eta13(b,1) + eta03(b,1))*( omega21temp(b) + omega22temp(b,0)*x1(0) + omega22temp(b,1)*x1(1) )   );
      for(int i = 0; i < N; i++){
        
        backVectemp += abs( 0.5*( ( eta13(b,0) - eta03(b,0) )*( omega11temp(b) + omega12temp(b,0)*x1(0) + omega12temp(b,1)*x1(1) + XiMatrix1(b,i) ) + ( eta13(b,1) - eta03(b,1) )*( omega21temp(b) + omega22temp(b,0)*x1(0) + omega22temp(b,1)*x1(1) + XiMatrix2(b,i) )  + (eta11temp(b) -  eta01temp(b)) + (eta12temp(b,0)- eta02temp(b,0))*x1(0) +  (eta12temp(b,1)- eta02temp(b,1))*x1(1)   ) );
          //backVec <- rowMeans(apply( XiArray, c(1,2), function(s){ abs(0.5*( sum((eta13 - eta03)*(c( omega11.temp + sum(omega12.temp*x1), omega21.temp + sum(omega22.temp*x1)  ) + s)  ) + (eta11.temp - eta01.temp) + sum((eta12.temp - eta02.temp)*x1)  ))   }  ))
        
      }
      backVectemp = backVectemp/(N*1.0);
      EYoptMatrix(b,(s-1)) = frontVectemp + backVectemp;
      
    }
    
    
  }
  
  return EYoptMatrix;
}

