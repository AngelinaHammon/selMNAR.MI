// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>

using namespace Rcpp;
using namespace arma;
using namespace std;




// Function for normal cdf
inline double N(double z) {
  if(z > 6.0) {return 1.0;};
  if(z < -6.0) {return 0.0;};
  
  double b1 = 0.31938153;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double p = 0.2316419;
  double c2 = 0.3989423;
  
  double a = fabs(z);
  double t = 1.0/(1.0+a*p);
  double b = c2*exp((-z)*(z/2.0));
  double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  n = 1.0-b*n;
  if(z < 0.0) n = 1.0 - n;
  return n;
}

// Function for bivariate normal cdf
inline double bivnor(double h1, double hk, double r)
{
  int NX=5;
  Rcpp::NumericVector X(NX);
  Rcpp::NumericVector W(NX);
  // data
  X[0]=.04691008;
  X[1]=.23076534;
  X[2]=.5;
  X[3]=.76923466;
  X[4]=.95308992;
  W[0]=.018854042;
  W[1]=.038088059;
  W[2]=.0452707394;
  W[3]=.038088059;
  W[4]=.018854042;
  // declarations
  double bv = 0;
  double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
  double cor_max = 0.7;
  double bv_fac1 = 0.13298076;
  double bv_fac2 = 0.053051647;
  
  // computation
  double h2 = hk;
  double h12 = (h1*h1+h2*h2)/2;
  double r_abs = std::abs(r);
  if (r_abs > cor_max){
    r2 = 1.0 - r*r;
    r3 = std::sqrt(r2);
    if (r<0){
      h2 = -h2;
    }
    h3 = h1*h2;
    h7 = std::exp( -h3 / 2.0);
    if ( r_abs < 1){
      h6 = std::abs(h1-h2);
      h5 = h6*h6 / 2.0;
      h6 = h6 / r3;
      aa = 0.5 - h3 / 8.0;
      ab = 3.0 - 2.0 * aa * h5;
      bv = bv_fac1*h6*ab*(1-N(h6))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
      for (int ii=0; ii<NX; ii++){
        r1 = r3*X[ii];
        rr = r1*r1;
        r2 = std::sqrt( 1.0 - rr);
        bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
      }
    }
    h11 = std::min(h1,h2);
    bv = bv*r3*h7 + N(h11);
    if (r < 0){
      bv = N(h1) - bv;
    }
    
  } else {
    h3=h1*h2;
    for (int ii=0; ii<NX; ii++){
      r1 = r*X[ii];
      rr2 = 1.0 - r1*r1;
      bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
    }
    bv = N(h1)*N(h2) + r*bv;
  }
  //--- OUTPUT
  return bv;
}

const double PI = M_PI;

// Function for normal pdf
inline double dn(double z) {
  return (1.0/sqrt(2.0*PI))*exp(-0.5*z*z);
}

// Function for bivariate normal pdf
inline double biv_pdf(double x, double y, double rho) {
  double u = x / 1 ;
  double v = y / 1 ;
  double c = 1 - rho*rho ;
  double p = (1 / (2 * M_PI * 1 * 1 * sqrt(c))) 
    * exp (-(u * u - 2 * rho * u * v + v * v) / (2 * c));
  return p;
}

// Function for group product
NumericVector groupProd(NumericVector v, NumericVector group){
  vec V = as<vec>(v);
  ivec G = as<ivec>(group);
  uword nGroup = G.n_elem - 1;
  vec z = zeros<vec> (nGroup);
  for(uword i=0; i< nGroup; i++) z[i] = prod(V.subvec(G[i], G[i+1]-1));
  return wrap(z);
}

// Function for group sum
NumericVector groupSum(NumericVector v, NumericVector group){
  vec V = as<vec>(v);
  ivec G = as<ivec>(group);
  uword nGroup = G.n_elem - 1;
  vec z = zeros<vec> (nGroup);
  for(uword i=0; i< nGroup; i++) z[i] = sum(V.subvec(G[i], G[i+1]-1));
  return wrap(z);
}



// From NumericMatrix to arma::mat
arma::mat trans(NumericMatrix x) {
  arma::mat y = as<arma::mat>(x) ;
  return(y) ;
}

// From NumericVector to arma::vec
arma::vec trans2(NumericVector x) {
  arma::vec y = as<arma::vec>(x) ;
  return(y) ;
}


// From arma::vec to NumericVector
NumericVector trans_back(arma::vec x) {
  NumericVector y = wrap(x) ;
  return(y) ;
}


// Function for Matrix multiplication with vector
NumericVector matmult(NumericMatrix x, NumericVector y) {
  arma::mat x1 = trans(x);
  arma::vec y1 = trans2(y);
  arma::vec mult = x1*y1;
  NumericVector res=trans_back(mult);
  return(res);
}

// Function which multiplies scalar with vector and the resulting matrix with a vector
NumericVector matmult_total(NumericMatrix x, NumericVector y, double z) {
  NumericMatrix res1 = x*z;
  NumericVector res2 = matmult(res1,y);
  return(res2);
}



// Equivalent function for "which()"
IntegerVector c_which(IntegerVector group, int elem){
  int n=group.size();
  
  std::vector<int> ind(n, 0); // Create a n element vector containing all 0's
  
  for(int i=0; i < n; ++i) { 
    ind[i] = group[i];
  };
  
  //std::vector<int> ind=group;
  std::vector<int> out(ind);
  std::vector<int>::iterator it;
  int j = 0;
  it = std::find(ind.begin(), ind.end(), elem);
  while(it++ != ind.end()){
    out[j++] = it - ind.begin();
    it = std::find(it, ind.end(), elem);
  }
  out.resize(j);
  
  int r=out.size();
  IntegerVector out2(r);
  
  for(int i=0; i < r; ++i) { 
    out2[i] = out[i];
  };
  
  return out2;
}


//' Internal functions written in C++ to speed up maximization process of the censored bivariate probit model with random intercepts
//'
//' Functions that calculate the modes of the random intercepts and that contain the loops for computing the log likelihood approximated by GHQ.
//' 
//' @keywords internal
// [[Rcpp::export]]
double mode(NumericVector par, int x, IntegerVector yS, IntegerVector yO, NumericVector linpredS, NumericVector linpredO, NumericVector group_aux, double rho, double tau, double sigma1, double sigma2) {
  
  double alpha1 = par[0];
  double alpha2 = par[1];
  
  int n = yS.size();
  
  NumericVector sum(n);
  
  if((round(tau*100)/100)==1){
    tau = 0.9999;
  };
  
  if((round(tau*100)/100)==-1){
    tau = -0.9999;
  };
  
  if((round(rho*100)/100)==1){
    rho = 0.9999;
  };
  
  if((round(rho*100)/100)==-1){
    rho = -0.9999;
  };
  
  if(sigma1 > 700){
    sigma1 = exp(300);
  };
  
  if(sigma2 > 700){
    sigma2 = exp(300);
  };
  
  if(sigma1 < -700){
    sigma1 = exp(-300);
  };
  
  if(sigma2 < -700){
    sigma2 = exp(-300);
  };
  
  
  if(sigma1 == 0){
    sigma1 =  1.1111e-300 ;
  };
  
  if(sigma2 == 0){
    sigma2 =  1.1111e-300 ;
  };
  
  for(int i=0; i < n; ++i) {  
    if (yS[i]==0) {
      sum[i] = N(-(linpredS[i]+alpha1));
    }
    else if (yS[i]==1 && yO[i]==1) {
      sum[i] = bivnor(linpredS[i]+alpha1,linpredO[i]+alpha2,rho);
    }
    else if (yS[i]==1 && yO[i]==0) {
      sum[i] = bivnor(linpredS[i]+alpha1,-(linpredO[i]+alpha2),-rho);
    };
    if (sum[i]==0) {
      sum[i] = 1.1111e-300;
    };
  };
  
  NumericVector prod = groupSum(log(sum), group_aux);
  
  double prod2 = prod[x-1];
  double sd1 = sqrt(sigma1);
  double sd2 = sqrt(sigma2);
  double temp5 = 1 - pow(tau,2);
  
  double zedd1 = alpha1/sd1;
  double zedd2 = alpha2/sd2;
  
  double pdf = exp(-log(2 * PI) - log(sd1) - log(sd2) - 0.5 * log(1-pow(tau,2)) -(0.5/temp5) * (pow(zedd1,2) + (-2 * tau * zedd1 + zedd2) * zedd2));
  
  
  if(pdf==0) {
    pdf = 1.1111e-300;
  };  
  
  double func = prod2+log(pdf);
  return (-func);
  
}

//' @rdname mode
//' @keywords internal
// [[Rcpp::export]]
double loop1(IntegerVector yS, IntegerVector yO, NumericVector linpredS, NumericVector linpredO, double rho, double tau, double sigma1, double sigma2, NumericVector w1, NumericVector a1, NumericVector group_aux, IntegerVector group, int M, NumericMatrix chol_vc_re) {
  
  int n = yS.size();
  
  NumericVector w2 = w1; 
  NumericVector a2 = a1;
  
  int la1 = a1.size();  
  int la2 = la1;
  
  NumericVector sum(n);
  
  NumericVector Li(M);
  
  for(int j=0; j < M; ++j) {
    Li[j] = 0;
  };  
  
  
  for(int p=0; p < la2; ++p) {
    
    for(int m=0; m < la1; ++m) {
      
      NumericVector a=NumericVector::create(a1[m],a2[p]);
      
      NumericVector a_st=matmult_total(chol_vc_re, a, sqrt(2));
      double a1_st = a_st[0];
      double a2_st = a_st[1];
      
      for(int i=0; i < n; ++i) {  
        if (yS[i]==0) {
          sum[i] = N(-(linpredS[i]+a1_st));
        }
        else if (yS[i]==1 && yO[i]==1) {
          sum[i] = bivnor(linpredS[i]+a1_st,linpredO[i]+a2_st,rho);
        }
        else if (yS[i]==1 && yO[i]==0) {
          sum[i] = bivnor(linpredS[i]+a1_st,-(linpredO[i]+a2_st),-rho);
        };
      };
      
      NumericVector prod = groupProd(sum,group_aux);
      Li += (1/PI*w1[m]*w2[p]*prod);
      
    };
    
  }; 
  
  for(int j=0; j < M; ++j) {
    if (Li[j]==0) {
      Li[j] = 1.1111e-300;
    };
  };  
  
  NumericVector Li_log = log(Li);
  
  double sum_Li = std::accumulate(Li_log.begin(), Li_log.end(), 0.000000000000000000000);
  
  double loglik = sum_Li;    
  return -loglik;
}


//' Internal functions written in C++ to speed up maximization process of the ordered probit model with sample selection and random intercepts
//'
//' Functions that calculate the modes of the random intercepts and that contain the loops for computing the log likelihood approximated by GHQ.
//' 
//' @keywords internal
// [[Rcpp::export]]
double mode_ord(NumericVector par, int x, IntegerVector yS, IntegerVector yO, NumericVector linpredS, NumericVector linpredO, NumericVector group_aux, double rho, double tau, NumericVector z, double sigma1, double sigma2, IntegerVector levels) {
  
  double alpha1 = par[0];
  double alpha2 = par[1];
  
  int n = yS.size();
  
  NumericVector sum(n);
  
  int lz = z.size();
  
  NumericVector z_aux(lz + 2);
  z_aux[0] = R_NegInf;
  z_aux[lz+1] = R_PosInf;
  
  for(int i=0; i < lz; ++i) {
    z_aux[i+1] = z[i];
  };
  
  
  if((round(tau*100)/100)==1){
    tau = 0.9999;
  };
  
  if((round(tau*100)/100)==-1){
    tau = -0.9999;
  };
  
  if((round(rho*100)/100)==1){
    rho = 0.9999;
  };
  
  if((round(rho*100)/100)==-1){
    rho = -0.9999;
  };
  
  if(sigma1 > 700){
    sigma1 = exp(300);
  };
  
  if(sigma2 > 700){
    sigma2 = exp(300);
  };
  
  if(sigma1 < -700){
    sigma1 = exp(-300);
  };
  
  if(sigma2 < -700){
    sigma2 = exp(-300);
  };
  
  
  if(sigma1 == 0){
    sigma1 =  1.1111e-300 ;
  };
  
  if(sigma2 == 0){
    sigma2 =  1.1111e-300 ;
  };
  
  
  for(int i=0; i < n; ++i) {  
    if (yS[i]==0) {
      sum[i] = N(-(linpredS[i]+alpha1));
    }
    else if(yS[i]==1 && yO[i]==levels[0]) {
      sum[i] = bivnor(linpredS[i]+alpha1,z_aux[1]-(linpredO[i]+alpha2),-rho);
    }
    else if(yS[i]==1 && yO[i]==levels[lz]) {
      sum[i] = bivnor(linpredS[i]+alpha1,100000000000000000 -(linpredO[i]+alpha2),-rho) -  bivnor(linpredS[i]+alpha1,z_aux[lz]-(linpredO[i]+alpha2),-rho);
    }
    else{
      for(int h=1; h < lz; ++h) {
        if (yS[i]==1 && yO[i]==levels[h]) {
          sum[i] = bivnor(linpredS[i]+alpha1,z_aux[h+1]-(linpredO[i]+alpha2),-rho) -  bivnor(linpredS[i]+alpha1,z_aux[h]-(linpredO[i]+alpha2),-rho);
        };
      };
    };
    if (sum[i]==0) {
      sum[i] = 1.1111e-300;
    };
    if(sum[i]<0){
      sum[i] = 1.1111e-30;
    };
  };
  
  
  NumericVector prod = groupSum(log(sum), group_aux);
  
  double prod2 = prod[x-1];
  double sd1 = sqrt(sigma1);
  double sd2 = sqrt(sigma2);
  double temp5 = 1 - pow(tau,2);
  
  double zedd1 = alpha1/sd1;
  double zedd2 = alpha2/sd2;
  
  double pdf = exp(-log(2 * PI) - log(sd1) - log(sd2) - 0.5 * log(1-pow(tau,2)) -(0.5/temp5) * (pow(zedd1,2) + (-2 * tau * zedd1 + zedd2) * zedd2));
  
  
  if(pdf==0) {
    pdf = 1.1111e-300;
  };  
  
  double func = prod2+log(pdf);
  return (-func);
  
}


//' @rdname mode_ord
//' @keywords internal
// [[Rcpp::export]]
double loop1_ord(IntegerVector yS, IntegerVector yO, NumericVector linpredS, NumericVector linpredO, double rho, double tau, double sigma1, double sigma2, NumericVector z, IntegerVector levels, NumericVector w1, NumericVector a1, NumericVector group_aux, IntegerVector group, int M, NumericMatrix chol_vc_re) {
  
  int n = yS.size();
  
  NumericVector w2 = w1; 
  NumericVector a2 = a1;
  
  int la1 = a1.size();  
  int la2 = la1;
  
  int lz = z.size();
  
  NumericVector sum(n);
  
  NumericVector Li(M);
  
  for(int j=0; j < M; ++j) {
    Li[j] = 0;
  };  
  
  NumericVector z_aux(lz + 2);
  z_aux[0] = R_NegInf;
  z_aux[lz+1] = R_PosInf;
  
  for(int i=0; i < lz; ++i) {
    z_aux[i+1] = z[i];
  };
  
  for(int p=0; p < la2; ++p) {
    
    for(int m=0; m < la1; ++m) {
      
      NumericVector a=NumericVector::create(a1[m],a2[p]);
      
      NumericVector a_st=matmult_total(chol_vc_re, a, sqrt(2));
      double a1_st = a_st[0];
      double a2_st = a_st[1];
      
      for(int i=0; i < n; ++i) {  
        if (yS[i]==0) {
          sum[i] = N(-(linpredS[i]+a1_st));
        }else if(yS[i]==1 && yO[i]==levels[0]) {
          sum[i] = bivnor(linpredS[i]+a1_st,z_aux[1]-(linpredO[i]+a2_st),-rho);
        }else if(yS[i]==1 && yO[i]==levels[lz]) {
          sum[i] = bivnor(linpredS[i]+a1_st,100000000000000000 -(linpredO[i]+a2_st),-rho) -  bivnor(linpredS[i]+a1_st,z_aux[lz]-(linpredO[i]+a2_st),-rho);
        }else{
          for(int h=1; h < lz; ++h) {
            if (yS[i]==1 && yO[i]==levels[h]) {
              sum[i] = bivnor(linpredS[i]+a1_st,z_aux[h+1]-(linpredO[i]+a2_st),-rho) -  bivnor(linpredS[i]+a1_st,z_aux[h]-(linpredO[i]+a2_st),-rho);
            };
          };
        };
      };
      
      
      NumericVector prod = groupProd(sum,group_aux);
      Li += (1/PI*w1[m]*w2[p]*prod);
      
    };
    
  }; 
  
  for(int j=0; j < M; ++j) {
    if (Li[j]==0) {
      Li[j] = 1.1111e-300;
    };
  };  
  
  NumericVector Li_log = log(Li);
  
  double sum_Li = std::accumulate(Li_log.begin(), Li_log.end(), 0.000000000000000000000);
  
  double loglik = sum_Li;    
  return -loglik;
}



