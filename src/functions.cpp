#include <Rcpp.h>
using namespace Rcpp;

NumericVector mids(NumericVector breaks){
  int len = breaks.length() ;
  NumericVector mids = NumericVector(len -1);
  for (int i = 0 ; i < len -1 ; i++ ) {
    mids[i] = (breaks[i +1] + breaks[i])/2;
  }
  return mids;
}
// [[Rcpp::export]]
NumericVector mybreaks_rcpp(NumericVector x, double nbr) {
  double* mx = std::min_element(x.begin(), x.end());
  double* Mx = std::max_element(x.begin(), x.end());
  NumericVector res = NumericVector(nbr);
  double offset = *mx;
  double h = (*Mx-*mx)/(nbr-1);
  for (int i=0; i<nbr; i++) {
    res[i] = offset;
    offset += h;
  }
  return res;
}

// [[Rcpp::export]]
List hist_rcpp(NumericVector x, NumericVector breaks) {
  
  /*x = x.sort();*/
  double len = breaks.length();
  double nb_br;
  NumericVector r_breaks;
  if (len == 1) {
    r_breaks = mybreaks_rcpp(x, breaks[0]);
    nb_br = breaks[0];
  }
  else {
    r_breaks = breaks;
    nb_br = len;
  }
  
  NumericVector counts = NumericVector(nb_br - 1);
  /*
   * Old code assuming x is sorted
   */
  /*int k=0;
   int j=1;
   for (int i=0; i<x.length(); i++) {
   if (x[i] >= r_breaks[j]) {
   counts[j-1] = k;
   k = 0;
   j = j + 1;
   }
   else {
   k = k + 1;
   }
   }*/
  /*
   * New code, not causing RStudio fatal error
   */
  for (int c=0; c<nb_br-1; c++) {
    counts[c] = 0;
  }
  for (int i=0; i<x.length(); i++) {
    for (int j=0; j<nb_br-1; j++) {
      if ((r_breaks[j] <= x[i]) && (x[i] <= r_breaks[j+1])) {
        counts[j] = counts[j] + 1;
        break;
      }
    }
  }
  NumericVector density = counts / sum(counts);
  density = density / diff(r_breaks);
  List res = List::create(
    Named("breaks")=r_breaks,
    Named("counts")=counts,
    Named("density")=density,
    Named("mids")=mids(r_breaks));
  
  return res;
}



// [[Rcpp::export]]
double riskhistRcpp(NumericVector obs01, double m) {
  
 
  double h = 1 / m;
  double n = obs01.length();
  NumericVector breaks = NumericVector(m+1);
  double s = 0;
  double interval = 1/(m);
  for (int i=0; i<m+1; i++) {
    breaks[i] = s;
    s = s+interval;
  }
  

  NumericVector p_hat = as<NumericVector>(hist_rcpp(obs01, breaks = breaks)["counts"]) / n;

  
  double res = 2 / h / (n - 1) - (n + 1) / (n - 1) / h * sum(pow(p_hat, 2));
  return res; /*(m * sum(p_hat^2)) */
}

// [[Rcpp::export]]
Rcpp::List broptRcpp(NumericVector x){
  int len = 5*floor(sqrt(x.length()));
  NumericVector Mgrid = NumericVector(len-1);

  for (int i = 2; i<=len; i++) {
    Mgrid[i-2] = i;
  }
  NumericVector J = NumericVector(Mgrid.length());
  //NumericVector lim1 = NumericVector::create(min(x)-0.5, max(x)+0.5);
  NumericVector xlim = range(x);
  xlim[0] -= 0.5;
  xlim[1] += 0.5;
  NumericVector obs01  = (x - xlim[0]) / (xlim[1] - xlim[0]);
  for(int m=0; m < Mgrid.length(); m++) {
    //Function riskhist("riskhist");
    J[m] = (riskhistRcpp(obs01, Mgrid[m]));
  }
  return Rcpp::List::create(Rcpp::Named("opt")=Mgrid[which_min(J)]);
}

Rcpp::List bropt_and_broptfp_Rcpp(NumericVector x) {
  int len = 5*floor(sqrt(x.length()));
  NumericVector Mgrid = NumericVector(len-1);
  
  for (int i = 2; i<=len; i++) {
    Mgrid[i-2] = i;
  }
  NumericVector J = NumericVector(Mgrid.length());
  NumericVector J_fp = NumericVector(Mgrid.length()); 
  //NumericVector lim1 = NumericVector::create(min(x)-0.5, max(x)+0.5);
  NumericVector xlim = range(x);
  xlim[0] -= 0.5;
  xlim[1] += 0.5;
  NumericVector obs01  = (x - xlim[0]) / (xlim[1] - xlim[0]);
  for(int m=0; m < Mgrid.length(); m++) {
    List vec = riskhist_and_fp_Rcpp(obs01,Mgrid[m])
    J[m] = vec[0]
    J_fp[m] = vec[1]
  }
  return Rcpp::List::create(Rcpp::Named("opt")=Mgrid[which_min(J)], Rcpp::Named("opt_fp")=max(5, Mgrid[which.min(J_fp)]));
}

// [[Rcpp::export]]
NumericVector riskhist_and_fp_Rcpp(NumericVector obs01, double m){
  
  double h      = 1 / m;
  int n      = obs01.length();
  NumericVector breaks = mybreaks_rcpp(NumericVector::create(0, 1), m + 1);
  
  NumericVector p_hat_fp  = hist_rcpp(obs01, breaks = breaks)["counts"];
  
  /*###riskhist*/
  NumericVector p_hat  = p_hat_fp / n;
  double res1 = sum(pow(p_hat, 2));
  double res = 2 / h / (n - 1) - (n + 1) / (n - 1) / h * res1;
  /*#######*/
  
  /*###fp####*/
  
  int len_p_hat = p_hat_fp.length();
  NumericVector vs1 = NumericVector(len_p_hat + 2);
  NumericVector vs2 = NumericVector(len_p_hat + 2);
  NumericVector vs3 = NumericVector(len_p_hat + 2);
  double p_i;
  for (int i=0; i<len_p_hat; i++) {
    p_i = p_hat_fp[i];
    vs3[i] = p_i;
    vs2[i+1] = -2 * p_i;
    vs1[i+2] = p_i;
  }
  
  double res1_fp = 0;
  for (int i=0; i<len_p_hat + 2; i++) {
    res1_fp = res1_fp + pow(vs1[i] + vs2[i] + vs3[i], 2);
  }
  double res_fp = 271 / (480 * n * h) + 49 / (2880 * pow(n, 2) * h) * res1_fp;
  /*#######*/
  
  return(NumericVector::create(res,res_fp));
}



