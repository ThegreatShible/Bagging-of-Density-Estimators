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
double riskhistRcpp(NumericVector obs, double m, NumericVector xlim) {
  NumericVector obs01  = (obs - xlim[1]) / (xlim[2] - xlim[1]);
  double h = 1 / m;
  double n = obs.length();
  NumericVector breaks = NumericVector(m+1);
  double s = 0;
  double interval = 1/(m+1);
  for (int i=1; i<=m+1; i++) {
    breaks[i] = s;
    s = s+interval;
  }
  //histogram p_hat = (hist_counts(obs01, Named("plot") = false, Named("breaks") = breaks, Named("warn.unused") = false));
  double p_hat = 3;
  double res = 2 / h / (n - 1) - (n + 1) / (n - 1) / h * pow(p_hat, 2);
  return res; /*(m * sum(p_hat^2)) */
}

// [[Rcpp::export]]
Rcpp::List broptRcpp(NumericVector x){
  int len = 5*floor(sqrt(x.length()));
  NumericVector Mgrid = NumericVector(len-1);
  for (int i=2; i<=len; i++) {
    Mgrid[i-1] = i;
  }
  NumericVector J = NumericVector(Mgrid.length());
  /*NumericVector lim = NumericVector::create(min(x)-0.5, max(x)+0.5);*/
  NumericVector lim = range(x);
  lim[1] -= 0.5;
  lim[2] += 0.5;
  for(int m=1; m<=Mgrid.length(); m++) {
    Function riskhist("riskhist");
    J[m] = (riskhistRcpp(x, Mgrid[m], lim));
  }
  return Rcpp::List::create(Rcpp::Named("opt")=Mgrid[which_min(J)]);
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
  
  Function f("gc");
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

