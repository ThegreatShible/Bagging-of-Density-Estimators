#include <Rcpp.h>
using namespace Rcpp;


// Obtain environment containing function
Rcpp::Environment base("package:stats"); 

// Make function callable from C++

// [[Rcpp::export]]
SEXP histRcpp(NumericVector x, double breaks) {
  Function f("hist");
  return f(x, Named("plot")=false, Named("breaks")=breaks, Named("warn.unused") = false);
}

void counts(NumericVector x, int breaks) {
  NumericVector x_sorted = x.sort();
  print(table(x));
  
}

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
List hist_rcpp(NumericVector x, double nb_breaks = 0, NumericVector v_breaks= NumericVector(0)) {
  
  Function f("gc");
  x = x.sort();
  
  NumericVector breaks;
  if (nb_breaks > 0) {
    breaks = mybreaks_rcpp(x, nb_breaks);
  }
  else if (v_breaks.length() > 0){
    breaks = v_breaks;
    nb_breaks = breaks.length();
  }
  
  /*if ((breaks[nb_breaks -1] < x[x.length()-1])  || (x[0] < breaks[0])){
    Rcout << breaks[nb_breaks - 1];
    Rcout << x[x.length()-1];
    Rcout << breaks[0];
    Rcout << x[0];
    ::Rf_error("error"); 
  }*/
    

  /*NumericVector mids = breaks;

  mids.erase(mids.length()-1);

  mids = mids + (diff(breaks) / 2);*/
  
  
  NumericVector counts = NumericVector(nb_breaks - 1);
  int k=0;
  int j=1;
  for (int i=0; i<x.length(); i++) {
    if (x[i] >= breaks[j]) {
      counts[j-1] = k;
      k = 0;
      j = j + 1;
    }
    else {
      k = k + 1;
    }
  }
  counts[j-1] = k;
  NumericVector density = counts / sum(counts);
  density = density / diff(breaks);
  //Rcout << mids<< "\n";
  List res = List::create(
    Named("breaks")=breaks,
    Named("counts")=counts,
    Named("density")=density,
    Named("mids")=mids(breaks));
  
  return res;
}

