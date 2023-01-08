#include <Rcpp.h>
using namespace Rcpp;

//' Function to sort an array using shellSort
//'  
//' @export
// [[Rcpp::export]]
NumericVector shellSort_c(NumericVector x) { 
  int n = x.size();
  NumericVector y = clone(x);
  
  // this shouldn't be necessary - clone(x) should do this but hey ho
  //for(int i=0; i<n; i++) y[i] = x[i];
  
  // Start with a big gap, then reduce the gap 
  for (int gap = n/2; gap > 0; gap /= 2) { 
    // Do a gapped insertion sort for this gap size. 
    // The first gap elements y[0..gap-1] are already in gapped order 
    // keep adding one more element until the entire array is 
    // gap sorted  
    for (int i = gap; i < n; i += 1) { 
      // add y[i] to the elements that have been gap sorted 
      // save y[i] in temp and make a hole at position i 
      double temp = y[i]; 
      
      // shift earlier gap-sorted elements up until the correct  
      // location for y[i] is found 
      int j;             
      for (j = i; j >= gap && y[j - gap] > temp; j -= gap) {
        y[j] = y[j - gap]; 
      }
      
      //  put temp (the original y[i]) in its correct location 
      y[j] = temp; 
    } 
  } 
  return y; 
} 

//' Function to return the sorting order of an array
//'  
//' @export
// [[Rcpp::export]]
IntegerVector order_c(NumericVector x) { 
  int n = x.size();
  NumericVector y = clone(x);
  IntegerVector result(n);

  for(int i=0; i<n; i++) result[i] = i+1;

  // Start with a big gap, then reduce the gap 
  for (int gap = n/2; gap > 0; gap /= 2) { 
    // Do a gapped insertion sort for this gap size. 
    // The first gap elements y[0..gap-1] are already in gapped order 
    // keep adding one more element until the entire array is 
    // gap sorted  
    for (int i = gap; i < n; i += 1) { 
      // add y[i] to the elements that have been gap sorted 
      // save y[i] in temp and make a hole at position i 
      double temp = y[i]; 
      int itemp = result[i];
      
      // shift earlier gap-sorted elements up until the correct  
      // location for y[i] is found 
      int j;             
      for (j = i; j >= gap && y[j - gap] > temp; j -= gap) {
        y[j] = y[j - gap];
        result[j] = result[j - gap];
      }
      
      //  put temp (the original y[i]) in its correct location 
      y[j] = temp; 
      result[j] = itemp;
    } 
  } 
  return result; 
} 

