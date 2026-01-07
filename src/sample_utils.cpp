#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector disc_cdf_from_quantiles(const NumericVector& quantiles,
                             const NumericVector& values,
                             const int k_max,
                             const int k_min) {

     int m = quantiles.size();

     int K = k_max - k_min + 1;

     NumericVector cdf(K, 0.0);   // discrete CDF F(k)

     // Build the CDF from quantile inversion on each linear segment
     for (int i = 0; i < m - 1; ++i) {
          double pL = quantiles[i];
          double pR = quantiles[i + 1];
          double vL = values[i];
          double vR = values[i + 1];

          if (vR <= vL) continue; // skip flat segments

          // integers in this quantile segment
          int k_start = std::max(k_min, static_cast<int>(std::ceil(vL)));
          int k_end   = std::min(k_max, static_cast<int>(std::floor(vR)));

          if (k_end < k_start) continue;

          double dv = vR - vL;
          double dp = pR - pL;

          // invert quantile function: find p such that Q(p) = k
          for (int k = k_start; k <= k_end; ++k) {
               double t = (static_cast<double>(k) - vL) / dv;
               double pk = pL + t * dp;
               int idx = k - k_min;

               // take sup over segments
               if (pk > cdf[idx])
                    cdf[idx] = pk;
          }
     }

     // Enforce monotonicity upwards
     double last = 0.0;
     for (int i = 0; i < K; ++i) {
          if (cdf[i] < last)
               cdf[i] = last;
          else
               last = cdf[i];
     }


     return cdf;
}

// [[Rcpp::export]]
IntegerVector rdisc_from_cdf(const IntegerVector& support,
                             const NumericVector& cdf,
                             double total_mass,
                             int n) {

  int K = cdf.size();
  IntegerVector out(n);

  for (int i = 0; i < n; ++i) {
    double u = R::runif(0.0, total_mass);

    int left = 0, right = K - 1;
    while (left < right) {
      int mid = (left + right) >> 1;
      if (cdf[mid] >= u)
        right = mid;
      else
        left = mid + 1;
    }
    out[i] = support[left];
  }

  return out;
}


// [[Rcpp::export]]
IntegerMatrix rdisc_from_quantile_matrix(const NumericMatrix& value_mat,
                                         const NumericVector& quantiles,
                                         const int n_samples) {
  const int n_days = value_mat.nrow();
  const int m      = value_mat.ncol(); // Number of quantiles 

  const int k_min = 0;
  IntegerMatrix out(n_days, n_samples);

  // raw pointer to quantiles (these *are* contiguous)
  const double* qptr = &quantiles[0];

  for (int d = 0; d < n_days; ++d) {

    // ---- k_max from FLOOR(max value) ----
    double v_max = value_mat(d, m - 1);   // value at highest quantile (0.999)
    int k_max = static_cast<int>(std::floor(v_max));

    if (k_max < k_min) {
      // degenerately all zeros for that day
      for (int j = 0; j < n_samples; ++j)
        out(d, j) = 0;
      continue;
    }

    int K = k_max - k_min + 1;

    // ---- CDF for this day ----
    NumericVector cdf(K);
    double* cptr = &cdf[0];

    // ---- Build the CDF ----
    for (int i = 0; i < m - 1; ++i) {

      double pL = qptr[i];
      double pR = qptr[i + 1];
      double vL = value_mat(d, i);
      double vR = value_mat(d, i + 1);

      if (vR <= vL) continue;

      int k_start = std::max(k_min, static_cast<int>(std::ceil(vL)));
      int k_end   = std::min(k_max, static_cast<int>(std::floor(vR)));
      if (k_end < k_start) continue;

      double dv = vR - vL;
      double dp = pR - pL;

      for (int k = k_start; k <= k_end; ++k) {
        double t  = (static_cast<double>(k) - vL) / dv;
        double pk = pL + t * dp;

        int idx = k - k_min;      // here idx == k since k_min = 0, but explicit

        if (pk > cptr[idx])
          cptr[idx] = pk;
      }
    }

    // ---- Enforce monotonicity ----
    double last = 0.0;
    for (int i = 0; i < K; ++i) {
      double ci = cptr[i];
      if (ci < last)
        cptr[i] = last;
      else
        last = ci;
    }

    double total_mass = cptr[K - 1];
    if (total_mass <= 0.0) {
      // fall back: all zeros for this day
      for (int j = 0; j < n_samples; ++j)
        out(d, j) = 0;
      continue;
    }

    // ---- Sampling using std::lower_bound ----
    for (int j = 0; j < n_samples; ++j) {
      double u = R::runif(0.0, total_mass);

      auto it = std::lower_bound(cptr, cptr + K, u);
      out(d, j) = static_cast<int>(it - cptr) + k_min;
    }
  }

  return out;
}

