#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace Rcpp;

// Harrell's C for every row of a score matrix, a larger score going with longer survival, with the
// conventions of survival::concordance(Surv(time, event) ~ x): a pair is comparable when the shorter
// time is an event, a censored time equal to an event time counting as the longer; two events at the
// same time are not comparable; a tie in the score counts one half. O(n log n) per row (Fenwick tree).

// [[Rcpp::export]]
NumericVector cindex_rows_cpp(NumericMatrix S, NumericVector time, IntegerVector event) {
  const int n = S.ncol(), m = S.nrow();
  if (time.size() != n || event.size() != n) stop("time and event must have one entry per column of S");
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  // time descending; at equal times censored before events, so that events see them as longer
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    if (time[a] != time[b]) return time[a] > time[b];
    return event[a] < event[b];
  });
  NumericVector out(m);
  std::vector<double> x(n);
  std::vector<int> byx(n), rk(n);
  std::vector<double> bit(n + 1);
  for (int r = 0; r < m; ++r) {
    bool bad = false;
    for (int i = 0; i < n; ++i) { x[i] = S(r, i); if (!std::isfinite(x[i])) bad = true; }
    if (bad) { out[r] = NA_REAL; continue; }
    std::iota(byx.begin(), byx.end(), 0);
    std::sort(byx.begin(), byx.end(), [&](int a, int b) { return x[a] < x[b]; });
    int u = 0;
    for (int k = 0; k < n; ++k) {
      if (k == 0 || x[byx[k]] != x[byx[k - 1]]) ++u;
      rk[byx[k]] = u;                       // dense rank, 1..u
    }
    std::fill(bit.begin(), bit.begin() + u + 1, 0.0);
    auto add = [&](int i) { for (; i <= u; i += i & (-i)) bit[i] += 1.0; };
    auto pre = [&](int i) { double s = 0; for (; i > 0; i -= i & (-i)) s += bit[i]; return s; };
    double conc = 0, disc = 0, tie = 0, ins = 0;
    int pos = 0;
    while (pos < n) {
      int end = pos;
      while (end < n && time[ord[end]] == time[ord[pos]]) ++end;
      for (int p = pos; p < end; ++p) if (event[ord[p]] == 0) { add(rk[ord[p]]); ins += 1; }
      for (int p = pos; p < end; ++p) if (event[ord[p]] == 1) {
        int ri = rk[ord[p]];
        double less = pre(ri - 1), leq = pre(ri);
        conc += ins - leq; disc += less; tie += leq - less;
      }
      for (int p = pos; p < end; ++p) if (event[ord[p]] == 1) { add(rk[ord[p]]); ins += 1; }
      pos = end;
    }
    double den = conc + disc + tie;
    out[r] = den > 0 ? (conc + 0.5 * tie) / den : NA_REAL;
  }
  return out;
}
