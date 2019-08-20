// A function that takes in two TH1's and returns the computation
// of their overlap.
#include <algorithm>
#include <stdio.h>
#include <math.h>

Double_t overlap(TH1 *h1, TH1 *h2) {

  // Assumes the same bin count for both.
  Int_t length = h1->GetXaxis()->GetNbins();

  Double_t overlap = 0.0;

  for (Int_t i = 0; i < length; i++) {
    Double_t ith_bin_1 = h1->GetBinContent(i);
    Double_t ith_bin_2 = h2->GetBinContent(i);
    Double_t product_i = sqrt(ith_bin_1 * ith_bin_2);
    overlap += product_i;
  }

  return overlap;
}
