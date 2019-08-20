// NPlotOverlay.C
// A macro that takes N root files that each store a TH1F and overlays the plots.
// USAGE:
// .x NPlotOverlay.C("file.root",

#include "overlap.h"
#include <algorithm>
#include <stdio.h>

void NPlotOverlay(const char *file, const char *title, const char *variable,
                  Double_t yheight, bool logscale = false) {

  // Set styling.
  gStyle->SetOptStat(0); // 0: disable legend 1: enable legend

  Int_t palette [6] = {3, 4, 6, 7, 2, 433}; // 46, 9 softer palette
  Int_t stipple [6] = {2, 3, 4, 5, 9, 1};
  Int_t style_i = 0;
  Int_t width = 2;
  Double_t font_size = 0.04;

  TCanvas *c = new TCanvas("c", variable, 800, 700);

  if (logscale) c->SetLogy();

  TH1F *h1 = new TH1F("", "", 20, 0., 1.);
  const char *first_hist_name = "";

  TFile *f = new TFile(file);
  TIter next(f->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    cout << "Found object: " << key->GetClassName() << endl;
    if (!cl->InheritsFrom("TH1F")) continue;
    h1 = (TH1F*)key->ReadObj();
    std::cout << key->GetTitle() << endl;
    first_hist_name = key->GetName();
    break;
  }

  h1->SetTitle(title);
  h1->GetXaxis()->SetTitleSize(font_size);
  h1->GetXaxis()->SetTitle(variable);
  h1->GetYaxis()->SetTitleSize(font_size);
  h1->GetYaxis()->SetTitle("a.u.");
  h1->SetMaximum(yheight);

  h1->SetLineColorAlpha(palette[0], 1);
  h1->SetLineWidth(width);
  h1->SetLineStyle(stipple[0]);
  h1->SetFillColorAlpha(palette[0], 0.0);

  h1->Draw("HIST");

  TLegend *legend = new TLegend(0.6, 0.75, 0.9, 0.90);
  legend->SetBorderSize(0.);
  legend->SetTextSize(font_size); // % of pad size
  legend->AddEntry(h1, first_hist_name, "f");

  TIter next2(f->GetListOfKeys());
  TKey *key2;
  while ((key2 = (TKey*)next2())) {
    TClass *cl2 = gROOT->GetClass(key2->GetClassName());
    cout << "Found object: " << key2->GetClassName() << endl;
    if (!cl2->InheritsFrom("TH1F")) continue;
    if (key2->GetName() == first_hist_name) continue;
    TH1F *this_hist = new TH1F("", "", 10, 0., 1.);

    this_hist = (TH1F*)key2->ReadObj();
    Double_t o = overlap(this_hist, h1);
    printf("Overlap with h1 is %f \n", o);
    style_i++;
    this_hist->SetLineColorAlpha(palette[style_i], 1);
    this_hist->SetFillColorAlpha(palette[style_i], 0.0);
    this_hist->SetLineWidth(width);
    this_hist->SetLineStyle(stipple[style_i]);

    this_hist->Draw("HIST SAME");

    legend->AddEntry(this_hist, key2->GetName(), "f");
  }

  legend->Draw();

}
