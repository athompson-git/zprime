// TopoPlotter - Adrian Thompson
// Last Updated April 29, 2018
//
// This is a macro that makes a topological plot of an event structure
// with N particles.
//
// It first reads a set of TH2's stored in a .root file that each hold
// 2D topological plots for a certain physics object, e.g the leading
// BTagged jet or the MissingET.
//
// Each TH2 should plot the particle ( pT * cos(dPhi) , pT * sin(dPhi) )
// where pT is the particle's transverse momentum and dPhi is the
// particle's azimuthal separation from a fixed object in the event;
// in my case I took everything relative to the leading hadronic Tau.
//
// From each TH2 we then draw the mean (x,y) as points and draw the
// standard deviations (x,y) as error ellipses in this topological
// graph.

{
   ///////////////////////////////////////////////////////////
   gROOT->Reset();
   c1 = new TCanvas("c1","gerrors2",200,10,700,500);
   c1->SetGrid();
   gStyle->SetOptStat(0); // 0: disable legend 1: enable legend

   // Set format parameters.
   Double_t kTextSize = 0.03;
   Int_t kMarkerStyle = 3;
   Int_t kLineWidth = 2;
   Double_t kScale = 0.1;

   // Create a 2D histogram to define the range
   TH2F *hist = new TH2F("hr","Z' (500 GeV) Ditau Topology",2,-150.,230.,2,-10.,10.);
   hist->SetXTitle("p_{T} cos(#Delta #phi_{#tau_{1}})");
   hist->SetYTitle("p_{T} sin(#Delta #phi_{#tau_{1}})");
   hist->Draw();
   c1->GetFrame()->SetBorderSize(12);

   // Prepare four arrays.
   vector<Double_t> x;
   vector<Double_t> y;
   vector<Double_t> s_x;
   vector<Double_t> s_y;
   vector<string> label;

   // 5 Particles:
   // Sample 1 (Dark): kRed-4, kBlue-4, kViolet-7, kOrange-7, kPink-7
   // Sample 2 (Light): kRed-7, kBlue-7, kViolet-4, kOrange-4, kPink-4
   Int_t palette [10] = {628, 596, 873, 793, 893, 625, 593, 876, 796, 896};


   ///////////////////////////////////////////////////////////
   // Read in the root file that contains all the plots.
   TFile *f = new TFile("hist_2D_event_topology.root");

   // Get the list of keys in the file.
   TKey *key;
   TIter next(f->GetListOfKeys());

   // Loop over each histogram in the file and collect the mean(x,y),
   // stddev(x,y) into organized arrays.
   while ((key = (TKey*)next())) {
     // Data quality check.
     TClass *cl = gROOT->GetClass(key->GetClassName());
     if (!cl->InheritsFrom("TH2")) continue;

     // Declare a dummy histo to read data from each TH2.
     TH2F *this_hist = new TH2F("", "", 2, -300.,300.,2,-300.,300.);
     this_hist = (TH2F*)key->ReadObj();

     // Read in each mean(x,y) and stddev(x,y) into 4 separate arrays.
     x.push_back(this_hist->GetMean(1));
     y.push_back(this_hist->GetMean(2));
     s_x.push_back(this_hist->GetStdDev(1));
     s_y.push_back(this_hist->GetStdDev(2));
     label.push_back(this_hist->GetName());
   }


   ///////////////////////////////////////////////////////////
   // Start by just drawing the FIRST point, line, and ellipse.
   TGraph *first = new TGraph();
   first->SetMarkerColor(palette[0]);
   first->SetMarkerStyle(kMarkerStyle);
   first->SetPoint(0,x[0],y[0]);

   // Draw lines to each point colored by sample type.
   TLine *line = new TLine(0, 0, x[0], y[0]);
   line->SetLineStyle(9);
   line->SetLineWidth(kLineWidth);
   line->SetLineColorAlpha(kRed,0.2);
   line->Draw();

   // Label the first point.
   first->GetPoint(0,x[0],y[0]);
   TLatex *latex = new TLatex(x[0]+1, y[0]+1, label[0].c_str());
   latex->SetTextSize(kTextSize);
   first->GetListOfFunctions()->Add(latex);

   // Draw the error ellipse for the first point.
   TEllipse *ell_first = new TEllipse(x[0], y[0], (2*s_x[0]*kScale),
                                      (s_y[0]*kScale));
   ell_first->SetFillColorAlpha(palette[0], 0.1);
   ell_first->Draw();

   first->Draw("P");

   TLegend *legend = new TLegend(0.91, 0.6, 1.1, 0.90);
   legend->SetBorderSize(0.);
   legend->SetTextSize(kTextSize); // % of pad size
   legend->AddEntry(first, label[0].c_str(), "p");


   ///////////////////////////////////////////////////////////
   // Now draw the rest of the points, lines, and ellipses.
   for (int i = 1; i<x.size(); i++) {
     Double_t xp = x[i];
     Double_t yp = y[i];

     TGraph *graph = new TGraph();
     graph->SetMarkerColor(palette[i]);
     graph->SetMarkerStyle(kMarkerStyle);
     graph->SetPoint(0,x[i],y[i]);

     // Label each marker.
     graph->GetPoint(0,x[i],y[i]);
     TLatex *l = new TLatex(x[i]+0.2*x[i], y[i]+0.2*y[i], label[i].c_str());
     l->SetTextSize(kTextSize);
     graph->GetListOfFunctions()->Add(l);

     // Draw lines to each particle colored by sample type.
     TLine *line = new TLine(0, 0, x[i], y[i]);
     if (i<5) {
       line->SetLineStyle(9);
       line->SetLineColorAlpha(kRed,0.5);
     } else {
       line->SetLineStyle(4);
       line->SetLineColorAlpha(kGreen,0.5);
     }
     line->SetLineWidth(kLineWidth);
     line->Draw();

     // Draw the error ellipses.
     TEllipse *ell = new TEllipse(x[i],y[i],(2*s_x[i]*kScale),(s_y[i]*kScale));
     ell->SetFillColorAlpha(palette[i], 0.1);
     ell->Draw();

     graph->Draw("P SAME");

     legend->AddEntry(graph, label[i].c_str(), "p");
   }


   legend->Draw();


}  // End Macro
