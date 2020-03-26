
TChain cA("toy"); 
TChain cB("toy"); 
TCut wA;
TCut wB;


void setAStyle(TH1 * h)
{
  h->SetLineWidth(3); 
  h->SetLineColor(46); 
  h->SetMarkerColor(46); 
}

void setBStyle(TH1 * h)
{
  h->SetLineWidth(3); 
  h->SetLineColor(38); 
  h->SetMarkerColor(38); 
}


void makePlots(double R = 15000, double h = 3000)
{
  // doh, I should have saved this lol 
  double volume = TMath::Pi() * h *R *R; 
  double area = 4*R*R+h*h; 
  double N = 1e8; 
  double factorA = volume/N; 
  double factorB = area/N; 

  wB = Form("pass*absorption_weight*interaction_weight * %f", factorB);
  wA = Form("pass*absorption_weight/lint * %f", factorA);

  cA.Add("toydisc_1*.root"); 
  cB.Add("toydisc_0*.root"); 
  cA.SetLineColor(2); 
  cB.SetLineColor(3); 
  cA.SetMarkerColor(2); 
  cB.SetMarkerColor(3); 
  TH1D * Aeff_A = new TH1D("Aeff_A", "Method A; log10(E [eV]); <A_{eff}> |_{4#pi} [m^{2}]", 5, 15.75, 18.25);
  TH1D * Aeff_B = new TH1D("Aeff_B", "Method B; log10(E [eV]); <A_{eff}> |_{4#pi} [m^{2}]", 5, 15.75, 18.25);

  setAStyle(Aeff_A);
  setBStyle(Aeff_B);


  cA.Draw("log10(E) >> Aeff_A",wA,"goff"); 
  cB.Draw("log10(E) >> Aeff_B",wB,"goff"); 


  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 
  TCanvas * c = new TCanvas; 
  c->Divide(1,2); 
  c->cd(1); 
  Aeff_A->Draw();
  Aeff_B->Draw("same");
  gPad->SetLogy(); 
  gPad->SetGridy(); 
  gPad->BuildLegend(0.1,0.7,0.3,0.9); 
  c->cd(2); 

  TH1D * ratio = (TH1D*) Aeff_A->Clone("ratio"); 
  ratio->SetTitle("Ratio; log10(E[eV]); A/B"); 
  ratio->Divide(Aeff_B); 
  ratio->Draw(); 

  TCanvas * c2 = new TCanvas; 

  TH1D * theta18_A = new TH1D("theta18A", "Method A (1e18); cos #theta_{#nu};  d <A_{eff}> |_{4#pi} / d cos #theta_{#nu} [m^{2}]", 60,-1,1);
  TH1D * theta18_B = new TH1D("theta18B", "Method B (1e18); cos #theta_{#nu};  d <A_{eff}> |_{4#pi} / d cos #theta_{#nu} [m^{2}]", 60,-1,1);
  TH1D * theta17_A = new TH1D("theta17A", "Method A (1e17); cos #theta_{#nu};  d <A_{eff}> |_{4#pi} / d cos #theta_{#nu} [m^{2}]", 60,-1,1);
  TH1D * theta17_B = new TH1D("theta17B", "Method B (1e17); cos #theta_{#nu};  d <A_{eff}> |_{4#pi} / d cos #theta_{#nu} [m^{2}]", 60,-1,1);

  setAStyle(theta18_A); 
  setAStyle(theta17_A); 
  setBStyle(theta18_B); 
  setBStyle(theta17_B); 
  theta17_A->SetLineColor(42);
  theta17_B->SetLineColor(33);


  cA.Draw("cos(direction.Theta()) >> theta18A", wA * TCut("E==1e18"),"goff"); 
  cA.Draw("cos(direction.Theta()) >> theta17A", wA * TCut("E==1e17"),"goff"); 
  cB.Draw("cos(direction.Theta()) >> theta18B", wB * TCut("E==1e18"),"goff"); 
  cB.Draw("cos(direction.Theta()) >> theta17B", wB * TCut("E==1e17"),"goff"); 

  theta18_A->Draw(); 
  theta17_A->Draw("same"); 
  theta18_B->Draw("same"); 
  theta17_B->Draw("same"); 
  gPad->BuildLegend(0.7,0.7,0.9,0.9); 

}
