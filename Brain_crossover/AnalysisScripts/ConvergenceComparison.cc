#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

void ConvergenceComparison()
{
  TFile *file_mutation = new TFile("c_dtime_generation_mutation_only.root");
  TFile *file_mutation_crossover = new TFile("c_dtime_generation_mutation_and_crossover.root");
  
  TH1F *h_dtime_generation_mutation = (TH1F*)((TCanvas*)file_mutation->Get("c_dtime_generation"))->GetPrimitive("h_dtime_generation");
  TH1F *h_dtime_generation_mutation_crossover = (TH1F*)((TCanvas*)file_mutation_crossover->Get("c_dtime_generation"))->GetPrimitive("h_dtime_generation");
  
  int rebin = 1;
  
  h_dtime_generation_mutation->Rebin(rebin);
  h_dtime_generation_mutation->Scale(1./rebin);
  h_dtime_generation_mutation->SetLineColor(kBlack);
  
  h_dtime_generation_mutation_crossover->Rebin(rebin);
  h_dtime_generation_mutation_crossover->Scale(1./rebin);
  h_dtime_generation_mutation_crossover->SetLineColor(kRed);
  
  h_dtime_generation_mutation->SetTitle(";generations; Average time to next generation");
  h_dtime_generation_mutation->GetYaxis()->SetTitleOffset(1.2);
  
  TCanvas *c_dtime_generation_comparison = new TCanvas("c_dtime_generation_comparison", "c_dtime_generation_comparison", 1000, 700);
  h_dtime_generation_mutation->SetMinimum(300);
  h_dtime_generation_mutation->GetXaxis()->SetRangeUser(100, 50000);
  c_dtime_generation_comparison->SetLogx();
  c_dtime_generation_comparison->SetLogy();
  h_dtime_generation_mutation->Draw("HIST E");
  h_dtime_generation_mutation_crossover->Draw("HIST E SAME");
  TLegend *leg = new TLegend(0.5, 0.75, 0.89, 0.89);
  leg->AddEntry(h_dtime_generation_mutation, "Mutation");
  leg->AddEntry(h_dtime_generation_mutation_crossover, "Crossover & mutation");
  leg->SetLineColor(0); leg->SetFillColor(0);
  leg->Draw();
  c_dtime_generation_comparison->SaveAs("c_dtime_generation_comparison.png");
}
