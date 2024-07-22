#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include <calobase/TowerInfov1.h>
#include <limits>
#include <string>
#include <vector>
#include <cmath> // Include this for fabs
//#include <omp.h>

R__LOAD_LIBRARY(libcaloTreeGen.so)


void ClusTrig(std::string outFile="./trigEmulator_run_49089_ClusTrig.root"){

  gROOT->SetBatch(kTRUE);
  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // TFile *out = new TFile(Form("%s", outFile.c_str()),"RECREATE");
  //  TFile *in = TFile::Open(Form("%s", inFile.c_str()));
  TFile *out = new TFile(Form("%s", outFile.c_str()), "RECREATE");

  // List of input files
  std::vector<std::string> fileList;
  for (int i = 1; i <= 40; ++i) {
    fileList.push_back(Form("/sphenix/tg/tg01/jets/jamesj3j3/output_48089_%d.root", i));
  }
  //gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  // TTree *ttree = (TTree*)in->Get("ttree");
  #pragma omp parallel for
  for (int i = 0; i < fileList.size(); ++i) {
    std::string fileName = fileList[i];
    TFile* in = TFile::Open(fileName.c_str());
    if (!in || in->IsZombie()) {
      std::cerr << "Failed to open file: " << fileName << std::endl;
      continue;
    }
    // Create a TChain to handle multiple input files
    TChain *ttree = new TChain("ttree");  // Replace "ttree" with the actual name of your tree

    // Add files to the chain
    for (int i = 1; i <= 40; ++i) {
      ttree->Add(Form("/sphenix/tg/tg01/jets/jamesj3j3/output_48089_%d.root", i));
    }  //T->Print();
  
  // Constants
  static const int n_hcal = 1536;
  static const int n_hcal_etabin = 24;
  static const int n_hcal_phibin = 64;
  static const int n_emcal = 24576;
  static const int n_emcal_etabin = 96;
  static const int n_emcal_phibin = 256;
  static const int n_hcaltrigger = 384;
  static const int n_hcaltrigger_etabin = 12;
  static const int n_hcaltrigger_phibin = 32;
  static const int n_emcaltrigger = 6144;
  static const int n_emcaltrigger_etabin = 48;
  static const int n_emcaltrigger_phibin = 128;
  static const int n_jettrigger = 288;
  static const int n_jettrigger_etabin = 9;
  static const int n_jettrigger_phibin = 32;
  const int jettrigger_min_etabin[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const int jettrigger_max_etabin[9] = {3, 4, 5, 6, 7, 8, 9, 10, 11};
  const int jettrigger_min_phibin[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
  const int jettrigger_max_phibin[32] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34}; // wrap around 0 and 31 in analysis part.
  const int n_trigger_vec = 64;
  
  //Tree Variables
  unsigned int b_trigger_sum_emcal[n_emcaltrigger_etabin][n_emcaltrigger_phibin];
  unsigned int b_trigger_sum_hcalin[n_hcaltrigger_etabin][n_hcaltrigger_phibin];
  unsigned int b_trigger_sum_hcalout[n_hcaltrigger_etabin][n_hcaltrigger_phibin];
  unsigned int b_trigger_sum_emcal_ll1[n_hcaltrigger_etabin][n_hcaltrigger_phibin];
  unsigned int b_trigger_sum_hcal_ll1[n_hcaltrigger_etabin][n_hcaltrigger_phibin];
  unsigned int b_trigger_sum_total_ll1[n_hcaltrigger_etabin][n_hcaltrigger_phibin];
  unsigned int b_trigger_sum_jet[n_jettrigger_etabin][n_jettrigger_phibin];
  float b_emcal_energy[n_emcal_etabin][n_emcal_phibin], b_emcal_time[n_emcal_etabin][n_emcal_phibin];
  float b_hcalin_energy[n_hcal_etabin][n_hcal_phibin], b_hcalin_time[n_hcal_etabin][n_hcal_phibin];
  float b_hcalout_energy[n_hcal_etabin][n_hcal_phibin], b_hcalout_time[n_hcal_etabin][n_hcal_phibin];

  
  //Cluster Info
  std::vector<float> *b_cluster_towere = 0;
  std::vector<float> *b_cluster_towereta = 0;
  std::vector<float> *b_cluster_towerphi = 0;

  //Vertex info
  float m_vertex{std::numeric_limits<float>::quiet_NaN()};

  //Trigger info
  //  std::vector<bool> *b_gl1_triggervec = 0;
  int b_gl1_triggervec[64];

  ttree -> SetCacheSize(100 * 1024 * 1024); // 100 MB cache size

  //  ttree -> SetBranchStatus("*", 0); // Disable all branches  
  ttree -> SetBranchAddress("emcal_tower_e", b_emcal_energy);
  ttree -> SetBranchAddress("ihcal_tower_e", b_hcalin_energy);
  ttree -> SetBranchAddress("ohcal_tower_e", b_hcalout_energy);
  ttree -> SetBranchAddress("emcal_tower_time", b_emcal_time);
  ttree -> SetBranchAddress("ihcal_tower_time", b_hcalin_time);
  ttree -> SetBranchAddress("ohcal_tower_time", b_hcalout_time);
  ttree -> SetBranchAddress("cluster_towere",&b_cluster_towere);
  ttree -> SetBranchAddress("cluster_towereta", &b_cluster_towereta); 
  ttree -> SetBranchAddress("cluster_towerphi", &b_cluster_towerphi);
  ttree -> SetBranchAddress("zvertex", &m_vertex);
  ttree -> SetBranchAddress("gl1_triggervec", &b_gl1_triggervec); //"gl1_triggervec[64]/I");
                  
  TCanvas *c1 = new TCanvas("Trigger Clusters 1","Trigger Clusters 1 Run 49089");
  TCanvas *c2 = new TCanvas("Trigger Clusters 2","Trigger Clusters 2 Run 49089");
  TCanvas *c3 = new TCanvas("Trigger Clusters 3","Trigger Clusters 3 Run 49089");
  TCanvas *c4 = new TCanvas("Trigger Clusters Phi 1","Trigger Clusters Phi 1 Run 49089");  
  TCanvas *c5 = new TCanvas("Trigger Clusters Phi 2","Trigger Clusters Phi 2 Run 49089");
  TCanvas *c6 = new TCanvas("Trigger Clusters Phi 3","Trigger Clusters Phi 3 Run 49089");
  TCanvas *c7 = new TCanvas("Trigger Clusters Eta 1","Trigger Clusters Eta 1 Run 49089");
  TCanvas *c8 = new TCanvas("Trigger Clusters Eta 2","Trigger Clusters Eta 2 Run 49089");
  TCanvas *c9 = new TCanvas("Trigger Clusters Eta 3","Trigger Clusters Eta 3 Run 49089");

// Defining the Histogram, setting number of bins for x and y axes                                                                                                                         
  //TH2D *ClusterPhotonTrigger26 = new TH2D("ClusterPhotonTrigger26","Tower Photon Trigger 26 Run 44224",50,0,90,50,0,250);                                                                  
  TH2D *TowerJetTrigger17 = new TH2D("TowerJetTrigger17","Tower Jet 2 (8 GeV) Trigger 17 Run 49089",96,0,96,256,0,256);
  TH2D *TowerJetTrigger18 = new TH2D("TowerJetTrigger18","Tower Jet 3 (10 GeV) Trigger 18 Run 49089",96,0,96,256,0,256);
  TH2D *TowerJetTrigger19 = new TH2D("TowerJetTrigger19","Tower Jet 4 (12 GeV) Trigger 19 Run 49089",96,0,96,256,0,256);
  TH1F *TowerJetTriggerPhi17 = new TH1F("TowerJetTriggerPhi17","Tower Jet 2 in Phi (8 GeV) Trigger 17 Run 49089",256,0,256);
  TH1F *TowerJetTriggerPhi18 = new TH1F("TowerJetTriggerPhi18","Tower Jet 3 in Phi (10 GeV) Trigger 18 Run 49089",256,0,256);
  TH1F *TowerJetTriggerPhi19 = new TH1F("TowerJetTriggerPhi19","Tower Jet 4 in Phi (12 GeV) Trigger 19 Run 49089",256,0,256);
  TH1F *TowerJetTriggerEta17 = new TH1F("TowerJetTriggerEta17","Tower Jet 2 in Eta (8 GeV) Trigger 17 Run 49089",96,0,96);
  TH1F *TowerJetTriggerEta18 = new TH1F("TowerJetTriggerEta18","Tower Jet 3 in Eta (10 GeV) Trigger 18 Run 49089",96,0,96);
  TH1F *TowerJetTriggerEta19 = new TH1F("TowerJetTriggerEta19","Tower Jet 4 in Eta (12 GeV) Trigger 19 Run 49089",96,0,96);
  
  int nEvents = ttree->GetEntries();
  //cout << "nEvents: " << nEvents << endl;
  for(int ev=0; ev<nEvents; ev++){

    //std::cout << "working on event " << ev << std::endl;
    ttree->GetEntry(ev);//Each event, this value is the same
    //cout << "This is the size of entries in the calo array" << " " << m_emcTowE->size() << "\n";

      /*Triggers*/
      /*
      3 -> ZDC Coincidence
      10 -> MBD N&S >=1 (scaled)
      16 -> Jet Threshold = 2 + MBD N&S >=1
      17 -> Jet Threshold = 3 + MBD N&S >=1
      18 -> Jet Threshold = 4 + MBD N&S >=1
      19 -> Jet Threshold = 5 + MBD N&S >=1
      25 -> Photon Threshold = 2 + MBD N&S >=1
      26 -> Photon Threshold = 3 + MBD N&S >=1
      27 -> Photon Threshold = 4 + MBD N&S >=1
      */
      

    /*Jet Triggers*/
   
    /*  if(m_triggerVector->at(17)){
      for(int i = 0; i < (*m_clusTowE).size(); i++){
	for(int j = 0; j < (*m_clusTowE)[0].size(); j++){
	  ClusterJetTrigger17->Fill((*m_clusTowEta)[i][j], (*m_clusTowPhi)[i][j], (*m_clusTowE)[i][j]); //Filling 2D cluster energy projection in eta and phi
	  //ClusterJetTrigger17->Fill(m_clusTowEta[1][i], m_clusTowPhi[1][i], m_clusTowE[1][i]); //Filling 2D cluster energy projection in eta and phi
	}
      }
    }
    */
    int nPhiBins17 = TowerJetTriggerPhi17->GetNbinsX();
    int nEtaBins17 = TowerJetTriggerEta17->GetNbinsX();

    // Print the number of bins
    //    std::cout << "Number of bins in eta 17: " << nEtaBins17 << std::endl;
    //  std::cout << "Number of bins in phi 17: " << nPhiBins17 << std::endl;

    const double epsilon = 1e-10; // Small value to prevent division by zero
    for (int i = 0; i < b_cluster_towere->size(); i++) {
      if (b_gl1_triggervec[17]) {
	double eta = b_cluster_towereta->at(i);
	double phi = b_cluster_towerphi->at(i);
	double energy = b_cluster_towere->at(i);

	TowerJetTrigger17->Fill(eta, phi, energy);
	if (eta > epsilon || eta < -epsilon) { // Ensure phi is not too small
	  TowerJetTriggerPhi17->Fill(phi, energy / nEtaBins17); // Weight by 1/eta
	}
	if (phi > epsilon || phi < -epsilon) { // Ensure eta is not too small
	  TowerJetTriggerEta17->Fill(eta, energy / nPhiBins17); // Weight by 1/phi
	}
      }
    }

    int nPhiBins18 = TowerJetTriggerPhi18->GetNbinsX();
    int nEtaBins18 = TowerJetTriggerEta18->GetNbinsX();
    // Print the number of bins
    //  std::cout << "Number of bins in eta 18: " << nEtaBins18 << std::endl;
    //  std::cout << "Number of bins in phi 18: " << nPhiBins18 << std::endl;
    for (int i = 0; i < b_cluster_towere->size(); i++) {
      if (b_gl1_triggervec[18]) {
        double eta = b_cluster_towereta->at(i);
        double phi = b_cluster_towerphi->at(i);
        double energy = b_cluster_towere->at(i);

        TowerJetTrigger18->Fill(eta, phi, energy);
        if (eta > epsilon || eta < -epsilon) { // Ensure phi is not too 
          TowerJetTriggerPhi18->Fill(phi, energy /nEtaBins18); // Weight by 1/
        }
        if (phi > epsilon || phi < -epsilon) { // Ensure eta is not too 
          TowerJetTriggerEta18->Fill(eta, energy / nPhiBins18); // Weight by 1/phi                                                                                                           
        }
      }
    }
    int nPhiBins19 = TowerJetTriggerPhi19->GetNbinsX();
    int nEtaBins19 = TowerJetTriggerEta19->GetNbinsX();
    // Print the number of bins
    //  std::cout << "Number of bins in eta: " << nEtaBins19 << std::endl;
    //  std::cout << "Number of bins in phi: " << nPhiBins19 << std::endl;
    for (int i = 0; i < b_cluster_towere->size(); i++) {
      if (b_gl1_triggervec[19]) {
        double eta = b_cluster_towereta->at(i);
        double phi = b_cluster_towerphi->at(i);
        double energy = b_cluster_towere->at(i);

        TowerJetTrigger19->Fill(eta, phi, energy);
        if (eta > epsilon || eta < -epsilon) { // Ensure phi is not too small                                           
          TowerJetTriggerPhi19->Fill(phi, energy / nEtaBins19); // Weight by 1/eta
        }
        if (phi > epsilon || phi < -epsilon) { // Ensure eta is not too small                                                                    
          TowerJetTriggerEta19->Fill(eta, energy / nPhiBins19); // Weight by 1/phi
        }
      }
    }
    /*
    for (int i = 0; i < b_cluster_towere->size(); i++) {
      if (b_gl1_triggervec[19]) {
        TowerJetTrigger19->Fill(b_cluster_towereta->at(i), b_cluster_towerphi->at(i), b_cluster_towere->at(i));
	TowerJetTriggerPhi19->Fill(b_cluster_towerphi->at(i), b_cluster_towere->at(i)/b_cluster_towereta->at(i));
        TowerJetTriggerEta19->Fill(b_cluster_towereta->at(i), b_cluster_towere->at(i)/b_cluster_towerphi->at(i));
      }
    }
    */
      
    b_cluster_towere->clear();
    b_cluster_towereta->clear();
    b_cluster_towerphi->clear();
    for (int i = 0; i < n_jettrigger_etabin; ++i) {
      for (int j = 0; j < n_jettrigger_phibin; ++j) {
	b_trigger_sum_jet[i][j] = 0;
      }
    }
    //    b_gl1_triggervec->clear();
    //m_vertex.clear();
    for (int ieta = 0; ieta < n_emcaltrigger_etabin; ++ieta) {
      for (int iphi = 0; iphi < n_emcaltrigger_phibin; ++iphi) {
	b_trigger_sum_emcal[ieta][iphi] = 0;
      }
    }
    for (int ieta = 0; ieta < n_hcaltrigger_etabin; ++ieta) {
      for (int iphi = 0; iphi < n_hcaltrigger_phibin; ++iphi) {
	b_trigger_sum_hcalin[ieta][iphi] = 0;
	b_trigger_sum_hcalout[ieta][iphi] = 0;
	b_trigger_sum_emcal_ll1[ieta][iphi] = 0;
	b_trigger_sum_hcal_ll1[ieta][iphi] = 0;
      }
    }
  }//end of event loop

  // out->cd();
  // Drawing the histogram and labelling axes
  
 
  c1->cd();
  c1->SetLogz();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);
  
  TowerJetTrigger17->Draw("COLZ");
  TowerJetTrigger17->GetXaxis()->SetTitle("#eta");
  TowerJetTrigger17->GetYaxis()->SetTitle("#phi");
  TowerJetTrigger17->GetZaxis()->SetTitle("Energy");
  
  c1->SaveAs("TowerJetTrigger17_49089.png");

  c2->cd();
  c2->SetLogz();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTrigger18->Draw("COLZ");
  TowerJetTrigger18->GetXaxis()->SetTitle("#eta");
  TowerJetTrigger18->GetYaxis()->SetTitle("#phi");
  TowerJetTrigger18->GetZaxis()->SetTitle("Energy");

  c2->SaveAs("TowerJetTrigger18_49089.png");

  c3->cd();
  c3->SetLogz();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);
   
  TowerJetTrigger19->Draw("COLZ");
  TowerJetTrigger19->GetXaxis()->SetTitle("#eta");
  TowerJetTrigger19->GetYaxis()->SetTitle("#phi");
  TowerJetTrigger19->GetZaxis()->SetTitle("Energy");
  
  c3->SaveAs("TowerJetTrigger19_49089.png");
  
  c4->cd();
  c4->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerPhi17->Draw("COLZ");
  TowerJetTriggerPhi17->GetXaxis()->SetTitle("#phi");
  TowerJetTriggerPhi17->GetYaxis()->SetTitle("Energy");
  
  // Save the canvas
  c4->SaveAs("TowerJetTriggerPhi17_49089w.png");
  
  c5->cd();
  c5->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerPhi18->Draw("COLZ");
  TowerJetTriggerPhi18->GetXaxis()->SetTitle("#phi");
  TowerJetTriggerPhi18->GetYaxis()->SetTitle("Energy");
  c5->SaveAs("TowerJetTriggerPhi18_49089w.png");
  
  c6->cd();
  c6->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerPhi19->Draw("COLZ");
  TowerJetTriggerPhi19->GetXaxis()->SetTitle("#phi");
  TowerJetTriggerPhi19->GetYaxis()->SetTitle("Energy");

  c6->SaveAs("TowerJetTriggerPhi19_49089w.png");

  c7->cd();
  c7->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerEta17->Draw("COLZ");
  TowerJetTriggerEta17->GetXaxis()->SetTitle("#eta");
  TowerJetTriggerEta17->GetYaxis()->SetTitle("Energy");

  c7->SaveAs("TowerJetTriggerEta17_49089w.png");
  
  c8->cd();
  c8->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerEta18->Draw("COLZ");
  TowerJetTriggerEta18->GetXaxis()->SetTitle("#eta");
  TowerJetTriggerEta18->GetYaxis()->SetTitle("Energy");
  c8->SaveAs("TowerJetTriggerEta18_49089w.png");
  
  c9->cd();
  c9->SetLogy();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

  TowerJetTriggerEta19->Draw("COLZ");
  TowerJetTriggerEta19->GetXaxis()->SetTitle("#eta");
  TowerJetTriggerEta19->GetYaxis()->SetTitle("Energy");

  c9->SaveAs("TowerJetTriggerEta19_49089w.png");
  

 // For creating text label for the Chi2 parameter: 
  // First draw the histogram, click on View and go to Toolbar, then click on TextPave or PaveText and click and drag to create a box on the canvas.
  // Then, left click to add text, edit borders and fill color of box
  // For this histogram, I chose text size 20, border size 0, fill color white, text font times, and text alignment 22 (middle, middle)
  out->Close();
  delete out;
  delete ttree;  
  }
  
}
