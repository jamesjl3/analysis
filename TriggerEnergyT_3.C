#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <TH2F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>
#include <TH2F.h>

void saveHistogramAsPNG(TH2F *histogram, const char *filename) {
  TCanvas *canvas = new TCanvas("canvas", "Canvas");
  histogram->Draw("COLZ");
  canvas->SaveAs(filename);
  delete canvas;
}

void TriggerEnergyT_3() {
  TChain chain_jet04("ttree");
  chain_jet04.Add("../macros/output_47089.root");

  // Constants
  static const int n_hcal_etabin = 24;
  static const int n_hcal_phibin = 64;
  static const int n_emcal_etabin = 96;
  static const int n_emcal_phibin = 256;
  static const int n_hcaltrigger_etabin = 12;
  static const int n_hcaltrigger_phibin = 32;
  static const int n_emcaltrigger_etabin = 48;
  static const int n_emcaltrigger_phibin = 128;
  static const int n_jettrigger_etabin = 9;
  static const int n_jettrigger_phibin = 32;

  const int jettrigger_min_etabin[9] = {0, 2, 4, 6, 8, 10, 12, 14, 16};
  const int jettrigger_max_etabin[9] = {7, 9, 11, 13, 15, 17, 19, 21, 23};
  const int jettrigger_min_phibin[32] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62};
  const int jettrigger_max_phibin[32] = {7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69};

  const float eta_map[n_hcal_etabin] = {-1.05417, -0.9625, -0.870833, -0.779167, -0.6875, -0.595833, -0.504167, -0.4125, -0.320833, -0.229167, -0.1375, -0.0458333, 0.0458333, 0.1375, 0.229167, 0.320833, 0.4125, 0.504167, 0.595833, 0.6875, 0.779167, 0.870833, 0.9625, 1.05417};

  // Tree Variables
  std::vector<float> *b_cluster_towere = 0;
  std::vector<float> *b_cluster_towereta = 0;
  std::vector<float> *b_cluster_towerphi = 0;
  std::vector<double> deltaR_values;
  std::vector<double> transverse_energy_values;
  int b_gl1_triggervec[64] = {0};
  float m_vertex = std::numeric_limits<float>::quiet_NaN();

  chain_jet04.SetBranchStatus("gl1_triggervec", 1);
  chain_jet04.SetBranchAddress("gl1_triggervec", b_gl1_triggervec);
  chain_jet04.SetBranchAddress("cluster_towere", &b_cluster_towere);
  chain_jet04.SetBranchAddress("cluster_towereta", &b_cluster_towereta);
  chain_jet04.SetBranchAddress("cluster_towerphi", &b_cluster_towerphi);
  chain_jet04.SetBranchAddress("zvertex", &m_vertex);
  
  // Create histograms for different energy ranges
  TH2F *hist_deltaR_transverseEnergy_0_2 = new TH2F("hist_deltaR_transverseEnergy_0_2", "Delta R vs Transverse Energy (0 < E < 2)", 50, 0, .5, 500, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_2_4 = new TH2F("hist_deltaR_transverseEnergy_2_4", "Delta R vs Transverse Energy (2 < E < 4)", 50, 0, .5, 500, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_4_6 = new TH2F("hist_deltaR_transverseEnergy_4_6", "Delta R vs Transverse Energy (4 < E < 6)", 50, 0, .5, 500, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_6_8 = new TH2F("hist_deltaR_transverseEnergy_6_8", "Delta R vs Transverse Energy (6 < E < 8)", 50, 0, .5, 500, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_8_plus = new TH2F("hist_deltaR_transverseEnergy_8_plus", "Delta R vs Transverse Energy (E > 8)", 50, 0, .5, 500, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_ALL = new TH2F ("hist_deltaR_transverseEnergy_ALL", "Delta R vs Transverse Energy", 50, 0 ,0.5, 500, 0, 5);
 
  int trigger_bit = 17; // Adjust this index as needed
  int nEntries = chain_jet04.GetEntries();

  double jet_radius = 0.4;  // Define jet size radius  


  for (int ev = 0; ev < nEntries; ++ev) {
    chain_jet04.GetEntry(ev);

    // Check trigger condition
    if (b_gl1_triggervec[trigger_bit] != 1) {
      continue; // Skip events that do not pass the trigger condition
    }

    double maxEnergy = 0;
    double maxEta = 0;
    double maxPhi = 0;

    // Loop over towers in the event
    for (int i = 0; i < b_cluster_towere->size(); ++i) {
      int eta_idx = static_cast<int>(b_cluster_towereta->at(i));
      int phi_idx = static_cast<int>(b_cluster_towerphi->at(i));
      double energy = b_cluster_towere->at(i);

      // Validate eta index
      if (eta_idx < 0 || eta_idx >= n_hcal_etabin) {
        continue; // Skip invalid eta index
      }

      // Map eta and phi indices to actual values
      double eta = eta_map[eta_idx];
      double phi = phi_idx * (2 * M_PI / 256);

      // Validate eta and phi values
      if (std::isnan(eta) || std::isnan(phi) || eta < -1.1 || eta > 1.1 || phi < 0 || phi > 2 * M_PI) {
        continue; // Skip invalid values
      }

      // Loop over trigger patch bins
      for (int ietabin = 0; ietabin < n_jettrigger_etabin; ++ietabin) {
        for (int iphibin = 0; iphibin < n_jettrigger_phibin; ++iphibin) {
          // Apply trigger patch eta-phi constraints
          if (eta >= eta_map[jettrigger_min_etabin[ietabin]] && eta <= eta_map[jettrigger_max_etabin[ietabin]] &&
              phi >= jettrigger_min_phibin[iphibin] * (2 * M_PI / 64) && phi <= jettrigger_max_phibin[iphibin] * (2 * M_PI / 64)) {

            // Find tower with maximum energy within the trigger patch
            if (energy > maxEnergy) {
              maxEnergy = energy;
              maxEta = eta;
              maxPhi = phi;
            }
          }
        }
      }
    

      // Skip filling histograms and calculations if maxEnergy is 0
      if (maxEnergy == 0) {
	continue;
      }
      // Debug: Print maxEnergy before determining histogram to fill
      //    std::cout << "maxEnergy: " << maxEnergy << std::endl;
      // Determine which histogram to fill based on maxEnergy
      
      TH2F *hist_to_fill = nullptr;
      // Fill the selected histogram
      if (hist_to_fill) {
	for (int i = 0; i < b_cluster_towere->size(); ++i) {
	  int eta_idx = static_cast<int>(b_cluster_towereta->at(i));
	  int phi_idx = static_cast<int>(b_cluster_towerphi->at(i));
	  double energy_i = b_cluster_towere->at(i);
	  
	  // Map eta and phi indices to actual values
	  double eta_i = eta_map[eta_idx];
	  double phi_i = phi_idx * (2 * M_PI / 256);
	  
	  // Validate eta and phi values
	  if (std::isnan(eta_i) || std::isnan(phi_i) || eta_i < -1.1 || eta_i > 1.1 || phi_i < 0 || phi_i > 2 * M_PI) {
	    continue; // Skip invalid values
	  }
	  
	  // Calculate delta R
	  double dPhi = maxPhi - phi_i;
	  if (dPhi > +M_PI) dPhi -= 2 * M_PI;
	  if (dPhi < -M_PI) dPhi += 2 * M_PI;
	  double dEta = maxEta - eta_i;
	  double deltaR = std::sqrt(dPhi * dPhi + dEta * dEta);
	  
	  //	std::cout << "dPhi: " << dPhi << std::endl;
	  //	std::cout << "dEta: " << dEta << std::endl;
	  if (deltaR < jet_radius) {
	    double transverseEnergy = energy_i * std::sin(std::atan2(std::sqrt(dEta * dEta + dPhi * dPhi), maxEta));
	  
	    // Store delta R and transverse energy values
	    deltaR_values.push_back(deltaR);
	    transverse_energy_values.push_back(transverseEnergy);

	    // Fill the appropriate histogram
            if (maxEnergy < 2) {
	      hist_deltaR_transverseEnergy_0_2->Fill(deltaR, transverseEnergy);
            } else if (maxEnergy < 4) {
	      hist_deltaR_transverseEnergy_2_4->Fill(deltaR, transverseEnergy);
            } else if (maxEnergy < 6) {
	      hist_deltaR_transverseEnergy_4_6->Fill(deltaR, transverseEnergy);
            } else if (maxEnergy < 8) {
	      hist_deltaR_transverseEnergy_6_8->Fill(deltaR, transverseEnergy);
            } else {
	      hist_deltaR_transverseEnergy_8_plus->Fill(deltaR, transverseEnergy);
            }
            hist_deltaR_transverseEnergy_ALL->Fill(deltaR, transverseEnergy);
	    
	    /*
	    // std::cout << "transverse_energy:" << transverseEnergy << std::endl;
	    
	    // Determine the Z-axis value (maxEnergy) for this bin
	    int binX = hist_to_fill->GetXaxis()->FindBin(deltaR);
	    int binY = hist_to_fill->GetYaxis()->FindBin(transverseEnergy);
	    double currentMaxEnergy = hist_to_fill->GetBinContent(binX, binY);
	    
	    // Update to maximum energy if current energy is higher
	    if (maxEnergy > currentMaxEnergy) {
	    hist_to_fill->SetBinContent(binX, binY, maxEnergy);
	    // hist_to_fill->Fill(deltaR, transverseEnergy, maxEnergy); // Uncomment this if using TH3 histograms
	    }
	    */
	    hist_to_fill->SetBinContent(deltaR, transverseEnergy);       
	  }
	}
      }
    }
  }
  hist_deltaR_transverseEnergy_0_2->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_0_2->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_0_2->SetZTitle("Trigger Patch Energy (GeV)");
  
  hist_deltaR_transverseEnergy_2_4->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_2_4->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_2_4->SetZTitle("Trigger Patch Energy (GeV)");

  hist_deltaR_transverseEnergy_4_6->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_4_6->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_4_6->SetZTitle("Trigger Patch Energy (GeV)");

  hist_deltaR_transverseEnergy_6_8->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_6_8->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_6_8->SetZTitle("Trigger Patch Energy (GeV)");

  hist_deltaR_transverseEnergy_8_plus->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_8_plus->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_8_plus->SetZTitle("Trigger Patch Energy (GeV)");

  hist_deltaR_transverseEnergy_ALL->SetXTitle("#Delta R");
  hist_deltaR_transverseEnergy_ALL->SetYTitle("Energy_{T} (GeV)");
  hist_deltaR_transverseEnergy_ALL->SetZTitle("Trigger Patch Energy (GeV)");

  hist_deltaR_transverseEnergy_ALL->Add(hist_deltaR_transverseEnergy_0_2);
  hist_deltaR_transverseEnergy_ALL->Add(hist_deltaR_transverseEnergy_2_4);
  hist_deltaR_transverseEnergy_ALL->Add(hist_deltaR_transverseEnergy_4_6);
  hist_deltaR_transverseEnergy_ALL->Add(hist_deltaR_transverseEnergy_6_8);
  hist_deltaR_transverseEnergy_ALL->Add(hist_deltaR_transverseEnergy_8_plus);

  // Save histograms to file
  TFile *output_file = new TFile("output_histograms.root", "RECREATE");
  hist_deltaR_transverseEnergy_0_2->Write();
  hist_deltaR_transverseEnergy_2_4->Write();
  hist_deltaR_transverseEnergy_4_6->Write();
  hist_deltaR_transverseEnergy_6_8->Write();
  hist_deltaR_transverseEnergy_8_plus->Write();
  hist_deltaR_transverseEnergy_ALL->Write();
  output_file->Close();

  // Save histograms as PNG files
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_0_2, "hist_deltaR_transverseEnergy_0_2.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_2_4, "hist_deltaR_transverseEnergy_2_4.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_4_6, "hist_deltaR_transverseEnergy_4_6.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_6_8, "hist_deltaR_transverseEnergy_6_8.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_8_plus, "hist_deltaR_transverseEnergy_8_plus.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_ALL, "hist_deltaR_transverseEnergy_ALL.png");

  // Clean up
  delete hist_deltaR_transverseEnergy_0_2;
  delete hist_deltaR_transverseEnergy_2_4;
  delete hist_deltaR_transverseEnergy_4_6;
  delete hist_deltaR_transverseEnergy_6_8;
  delete hist_deltaR_transverseEnergy_8_plus;
  delete hist_deltaR_transverseEnergy_ALL;
  delete output_file;
}
