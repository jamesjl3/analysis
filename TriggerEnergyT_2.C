#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <TH2F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>

void TriggerEnergyT_2() {
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
  TH2F *hist_deltaR_transverseEnergy_0_2 = new TH2F("hist_deltaR_transverseEnergy_0_2", "Delta R vs Transverse Energy (0 < E < 2)", 100, 0, .5, 100, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_2_4 = new TH2F("hist_deltaR_transverseEnergy_2_4", "Delta R vs Transverse Energy (2 < E < 4)", 100, 0, .5, 100, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_4_6 = new TH2F("hist_deltaR_transverseEnergy_4_6", "Delta R vs Transverse Energy (4 < E < 6)", 100, 0, .5, 100, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_6_8 = new TH2F("hist_deltaR_transverseEnergy_6_8", "Delta R vs Transverse Energy (6 < E < 8)", 100, 0, .5, 100, 0, 5);
  TH2F *hist_deltaR_transverseEnergy_8_plus = new TH2F("hist_deltaR_transverseEnergy_8_plus", "Delta R vs Transverse Energy (E > 8)", 100, 0, .5, 100, 0, 5);

  int trigger_bit = 17; // Adjust this index as needed
  int nEntries = chain_jet04.GetEntries();

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

      // Print cluster information for debugging
      //      std::cout << "Cluster " << i << ": eta=" << eta << ", phi=" << phi << ", energy=" << energy << std::endl;

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
  

    // Debugging: Print maxEnergy, maxEta, and maxPhi for each event
    // std::cout << "Event: " << ev + 1 << std::endl;
    // std::cout << "maxEnergy: " << maxEnergy << std::endl;
    // std::cout << "maxEta: " << maxEta << std::endl;
    // std::cout << "maxPhi: " << maxPhi << std::endl;
    // std::cout << "---------------------------" << std::endl;
      if (maxEnergy == 0) {
	return;  // Skip filling histograms and calculations if maxEnergy is 0
      }
      // Determine which histogram to fill based on maxEnergy
      TH2F *hist_to_fill = nullptr;
      if (maxEnergy > 0 && maxEnergy <= 2) {
	hist_to_fill = hist_deltaR_transverseEnergy_0_2;
      } else if (maxEnergy > 2 && maxEnergy <= 4) {
	hist_to_fill = hist_deltaR_transverseEnergy_2_4;
      } else if (maxEnergy > 4 && maxEnergy <= 6) {
	hist_to_fill = hist_deltaR_transverseEnergy_4_6;
      } else if (maxEnergy > 6 && maxEnergy <= 8) {
	hist_to_fill = hist_deltaR_transverseEnergy_6_8;
      } else if (maxEnergy > 8) {
	hist_to_fill = hist_deltaR_transverseEnergy_8_plus;
      }
      
      if (hist_to_fill) {
	// Calculate delta R and transverse energy for all towers relative to the highest energy tower
	for (int i = 0; i < b_cluster_towere->size(); ++i) {
	  int eta_idx = static_cast<int>(b_cluster_towereta->at(i));
	  int phi_idx = static_cast<int>(b_cluster_towerphi->at(i));
	  double energy_i = b_cluster_towere->at(i);
	  
	  // Validate eta index
	  if (eta_idx < 0 || eta_idx >= n_hcal_etabin) {
	    continue; // Skip invalid eta index
	  }
	  
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
	  
	  // Calculate transverse energy
	  double transverseEnergy = energy_i * std::sin(std::atan2(std::sqrt(dEta * dEta + dPhi * dPhi), maxEta));
	  
	  // Store delta R and transverse energy values
	  deltaR_values.push_back(deltaR);
	  transverse_energy_values.push_back(transverseEnergy);
	}
      }
    }
  }   
  // Debugging: Print the sizes of the vectors
  std::cout << "Number of deltaR values: " << deltaR_values.size() << std::endl;
  std::cout << "Number of transverse energy values: " << transverse_energy_values.size() << std::endl;
  
  // Debugging: Print some of the values
  for (size_t i = 0; i < std::min(deltaR_values.size(), size_t(10)); ++i) {
    std::cout << "deltaR: " << deltaR_values[i] << ", transverseEnergy: " << transverse_energy_values[i] << std::endl;
  }
  // Save histograms
  TCanvas *c1 = new TCanvas("c1", "Delta R vs Transverse Energy (0 < E < 2)", 800, 600);
  hist_deltaR_transverseEnergy_0_2->Draw("COLZ");
  c1->SaveAs("deltaR_vs_transverseEnergy_0_2.png");
  // delete c1;
  
  TCanvas *c2 = new TCanvas("c2", "Delta R vs Transverse Energy (2 < E < 4)", 800, 600);
  hist_deltaR_transverseEnergy_2_4->Draw("COLZ");
  c2->SaveAs("deltaR_vs_transverseEnergy_2_4.png");
  // delete c2;
  
  TCanvas *c3 = new TCanvas("c3", "Delta R vs Transverse Energy (4 < E < 6)", 800, 600);
  hist_deltaR_transverseEnergy_4_6->Draw("COLZ");
  c3->SaveAs("deltaR_vs_transverseEnergy_4_6.png");
  // delete c3;
  
  TCanvas *c4 = new TCanvas("c4", "Delta R vs Transverse Energy (6 < E < 8)", 800, 600);
  hist_deltaR_transverseEnergy_6_8->Draw("COLZ");
  c4->SaveAs("deltaR_vs_transverseEnergy_6_8.png");
  // delete c4;
  
  TCanvas *c5 = new TCanvas("c5", "Delta R vs Transverse Energy (E > 8)", 800, 600);
  hist_deltaR_transverseEnergy_8_plus->Draw("COLZ");
  c5->SaveAs("deltaR_vs_transverseEnergy_8_plus.png");
  //  delete c5;
}
