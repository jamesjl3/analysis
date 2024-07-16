#include "sPhenixStyle.h"
#include "sPhenixStyle.C"
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <TH2F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>

void saveHistogramAsPNG(TH2F *histogram, const char *filename) {
  //  SetsPhenixStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 700);
  histogram->GetXaxis()->SetTitle("#Delta R");
  histogram->GetYaxis()->SetTitle("Transverse Energy (GeV)");
  canvas->SetLogz();
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.2);

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
  int b_gl1_triggervec[64] = {0};
  float m_vertex = std::numeric_limits<float>::quiet_NaN();
  int trigger_jet_energy[n_jettrigger_etabin][n_jettrigger_phibin] = {{0}};
  int trigger_jet_patch[n_jettrigger_etabin][n_jettrigger_phibin];

  chain_jet04.SetBranchStatus("gl1_triggervec", 1);
  chain_jet04.SetBranchAddress("gl1_triggervec", b_gl1_triggervec);
  chain_jet04.SetBranchAddress("cluster_towere", &b_cluster_towere);
  chain_jet04.SetBranchAddress("cluster_towereta", &b_cluster_towereta);
  chain_jet04.SetBranchAddress("cluster_towerphi", &b_cluster_towerphi);
  chain_jet04.SetBranchAddress("zvertex", &m_vertex);
  chain_jet04.SetBranchAddress("trigger_sum_jet", trigger_jet_patch);

  // Create histograms for different energy ranges
  TH2F *hist_deltaR_transverseEnergy_0_2 = new TH2F("hist_deltaR_transverseEnergy_0_2", "Delta R vs Transverse Energy (0 < E < 2 GeV)", 360, 0, 3.60, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_2_4 = new TH2F("hist_deltaR_transverseEnergy_2_4", "Delta R vs Transverse Energy (2 < E < 4 GeV)", 360, 0, 3.60, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_4_6 = new TH2F("hist_deltaR_transverseEnergy_4_6", "Delta R vs Transverse Energy (4 < E < 6 GeV)", 360, 0, 3.60, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_6_8 = new TH2F("hist_deltaR_transverseEnergy_6_8", "Delta R vs Transverse Energy (6 < E < 8 GeV)", 360, 0, 3.60, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_8_plus = new TH2F("hist_deltaR_transverseEnergy_8_plus", "Delta R vs Transverse Energy (E > 8 GeV)", 360, 0, 3.60, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_ALL = new TH2F("hist_deltaR_transverseEnergy_ALL", "Delta R vs Transverse Energy (All Energy)", 360, 0, 3.60, 100, 0, 10);

  int trigger_bit = 17; // Adjust this index as needed
  int nEntries = chain_jet04.GetEntries();

  for (int ev = 0; ev < nEntries; ++ev) {
    chain_jet04.GetEntry(ev);

    // Check trigger condition
    if (b_gl1_triggervec[trigger_bit] != 1) {
      continue; // Skip events that do not pass the trigger condition
    }

    double maxEnergyPatch = 0;
    double maxEtaPatch = 0;
    double maxPhiPatch = 0;
    bool foundTriggerPatch = false; // Flag to track if the trigger patch is found

    for (int ietabin = 0; ietabin < n_jettrigger_etabin; ++ietabin) {
      for (int iphibin = 0; iphibin < n_jettrigger_phibin; ++iphibin) {
        // Check if the current patch fired the trigger
        if (trigger_jet_patch[ietabin][iphibin] == 1) {
	  // Now check if it has higher energy than the current maxEnergyPatch
	  if (trigger_jet_energy[ietabin][iphibin] > maxEnergyPatch) {
	    maxEnergyPatch = trigger_jet_energy[ietabin][iphibin];
	    maxEtaPatch = eta_map[jettrigger_min_etabin[ietabin]];
	    maxPhiPatch = jettrigger_min_phibin[iphibin] * (2 * M_PI / 256);
	    foundTriggerPatch = true; // Set flag when the trigger patch is found
	  }
        }
      }
      if (foundTriggerPatch) {
        break; // Exit outer loop once trigger patch is found
      }
    }

    // Ensure that the trigger patch was found before proceeding
    if (!foundTriggerPatch) {
      continue; // Skip the rest of the processing if no trigger patch was found
    }
      
    // Calculate delta R and transverse energy for sub-towers within the event
    for (int j = 0; j < b_cluster_towere->size(); ++j) {
      double sub_tower_eta = eta_map[static_cast<int>(b_cluster_towereta->at(j))];
      double sub_tower_phi = b_cluster_towerphi->at(j) * (2 * M_PI / 256);
      double sub_tower_energy = b_cluster_towere->at(j);
      
      // Validate sub-tower eta and phi values
      if (std::isnan(sub_tower_eta) || std::isnan(sub_tower_phi) || sub_tower_eta < -1.1 || sub_tower_eta > 1.1 || sub_tower_phi < 0 || sub_tower_phi > 2 * M_PI) {
	continue; // Skip invalid values
      }
      
      // Calculate delta R between maxEtaPatch, maxPhiPatch and sub-tower
      double delta_eta = (maxEtaPatch - sub_tower_eta);
      double delta_phi = (maxPhiPatch - sub_tower_phi);

      //      std::cout << "maxEtaPatch: " << maxEtaPatch << std::endl;
      //      std::cout << "maxPhiPatch: " << maxPhiPatch << std::endl;
      //      std::cout << "delta_eta: " << delta_eta << std::endl;     

      // Wrap delta_eta within [-2.2, 2.2)  
      if (delta_eta > 2.2) {
	delta_eta -= 4.4;
      } else if (delta_eta < -2.2) {
	delta_eta += 4.4;
      }
      // Wrap delta_phi within [-π, π)
      if (delta_phi > M_PI) {
	delta_phi -= 2 * M_PI;
      } else if (delta_phi < -M_PI) {
	delta_phi += 2 * M_PI;
      }
      double deltaR = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
      
      // Calculate transverse energy for the sub-tower
      double transverse_energy = sub_tower_energy / std::cosh(sub_tower_eta);
         
      // Fill the appropriate histogram based on maxEnergyPatch
      TH2F *hist_to_fill = nullptr;
      if (maxEnergyPatch < 2) {
	hist_to_fill = hist_deltaR_transverseEnergy_0_2;
      } else if (maxEnergyPatch < 4) {
	hist_to_fill = hist_deltaR_transverseEnergy_2_4;
      } else if (maxEnergyPatch < 6) {
	hist_to_fill = hist_deltaR_transverseEnergy_4_6;
      } else if (maxEnergyPatch < 8) {
	hist_to_fill = hist_deltaR_transverseEnergy_6_8;
      } else {
	hist_to_fill = hist_deltaR_transverseEnergy_8_plus;
      }
      
      // Fill histograms with delta R and transverse energy values
      hist_to_fill->Fill(deltaR, transverse_energy);
      hist_deltaR_transverseEnergy_ALL->Fill(deltaR, transverse_energy);
    }
  }
    // Save histograms as PNG files
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_0_2, "hist_deltaR_transverseEnergy_0_2.png");
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_2_4, "hist_deltaR_transverseEnergy_2_4.png");
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_4_6, "hist_deltaR_transverseEnergy_4_6.png");
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_6_8, "hist_deltaR_transverseEnergy_6_8.png");
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_8_plus, "hist_deltaR_transverseEnergy_8_plus.png");
    saveHistogramAsPNG(hist_deltaR_transverseEnergy_ALL, "hist_deltaR_transverseEnergy_ALL.png");

    std::cout << "Done!" << std::endl;
}
