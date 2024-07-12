#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <TH2F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TMath.h>

void TriggerEnergyT_3() {
  TChain chain_jet04("ttree");
  chain_jet04.Add("../macros/output_47089.root");

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
  const int jettrigger_max_phibin[32] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34}; // wrap around 0 and 31

  // Tree Variables
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
  chain_jet04.SetBranchAddress("trigger_sum_jet", b_trigger_sum_jet);

  TH2D *TowerJetTrigger17 = new TH2D("TowerJetTrigger17", "Tower Jet 2 (8 GeV) Trigger 17 Run 47089", 96, 0, 96, 256, 0, 256);
  TH2D *TowerJetTrigger18 = new TH2D("TowerJetTrigger18", "Tower Jet 3 (10 GeV) Trigger 18 Run 47089", 96, 0, 96, 256, 0, 256);
  TH2D *TowerJetTrigger19 = new TH2D("TowerJetTrigger19", "Tower Jet 4 (12 GeV) Trigger 19 Run 47089", 96, 0, 96, 256, 0, 256);

    // Variables to store highest energy trigger tower information
    double highest_trigger_et = 0;
    double highest_trigger_emcal_et = 0;
    double highest_trigger_ihcal_et = 0;
    double highest_trigger_ohcal_et = 0;
    double highest_trigger_e = 0;
    double highest_trigger_emcal_e = 0;
    double highest_trigger_ihcal_e = 0;
    double highest_trigger_ohcal_e = 0;
    int highest_trigger_eta = 0;
    int highest_trigger_phi = 0;

    float emcal_tower_e[n_emcal_etabin][n_emcal_phibin];
    float ihcal_tower_e[n_hcal_etabin][n_hcal_phibin];
    float ohcal_tower_e[n_hcal_etabin][n_hcal_phibin];
    float eta_map[n_jettrigger_etabin];
    // Assuming you have defined n_jettrigger_etabin and n_jettrigger_phibin appropriately
    int trigger_bit = 17; // Adjust this index as needed
    int nEntries = chain_jet04.GetEntries();

    for (int ev = 0; ev < nEntries; ++ev) {
      chain_jet04.GetEntry(ev);

      // Check trigger condition
      if (b_gl1_triggervec[trigger_bit] != 1) {
	continue; // Skip events that do not pass the trigger condition
      }

      // Loop over towers in the event
      for (int i = 0; i < b_cluster_towere->size(); ++i) {
	double eta = b_cluster_towereta->at(i);
	double phi = b_cluster_towerphi->at(i);
	double energy = b_cluster_towere->at(i);

	// Initialize variables to track maximum energy tower within trigger patch
	double maxEnergy = 0;
	double maxEta = 0;
	double maxPhi = 0;

	// Loop over trigger patch bins
	for (int ietabin = 0; ietabin < n_jettrigger_etabin; ++ietabin) {
	  for (int iphibin = 0; iphibin < n_jettrigger_phibin; ++iphibin) {

	    // Apply trigger patch eta-phi constraints
	    if (eta >= jettrigger_min_etabin[ietabin] && eta <= jettrigger_max_etabin[ietabin] &&
		phi >= jettrigger_min_phibin[iphibin] && phi <= jettrigger_max_phibin[iphibin]) {

	      // Find tower with maximum energy within the trigger patch
	      if (energy > maxEnergy) {
		maxEnergy = energy;
		maxEta = eta;
		maxPhi = phi;
	      }
	    }
	  }
	}

	// Print maxEnergy, maxEta, and maxPhi for each event
	std::cout << "Event: " << ev + 1 << std::endl;
	std::cout << "maxEnergy: " << maxEnergy << std::endl;
	std::cout << "maxEta: " << maxEta << std::endl;
	std::cout << "maxPhi: " << maxPhi << std::endl;
	std::cout << "---------------------------" << std::endl;
   
    
	// Transverse energy calculations
	for (int i = 0; i < n_jettrigger_etabin; ++i) {
	  for (int j = 0; j < n_jettrigger_phibin; ++j) {
	    double trigger_emcal_e[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_ihcal_e[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_ohcal_e[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_e[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_emcal_et[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_ihcal_et[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_ohcal_et[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    double trigger_et[n_jettrigger_etabin][n_jettrigger_phibin] = {0};
	    
	    for (int k = jettrigger_min_etabin[i]; k <= jettrigger_max_etabin[i]; ++k) {
	      for (int l = jettrigger_min_phibin[j]; l <= jettrigger_max_phibin[j]; ++l) {
		if (l >= 64) {
		  trigger_emcal_e[i][j] += emcal_tower_e[k][l-64];
		  trigger_ihcal_e[i][j] += ihcal_tower_e[k][l-64];
		  trigger_ohcal_e[i][j] += ohcal_tower_e[k][l-64];
		  trigger_e[i][j] += emcal_tower_e[k][l-64] + ihcal_tower_e[k][l-64] + ohcal_tower_e[k][l-64];
		  trigger_emcal_et[i][j] += emcal_tower_e[k][l-64] / (float)(TMath::CosH(eta_map[k]));
		  trigger_ihcal_et[i][j] += ihcal_tower_e[k][l-64] / (float)(TMath::CosH(eta_map[k]));
		  trigger_ohcal_et[i][j] += ohcal_tower_e[k][l-64] / (float)(TMath::CosH(eta_map[k]));
		  trigger_et[i][j] += (emcal_tower_e[k][l-64] + ihcal_tower_e[k][l-64] + ohcal_tower_e[k][l-64]) / (float)(TMath::CosH(eta_map[k]));
		} else {
		  trigger_emcal_e[i][j] += emcal_tower_e[k][l];
		  trigger_ihcal_e[i][j] += ihcal_tower_e[k][l];
		  trigger_ohcal_e[i][j] += ohcal_tower_e[k][l];
		  trigger_e[i][j] += emcal_tower_e[k][l] + ihcal_tower_e[k][l] + ohcal_tower_e[k][l];
		  trigger_emcal_et[i][j] += emcal_tower_e[k][l] / (float)(TMath::CosH(eta_map[k]));
		  trigger_ihcal_et[i][j] += ihcal_tower_e[k][l] / (float)(TMath::CosH(eta_map[k]));
		  trigger_ohcal_et[i][j] += ohcal_tower_e[k][l] / (float)(TMath::CosH(eta_map[k]));
		  trigger_et[i][j] += (emcal_tower_e[k][l] + ihcal_tower_e[k][l] + ohcal_tower_e[k][l]) / (float)(TMath::CosH(eta_map[k]));
		}
	      }
	    }
	    if (trigger_et[i][j] > highest_trigger_et) {
	      highest_trigger_et = trigger_et[i][j];
	      highest_trigger_emcal_et = trigger_emcal_et[i][j];
	      highest_trigger_ihcal_et = trigger_ihcal_et[i][j];
	      highest_trigger_ohcal_et = trigger_ohcal_et[i][j];
	      highest_trigger_e = trigger_e[i][j];
	      highest_trigger_emcal_e = trigger_emcal_e[i][j];
	      highest_trigger_ihcal_e = trigger_ihcal_e[i][j];
	      highest_trigger_ohcal_e = trigger_ohcal_e[i][j];
	      highest_trigger_eta = i;
	      highest_trigger_phi = j;
	    }
	  }
	}
	
	// Calculate delta R and transverse energy for all towers relative to the highest energy tower
	for (int i = 0; i < b_cluster_towere->size(); ++i) {
	  double eta_i = b_cluster_towereta->at(i);
	  double phi_i = b_cluster_towerphi->at(i);
	  double energy_i = b_cluster_towere->at(i);
	  
	  // Calculate delta R
	  double dPhi = maxPhi - phi_i;
	  if (dPhi > +3.14159) dPhi -= 2 * 3.14159;
	  if (dPhi < -3.14159) dPhi += 2 * 3.14159;
	  double dEta = maxEta - eta_i;
	  double deltaR = std::sqrt(dPhi * dPhi + dEta * dEta);
	  
	  // Calculate transverse energy
	  double transverseEnergy = energy_i * std::sin(std::atan2(std::sqrt(dEta * dEta + dPhi * dPhi), highest_trigger_eta));
	  
	  // Store delta R and transverse energy values
	  deltaR_values.push_back(deltaR);
	  transverse_energy_values.push_back(transverseEnergy);
	}
      }
      
      // Debugging: Print the sizes of the vectors
      std::cout << "Number of deltaR values: " << deltaR_values.size() << std::endl;
      std::cout << "Number of transverse energy values: " << transverse_energy_values.size() << std::endl;
      
      // Debugging: Print some of the values
      for (size_t i = 0; i < std::min(deltaR_values.size(), size_t(10)); ++i) {
	std::cout << "deltaR: " << deltaR_values[i] << ", transverseEnergy: " << transverse_energy_values[i] << std::endl;
      }
    }
    // Create a histogram to plot delta R vs transverse energy
    TH2F *hist_deltaR_transverseEnergy = new TH2F("hist_deltaR_transverseEnergy", "Delta R vs Transverse Energy", 100, 0, 5, 100, 0, 50);
    for (size_t i = 0; i < deltaR_values.size(); ++i) {
      hist_deltaR_transverseEnergy->Fill(deltaR_values[i], transverse_energy_values[i]);
    }
    
    TCanvas *c1 = new TCanvas("c1", "Delta R vs Transverse Energy", 800, 600);
    hist_deltaR_transverseEnergy->Draw("COLZ");
    c1->SaveAs("deltaR_vs_transverseEnergy.png");
}
