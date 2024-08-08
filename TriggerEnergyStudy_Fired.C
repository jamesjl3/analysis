#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TROOT.h>
//#include <omp.h>

void saveHistogramAsPNG(TH2F *histogram, const char *filename) {
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

void TriggerEnergyStudy_Fired() {
  gROOT->SetBatch(kTRUE);
  // SetsPhenixStyle(); // Uncomment if you need this

  // Create output file
  std::string outFile = "output_histograms_48089_fired.root";
  TFile *out = new TFile(outFile.c_str(), "RECREATE");

  // List of input files
  std::vector<std::string> fileList;
  for (int i = 1; i <= 40; ++i) {
    fileList.push_back(Form("/sphenix/tg/tg01/jets/jamesj3j3/run48089/extracted_tree_48089_%d.root", i));
  }

  // Create histograms
  TH2F *hist_deltaR_transverseEnergy_0_1_fired = new TH2F("hist_deltaR_transverseEnergy_0_1_fired", "Delta R vs Transverse Energy (0 < E < 1 GeV)", 100, 0, 1, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_1_2_fired = new TH2F("hist_deltaR_transverseEnergy_1_2_fired", "Delta R vs Transverse Energy (1 < E < 2 GeV)", 100, 0, 1, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_2_3_fired = new TH2F("hist_deltaR_transverseEnergy_2_3_fired", "Delta R vs Transverse Energy (2 < E < 3 GeV)", 100, 0, 1, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_3_4_fired = new TH2F("hist_deltaR_transverseEnergy_3_4_fired", "Delta R vs Transverse Energy (3 < E < 4 GeV)", 100, 0, 1, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_4_plus_fired = new TH2F("hist_deltaR_transverseEnergy_4_plus_fired", "Delta R vs Transverse Energy (E > 4 GeV)", 100, 0, 1, 100, 0, 10);
  TH2F *hist_deltaR_transverseEnergy_ALL_fired = new TH2F("hist_deltaR_transverseEnergy_ALL_fired", "Delta R vs Transverse Energy (All Energy)", 100, 0, 1, 100, 0, 10);

  int trigger_bit = 17; // Adjust this index as needed

  #pragma omp parallel
  {
    // Iterate over files in parallel
    #pragma omp for
    for (int i = 0; i < fileList.size(); ++i) {
      std::string fileName = fileList[i];
      TFile *in = TFile::Open(fileName.c_str());
      if (!in || in->IsZombie()) {
	std::cerr << "Failed to open file: " << fileName << std::endl;
        continue;
      }

      // Create a TChain and add the file
      TChain chain_jet04("ttree_48089");
      chain_jet04.Add(fileName.c_str());

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

      int trigger_jet_patch[n_jettrigger_etabin][n_jettrigger_phibin];
      chain_jet04.SetBranchAddress("trigger_sum_jet", trigger_jet_patch);

      int nEntries = chain_jet04.GetEntries();

      for (int ev = 0; ev < nEntries; ++ev) {
        chain_jet04.GetEntry(ev);

        // Check trigger condition
        if (b_gl1_triggervec[trigger_bit] != 1) {
          continue; // Skip events that do not pass the trigger condition
        }

        double maxEnergy = -1;
        double maxEta = -1;
        double maxPhi = -1;

        int trigger_patch_eta_idx = -1;
        int trigger_patch_phi_idx = -1;
        int max_jet_patch = -1;

	// Determine the maximum energy sum patch
	for (int ietabin = 0; ietabin < n_jettrigger_etabin; ++ietabin) {
	  for (int iphibin = 0; iphibin < n_jettrigger_phibin; ++iphibin) {
            if (trigger_jet_patch[ietabin][iphibin] > max_jet_patch) {
	      max_jet_patch = trigger_jet_patch[ietabin][iphibin];
	      trigger_patch_eta_idx = ietabin;
	      trigger_patch_phi_idx = iphibin;
            }
	  }
	}

	if (trigger_patch_eta_idx == -1 || trigger_patch_phi_idx == -1) {
	  std::cerr << "Warning: No valid jet patch found for this event." << std::endl;
	  continue;  // Skip events with no valid jet patch
	}
	
	// Loop over towers in the patch for the event
	for (int i = 0; i < b_cluster_towere->size(); ++i) {
	  int eta_idx = static_cast<int>(b_cluster_towereta->at(i));
	  int phi_idx = static_cast<int>(b_cluster_towerphi->at(i));
	  double energy = b_cluster_towere->at(i);
	  
	  if (eta_idx < 0 || eta_idx >= n_hcal_etabin) {
            continue;
	  }
	  
	  double eta = eta_map[eta_idx];
	  double phi = phi_idx * (2 * M_PI / 256);
	  
	  if (std::isnan(eta) || std::isnan(phi) || eta < -1.1 || eta > 1.1 || phi < 0 || phi > 2 * M_PI) {
            continue;
	  }
	  
	  // Ensure we are in the valid patch area
	  if (eta >= eta_map[jettrigger_min_etabin[trigger_patch_eta_idx]] &&
	      eta <= eta_map[jettrigger_max_etabin[trigger_patch_eta_idx]] &&
	      phi >= jettrigger_min_phibin[trigger_patch_phi_idx] * (2 * M_PI / 64) &&
	      phi <= jettrigger_max_phibin[trigger_patch_phi_idx] * (2 * M_PI / 64)) {
	    
	    if (energy > maxEnergy) {
	      maxEnergy = energy;
	      maxEta = eta;
	      maxPhi = phi;
	    }
	  }
	}

	if (maxEnergy <= 0) {
	  std::cerr << "Skipping event due to no valid maxEnergy found." << std::endl;
	  continue;  // Skip to the next event
	}
	
	std::cout << "maxEnergy: " << maxEnergy << std::endl;
	
	// Check for inconsistencies with tower energies
	for (int j = 0; j < b_cluster_towere->size(); ++j) {
	  double sub_tower_eta = eta_map[static_cast<int>(b_cluster_towereta->at(j))];
	  double sub_tower_phi = b_cluster_towerphi->at(j) * (2 * M_PI / 256);
	  double sub_tower_energy = b_cluster_towere->at(j);
	  
	  if (sub_tower_energy > maxEnergy) {
	    std::cerr << "Inconsistency detected: sub_tower_energy = " << sub_tower_energy << " > maxEnergy = " << maxEnergy << std::endl;
	  }
	  
	  
	  // Ensure the sub-tower is within the eta and phi ranges of the fired patch
	  if (sub_tower_eta >= eta_map[jettrigger_min_etabin[trigger_patch_eta_idx]] &&
	      sub_tower_eta <= eta_map[jettrigger_max_etabin[trigger_patch_eta_idx]] &&
	      sub_tower_phi >= jettrigger_min_phibin[trigger_patch_phi_idx] * (2 * M_PI / 64) &&
	      sub_tower_phi <= jettrigger_max_phibin[trigger_patch_phi_idx] * (2 * M_PI / 64)) {
	    
	    // Skip invalid eta or phi values
	    if (std::isnan(sub_tower_eta) || std::isnan(sub_tower_phi) || sub_tower_eta < -1.1 || sub_tower_eta > 1.1 || sub_tower_phi < 0 || sub_tower_phi > 2 * M_PI) {
	      continue;
	    }
	    
	    // Skip if the tower is the leading tower itself
	    if (sub_tower_eta == maxEta && sub_tower_phi == maxPhi) {
	      continue;
	    }
	    
	    // Calculate delta R
	    double delta_eta = maxEta - sub_tower_eta;
	    double delta_phi = maxPhi - sub_tower_phi;
	    // Handle the phi wrap-around
	    if (delta_phi > M_PI) {
	      delta_phi -= 2 * M_PI;
	    } else if (delta_phi < -M_PI) {
	      delta_phi += 2 * M_PI;
	    }
	    double deltaR = TMath::Sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
	    
	    // Calculate transverse energy
	    double transverseEnergy = sub_tower_energy/TMath::CosH(sub_tower_eta);
	    
	    TH2F *hist_to_fill = nullptr;
	    if (maxEnergy > 0 && maxEnergy < 1) {
	      hist_to_fill = hist_deltaR_transverseEnergy_0_1_fired;
	    } else if (maxEnergy < 2) {
	      hist_to_fill = hist_deltaR_transverseEnergy_1_2_fired;
	    } else if (maxEnergy < 3) {
	      hist_to_fill = hist_deltaR_transverseEnergy_2_3_fired;
	    } else if (maxEnergy < 4) {
	      hist_to_fill = hist_deltaR_transverseEnergy_3_4_fired;
	    } else {
	      hist_to_fill = hist_deltaR_transverseEnergy_4_plus_fired;
	    }
	    
	    if(hist_to_fill){
	      hist_to_fill->Fill(deltaR, transverseEnergy);
	      hist_deltaR_transverseEnergy_ALL_fired->Fill(deltaR, transverseEnergy);
	    }
	  }
	}
      }
      // Close file and clean up                                                                                                                                              
      in->Close();
      delete in;
    }
  }
  
  // Save histograms to output file
  out->cd();
  hist_deltaR_transverseEnergy_0_1_fired->Write();
  hist_deltaR_transverseEnergy_1_2_fired->Write();
  hist_deltaR_transverseEnergy_2_3_fired->Write();
  hist_deltaR_transverseEnergy_3_4_fired->Write();
  hist_deltaR_transverseEnergy_4_plus_fired->Write();
  hist_deltaR_transverseEnergy_ALL_fired->Write();
  
  // Save histograms as PNG files
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_0_1_fired, "48089_hist_deltaR_transverseEnergy_0_1_fired.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_1_2_fired, "48089_hist_deltaR_transverseEnergy_1_2_fired.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_2_3_fired, "48089_hist_deltaR_transverseEnergy_2_3_fired.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_3_4_fired, "48089_hist_deltaR_transverseEnergy_3_4_fired.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_4_plus_fired, "48089_hist_deltaR_transverseEnergy_4_plus_fired.png");
  saveHistogramAsPNG(hist_deltaR_transverseEnergy_ALL_fired, "48089_hist_deltaR_transverseEnergy_ALL_fired.png");
  
  std::cout << "Done!" << std::endl;
  out->Close();
}
 
