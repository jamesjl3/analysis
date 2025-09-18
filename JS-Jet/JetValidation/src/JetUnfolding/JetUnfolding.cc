//module for producing a TTree with jet information for doing jet validation studies
// for questions/bugs please contact Jennifer James  jennifer.l.james@vanderbilt.edu
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <jetbase/Jetv1.h>
#include <jetbase/Jet.h>
#include <jetbase/JetAlgo.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetInput.h>
#include <g4jets/TruthJetInput.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Hit.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetMapv1.h>
#include <jetbase/JetContainerv1.h>
#include <centrality/CentralityInfo.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <JetUnfolding.h>
#include <mbd/MbdOut.h>
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/JetDefinition.hh>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh" // In external code, this should be fastjet/contrib/SoftDrop.hh              
#include <map>
#include "TProfile.h"
#include <utility>
#include <cstdlib>  // for exit             
#include <memory>  // for allocator_traits<>::value_type             
#include <Pythia8/Pythia.h> // Include the Pythia header                                                                     
#include <jetbackground/TowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
// tower rho
#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>
//#include <cstdlib>  // for exit               
#include <fastjet/PseudoJet.hh>
#include <iostream>
//#include <sstream>
//#include <iomanip>
#include <numeric>
#include <cmath>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <TTree.h>
#include <TGraph.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenCrossSection.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHGenIntegral.h>
#include <phhepmc/PHGenIntegralv1.h>

using namespace fastjet;
//____________________________________________________________________________..
JetUnfolding::JetUnfolding(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename):
  SubsysReco("JetUnfolding_" + recojetname + "_" + truthjetname)
  , hDeltaR(nullptr)
  , hRecoJetPt(nullptr)
  , hTruthJetPt(nullptr)
  , hFakeRecoJetPt(nullptr)
  , hFakeTruthJetPt(nullptr)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1.1, 1.1)
  , m_ptRange(5, 100)
  , m_doTruthJets(0)
  , m_event(-1)
  , m_nJet(-1)
  , m_nTruthJet(-1)
  , m_id()
  , m_nComponent()
  , m_eta()
  , m_phi()
  , m_e()
  , m_pt()
  , m_id_truth()
  , m_nComponent_truth()
  , m_eta_truth()
  , m_phi_truth()
  , m_e_truth()
  , m_pt_truth()
{
  std::cout << "JetUnfolding::JetUnfolding(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
JetUnfolding::~JetUnfolding()
{
  std::cout << "JetUnfolding::~JetUnfolding() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
//____________________________________________________________________________..
// Define the member function                                                                                                                                                                                                                                                   
float JetUnfolding::calculateDeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if (dphi > M_PI) dphi -= 2 * M_PI;
  if (dphi < -M_PI) dphi += 2 * M_PI;
  return std::sqrt(deta * deta + dphi * dphi);
}
//____________________________________________________________________________..                                                                                                                                                                                                
int JetUnfolding::Init(PHCompositeNode *topNode)
{
  std::cout << "JetUnfolding::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "JetUnfolding::Init - Output to " << m_outputFileName << std::endl;
  // Initialize histograms                                                                                                                                                                                                                                                      
  hDeltaR = new TH1F("hDeltaR", "Delta R Distribution;#DeltaR;Counts", 50, 0, 1.0);
  hRecoJetPtTotal = new TH1F("hRecoJetPtTotal", "Reco Jet pT; p_{T} (GeV/c);Counts", 11, 5, 60);
  hTruthJetPtTotal = new TH1F("hTruthJetPtTotal", "Truth Jet pT; p_{T} (GeV/c);Counts", 11, 5, 60);
  hRecoJetPt = new TH1F("hRecoJetPt", "Reco Jet pT; p_{T} (GeV/c);Counts", 11, 5, 60);
  hTruthJetPt = new TH1F("hTruthJetPt", "Truth Jet pT; p_{T} (GeV/c);Counts", 11, 5, 60);
  hFakeRecoJetPt = new TH1F("hFakeRecoJetPt", "Fake Reco Jet pT; p_{T} (GeV/c);Counts", 11, 5 ,60);
  hFakeTruthJetPt = new TH1F("hFakeTruthJetPt", "Fake Truth Jet pT; p_{T} (GeV/c);Counts", 11, 5, 60);
  PurityRatio = new TH1F("PurityRatio", "Purity Ratio;Reco Jet p_{T};Purity", 11, 5, 60);
  EfficiencyRatio = new TH1F("EfficiencyRatio", "Efficiency Ratio;Truth Jet p_{T};Efficiency", 11, 5, 60);
  spectrum = new TH1F("spectrum","Jet Spectrum;p_{T} [GeV];1/N_{jets}", 11, 5, 60);
  energy_spectrum = new TH1F("energy_spectrum", "Energy Sepctrum in CEMC;E [GeV];1/N_{jets}", 11, 5, 60);

  eta_pt = new TH2D("eta_pt",";p_{T} [GeV];#eta",9,5,50,23,-1.1,1.1);
  phi_pt = new TH2D("phi_pt",";p_{T} [GeV];#phi",9,5,50,23,-M_PI,M_PI);
  M_pt = new TH2D("M_pt",";p_{T} [GeV];M [GeV]",9,5,50,23,0,50);

  spectrum_truth = new TH1F("spectrum_truth","Jet Spectrum;p_{T} [GeV];1/N_{jets}", 11, 5, 60);
  energy_spectrum_truth = new TH1F("energy_spectrum_truth", "Energy Sepctrum in CEMC;E [GeV];1/N_{jets}", 11, 5, 60);

  eta_pt_truth = new TH2D("eta_pt_truth",";p_{T} [GeV];#eta",9,5,50,23,-1.1,1.1);
  phi_pt_truth = new TH2D("phi_pt_truth",";p_{T} [GeV];#phi",9,5,50,23,-M_PI,M_PI);
  M_pt_truth = new TH2D("M_pt_truth",";p_{T} [GeV];M [GeV]",9,5,50,23,0,50);
  
  jet_AN = new TH1F("jet_AN",";AN",10,0,1);

  //avg_ratio = new TH1D("avg_ratio", "<Det p_{T}/Truth p_{T}>;Det. p_{T, unc.} (GeV/c);<Det p_{T}/Truth p_{T}>", 11, 5, 60);
  // Define 2D histograms to store ratios
  h2_ratio_RecoJetPt = new TH2D("h2_ratio_RecoJetPt", "Reco Jet pT Correction", 11, 5, 60, 40, 0, 2);
  h2_ratio_CorrectedRecoJetPt = new TH2D("h2_ratio_CorrectedRecoJetPt", "Fake Reco Jet pT Correction", 11, 5, 60, 40, 0, 2);
  h2_ratio_RecoJetPtTotal = new TH2D("h2_ratio_RecoJetPtTotal", "Total Reco Jet pT Correction", 11, 5, 60, 40, 0, 2);
  //  ratioGraph = new TGraph;
  corrected_jet_spectra_total = new TH1F("corrected_jet_spectra_total","Corrected Jet Spectra ALL jets; Det. p_{T, cor.} (GeV/c);1/N_{jets}",11, 5, 60);
  corrected_fake_jet_spectra = new TH1F("corrected_fake_jet_spectra", "Corrected Fake Jet Spectra; Det. p_{T, cor.} (GeV/x);1/N_{jets}", 11, 5, 60);
  corrected_jet_spectra = new TH1F("corrected_jet_spectra","Corrected Jet Spectra; Det. p_{T, cor.} (GeV/c);1/N_{jets}",11, 5, 60);
  uncorrected_jet_spectra = new TH1F("uncorrected_jet_spectra","Uncorrected Jet Spectra; Det. p_{T, cor.} (GeV/c);1/N_{jets}", 11, 5, 60);
  hResponseMatrix = new TH2F("hResponseMatrix", "Response Matrix; Reco Jet p_{T} (GeV); Truth Jet p_{T} (GeV)", 11, 5, 60, 11, 5, 60);
  hResponseMatrix_Corrected = new TH2F("hResponseMatrix_Corrected", "Response Matrix; Truth Jet p_{T,corr} (GeV); Reco Jet p_{T,corr} (GeV)", 11, 5, 60, 11, 5, 60);

  truth_jet_spectra = new TH1F("truth_jet_spectra", "; [GeV];1/N_{jets}", 11, 5, 60);
  hCrossSection = new TH1F("hCrossSection", "Cross Section per Event;Event;Cross Section [mb]", 1000, 0, 1000);
  hCrossSectionDist = new TH1F("hCrossSectionDist", "Cross Section Distribution;Cross Section [mb];Counts", 100, 0, 100);
  // Histograms for differential jet cross sections
  hTruthJetXsec = new TH1F("hTruthJetXsec", "Truth Jet Cross Section; p_{T} [GeV]; d#sigma/dp_{T} [mb/GeV]", 11, 5, 60);
  hRecoJetXsec = new TH1F("hRecoJetXsec", "Reco Jet Cross Section; p_{T} [GeV]; d#sigma/dp_{T} [mb/GeV]", 11, 5, 60);
  hRecoJetXsec_Corrected = new TH1F("hRecoJetXsec_Corrected", "Corrected Reco Jet Cross Section; p_{T} [GeV]; d#sigma/dp_{T} [mb/GeV]", 11, 5, 60);
 
  //m_event_tree = new TTree("event_info", "EventInfo");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetUnfolding::InitRun(PHCompositeNode* topNode)
{
  std::cout << "JetUnfolding::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  PHGenIntegral *genintegral = findNode::getClass<PHGenIntegral>(topNode, "PHGenIntegral");

  if (genintegral)
  {
    float cross_sec_pb = genintegral->get_CrossSection_Processed_Event(); // in picobarns
    float lumi_pb = genintegral->get_Integrated_Lumi();                   // in pb^-1
    int ngen = genintegral->get_N_Processed_Event();                      // total generated events

    float cross_sec_mb = cross_sec_pb * 1e-9; // convert to millibarns (mb)

    std::cout << ">> PHGenIntegral Info <<\n";
    std::cout << "  Events Processed: " << ngen << "\n";
    std::cout << "  Cross Section: " << cross_sec_pb << " pb (" << cross_sec_mb << " mb)\n";
    std::cout << "  Integrated Lumi: " << lumi_pb << " pb^-1\n";

    // Store for use in End()
    m_totalCrossSectionMB = cross_sec_mb;
    m_totalEvents = ngen;
  }
  else
  {
    std::cerr << "PHGenIntegral node not found! Can't retrieve cross section.\n";
    m_totalCrossSectionMB = 0;
    m_totalEvents = 0;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetUnfolding::process_event(PHCompositeNode *topNode)
{
  std::cout << "JetUnfolding::process_event(PHCompositeNode *topNode) Processing Event " << m_event+1 << std::endl;
  ++m_event;
  
  //particle map 
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
    {
      std::cout
        << "MyJetAnalysis:;process_event - Error can not find DST G4TruthInfo node"
        << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  // Interface to reco jets
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  if (!jets) {
    std::cerr << "MyJetAnalysis::process_event - Error: cannot find DST Reco JetContainer node " << m_recoJetName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // JetContainer for truth jets
  JetContainer* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  if (!jetsMC && m_doTruthJets) {
    std::cerr << "MyJetAnalysis::process_event - Error: cannot find DST Truth JetMap node " << m_truthJetName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // Change exit to return appropriate code
  }
  
  //centrality
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent_node)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find centrality node "
        << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  //zvertex
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find global vertex  node "
        << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  if (vertexmap->empty())
    {
      std::cout
        << "MyJetAnalysis::process_event - global vertex node is empty "
        << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  else
    {
      GlobalVertex *vtx = vertexmap->begin()->second;
      m_zvtx = vtx->get_z();
      CLHEP::Hep3Vector vertex;
      vertex.set(vtx->get_x(),vtx->get_y(),vtx->get_z());
    }
  if (fabs(m_zvtx) > 30.0)
    {
      std::cout << "JetUnfolding::process_event - bad vertex" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  //calorimeter towers
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if(!towersEM3 || !towersIH3 || !towersOH3){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower node "
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(!tower_geom || !tower_geomOH || !tower_geomEM){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower geometry "
      << std::endl;
    exit(-1);
  }
  
  std::vector<PseudoJet> particles_reco;
  particles_reco.clear();
  float totalPx_reco = 0;
  float totalPy_reco = 0;
  float totalPz_reco = 0;
  float totalE_reco = 0;
  m_centrality = (int)(100*cent_node->get_centile(CentralityInfo::PROP::mbd_NS));
  m_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);

  //-------------Jet Matching--------------------------------------//
  // Vectors to hold matched jets
  std::vector<const Jet*> matchedTruthJets;
  std::vector<const Jet*> matchedRecoJets;

  // Counters for selected truth jets
  //  int numTruthJetsSelected = 0;

  // Create sorted indices of truth and reco jets based on pT
  std::vector<size_t> truthIndices(jetsMC->size());
  std::vector<size_t> recoIndices(jets->size());
  std::vector<int> recoToTruthMatchIndex(jets->size(), -1);  // Initialize with -1 (invalid)
  
  std::iota(truthIndices.begin(), truthIndices.end(), 0);
  std::iota(recoIndices.begin(), recoIndices.end(), 0);
  
  std::sort(truthIndices.begin(), truthIndices.end(), [&](size_t i, size_t j) {
    return (*jetsMC)[i]->get_pt() > (*jetsMC)[j]->get_pt();
  });
  
  std::sort(recoIndices.begin(), recoIndices.end(), [&](size_t i, size_t j) {
    return (*jets)[i]->get_pt() > (*jets)[j]->get_pt();
  });
  
  // Keep track of which reco jets have been matched
  std::vector<bool> recoJetMatched(jets->size(), false);
  
  // Define the range of truth jet pT bins
  std::vector<std::pair<float, float>> truthPtBins = {
    {30.0, 35.0}, {35.0, 40.0}, {40.0, 45.0}, {45.0, 50.0}, {50.0, 55.0}, {55.0, 60.0}
  };
  std::vector<std::pair<float, float>> recoPtBins = {
    {5.0, 10.0}, {10.0, 15.0}, {15.0, 20.0}, {20.0, 25.0}, {25.0, 30.0}, {30.0, 35.0}, {35.0, 40.0}, {40.0, 45.0}, {45.0, 50.0}, {50.0, 55.0}, {55.0, 60.0}
  };
  struct MatchedJetPair {
    const Jet* recoJet;
    float truthPt;
    size_t recoPtBin;
    size_t truthPtBin;
  };
  /*
  // Struct to store matched pairs and the bin index
  struct MatchedJetPair {
    const Jet* recoJet;
    float truthPt;
    size_t truthPtBin;
  };
  */
  std::vector<MatchedJetPair> matchedJets;

  // Vectors to store the sums of reco jet pt and truth jet pt for each truth pt bin
  std::vector<float> recoPtSum(recoPtBins.size(), 0.0);
  std::vector<float> truthPtSum(recoPtBins.size(), 0.0);
  // 1. Loop over truth jets
  for (size_t truthIdx = 0; truthIdx < jetsMC->size(); ++truthIdx) {
    const auto& truthJet = (*jetsMC)[truthIdx];
    float truthJetPt = truthJet->get_pt();
        std::cout << "Truth Jet pt: " << truthJetPt << std::endl;
    float truthJetEta = truthJet->get_eta();
    std::cout << "Truth jet eta: " << truthJetEta << std::endl;
    float truthJetPhi = truthJet->get_phi();
    std::cout << "Truth jet phi: " << truthJetPhi << std::endl;
    
    hTruthJetPtTotal->Fill(truthJetPt);
    
    if (truthJetPt < 30.0 || std::abs(truthJetEta) >= 1.1) continue;
    // Find which reco pt bin the *recoJet->pt()* falls into
    int recoPtBin = -1;
    // Commented out: We now bin by reco pT, not truth pT
    /*
      int truthPtBin = -1;
      for (size_t binIdx = 0; binIdx < truthPtBins.size(); ++binIdx) {
      if (truthJetPt >= truthPtBins[binIdx].first && truthJetPt < truthPtBins[binIdx].second) {
      truthPtBin = binIdx;
      break;
      }
      }
      
      if (truthPtBin == -1) continue; // Skip if no matching bin found
    */
    
    bool matched = false;
    const Jet* bestRecoJet = nullptr;
    size_t bestRecoIdx = 0;
    float minDeltaR = 0.2;
    
    for (size_t recoIdx = 0; recoIdx < jets->size(); ++recoIdx) {
      const auto& recoJet = (*jets)[recoIdx];
      if (recoJetMatched[recoIdx]) continue;
      float deltaR = calculateDeltaR(truthJetEta, truthJetPhi, recoJet->get_eta(), recoJet->get_phi());
      
      if (deltaR < minDeltaR && recoJet->get_pt() > 5.0 && recoJet->get_eta() <= 1.1) {
	minDeltaR = deltaR;
	bestRecoJet = recoJet;
	bestRecoIdx = recoIdx;
	matched = true;
      }
    }
    
    if (matched) {
      int matchedTruthIdx = truthIdx;
      matchedTruthJets.push_back(truthJet);
      matchedRecoJets.push_back(bestRecoJet);
      hDeltaR->Fill(minDeltaR);
      hRecoJetPt->Fill(bestRecoJet->get_pt());
      hTruthJetPt->Fill(truthJetPt);
      recoJetMatched[bestRecoIdx] = true;
      recoToTruthMatchIndex[bestRecoIdx] = matchedTruthIdx;
      
      // NEW: Bin based on reco jet pT
      float recoPt = bestRecoJet->get_pt();
      for (size_t binIdx = 0; binIdx < recoPtBins.size(); ++binIdx) {
	if (recoPt >= recoPtBins[binIdx].first && recoPt < recoPtBins[binIdx].second) {
	  recoPtBin = binIdx;
	  break;
	}
      }
      
      if (recoPtBin == -1) continue; // Skip if reco pT not in any bin
      
      // Save pt sums
      recoPtSum[recoPtBin] += recoPt;
      truthPtSum[recoPtBin] += truthJetPt;
      
      // Store for later correction (now using recoPtBin)
      //matchedJets.push_back({bestRecoJet, truthJetPt, recoPtBin});
      matchedJets.push_back({bestRecoJet, truthJetPt, static_cast<size_t>(recoPtBin)});
    } else {
      hFakeTruthJetPt->Fill(truthJetPt);
    }
  }  
  // Calculate bin ratios (reco-based binning)
  std::vector<double> ptCorrectionRatio(recoPtBins.size(), 1.0);  // Default to 1.0
  for (size_t binIdx = 0; binIdx < recoPtBins.size(); ++binIdx) {
    if (truthPtSum[binIdx] > 0.0) {
      ptCorrectionRatio[binIdx] = recoPtSum[binIdx] / truthPtSum[binIdx];
      std::cout << "Reco Bin (" << recoPtBins[binIdx].first << ", " << recoPtBins[binIdx].second << ") "
		<< "Ratio: " << ptCorrectionRatio[binIdx] << std::endl;
      double binCenter = 0.5 * (truthPtBins[binIdx].first + truthPtBins[binIdx].second);
      h2_ratio_RecoJetPt->Fill(binCenter, ptCorrectionRatio[binIdx]);
      //h2_ratio_RecoJetPt->Fill(truthPtBins[binIdx].first, ptCorrectionRatio[binIdx]);
    }
  }
  std::vector<double> correctedRecoPtSum(truthPtBins.size(), 0.0);
  std::vector<double> truthPtSum_truthBinned(truthPtBins.size(), 0.0);
  std::vector<int> truthBinCounts(truthPtBins.size(), 0);
  
  // Apply correction and track per truth pT bin
  for (const auto& match : matchedJets)
    {
      double ratio = ptCorrectionRatio[match.recoPtBin];
      double correctedPt = match.recoJet->get_pt();
      float truthPt = match.truthPt;
      
      if (ratio > 0.0)
	{
	  correctedPt /= ratio;
	}
      
      // Fill corrected spectra (no change)
      corrected_jet_spectra->Fill(correctedPt);

      hResponseMatrix->Fill(match.recoJet->get_pt(), truthPt);        // raw reco pt
      hResponseMatrix_Corrected->Fill(correctedPt, truthPt);
      
      // Find truth pt bin (now inline)
      int truthBin = -1;
      for (size_t binIdx = 0; binIdx < truthPtBins.size(); ++binIdx)
	{
	  if (truthPt >= truthPtBins[binIdx].first && truthPt < truthPtBins[binIdx].second)
	    {
	      truthBin = binIdx;
	      break;
	    }
	}

      if (truthBin >= 0)
	{
	  correctedRecoPtSum[truthBin] += correctedPt;
	  truthPtSum_truthBinned[truthBin] += truthPt;
	  truthBinCounts[truthBin]++;
	}
    }
  for (size_t binIdx = 0; binIdx < truthPtBins.size(); ++binIdx)
    {
      if (truthBinCounts[binIdx] > 0)
  {
    double meanCorrected = correctedRecoPtSum[binIdx] / truthBinCounts[binIdx];
    double meanTruth = truthPtSum_truthBinned[binIdx] / truthBinCounts[binIdx];
    double correctionRatio = meanCorrected / meanTruth;
    
    double binCenter = 0.5 * (truthPtBins[binIdx].first + truthPtBins[binIdx].second);
    h2_ratio_CorrectedRecoJetPt->Fill(binCenter, correctionRatio);
    
    std::cout << "Corrected Ratio Bin (" << binCenter << "): " << correctionRatio << std::endl;
  }
    }
  
  for (size_t recoIdx = 0; recoIdx < jets->size(); ++recoIdx) {
    if (recoJetMatched[recoIdx]) continue;
    const auto& recoJet = (*jets)[recoIdx];
    if (recoJet->get_pt() > 5.0 && fabs(recoJet->get_eta()) <= 1.1) {
      hFakeRecoJetPt->Fill(recoJet->get_pt());
    }
  }
  
  for (size_t recoIdx = 0; recoIdx < jets->size(); ++recoIdx) {
    const auto& recoJet = (*jets)[recoIdx];
    float recoJetPt = recoJet->get_pt();
    // Skip reco jets below 5 GeV
    if (recoJetPt < 5) continue;    
    // Fill total reco jet pT spectrum
    hRecoJetPtTotal->Fill(recoJetPt);
  }

  // Loop over all reco jets to fill corrected total spectrum
  for (size_t recoIdx = 0; recoIdx < jets->size(); ++recoIdx) {
    const auto& recoJet = (*jets)[recoIdx];
    double recoJetPt = recoJet->get_pt();
    
    // Apply kinematic cuts
    if (recoJetPt < 5.0 || fabs(recoJet->get_eta()) > 1.1) continue;
    
    // Determine reco pt bin
    int recoPtBin = -1;
    for (size_t binIdx = 0; binIdx < recoPtBins.size(); ++binIdx) {
      if (recoJetPt >= recoPtBins[binIdx].first && recoJetPt < recoPtBins[binIdx].second) {
	recoPtBin = binIdx;
	break;
      }
    }
    
    // If valid bin and nonzero correction ratio, apply correction and fill
    if (recoPtBin != -1 && ptCorrectionRatio[recoPtBin] > 0.0) {
      double correctedPt = recoJetPt / ptCorrectionRatio[recoPtBin];
      corrected_jet_spectra_total->Fill(correctedPt);
    }
  }
  
  // Loop again for fake jets (unmatched)
  for (size_t recoIdx = 0; recoIdx < jets->size(); ++recoIdx) {
    if (recoJetMatched[recoIdx]) continue;
    
    const auto& recoJet = (*jets)[recoIdx];
    double recoJetPt = recoJet->get_pt();
    
    if (recoJetPt < 5.0 || fabs(recoJet->get_eta()) > 1.1) continue;
    
    int recoPtBin = -1;
    for (size_t binIdx = 0; binIdx < recoPtBins.size(); ++binIdx) {
    if (recoJetPt >= recoPtBins[binIdx].first && recoJetPt < recoPtBins[binIdx].second) {
      recoPtBin = binIdx;
      break;
    }
    }
    
    if (recoPtBin != -1 && ptCorrectionRatio[recoPtBin] > 0.0) {
      double correctedPt = recoJetPt / ptCorrectionRatio[recoPtBin];
      corrected_fake_jet_spectra->Fill(correctedPt);
    }
  }
    
  //get reco jets
  m_nJet = 0;
  
  double leading_pT = -999;
  double subleading_pT = -999;
  for (size_t i = 0; i < matchedRecoJets.size(); ++i) {
    const auto& recoJet = matchedRecoJets[i];
    //    const auto& truthJet = matchedTruthJets[i];
    m_nJet++;
    double recoJetPt = recoJet->get_pt();
    if(recoJet->get_pt() > leading_pT)
      {
        subleading_pT = leading_pT;
        leading_pT = recoJet->get_pt();
      }
    else if(recoJet->get_pt() > subleading_pT)
      {
        subleading_pT = recoJet->get_pt();
      }
    
    m_id.push_back(recoJet->get_id());
    m_nComponent.push_back(recoJet->size_comp());
    m_eta.push_back(recoJet->get_eta());
    m_phi.push_back(recoJet->get_phi());
    if(m_phi[m_phi.size()-1] < -M_PI)
      {
        m_phi[m_phi.size()-1] += 2*M_PI;
      }
    if(m_phi[m_phi.size()-1] > M_PI)
      {
        m_phi[m_phi.size()-1] -= 2*M_PI;
      }
    m_e.push_back(recoJet->get_e());
    energy_spectrum->Fill(recoJet->get_e());
    m_pt.push_back(recoJet->get_pt());
    spectrum->Fill(recoJetPt);
    eta_pt->Fill(recoJetPt, recoJet->get_eta());
    phi_pt->Fill(recoJetPt, recoJet->get_phi());
    M_pt->Fill(recoJetPt,recoJet->get_mass());    
    if(subleading_pT > 0 && leading_pT > 0) 
      {
	jet_AN->Fill(subleading_pT/leading_pT);
      }
    int nconst = 0;
    
    auto* nonConstJet = const_cast<Jet*>(recoJet);
    for (auto comp: nonConstJet->get_comp_vec())
      {
	// std::cout << "jet " << m_nJet << " consituent " << nconst << std::endl;
        TowerInfo *tower;
        nconst++;
        unsigned int channel = comp.second;
	//    std::cout << "Number of components in jet: "
	//                  << nonConstJet->get_comp_vec().size() << std::endl;
	//  std::cout << "Component type: " << comp.first << std::endl;
        if (comp.first == 30) //HCALIN: Unsub componenet = 26, sub component = 30
          {
	    //  std::cout << "Processing component ihcal" << std::endl;
            tower = towersIH3->get_tower_at_channel(channel);
            if(!tower || !tower_geom){
              continue;
            }
	    unsigned int calokey = towersIH3->encode_key(channel);
            int ieta = towersIH3->getTowerEtaBin(calokey);
            int iphi = towersIH3->getTowerPhiBin(calokey);
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
            // float UE = background->get_UE(1).at(ieta);
            float tower_phi = tower_geom->get_tower_geometry(key)->get_phi();
            float tower_eta = tower_geom->get_tower_geometry(key)->get_eta();
	    
            // UE = UE * (1 + 2 * background_v2 * cos(2 * (tower_phi - background_Psi2)));
            totalE_reco += tower->get_energy()/* + UE*/;
            double pt = (tower->get_energy() /*UE*/) / cosh(tower_eta);
            totalPx_reco += pt * cos(tower_phi);
            totalPy_reco += pt * sin(tower_phi);
            totalPz_reco += pt * sinh(tower_eta);
          }
	else if (comp.first == 31) //HCALOUT: Unsub component = 27, sub component = 31
          {
	    //            std::cout << "Processing component ohcal" << std::endl;
            tower = towersOH3->get_tower_at_channel(channel);
            if(!tower || !tower_geomOH)
              {
                continue;
              }
            unsigned int calokey = towersOH3->encode_key(channel);
            int ieta = towersOH3->getTowerEtaBin(calokey);
            int iphi = towersOH3->getTowerPhiBin(calokey);
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
            // float /*UE*/ = background->get_/*UE*/(2).at(ieta);
            float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
            float tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();
            // /*UE*/ = /*UE*/ * (1 + 2 * background_v2 * cos(2 * (tower_phi - background_Psi2)));
            totalE_reco +=tower->get_energy() /*+UE*/;
            double pt = (tower->get_energy() /*+UE*/) / cosh(tower_eta);
            totalPx_reco += pt * cos(tower_phi);
            totalPy_reco += pt * sin(tower_phi);
            totalPz_reco += pt * sinh(tower_eta);
          }
	else if (comp.first == 29) //EMCAL unsub componenet = 25, sub compenent = 29
          {
	    // std::cout << "Processing component emcal" << std::endl;
            tower = towersEM3->get_tower_at_channel(channel);
            if(!tower || !tower_geomEM)
              {
		// std::cout << "skipping event if needed" << std::endl;
                continue;
              }
	    
            unsigned int calokey = towersEM3->encode_key(channel);
            int ieta = towersEM3->getTowerEtaBin(calokey);
            int iphi = towersEM3->getTowerPhiBin(calokey);
	    // std::cout << "calokey: " << calokey << ", ieta: " << ieta << ", iphi: " << iphi << ", channel: " << channel << std::endl;
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
            // float /*UE*/ = background->get_/*UE*/(0).at(ieta);
            auto tower_geometry_object = tower_geom->get_tower_geometry(key);
            if (!tower_geometry_object) {
              std::cerr << "Error: Tower geometry object is null for key " << key << std::endl;
              continue; // Skip this component and go to the next one
            }
	    float tower_phi = tower_geometry_object->get_phi();
	    //        std::cout << "tower_phi: " << tower_phi << std::endl;
            float tower_eta = tower_geometry_object->get_eta();
	    //  std::cout << "tower_eta: " << tower_eta << std::endl;
            double ETmp = tower->get_energy() /*+ UE*/;
	    //  std::cout << "Etmp: " << ETmp << std::endl;
            double pt = (tower->get_energy() /*+ UE*/) / cosh(tower_eta);
	    
            double pxTmp = pt * cos(tower_phi);
            double pyTmp = pt * sin(tower_phi);
            double pzTmp = pt * sinh(tower_eta);
	    
            totalE_reco += ETmp;
            totalPx_reco += pxTmp;
            totalPy_reco += pyTmp;
            totalPz_reco += pzTmp;
            if (std::isnan(pxTmp) || std::isnan(pyTmp) || std::isnan(pzTmp) || std::isnan(totalE_reco)) {
              std::cerr << "Error: Invalid PseudoJet parameters." << std::endl;
              continue; // Skip this component and go to the next one
            }
            particles_reco.push_back( PseudoJet( pxTmp, pyTmp, pzTmp, totalE_reco) );
	    /* std::cout << "Pushing back PseudoJet with px: " << pxTmp
		      << ", py: " << pyTmp
		      << ", pz: " << pzTmp
		      << ", E: " << totalE_reco << std::endl;
	    */
          }
      }
    
    
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    // get truth reco 
    m_nTruthJet = 0;
    
    double leadingTruth_pT = -999;
    double subleadingTruth_pT = -999;
    for (size_t i = 0; i < matchedTruthJets.size(); ++i) {
      //    const auto& recoJet = matchedRecoJets[i];
      const auto& truthJet = matchedTruthJets[i];
      
      m_nTruthJet++;
      
      double truthJetPt = truthJet->get_pt();
      if (truthJetPt > leadingTruth_pT)
	{
	  subleadingTruth_pT = leadingTruth_pT;
	  leadingTruth_pT = truthJetPt;
	}
      else if (truthJetPt > subleadingTruth_pT)
	{
	  subleadingTruth_pT = truthJetPt;
	}
      
      // Store truth jet properties
      m_id_truth.push_back(truthJet->get_id());
      m_nComponent_truth.push_back(truthJet->size_comp());
      m_eta_truth.push_back(truthJet->get_eta());
      m_phi_truth.push_back(truthJet->get_phi());
      m_e_truth.push_back(truthJet->get_e());
      m_pt_truth.push_back(truthJetPt);
      
      energy_spectrum_truth->Fill(truthJet->get_e());
      spectrum_truth->Fill(truthJetPt);
      eta_pt_truth->Fill(truthJetPt, truthJet->get_eta());
      phi_pt_truth->Fill(truthJetPt, truthJet->get_phi());
      M_pt_truth->Fill(truthJetPt, truthJet->get_mass());
      /*    
      // Inside event loop                                                                                                                        
      double ratio = recoJetPt / truthJetPt;
      h2_ratio->Fill(truthJetPt, ratio);  // X = truthJetPt, Y = ratio                                                                            
      // After event loop: Profile X to get mean corrections                                                                                      
      TH1D* mean_ratio = (TH1D*)h2_ratio->ProfileX("mean_ratio")->Clone("mean_ratio_TH1D");
      // Loop over bins of the mean ratio and print correction values for truth jet pT                                                            
      for (int i = 1; i <= mean_ratio->GetNbinsX(); i++) {
        double binCenter = mean_ratio->GetBinCenter(i);
        double meanValue = mean_ratio->GetBinContent(i);
        std::cout << "Truth Jet pT: " << binCenter << " GeV -> Mean Correction: " << meanValue << std::endl;
      }
      
      // Inside second event loop (applying correction)                                                                                           
      double correction = mean_ratio->GetBinContent(mean_ratio->FindBin(recoJetPt));
      
      // Check if correction factor is valid and apply the correction                                                                             
      double correctedRecoJetPt = recoJetPt;  // Default to recoJetPt if correction is invalid                                               
      if (correction > 0) {
        correctedRecoJetPt = recoJetPt / correction;
      }
      // Fill histograms with the corrected and uncorrected jet pT                                                                                
      corrected_jet_spectra->Fill(correctedRecoJetPt);
      */     
      uncorrected_jet_spectra->Fill(recoJetPt);
      
      // Fill truth jet spectrum
      truth_jet_spectra->Fill(truthJetPt);
      
      // Loop over bins to calculate purity (5 GeV bins)
      for (int bin = 1; bin <= PurityRatio->GetNbinsX(); ++bin) {
	// Get the lower and upper edges of the current bin (5 GeV increments)
	double reco_pt_low = PurityRatio->GetXaxis()->GetBinLowEdge(bin);
	double reco_pt_high = reco_pt_low + 2.5;  // 5 GeV bin width
	
	// Get the total number of reco jets in this 5 GeV range
	double N_total_reco = corrected_jet_spectra_total->Integral(
							corrected_jet_spectra_total->FindBin(reco_pt_low),
							corrected_jet_spectra_total->FindBin(reco_pt_high)
							);
	
	
	// Get the number of unmatched (fake) reco jets in this 5 GeV range
	double N_fake_reco = corrected_fake_jet_spectra->Integral(
						      corrected_fake_jet_spectra->FindBin(reco_pt_low),
						      corrected_fake_jet_spectra->FindBin(reco_pt_high)
						      );
	
	// Compute purity for the 5 GeV bin
	double purity = (N_total_reco > 0) ? (1.0 - (N_fake_reco / N_total_reco)) : 0.0;
	
	// Set the purity value in the histogram
	PurityRatio->SetBinContent(bin, purity);
      }
      // Loop over bins to calculate efficiency (5 GeV bins)
      for (int bin = 1; bin <= EfficiencyRatio->GetNbinsX(); ++bin) {
	// Get the lower and upper edges of the current bin (5 GeV increments)
	double truth_pt_low = EfficiencyRatio->GetXaxis()->GetBinLowEdge(bin);
	double truth_pt_high = truth_pt_low + 2.5;  // 5 GeV bin width
	
	// Get the total number of truth jets in this 5 GeV range
	double N_total_truth = hTruthJetPtTotal->Integral(
							  hTruthJetPtTotal->FindBin(truth_pt_low),
							  hTruthJetPtTotal->FindBin(truth_pt_high)
							  );
	// Get the number of unmatched (fake) truth jets in this 5 GeV range
	double N_fake_truth = hFakeTruthJetPt->Integral(
							hFakeTruthJetPt->FindBin(truth_pt_low),
							hFakeTruthJetPt->FindBin(truth_pt_high)
							);	
	// Compute efficiency for the 5 GeV bin
	double efficiency = (N_total_truth > 0) ? (1.0 - (N_fake_truth / N_total_truth)) : 0.0;
	
	// Set the efficiency value in the histogram
	EfficiencyRatio->SetBinContent(bin, efficiency);
      }
    }   
    // fill event tree
    //m_event_tree->Fill();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetUnfolding::ResetEvent(PHCompositeNode* /*topNode*/)
{
  //std::cout << "JetUnfolding::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  m_id.clear();
  m_nComponent.clear();
  m_eta.clear();
  m_phi.clear();
  m_e.clear();
  m_pt.clear();

  m_id_truth.clear();
  m_nComponent_truth.clear();
  m_eta_truth.clear();
  m_phi_truth.clear();
  m_e_truth.clear();
  m_pt_truth.clear();
    
  m_triggerVector.clear();
    
  return Fun4AllReturnCodes::EVENT_OK;
}
  
//____________________________________________________________________________..
int JetUnfolding::EndRun(const int runnumber)
{
  std::cout << "JetUnfolding::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetUnfolding::End(PHCompositeNode* /*topNode*/)
{
  if (m_totalEvents > 0 && m_totalCrossSectionMB > 0)
    {
      float norm_factor = m_totalCrossSectionMB / m_totalEvents;
      hTruthJetXsec->Add(hTruthJetPt);
      hTruthJetXsec->Scale(norm_factor);

      hRecoJetXsec->Add(hRecoJetPt);
      hRecoJetXsec->Scale(norm_factor);
      
      hRecoJetXsec_Corrected->Add(corrected_jet_spectra);
      hRecoJetXsec_Corrected->Scale(norm_factor);
    }
  
  std::cout << "JetUnfolding::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  hDeltaR->Write();
  hRecoJetPt->Write();
  hTruthJetPt->Write();
  hRecoJetPtTotal->Write();
  hTruthJetPtTotal->Write();
  hFakeRecoJetPt->Write();
  hFakeTruthJetPt->Write();
  PurityRatio->Write();
  EfficiencyRatio->Write();
  
  spectrum->Write();
  energy_spectrum->Write();
    
  eta_pt->Write();
  phi_pt->Write();
  M_pt->Write();
 
  spectrum_truth->Write();
  energy_spectrum_truth->Write();

  eta_pt_truth->Write();
  phi_pt_truth->Write();
  M_pt_truth->Write();
   
  jet_AN->Write();

  h2_ratio_RecoJetPt->Write();
  h2_ratio_CorrectedRecoJetPt->Write();
  h2_ratio_RecoJetPtTotal->Write();
  //ratioGraph->Write();
  corrected_jet_spectra_total->Write();
  corrected_fake_jet_spectra->Write();
  corrected_jet_spectra->Write();
  uncorrected_jet_spectra->Write();
  hResponseMatrix->Write();
  hResponseMatrix_Corrected->Write();
  truth_jet_spectra->Write();
  hCrossSection->Write();
  hCrossSectionDist->Write();
  hTruthJetXsec->Write();
  hRecoJetXsec->Write();
  hRecoJetXsec_Corrected->Write();

  std::cout << "JetUnfolding::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetUnfolding::Reset(PHCompositeNode* /*topNode*/)
{
  std::cout << "JetUnfolding::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetUnfolding::Print(const std::string &what) const
{
  std::cout << "JetUnfolding::Print(const std::string &what) const Printing info for " << what << std::endl;
}
