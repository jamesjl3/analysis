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
#include <JetResoStudy.h>

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
#include <utility>
#include <cstdlib>  // for exit             
#include <memory>  // for allocator_traits<>::value_type             
#include <Pythia8/Pythia.h> // Include the Pythia header                                                                     
#include <jetbackground/TowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>


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

using namespace fastjet;
//____________________________________________________________________________..
JetResoStudy::JetResoStudy(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename):
  SubsysReco("JetResoStudy_" + recojetname + "_" + truthjetname)
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
  std::cout << "JetResoStudy::JetResoStudy(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
JetResoStudy::~JetResoStudy()
{
  std::cout << "JetResoStudy::~JetResoStudy() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
//____________________________________________________________________________..
// Define the member function                                                                                                                                                                                                                                                   
float JetResoStudy::calculateDeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if (dphi > M_PI) dphi -= 2 * M_PI;
  if (dphi < -M_PI) dphi += 2 * M_PI;
  return std::sqrt(deta * deta + dphi * dphi);
}
//____________________________________________________________________________..                                                                                                                                                                                                
int JetResoStudy::Init(PHCompositeNode *topNode)
{
  std::cout << "JetResoStudy::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "JetResoStudy::Init - Output to " << m_outputFileName << std::endl;
  // Initialize histograms                                                                                                                                                                                                                                                      
  hDeltaR = new TH1F("hDeltaR", "Delta R Distribution;#DeltaR;Counts", 100, 0, 1.0);
  hRecoJetPt = new TH1F("hRecoJetPt", "Reco Jet pT; p_{T} (GeV/c);Counts", 100, 0, 100);
  hTruthJetPt = new TH1F("hTruthJetPt", "Truth Jet pT; p_{T} (GeV/c);Counts", 100, 0, 100);
  hFakeRecoJetPt = new TH1F("hFakeRecoJetPt", "Fake Reco Jet pT; p_{T} (GeV/c);Counts", 100, 0 ,100);
  hFakeTruthJetPt = new TH1F("hFakeTruthJetPt", "Fake Truth Jet pT; p_{T} (GeV/c);Counts", 100, 0, 100);
  hEfficiency = new TH1F("hEfficiency", "Jet Matching Efficiency; p_{T} bin; Efficiency", 100, 0, 100);
  hPurity = new TH1F("hPurity", "Jet Purity per Event;Purity;Number of Events", 100, 0, 1);
  hFakeRate = new TH1F("hFakeRate","Jet Fake Rate per Event;Fake Rate;Number of Events",100, 0, 1);
  hMomentumResolution = new TH1F("hMomentumResolution", "Momentum Resolution;Resolution;Counts", 100, -2, 2);
  hScaleFactor = new TH1F("hScaleFactor", "Scale Factor;Scale Factor;Counts", 100, 0, 2);

  spectrum = new TH1F("spectrum","Jet Spectrum;p_{T} [GeV];1/N_{jets}",9,5,100);
  energy_spectrum = new TH1F("energy_spectrum", "Energy Sepctrum in CEMC;E [GeV];1/N_{jets}", 20, 0, 100);

  eta_pt = new TH2D("eta_pt",";p_{T} [GeV];#eta",9,5,50,23,-1.1,1.1);
  phi_pt = new TH2D("phi_pt",";p_{T} [GeV];#phi",9,5,50,23,-M_PI,M_PI);
  M_pt = new TH2D("M_pt",";p_{T} [GeV];M [GeV]",9,5,50,23,0,50);

  spectrum_truth = new TH1F("spectrum_truth","Jet Spectrum;p_{T} [GeV];1/N_{jets}",9,5,100);
  energy_spectrum_truth = new TH1F("energy_spectrum_truth", "Energy Sepctrum in CEMC;E [GeV];1/N_{jets}", 20, 0, 100);

  eta_pt_truth = new TH2D("eta_pt_truth",";p_{T} [GeV];#eta",9,5,50,23,-1.1,1.1);
  phi_pt_truth = new TH2D("phi_pt_truth",";p_{T} [GeV];#phi",9,5,50,23,-M_PI,M_PI);
  M_pt_truth = new TH2D("M_pt_truth",";p_{T} [GeV];M [GeV]",9,5,50,23,0,50);
  
  jet_AN = new TH1F("jet_AN",";AN",10,0,1);

  //  Double_t bins[33] = {0.0000000, 0.0014720269, 0.0031607401, 0.0050980365, 0.0073205081, 0.0098701335, 0.012795071, 0.016150566, 0.020000000, 0.024416081, 0.029482220, 0.035294109, 0.041961524, 0.049610400, 0.058385212, 0.068451699, 0.080000000, 0.093248242, 0.10844666, 0.12588233, 0.14588457, 0.16883120, 0.19515564, 0.22535510, 0.26000000, 0.29974473, 0.34533998, 0.39764699, 0.45765372, 0.52649360, 0.60546691, 0.69606529, 0.80000000};

  ////// data //////////////////// Substructure hists
  _h_R04_z_sj_30_35_data= new TH1F("R04_z_sj_30_35_data","z_sj in subjets 1 & 2", 10, 0, 0.5);                 //jet pT 5 -> 10 GeV                                                  
  _h_R04_theta_sj_30_35_data= new TH1F("R04_theta_sj_30_35_data","theta_sj in subjets 1 & 2", 10, 0, 0.5);     //jet pT 5 -> 10 GeV                    
  _h_R04_z_sj_35_40_data= new TH1F("R04_z_sj_35_40_data","z_sj in subjets 1 & 2", 10, 0, 0.5);               //jet pT 10 -> 20 GeV     
  _h_R04_theta_sj_35_40_data= new TH1F("R04_theta_sj_35_40_data","theta_sj in subjets 1 & 2", 10, 0, 0.5);   //jet pT 10 -> 20 GeV
  _h_R04_z_sj_40_45_data= new TH1F("R04_z_sj_40_45_data","z_sj in subjets 1 & 2", 10, 0, 0.5);               //jet pT 20 -> 30 GeV 
  _h_R04_theta_sj_40_45_data= new TH1F("R04_theta_sj_40_45_data","theta_sj in subjets 1 & 2", 10, 0, 0.5);   //jet pT 20 -> 30 GeV                                   
  _h_R04_z_sj_45_data= new TH1F("R04_z_sj_45_data","z_sj in subjets 1 & 2", 10, 0, 0.5);                     //jet pT > 30 GeV  
  _h_R04_theta_sj_45_data= new TH1F("R04_theta_sj_45_data","theta_sj in subjets 1 & 2", 10, 0, 0.5);         //jet pT > 30 GeV
 
  ////// truth //////////////////// Substructure hists
  _h_R04_z_sj_30_35_truth= new TH1F("R04_z_sj_30_35_truth","z_sj in subjets 1 & 2", 10, 0, 0.5);                 //jet pT 5 -> 10 GeV                                                  
  _h_R04_theta_sj_30_35_truth= new TH1F("R04_theta_sj_30_35_truth","theta_sj in subjets 1 & 2", 10, 0, 0.5);     //jet pT 5 -> 10 GeV                    
  _h_R04_z_sj_35_40_truth= new TH1F("R04_z_sj_35_40_truth","z_sj in subjets 1 & 2", 10, 0, 0.5);               //jet pT 10 -> 20 GeV     
  _h_R04_theta_sj_35_40_truth= new TH1F("R04_theta_sj_35_40_truth","theta_sj in subjets 1 & 2", 10, 0, 0.5);   //jet pT 10 -> 20 GeV
  _h_R04_z_sj_40_45_truth= new TH1F("R04_z_sj_40_45_truth","z_sj in subjets 1 & 2", 10, 0, 0.5);               //jet pT 20 -> 30 GeV 
  _h_R04_theta_sj_40_45_truth= new TH1F("R04_theta_sj_40_45_truth","theta_sj in subjets 1 & 2", 10, 0, 0.5);   //jet pT 20 -> 30 GeV                                   
  _h_R04_z_sj_45_truth= new TH1F("R04_z_sj_45_truth","z_sj in subjets 1 & 2", 10, 0, 0.5);                     //jet pT > 30 GeV  
  _h_R04_theta_sj_45_truth= new TH1F("R04_theta_sj_45_truth","theta_sj in subjets 1 & 2", 10, 0, 0.5);         //jet pT > 30 GeV
      
  // data softdrop hists zcut = 0.1                                                                                                                                               
  _h_R04_z_g_30_35_data= new TH1F("R04_z_g_30_35_data","z_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);                    //jet pT 5 -> 10 GeV                              
  _h_R04_theta_g_30_35_data= new TH1F("R04_theta_g_30_35_data","theta_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);        //jet pT 5 -> 10 GeV                               
  _h_R04_z_g_35_40_data= new TH1F("R04_z_g_35_40_data","z_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);                  //jet pT 10 -> 20 GeV
  _h_R04_theta_g_35_40_data= new TH1F("R04_theta_g_35_40_data","theta_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);      //jet pT 10 -> 20 GeV                
  _h_R04_z_g_40_45_data= new TH1F("R04_z_g_40_45_data","z_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);                  //jet pT 20 -> 30 GeV
  _h_R04_theta_g_40_45_data= new TH1F("R04_theta_g_40_45_data","theta_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);      //jet pT 20 -> 30 GeV
  _h_R04_z_g_45_data= new TH1F("R04_z_g_45_data","z_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);                        //jet pT > 30 GeV
  _h_R04_theta_g_45_data= new TH1F("R04_theta_g_45_data","theta_g in subjets 1 & 2 w/ z_{cut} = 0.1 & #beta = 0.0", 10, 0, 0.5);            //jet pT > 30 GeV                    
  // truth softdrop hists zcut = 0.1
  _h_R04_z_g_30_35_truth= new TH1F("R04_z_g_30_35_truth","z_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);                    //jet pT 5 -> 10 GeV                       
  _h_R04_theta_g_30_35_truth= new TH1F("R04_theta_g_30_35_truth","theta_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);        //jet pT 5 -> 10 GeV
  _h_R04_z_g_35_40_truth= new TH1F("R04_z_g_35_40_truth","z_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);                  //jet pT 10 -> 20 GeV
  _h_R04_theta_g_35_40_truth= new TH1F("R04_theta_g_35_40_truth","theta_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);      //jet pT 10 -> 20 GeV                
  _h_R04_z_g_40_45_truth= new TH1F("R04_z_g_40_45_truth","z_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);                  //jet pT 20 -> 30 GeV
  _h_R04_theta_g_40_45_truth= new TH1F("R04_theta_g_40_45_truth","theta_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);      //jet pT 20 -> 30 GeV
  _h_R04_z_g_45_truth= new TH1F("R04_z_g_45_truth","z_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);                        //jet pT > 30 GeV        
  _h_R04_theta_g_45_truth= new TH1F("R04_theta_g_45_truth","theta_g in subjets 1 & 2 w/ z_{cut} = 0.2 & #beta = 0.0", 10, 0, 0.5);            //jet pT > 30 GeV
 


  // corellation hists
  correlation_theta_30_35_data = new TH2D("correlation_theta_30_35_data", "Correlation Plot 5 < p_{T} < 10 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_35_40_data = new TH2D("correlation_theta_35_40_data", "Correlation Plot 10 < p_{T} < 20 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_40_45_data = new TH2D("correlation_theta_40_45_data", "Correlation Plot 20 < p_{T} < 30 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_45_data = new TH2D("correlation_theta_45_data", "Correlation Plot p_{T} > 30 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_30_35_truth = new TH2D("correlation_theta_30_35_truth", "Correlation Plot 5 < p_{T} < 10 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_35_40_truth = new TH2D("correlation_theta_35_40_truth", "Correlation Plot 10 < p_{T} < 20 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_40_45_truth = new TH2D("correlation_theta_40_45_truth", "Correlation Plot 20 < p_{T} < 30 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_theta_45_truth = new TH2D("correlation_theta_45_truth", "Correlation Plot p_{T} > 30 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; R_{g}; #theta_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_30_35_data = new TH2D("correlation_z_5_10_data", "Correlation Plot 5 < p_{T} < 10 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_35_40_data = new TH2D("correlation_z_35_40_data", "Correlation Plot 10 < p_{T} < 20 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_40_45_data = new TH2D("correlation_z_20_30_data", "Correlation Plot 20 < p_{T} < 30 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_45_data = new TH2D("correlation_z_45_data", "Correlation Plot p_{T} > 30 [GeV/c] w/ z_{cut} = 0.1 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_30_35_truth = new TH2D("correlation_z_5_10_truth", "Correlation Plot 5 < p_{T} < 10 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_35_40_truth = new TH2D("correlation_z_35_40_truth", "Correlation Plot 10 < p_{T} < 20 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_40_45_truth = new TH2D("correlation_z_20_30_truth", "Correlation Plot 20 < p_{T} < 30 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);
  correlation_z_45_truth = new TH2D("correlation_z_45_truth", "Correlation Plot p_{T} > 30 [GeV/c] w/ z_{cut} = 0.2 & #beta = 0.0; z_{g}; z_{sj}", 20, 0, 0.5, 20, 0, 0.5);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetResoStudy::InitRun(PHCompositeNode* /*topNode*/)
{
  std::cout << "JetResoStudy::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetResoStudy::process_event(PHCompositeNode *topNode)
{
  std::cout << "JetResoStudy::process_event(PHCompositeNode *topNode) Processing Event " << m_event+1 << std::endl;
  ++m_event;
  /*  
      TruthJetInput truthJetInput(Jet::PARTICLE);

      // Configure the TruthJetInput
      truthJetInput.set_eta_range(m_etaRange.first, m_etaRange.second); // Set eta range

      // Retrieve the truth jets
      std::vector<Jet*> truthJets = truthJetInput.get_input(topNode);
  */
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
  // Initialize counters                                                                                                                                                                           
  int numTruthJetsBefore = jetsMC->size();
  int numRecoJetsBefore = jets->size();
  int numMatchedPairs = 0;

  // Vectors to hold matched jets                                                                                                                                                                  
  std::vector<const Jet*> matchedTruthJets;
  std::vector<const Jet*> matchedRecoJets;

  // Boolean vector to track matched truth jets                                                                                                                                                    
  std::vector<bool> truthJetMatched(jetsMC->size(), false);

  // Boolean vector to track matched reco jets                                                                                                                                                    
  std::vector<bool> recoJetMatched(jets->size(), false); // Add this line to track matched reco jets

  // Counter for selected truth jets
  int numTruthJetsSelected = 0; // Initialize counter for selected truth jets

  // Vector to store delta p_T values
  std::vector<double> deltaPTs;

  // Count matched pairs                                                                                                                                                                          
  for (size_t truthIndex = 0; truthIndex < jetsMC->size(); ++truthIndex) {
    const auto& truthJet = (*jetsMC)[truthIndex];

    // Apply pT cut on truth jets (30 GeV)
    float truthJetPt = truthJet->get_pt(); // Get the pT of the truth jet                                                                                                                          
    float truthJetEta = truthJet->get_eta(); // Get the eta of the truth jet                                                                                                                       
    // Calculate efficiency: Matched truth jets / Selected truth jets                                                                                                                             
    if (truthJetPt >= 30.0 && std::abs(truthJetEta) < 1.1) {
      ++numTruthJetsSelected; // Increment counter for selected truth jets

      bool foundMatchForTruthJet = false; // Track if a match is found for this truth jet

      for (size_t recoIndex = 0; recoIndex < jets->size(); ++recoIndex) {
	const auto& recoJet = (*jets)[recoIndex];

	// Check if this reco jet has already been matched
	if (recoJetMatched[recoIndex]) continue; // Skip if already matched

	float deltaR = calculateDeltaR(truthJet->get_eta(), truthJet->get_phi(),
				       recoJet->get_eta(), recoJet->get_phi());
	if (deltaR < 0.2 && recoJet->get_pt() > 5.0 && truthJet->get_pt() > 5.0) {
	  // Print the delta R value                                                                                                                  
	  hDeltaR->Fill(deltaR);                                               
	  if (!foundMatchForTruthJet) {   
	    //	    std::cout << "Delta R = " << deltaR << std::endl;

	    matchedTruthJets.push_back(truthJet);
	    matchedRecoJets.push_back(recoJet);
	    hRecoJetPt->Fill(recoJet->get_pt());
	    hTruthJetPt->Fill(truthJet->get_pt());
	    /*	    std::cout << "Matched Truth Jet: px = " << truthJet->get_px()
		      << ", py = " << truthJet->get_py()
		      << ", pz = " << truthJet->get_pz()
		      << ", E = " << truthJet->get_e()
		      << ", eta = " << truthJet->get_eta()
		      << ", phi = " << truthJet->get_phi() << std::endl;

	    std::cout << "Matched Reco Jet: px = " << recoJet->get_px()
		      << ", py = " << recoJet->get_py()
		      << ", pz = " << recoJet->get_pz()
		      << ", E = " << recoJet->get_e()
		      << ", eta = " << recoJet->get_eta()
		      << ", phi = " << recoJet->get_phi() << std::endl;
	    */  ++numMatchedPairs; // Increment matched pairs count                                                                                                                                      
	    truthJetMatched[truthIndex] = true; // Mark this truth jet as matched 
	    recoJetMatched[recoIndex] = true; // Mark this reco jet as matched                                                                                                                    
	    foundMatchForTruthJet = true; // A match was found for this truth jet

	    // Calculate and store delta p_T
	    double recoPt = recoJet->get_pt(); // Get the matched reco jet p_T
	    double deltaPT = recoPt - truthJetPt; // Calculate delta p_T
	    deltaPTs.push_back(deltaPT); // Store delta p_T

	    // Calculate and fill the resolution
	    if (truthJetPt > 0) { // Avoid division by zero
	      double resolution = deltaPT / truthJetPt;
	      hMomentumResolution->Fill(resolution); // Fill histogram with resolution
	    }

	    // Calculate efficiency for this truth jet
	    int ptBin = static_cast<int>(truthJetPt / 10);  // Example: 10 GeV bins
	    double efficiency = (numTruthJetsSelected > 0) ?
	      static_cast<double>(matchedTruthJets.size()) / numTruthJetsSelected : 0.0;
	    hEfficiency->Fill(ptBin, efficiency);

	    break; // Exit loop after matching one reco jet to one truth jet                                                               
	  }
	}
      }
    }
  }
  // Declare purity variable outside the conditional block to ensure it's in scope
  double purity = 0.0;
  // Log unmatched events but continue processing
  if (matchedTruthJets.empty() || matchedRecoJets.empty()) {
    // std::cout << "Event " << m_event << " has no matched jets but will still be processed." << std::endl;
  } else {
    // Proceed only if matched jets exist to avoid segmentation faults
    // Calculate the standard deviation of delta p_T
    if (!deltaPTs.empty()) {  
      double mean = std::accumulate(deltaPTs.begin(), deltaPTs.end(), 0.0) / deltaPTs.size();
      double sq_sum = std::inner_product(deltaPTs.begin(), deltaPTs.end(), deltaPTs.begin(), 0.0);
      double stdev = std::sqrt(sq_sum / deltaPTs.size() - mean * mean);

      std::cout << "Momentum Resolution (standard deviation of delta p_T): " << stdev << std::endl;
    } else {
      std::cout << "No delta p_T values to calculate momentum resolution." << std::endl;
    }

    // Calculate per-event purity
    double purity = (numRecoJetsBefore > 0) ? 
      static_cast<double>(numMatchedPairs) / numRecoJetsBefore : 0.0;
    hPurity->Fill(purity);
    //  std::cout << "Purity for this event: " << purity << std::endl;
  }

  // Count unmatched reco jets (only consider reco jets with pT > 5 GeV)                                                     
  std::vector<const Jet*> unmatchedRecoJets;
  int totalRecoJets = 0; // Initialize total reco jet counter  
  double totalRecoPt = 0.0; // Initialize total p_T for reco jets
  double totalTruthPt = 0.0; // Initialize total p_T for truth jets
  for (size_t recoIndex = 0; recoIndex < jets->size(); ++recoIndex) {
    const auto& recoJet = (*jets)[recoIndex];
    // Increment the total reco jet counter
    totalRecoJets++;
    if (recoJet->get_pt() > 5.0 && !recoJetMatched[recoIndex]) {
      unmatchedRecoJets.push_back(recoJet);
      hFakeRecoJetPt->Fill(recoJet->get_pt());
      /*std::cout << "Unmatched Reco Jet: px = " << recoJet->get_px()
		<< ", py = " << recoJet->get_py()
		<< ", pz = " << recoJet->get_pz()
		<< ", E = " << recoJet->get_e()
		<< ", eta = " << recoJet->get_eta()
		<< ", phi = " << recoJet->get_phi() << std::endl;
    */
    } else {
      // Only accumulate p_T for matched reco jets to calculate scale factor later
      if (recoJetMatched[recoIndex]) {
	totalRecoPt += recoJet->get_pt();
      }
    }
  }

  // Calculate the number of unmatched reco jets
  int numUnmatchedRecoJets = unmatchedRecoJets.size(); // Number of unmatched reco jets
  // Calculate fake rate
  double fakeRate = (totalRecoJets > 0) ? 
    static_cast<double>(numUnmatchedRecoJets) / totalRecoJets : 0.0;
  // Fill the fake rate histogram (assuming you have hFakeRate defined)
  hFakeRate->Fill(fakeRate);
  // Print fake rate for debugging
  // std::cout << "Fake Rate for this event: " << fakeRate << std::endl;
  
  // Count unmatched truth jets (only consider truth jets with pT > 30 GeV)
  std::vector<const Jet*> unmatchedTruthJets;
  int totalTruthJets = 0; // Initialize total truth jet counter

  truthJetMatched.resize(jets->size(), false);

  for (size_t truthIndex = 0; truthIndex < jets->size(); ++truthIndex) {
    const auto& truthJet = (*jets)[truthIndex];
       // Increment the total truth jet counter
    totalTruthJets++;
    /*   // Check if the truth jet is valid
    if (!truthJet) {
      //      std::cerr << "Error: Found an invalid truthJet at index " << truthIndex << std::endl;
      continue;
    }
    //  std::cout << "Checking if jets is valid..." << std::endl;
    if (!jets) {
      //   std::cerr << "Error: jets pointer is null!" << std::endl;
      return -1; // Or handle accordingly
    }

    // std::cout << "Total jets: " << jets->size() << std::endl;
    if (jets->size() == 0) {
      //  std::cerr << "Error: No jets available!" << std::endl;
      return -1; // Or handle accordingly
    }

    // Check if truthJetMatched is correctly sized
    if (truthJetMatched.size() < jets->size()) {
      std::cerr << "Error: truthJetMatched size is less than jets size!" << std::endl;
      return -1; // Or handle accordingly
    }
*/
     // Check for unmatched truth jets with pT > 30 GeV
    // Check if the truth jet is valid
    if (!truthJet) {
      std::cerr << "Error: Found an invalid truthJet at index " << truthIndex << std::endl;
      continue; // Skip to the next iteration
    }
  
    if (truthJet->get_pt() > 30.0 && !truthJetMatched[truthIndex]) {
      std::cout << "442"<< std::endl;
      unmatchedTruthJets.push_back(truthJet);
      std::cout << "444"<< std::endl;
      hFakeTruthJetPt->Fill(truthJet->get_pt()); // Fill the histogram with unmatched truth jet pT
      std::cout << "446"<< std::endl;
      // Optional: Print debug information
      std::cout << "Unmatched Truth Jet: pt = " << truthJet->get_pt()
		<< ", eta = " << truthJet->get_eta() 
		<< ", phi = " << truthJet->get_phi() << std::endl;
    } else if (truthJetMatched[truthIndex]) {
      std::cout << "452"<< std::endl;
      // Accumulate matched truth jet p_T for later scale factor calculation
      totalTruthPt += truthJet->get_pt();
      std::cout << "455"<< std::endl;
    }
  }

  // Debugging output for total truth jets
  std::cout << "Total Truth Jets processed: " << totalTruthJets << std::endl;

  // Check for no available truth jets to avoid division by zero
  if (totalTruthJets == 0) {
    std::cerr << "Error: No truth jets available for event " << m_event << std::endl;
    return -1; // Or skip processing this event
  }

  double averageRecoPt = (totalRecoJets > 0) ? totalRecoPt / totalRecoJets : 0.0;
  double averageTruthPt = (totalTruthJets > 0) ? totalTruthPt / totalTruthJets : 0.0;

  // Debugging output
  std::cout << "Average Reco p_T: " << averageRecoPt 
	    << ", Average Truth p_T: " << averageTruthPt << std::endl;

  // Calculate the scale factor only if averageTruthPt is non-zero
  double scaleFactor = 0.0;
  if (averageTruthPt != 0.0) {
    scaleFactor = averageRecoPt / averageTruthPt;
  } else {
    std::cerr << "Warning: Average Truth p_T is zero for event " << m_event << std::endl;
  }
  if (!hScaleFactor) {
    std::cerr << "Error: hScaleFactor histogram is not initialized!" << std::endl;
    return -1; // Or handle the error appropriately
  }
  hScaleFactor->Fill(scaleFactor);

  // Print average p_T and scale factor for debugging
  std::cout << "Average p_T for Reco Jets: " << averageRecoPt << std::endl;
  std::cout << "Average p_T for Truth Jets: " << averageTruthPt << std::endl;
  std::cout << "Scale Factor: " << scaleFactor << std::endl;
 
  // After processing matchedTruthJets and matchedRecoJets
  int numMatchedTruthJets = matchedTruthJets.size();
  int numMatchedRecoJets = matchedRecoJets.size();

  // Print purity for debugging
  std::cout << "Purity for this event: " << purity << std::endl;
  // Print summary                                                                                                                                                                                
  std::cout << "Number of matched Truth Jets: " << numMatchedTruthJets << std::endl;
  std::cout << "Number of matched Reco Jets: " << numMatchedRecoJets << std::endl;
 
  std::cout << "Number of unmatched Truth Jets before matching: " << numTruthJetsBefore << std::endl;
  std::cout << "Number of unmatched Reco Jets before matching: " << numRecoJetsBefore << std::endl;
  std::cout << "Number of selected Truth Jets: " << numTruthJetsSelected << std::endl; // Print selected truth jets
 
  
  /*
  // Match truth jets with reconstructed jets
  std::vector<std::pair<Jet*, Jet*>> matchedJetPairs;
  for (Jet* truthJet : *truthJets) {
  for (Jet* recoJet : *jets) {
  float deltaR = calculateDeltaR(truthJet->get_eta(), truthJet->get_phi(),
  recoJet->get_eta(), recoJet->get_phi());
  if (deltaR < 0.4) {
  // Add the matched pair (recoJet, truthJet) to the vector
  matchedJetPairs.push_back(std::make_pair(recoJet, truthJet));
  break; // Assuming one-to-one matching
  }
  }
  }
  */
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
      std::cout << "JetResoStudy::process_event - bad vertex" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  /*
    if(!m_doSim)
    {
    //trigger check
    //grab the gl1 data                                                                                                                                                                             
    Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1PacketInfo)
    {
    std::cout << PHWHERE << "caloTreeGen::process_event: GL1Packet node is missing. Output related to this node will be empty" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
    }
    
    if (gl1PacketInfo)
    {
    boost::dynamic_bitset<> triggers(32, gl1PacketInfo->getTriggerInput());
    bool hasTrigger = false;
      
    for(uint32_t iTrg = 0; iTrg < 32; iTrg++)
    {
    if (iTrg != 16 && iTrg != 17 && iTrg != 18 && iTrg != 19)
    {
    continue;
    }
	
    if (triggers.test(iTrg))
    {
    hasTrigger = true;
    break;
    }
    }   
    if (!hasTrigger) {
    std::cout << "jet trigger not found" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
    }
    } 
    }
  */
  //calorimeter towers
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  /*
    TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
    TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
    RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  */
  if(!towersEM3 || !towersIH3 || !towersOH3){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower node "
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(!tower_geomEM || !tower_geom || !tower_geomOH){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower geometry "
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::vector<PseudoJet> particles_reco;
  particles_reco.clear();
  float totalPx_reco = 0;
  float totalPy_reco = 0;
  float totalPz_reco = 0;
  float totalE_reco = 0;

  /*  
  std::vector<PseudoJet> particles_truth;
  particles_truth.clear();
  float pxpart = 0;
  float pypart = 0;
  float pzpart = 0;
  float totalE_truth = 0;
  */

  //get the event centrality/impact parameter from HIJING
  //m_centrality =  cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
  m_centrality = (int)(100*cent_node->get_centile(CentralityInfo::PROP::mbd_NS));
  m_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);

  //get reco jets
  m_nJet = 0;
  
  double leading_pT = -999;
  double subleading_pT = -999;
  
  //  for (const auto& recoJet : matchedRecoJets) 
  // {
  //  if(recoJet->get_pt() < 5 || recoJet->size_comp() < 10 || recoJet->get_eta() < m_etaRange.first || recoJet->get_eta() > m_etaRange.second) continue; // to remove noise jets
  //std::cout << "working on jet " << m_nJet << " which has " << recoJet->size_comp() << " components" << std::endl;
  for (size_t i = 0; i < matchedRecoJets.size(); ++i) {
    const auto& recoJet = matchedRecoJets[i];
    //    const auto& truthJet = matchedTruthJets[i];  
    m_nJet++;
    
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
    energy_spectrum->Fill(m_e[m_e.size()-1]);
    m_pt.push_back(recoJet->get_pt());
    spectrum->Fill(m_pt[m_pt.size()-1]);
    eta_pt->Fill(m_pt[m_pt.size()-1],m_eta[m_eta.size()-1]);
    phi_pt->Fill(m_pt[m_pt.size()-1],m_phi[m_phi.size()-1]);
    M_pt->Fill(m_pt[m_pt.size()-1],recoJet->get_mass());
    
    int nconst = 0;
    // Use const_cast to call the non-const version of get_comp_vec()
    for (auto comp : const_cast<Jet*>(recoJet)->get_comp_vec())
      //      for (auto comp: recoJet->get_comp_vec())
      {
	//std::cout << "jet " << m_nJet << " consituent " << nconst << std::endl;
	TowerInfo *tower;
	nconst++;
	unsigned int channel = comp.second;
	// std::cout << "Number of components in jet: " << recoJet->get_comp_vec().size() << std::endl;
	//  std::cout << "Component type: " << comp.first << std::endl;
	if (comp.first == 26)
	  {
	    //	      std::cout << "Processing component ihcal" << std::endl;
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
	else if (comp.first == 27)
	  {
	    //	      std::cout << "Processing component ohcal" << std::endl;
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
	else if (comp.first == 25)
	  {
	    //	      std::cout << "Processing component emcal" << std::endl;
	    tower = towersEM3->get_tower_at_channel(channel);
	    if(!tower || !tower_geom)
	      {
		//  std::cout << "skipping event if needed" << std::endl;
		continue;
	      }
	      
	    unsigned int calokey = towersEM3->encode_key(channel);
	    int ieta = towersEM3->getTowerEtaBin(calokey);
	    int iphi = towersEM3->getTowerPhiBin(calokey);
	    //  std::cout << "calokey: " << calokey << ", ieta: " << ieta << ", iphi: " << iphi << ", channel: " << channel << std::endl;
	    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
	    // float /*UE*/ = background->get_/*UE*/(0).at(ieta);
	    auto tower_geometry_object = tower_geom->get_tower_geometry(key);
	    if (!tower_geometry_object) {
	      //		std::cerr << "Error: Tower geometry object is null for key " << key << std::endl;
	      continue; // Skip this component and go to the next one
	    }
	    /*	      float tower_phi = tower_geom->get_tower_geometry->get_phi();
		      std::cout << "tower_phi: " << tower_phi << std::endl;
		      float tower_eta = tower_geom->get_tower_geometry->get_eta();
		      std::cout << "tower_eta: " << tower_eta << std::endl;
	    */
	    // /*UE*/ = /*UE*/ * (1 + 2 * background_v2 * cos(2 * (tower_phi - background_Psi2)));
	    float tower_phi = tower_geometry_object->get_phi();
	    //  std::cout << "tower_phi: " << tower_phi << std::endl;
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
	      //	std::cerr << "Error: Invalid PseudoJet parameters." << std::endl;
	      continue; // Skip this component and go to the next one
	    }	      
	    particles_reco.push_back( PseudoJet( pxTmp, pyTmp, pzTmp, totalE_reco) );
	    // std::cout << "Pushing back PseudoJet with px: " << pxTmp
	    //	<< ", py: " << pyTmp
	    //	<< ", pz: " << pzTmp
	    //	<< ", E: " << totalE_reco << std::endl;

	  }
      }
    //Insert jenn's analysis                                                                                       
    double radius[5] = {0.05, 0.1, 0.2, 0.4, 0.6}; // jet radius                                                                                     
    //	double pseudorapidity = -999.; // pseudorapidity
    double theta_sj = -1.; // delta radius (value describes an unachievable value                                              
    double z_sj = -1.; // delta radius (value describes an unachievable value)
      
    JetDefinition jetDefAKT_R01( antikt_algorithm, radius[1]);
    JetDefinition jetDefAKT_R04( antikt_algorithm, radius[3]);
    JetDefinition jetDefCA(cambridge_algorithm, radius[3]);
      
    ClusterSequence clustSeq_R04( particles_reco, jetDefAKT_R04 );
    std::vector<PseudoJet> sortedJets_R04 = sorted_by_pt( clustSeq_R04.inclusive_jets() );
      
    for (int j = 0; j < (int)sortedJets_R04.size(); j++) {
      PseudoJet jet_reco = sortedJets_R04.at(j);
      if(fabs(jet_reco.eta()) > 0.7) //0.6 originally
	continue;
      /*	
      // Get and print constituents of jet_truth
      std::vector<PseudoJet> constituents = jet_truth.constituents();
      std::cout << “Constituents of jet ” << j << “: “;
      for (size_t i = 0; i < constituents.size(); ++i) {
      std::cout << “Pt: ” << constituents[i].pt() << “, Eta: ” << constituents[i].eta()
      << “, Phi: ” << constituents[i].phi() << “; “;
      }
      std::cout << std::endl;
      */
      ClusterSequence clustSeq_R01_con(jet_reco.constituents() , jetDefAKT_R01 );
      std:: vector<PseudoJet> sortedJets_R01_con = sorted_by_pt( clustSeq_R01_con.inclusive_jets() );
      if (sortedJets_R01_con.size() < 2){
	continue;
      }
      PseudoJet sj1 = sortedJets_R01_con.at(0);
      PseudoJet sj2 = sortedJets_R01_con.at(1);
      if (sj1.pt() < 0 || sj2.pt() < 0 )
	continue;
      theta_sj = sj1.delta_R(sj2);
      z_sj = sj2.pt()/(sj2.pt()+sj1.pt());
      double z_cut = 0.1;
      double beta = 0.0;
	
      // 30 to 35 pT
      if (recoJet->get_pt() > 30 && recoJet->get_pt() < 35 ){
	ClusterSequence clustSeqCA(jet_reco.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	// SoftDrop parameters
	contrib::SoftDrop sd(beta, z_cut);
	//! get subjets
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_30_35_data->Fill(z_sj);
	  _h_R04_theta_sj_30_35_data->Fill(theta_sj);
	}
	// Apply SoftDrop to the jet
	PseudoJet sd_jet = sd(jet_reco);
	// std::cout << "ran sd" << std::endl;
	if (sd_jet == 0)
	  continue;
	double delta_R_subjets = sd_jet.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_jet.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_30_35_data->Fill(z_subjets);
	_h_R04_theta_g_30_35_data->Fill(delta_R_subjets);
	correlation_theta_30_35_data->Fill(delta_R_subjets, theta_sj);                                                                                                                         
	correlation_z_30_35_data->Fill(z_subjets, z_sj); 
	// e.g., skip this jet or perform alternative analysis      
      }
      ////////////////////////////////////////////////////////////////////////////////// End loop over 30 -> 35 GeV
      // 35 to 40 pT                                                                                                                                    

      if (recoJet->get_pt() > 35 && recoJet->get_pt() < 40 ){
	  
	ClusterSequence clustSeqCA(jet_reco.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	// SoftDrop parameters                                                                                                                                                                
	contrib::SoftDrop sd2(beta, z_cut);
	//! get subjets                                                                                                                                                                       
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_35_40_data->Fill(z_sj);
	  _h_R04_theta_sj_35_40_data->Fill(theta_sj);
	}
	// Apply SoftDrop to the jet                                                                                                                                                          
	PseudoJet sd_jet2 = sd2(jet_reco);
	if (sd_jet2 == 0)
	  continue;
	double delta_R_subjets = sd_jet2.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_jet2.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_35_40_data->Fill(z_subjets);
	_h_R04_theta_g_35_40_data->Fill(delta_R_subjets);
	correlation_theta_35_40_data->Fill(delta_R_subjets, theta_sj);
	correlation_z_35_40_data->Fill(z_subjets, z_sj);
      }
      ////////////////////////////////////////////////////////////////////// End loop over 35 -> 40 GeV 
      // 40 to 45 pT                                                                                                                                         

      if (recoJet->get_pt() > 40 && recoJet->get_pt() < 45 ){
	  
	ClusterSequence clustSeqCA(jet_reco.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	contrib::SoftDrop sd4(beta, z_cut);
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_40_45_data->Fill(z_sj);
	  _h_R04_theta_sj_40_45_data->Fill(theta_sj);
	}
	PseudoJet sd_jet4 = sd4(jet_reco);
	if (sd_jet4 == 0)
	  continue;
	double delta_R_subjets = sd_jet4.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_jet4.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_40_45_data->Fill(z_subjets);
	_h_R04_theta_g_40_45_data->Fill(delta_R_subjets);
	correlation_theta_40_45_data->Fill(delta_R_subjets, theta_sj);
	correlation_z_40_45_data->Fill(z_subjets, z_sj);
      }
      //////////////////////////////////////////////////////////////////////// End loop over 40 -> 45 GeV
      // 45 pT or more                                                                                                                                                                       
      if (recoJet->get_pt() > 45 ){
	  
	ClusterSequence clustSeqCA(jet_reco.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	contrib::SoftDrop sd6(beta, z_cut);
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_45_data->Fill(z_sj);
	  _h_R04_theta_sj_45_data->Fill(theta_sj);
	}
	PseudoJet sd_jet6 = sd6(jet_reco);
	if (sd_jet6 == 0)
	  continue;
	double delta_R_subjets = sd_jet6.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_jet6.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_45_data->Fill(z_subjets);
	_h_R04_theta_g_45_data->Fill(delta_R_subjets);
	correlation_theta_45_data->Fill(delta_R_subjets, theta_sj);
	correlation_z_45_data->Fill(z_subjets, z_sj);
      }
      /////////////////////////////////////////////////////////////////////// End loop over 45+ GeV
    }
      
    if(subleading_pT > 0 && leading_pT > 0) 
      {
	jet_AN->Fill(subleading_pT/leading_pT);
      }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // get truth reco 
  m_nTruthJet = 0;
  
  double leadingTruth_pT = -999;
  double subleadingTruth_pT = -999;
  for (size_t i = 0; i < matchedTruthJets.size(); ++i) {
    //    const auto& recoJet = matchedRecoJets[i];
    const auto& truthJet = matchedTruthJets[i];
    
    //  for (const auto& truthJet : matchedTruthJets)
    // {
    // if(truthJet->get_pt() < 5 /* || truthJet->size_comp() < 2*/ || truthJet->get_eta() < m_etaRange.first || truthJet->get_eta() > m_etaRange.second)
    //	{
    //	  //	std::cout << "noise truthJets removed" << std::endl;
    //	  continue; // to remove noise truthJet	
    //	}
      
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
    
    // Correct phi to be within [-pi, pi]
    float& phi = m_phi_truth.back();
    //  while (phi < -M_PI) phi += 2 * M_PI;
    //  while (phi > M_PI) phi -= 2 * M_PI;
      
    // Fill histograms
    energy_spectrum_truth->Fill(truthJet->get_e());
    spectrum_truth->Fill(truthJetPt);
    eta_pt_truth->Fill(truthJetPt, truthJet->get_eta());
    phi_pt_truth->Fill(truthJetPt, phi);
    M_pt_truth->Fill(truthJetPt, truthJet->get_mass());
    
    // Create pseudojets for substructure analysis
    float tower_phi2 = phi;
    float tower_eta2 = truthJet->get_eta();
    double ETmp2 = truthJet->get_e();
    double pt2 = ETmp2 / cosh(tower_eta2);
    double pxTmp2 = pt2 * cos(tower_phi2);
    double pyTmp2 = pt2 * sin(tower_phi2);
    double pzTmp2 = pt2 * sinh(tower_eta2);
    
    // Check for NaN values before creating the PseudoJet
    if (std::isnan(pxTmp2) || std::isnan(pyTmp2) || std::isnan(pzTmp2) || std::isnan(ETmp2)) {
      //  std::cerr << "Error: Invalid PseudoJet parameters." << std::endl;
      continue;
    }
    std::vector<PseudoJet> particles_truth;
    particles_truth.clear();
    for (auto comp : const_cast<Jet*>(truthJet)->get_comp_vec()) {
      PHG4TruthInfoContainer::Map particlemap = truthinfo->GetMap();

      for (auto iter = particlemap.begin(); iter != particlemap.end(); ++iter) {
        if (static_cast<unsigned int>(iter->first) == comp.second) {
	  PHG4Particle* part = iter->second;

	  if (part == nullptr) {
	    std::cout << "Error: Particle not found for component " << comp.second << std::endl;
	    continue;
	  }

	  float Epart = part->get_e();
	  float pxpart = part->get_px();
	  float pypart = part->get_py();
	  float pzpart = part->get_pz();

	  particles_truth.push_back(PseudoJet(pxpart, pypart, pzpart, Epart));
        }
      }
    }
    //Insert jenn's analysis                                                                                       
    double radius[5] = {0.05, 0.1, 0.2, 0.4, 0.6}; // jet radius                                                                                     
    //	double pseudorapidity = -999.; // pseudorapidity
    double theta_sj = -1.; // delta radius (value describes an unachievable value                                              
    double z_sj = -1.; // delta radius (value describes an unachievable value)
      
    JetDefinition jetDefAKT_R01( antikt_algorithm, radius[1]);
    JetDefinition jetDefAKT_R04( antikt_algorithm, radius[3]);
    JetDefinition jetDefCA(cambridge_algorithm, radius[3]);
    
    ClusterSequence clustSeq_R04_truth( particles_truth, jetDefAKT_R04 );
    std::vector<PseudoJet> sortedJets_R04_truth = sorted_by_pt( clustSeq_R04_truth.inclusive_jets() );
      
    for (int j = 0; j < (int)sortedJets_R04_truth.size(); j++) {
      PseudoJet jet_truth = sortedJets_R04_truth.at(j);
      if(fabs(jet_truth.eta()) > 0.7) //0.6 originally                                                                                                                                                     
	continue;
	
      // Get and print constituents of jet_truth
      std::vector<PseudoJet> constituents = jet_truth.constituents();
      // std::cout << “Constituents of jet ” << j << “: “;
      for (size_t i = 0; i < constituents.size(); ++i) {
	std::cout << "pt:" << constituents[i].pt() << " Eta: " << constituents[i].eta()
		  << ", Phi: " << constituents[i].phi() << "; ";
      }
      std::cout << std::endl;

      ClusterSequence clustSeq_R01_con(constituents, jetDefAKT_R01 );
      std:: vector<PseudoJet> sortedJets_R01_con = sorted_by_pt( clustSeq_R01_con.inclusive_jets() );
      if (sortedJets_R01_con.size() < 2){
	continue;
      }
      PseudoJet sj1 = sortedJets_R01_con.at(0);
      PseudoJet sj2 = sortedJets_R01_con.at(1);
      if (sj1.pt() < 0 || sj2.pt() < 0 )
	continue;
      theta_sj = sj1.delta_R(sj2);
      z_sj = sj2.pt()/(sj2.pt()+sj1.pt());
      double z_cut = 0.1;
      double beta = 0.0;

      // 30 to 35 pT                                                                                                                                              
                                       
      if (truthJet->get_pt() > 30 && truthJet->get_pt() < 35 ){
	ClusterSequence clustSeqCA(jet_truth.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	// SoftDrop parameters                                                                                                                                  

	contrib::SoftDrop sd1(beta, z_cut);
	//! get subjets                                                                                                                                      
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_30_35_truth->Fill(z_sj);
	  _h_R04_theta_sj_30_35_truth->Fill(theta_sj);
	}
	// Apply SoftDrop to the jet                                                                                  
	PseudoJet sd_truth1 = sd1(jet_truth);
	// std::cout << "ran sd" << std::endl;
	if (sd_truth1 == 0)
	  continue;
	double delta_R_subjets = sd_truth1.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_truth1.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_30_35_truth->Fill(z_subjets);
	_h_R04_theta_g_30_35_truth->Fill(delta_R_subjets);
	correlation_theta_30_35_truth->Fill(delta_R_subjets, theta_sj);
	correlation_z_30_35_truth->Fill(z_subjets, z_sj);
	// e.g., skip this jet or perform alternative analysis 
      }
      ////////////////////////////////////////////////////////////////////////////////// End loop pver 30 -> 35 GeV
      // 35 to 40 pT                                                                                                                                               

      if (truthJet->get_pt() > 35 && truthJet->get_pt() < 40 ){
	
	ClusterSequence clustSeqCA(jet_truth.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	// SoftDrop parameters                                                                                              
	contrib::SoftDrop sd3(beta, z_cut);
	//! get subjets                                                                                             
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_35_40_truth->Fill(z_sj);
	  _h_R04_theta_sj_35_40_truth->Fill(theta_sj);
	}
	// Apply SoftDrop to the jet                                           
	PseudoJet sd_truth3 = sd3(jet_truth);
	if (sd_truth3 == 0)
	  continue;
	double delta_R_subjets = sd_truth3.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_truth3.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_35_40_truth->Fill(z_subjets);
	_h_R04_theta_g_35_40_truth->Fill(delta_R_subjets);
	correlation_theta_35_40_truth->Fill(delta_R_subjets, theta_sj);
	correlation_z_35_40_truth->Fill(z_subjets, z_sj);
      }
      ////////////////////////////////////////////////////////////////////////////////// End loop over 35 -> 40 GeV 
      // 40 to 45 pT                                                                               
      if (truthJet->get_pt() > 40 && truthJet->get_pt() < 45){
	
	ClusterSequence clustSeqCA(jet_truth.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	contrib::SoftDrop sd5(beta, z_cut);
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_40_45_truth->Fill(z_sj);
	  _h_R04_theta_sj_40_45_truth->Fill(theta_sj);
	}
	PseudoJet sd_truth5 = sd5(jet_truth);
	if (sd_truth5 == 0)
	  continue;
	double delta_R_subjets = sd_truth5.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_truth5.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_40_45_truth->Fill(z_subjets);
	_h_R04_theta_g_40_45_truth->Fill(delta_R_subjets);
	correlation_theta_40_45_truth->Fill(delta_R_subjets, theta_sj);
	correlation_z_40_45_truth->Fill(z_subjets, z_sj);
      }
      ///////////////////////////////////////////////////////////////////////////////////// End loop over 40 -> 45 GeV
      // 45 pT or more                                                                                                                     
      if (truthJet->get_pt() > 45){
	
	ClusterSequence clustSeqCA(jet_truth.constituents(), jetDefCA);
	std::vector<PseudoJet> cambridgeJets = sorted_by_pt(clustSeqCA.inclusive_jets());
	contrib::SoftDrop sd7(beta, z_cut);
	if (!isnan(theta_sj) && !isnan(z_sj) && !isinf(theta_sj) && !isinf(z_sj)){
	  _h_R04_z_sj_45_truth->Fill(z_sj);
	  _h_R04_theta_sj_45_truth->Fill(theta_sj);
	}
	PseudoJet sd_truth7 = sd7(jet_truth);
	if (sd_truth7 == 0)
	  continue;
	double delta_R_subjets = sd_truth7.structure_of<contrib::SoftDrop>().delta_R();
	double z_subjets = sd_truth7.structure_of<contrib::SoftDrop>().symmetry();
	_h_R04_z_g_45_truth->Fill(z_subjets);
	_h_R04_theta_g_45_truth->Fill(delta_R_subjets);
	correlation_theta_45_truth->Fill(delta_R_subjets, theta_sj);
	correlation_z_45_truth->Fill(z_subjets, z_sj);
      }
      ///////////////////////////////////////////////////////////////////////////////////// End loop over 45+ GeV
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int JetResoStudy::ResetEvent(PHCompositeNode* /*topNode*/)
{
  //std::cout << "JetResoStudy::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
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
int JetResoStudy::EndRun(const int runnumber)
{
  std::cout << "JetResoStudy::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetResoStudy::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "JetResoStudy::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  hDeltaR->Write();
  hRecoJetPt->Write();
  hTruthJetPt->Write();
  hFakeRecoJetPt->Write();
  hFakeTruthJetPt->Write();
  hEfficiency->Write();    
  hPurity->Write();
  hFakeRate->Write();
  hMomentumResolution->Write();
  hScaleFactor->Write();

  double averagePurity = hPurity->GetMean();
  std::cout << "Average Purity over all events: " << averagePurity << std::endl;

  double averageFakeRate = hFakeRate->GetMean();
  std::cout << "Average Fake Rate over all events: " << averageFakeRate << std::endl;

  double averageMomentumReso = hMomentumResolution->GetMean();
  std::cout << "Average Momentum Resolution over all event: " << averageMomentumReso << std::endl;
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
    
  //std::cout << "number of jets:" << _h_R04_z_sj_35_40->Integral() << std::endl;
  //  5-10 hists data                                                                                                                                                                 
  _h_R04_z_sj_30_35_data->Write();
  _h_R04_theta_sj_30_35_data->Write();
  //  10-20 hists                                                                                                                                                                                
  _h_R04_z_sj_35_40_data->Write();
  _h_R04_theta_sj_35_40_data->Write();
  //  20-30 hists                                                                                                                                                                                
  _h_R04_z_sj_40_45_data->Write();
  _h_R04_theta_sj_40_45_data->Write();
  //  30 hists                                                                                                                                                                                
  _h_R04_z_sj_45_data->Write();
  _h_R04_theta_sj_45_data->Write();
  //truth anti-kT
  //  5-10 hists truth                                                                              
  _h_R04_z_sj_30_35_truth->Write();
  _h_R04_theta_sj_30_35_truth->Write();
  //  10-20 hists                                                                                                                                                                                 
  _h_R04_z_sj_35_40_truth->Write();
  _h_R04_theta_sj_35_40_truth->Write();
  //  20-30 hists                                                                                                                                                                                 
  _h_R04_z_sj_40_45_truth->Write();
  _h_R04_theta_sj_40_45_truth->Write();
  //  30 hists                                                                                                                                                                                    
  _h_R04_z_sj_45_truth->Write();
  _h_R04_theta_sj_45_truth->Write();
  //SoftDrop  zcut=0.1 data
  //  5-10 hists                                                                                                                                                                                 
  _h_R04_z_g_30_35_data->Write();
  _h_R04_theta_g_30_35_data->Write();                                                                                                                                                                       
  //  10-20 hists                                                                                                                                                                                
  _h_R04_z_g_35_40_data->Write();
  _h_R04_theta_g_35_40_data->Write();
  //  20-30 hists                                                                                                                                                                                 
  _h_R04_z_g_40_45_data->Write();
  _h_R04_theta_g_40_45_data->Write();
  //  30 hists                                                                                                                                                                                 
  _h_R04_z_g_45_data->Write();
  _h_R04_theta_g_45_data->Write();
  //SoftDrop  zcut=0.1 truth                                                                                                                                                                        
  //  5-10 hists                                                                                                                                                                                  
  _h_R04_z_g_30_35_truth->Write();
  _h_R04_theta_g_30_35_truth->Write();
  //  10-20 hists                                                                                                                                                                                 
  _h_R04_z_g_35_40_truth->Write();
  _h_R04_theta_g_35_40_truth->Write();
  //  20-30 hists                                                                                                                                                                                 
  _h_R04_z_g_40_45_truth->Write();
  _h_R04_theta_g_40_45_truth->Write();
  //  30 hists                                                                                                                                                                                    
  _h_R04_z_g_45_truth->Write();
  _h_R04_theta_g_45_truth->Write();
  // Correlation plots
  correlation_theta_30_35_data->Write();
  correlation_theta_35_40_data->Write();
  correlation_theta_40_45_data->Write();
  correlation_theta_45_data->Write();
  correlation_theta_30_35_truth->Write();
  correlation_theta_35_40_truth->Write();
  correlation_theta_40_45_truth->Write();
  correlation_theta_45_truth->Write();

  correlation_z_30_35_data->Write();
  correlation_z_35_40_data->Write();
  correlation_z_40_45_data->Write();
  correlation_z_45_data->Write();
  correlation_z_30_35_truth->Write();
  correlation_z_35_40_truth->Write();
  correlation_z_40_45_truth->Write();
  correlation_z_45_truth->Write();
 
  std::cout << "JetResoStudy::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetResoStudy::Reset(PHCompositeNode* /*topNode*/)
{
  std::cout << "JetResoStudy::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetResoStudy::Print(const std::string &what) const
{
  std::cout << "JetResoStudy::Print(const std::string &what) const Printing info for " << what << std::endl;
}
