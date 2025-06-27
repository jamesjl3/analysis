// Rewritten JetMatchingSubjets.cc to include jet matching and apply matched jets to subjet analysis with TTree support
#include "JetMatchingSubjets.h"
#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetMap.h>
#include <jetbase/Jet.h>
#include <centrality/CentralityInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <jetbackground/TowerBackground.h>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h> 
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <map>

using namespace std;
using namespace fastjet;

JetMatchingSubjets::JetMatchingSubjets(const std::string& recojetname,
                                       const std::string& truthjetname,
                                       const std::string& outputfilename)
  : SubsysReco("JetMatchingSubjets_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1.1, 1.1)
  , m_ptRange(5, 100)
  , m_doTruthJets(0)
{}

JetMatchingSubjets::~JetMatchingSubjets() {}

TTree* m_T = nullptr;
int m_event = -1;
float m_centrality = 0;
float m_impactparam = 0;
std::vector<float> m_pt, m_eta, m_phi;
std::vector<float> m_pt_truth, m_eta_truth, m_phi_truth;

TH1F* _h_R04_z_sj_10 = nullptr;
TH1F* _h_R04_theta_sj_10 = nullptr;
TH1F* _h_R04_z_g_10_01 = nullptr;
TH1F* _h_R04_theta_g_10_01 = nullptr;

TH1F* _h_R04_truth_z_sj_10 = nullptr;
TH1F* _h_R04_truth_theta_sj_10 = nullptr;
TH1F* _h_R04_truth_z_g_10_01 = nullptr;
TH1F* _h_R04_truth_theta_g_10_01 = nullptr;

float deltaR(Jet* a, Jet* b) {
  float deta = a->get_eta() - b->get_eta();
    float dphi = TVector2::Phi_mpi_pi(a->get_phi() - b->get_phi());
    return sqrt(deta * deta + dphi * dphi);
}
void JetMatchingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax)
{
  recoToTruth.clear();
  truthToReco.clear();

  std::vector<std::tuple<float, Jet*, Jet*>> pairs;

  for (auto reco : *recoJets)
  {
    // Apply reco cuts BEFORE matching:
    if (reco->get_pt() < 5.0) continue;  // Use your desired reco cut (30.0 here as example)
    if (fabs(reco->get_eta()) > 0.7) continue;  // Example eta cut

    for (auto truth : *truthJets)
    {
      // Apply truth cuts BEFORE matching:
      if (truth->get_pt() < 10.0) continue;  // Same cut
      if (fabs(truth->get_eta()) > 0.7) continue;

      float dR = deltaR(reco, truth);
      if (dR < dRMax)
      {
        pairs.emplace_back(dR, reco, truth);
      }
    }
  }

  std::sort(pairs.begin(), pairs.end());

  std::set<Jet*> matchedReco, matchedTruth;
  for (const auto& tup : pairs)
  {
    Jet* reco = std::get<1>(tup);
    Jet* truth = std::get<2>(tup);
    if (matchedReco.count(reco) == 0 && matchedTruth.count(truth) == 0)
    {
      recoToTruth[reco] = truth;
      truthToReco[truth] = reco;
      matchedReco.insert(reco);
      matchedTruth.insert(truth);
    }
  }
}
std::vector<PseudoJet> BuildTruthPseudoJets(Jet* truthJet, PHG4TruthInfoContainer* truthInfo)
{
    std::vector<PseudoJet> particles;
    for (auto comp : truthJet->get_comp_vec())
    {
        unsigned int truth_id = comp.second;

        PHG4Particle *particle = truthInfo->GetParticle(truth_id);
        if (!particle) continue;

        float px = particle->get_px();
        float py = particle->get_py();
        float pz = particle->get_pz();
        float e  = particle->get_e();

        particles.emplace_back(px, py, pz, e);
    }

    std::cout << "DEBUG: Built " << particles.size() << " PseudoJets from truthJet with pt = "
              << truthJet->get_pt() << std::endl;

    return particles;
}

std::vector<PseudoJet> BuildPseudoJets(Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh, RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH, TowerBackground* bg, float v2, float psi2, bool doUnsub) {
  std::vector<PseudoJet> particles;
  for (auto comp : jet->get_comp_vec()) {
    TowerInfo* tower = nullptr;
    unsigned int ch = comp.second;
    float eta = 0, phi = 0, UE = 0;

    if (comp.first == 14 || comp.first == 29) {
      tower = em->get_tower_at_channel(ch);
      if (!tower || !geomEM) continue;
      auto calokey = em->encode_key(ch);
      int ieta = em->getTowerEtaBin(calokey);
      int iphi = em->getTowerPhiBin(calokey);
      auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALIN, ieta, iphi);
      eta = geomEM->get_tower_geometry(key)->get_eta();
      phi = geomEM->get_tower_geometry(key)->get_phi();
      UE = bg->get_UE(0).at(ieta);
    } else if (comp.first == 15 || comp.first == 30) {
      tower = ih->get_tower_at_channel(ch);
      if (!tower || !geomEM) continue;
      auto calokey = ih->encode_key(ch);
      int ieta = ih->getTowerEtaBin(calokey);
      int iphi = ih->getTowerPhiBin(calokey);
      auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALIN, ieta, iphi);
      eta = geomEM->get_tower_geometry(key)->get_eta();
      phi = geomEM->get_tower_geometry(key)->get_phi();
      UE = bg->get_UE(1).at(ieta);
    } else if (comp.first == 16 || comp.first == 31) {
      tower = oh->get_tower_at_channel(ch);
      if (!tower || !geomOH) continue;
      auto calokey = oh->encode_key(ch);
      int ieta = oh->getTowerEtaBin(calokey);
      int iphi = oh->getTowerPhiBin(calokey);
      auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALOUT, ieta, iphi);
      eta = geomOH->get_tower_geometry(key)->get_eta();
      phi = geomOH->get_tower_geometry(key)->get_phi();
      UE = bg->get_UE(2).at(ieta);
    } else continue;

    UE *= (1 + 2 * v2 * cos(2 * (phi - psi2)));
    float energy = tower->get_energy() + UE;
    if (doUnsub) energy += UE;
    float pt = energy / cosh(eta);
    float px = pt * cos(phi);
    float py = pt * sin(phi);
    float pz = pt * sinh(eta);
    particles.emplace_back(px, py, pz, energy);
  }
  return particles;
}
void JetMatchingSubjets::AnalyzeMatchedJets(JetContainer* recoJets,
                                            TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                                            RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
                                            TowerBackground* bg, float v2, float psi2)
{
  JetDefinition jetDefAKT_R04(antikt_algorithm, 0.4);
  JetDefinition jetDefAKT_R01(antikt_algorithm, 0.1);

  for (auto& pair : recoToTruth)
  {
    Jet* recoJet = pair.first;
    //    Jet* truthJet = pair.second;

    if (recoJet->get_pt() < 5.0) continue;

    m_pt.push_back(recoJet->get_pt());
    m_eta.push_back(recoJet->get_eta());
    m_phi.push_back(recoJet->get_phi());

    auto particles = BuildPseudoJets(recoJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false);
    ClusterSequence clustSeq(particles, jetDefAKT_R04);
    auto jets = sorted_by_pt(clustSeq.inclusive_jets());
    if (jets.empty()) continue;

    PseudoJet leading = jets[0];
    if (fabs(leading.eta()) > 0.7) continue;

    if (leading.pt() > 10) {
      // ---- Subjet clustering (R=0.1) ----
      ClusterSequence subClust(leading.constituents(), jetDefAKT_R01);
      auto subjets = sorted_by_pt(subClust.inclusive_jets());
      if (subjets.size() >= 2) {
        PseudoJet sj1 = subjets[0], sj2 = subjets[1];
        if (sj1.pt() >= 3 && sj2.pt() >= 3) {
          double theta_sj = sj1.delta_R(sj2);
          double z_sj = sj2.pt() / (sj1.pt() + sj2.pt());

          if (!std::isnan(theta_sj) && !std::isnan(z_sj)) {
            _h_R04_z_sj_10->Fill(z_sj);
            _h_R04_theta_sj_10->Fill(theta_sj);
          }
        }
      }

      // ---- SoftDrop grooming ----
      contrib::SoftDrop sd(0.0, 0.1);
      PseudoJet sd_jet = sd(leading);
      if (sd_jet != 0) {
        double zg = sd_jet.structure_of<contrib::SoftDrop>().symmetry();
        double rg = sd_jet.structure_of<contrib::SoftDrop>().delta_R();
        _h_R04_z_g_10_01->Fill(zg);
        _h_R04_theta_g_10_01->Fill(rg);
      }
    }
  }
}
/*
  void JetMatchingSubjets::AnalyzeTruthJets(JetContainer* truthJets,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
  TowerBackground* bg, float v2, float psi2)
  */
void JetMatchingSubjets::AnalyzeTruthJets(JetContainer* truthJets,
                                          PHG4TruthInfoContainer* truthInfo)

{
  JetDefinition jetDefAKT_R04(antikt_algorithm, 0.4);
  JetDefinition jetDefAKT_R01(antikt_algorithm, 0.1);

  // Loop over all matched truth jets (truthToReco map):
  for (auto& pair : truthToReco)
  {
    Jet* truthJet = pair.first;
    // Jet* recoJet  = pair.second;
    if (truthJet->get_pt() < 10.0) continue;

    m_pt_truth.push_back(truthJet->get_pt());
    m_eta_truth.push_back(truthJet->get_eta());
    m_phi_truth.push_back(truthJet->get_phi());

    // Build pseudo-jets from truth jet constituents:
    //auto particles = BuildPseudoJets(truthJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false);
    auto particles = BuildTruthPseudoJets(truthJet, truthInfo);
 
    std::cout << "DEBUG: truthJet pt = " << truthJet->get_pt()
	      << ", eta = " << truthJet->get_eta()
	      << ", N constituents = " << truthJet->get_comp_vec().size() << std::endl;
    
    ClusterSequence clustSeq(particles, jetDefAKT_R04);
    auto jets = sorted_by_pt(clustSeq.inclusive_jets());
    if (jets.empty()) continue;

    PseudoJet leading = jets[0];
    if (fabs(leading.eta()) > 0.7) continue;

    if (leading.pt() > 10)
    {
      // ---- Subjet clustering (R=0.1) ----
      ClusterSequence subClust(leading.constituents(), jetDefAKT_R01);
      auto subjets = sorted_by_pt(subClust.inclusive_jets());
      if (subjets.size() >= 2)
      {
        PseudoJet sj1 = subjets[0], sj2 = subjets[1];
        if (sj1.pt() >= 3 && sj2.pt() >= 3)
        {
          double theta_sj = sj1.delta_R(sj2);
          double z_sj = sj2.pt() / (sj1.pt() + sj2.pt());

          if (!std::isnan(theta_sj) && !std::isnan(z_sj))
          {
            _h_R04_truth_z_sj_10->Fill(z_sj);
            _h_R04_truth_theta_sj_10->Fill(theta_sj);
          }
        }
      }

      // ---- SoftDrop grooming ----
      contrib::SoftDrop sd(0.0, 0.1);
      PseudoJet sd_jet = sd(leading);
      if (sd_jet != 0)
      {
        double zg = sd_jet.structure_of<contrib::SoftDrop>().symmetry();
        double rg = sd_jet.structure_of<contrib::SoftDrop>().delta_R();
        _h_R04_truth_z_g_10_01->Fill(zg);
        _h_R04_truth_theta_g_10_01->Fill(rg);
      }
    }
  }
}
int JetMatchingSubjets::process_event(PHCompositeNode* topNode) {
  ++m_event;
  topNode->print(); // or printTree(topNode);

  m_pt.clear(); m_eta.clear(); m_phi.clear();
  m_pt_truth.clear(); m_eta_truth.clear(); m_phi_truth.clear();
  
  recoToTruth.clear();
  truthToReco.clear();

  std::cout << "DEBUG: Attempting to get truth jets from node: " << m_truthJetName << std::endl;

  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  JetContainer* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  TowerBackground* bg = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

  if (!jets || !jetsMC || !towersEM3 || !towersIH3 || !towersOH3 || !geomEM || !geomOH || !bg || !cent_node){
    std::cerr << "Missing required node(s):" << std::endl;
    if (!jets) std::cerr << " - JetContainer " << m_recoJetName << std::endl;
    if (!jetsMC) std::cerr << " - JetContainer " << m_truthJetName << std::endl;
    if (!towersEM3) std::cerr << " - TOWERINFO_CALIB_CEMC_RETOWER_SUB1" << std::endl;
    if (!towersIH3) std::cerr << " - TOWERINFO_CALIB_HCALIN_SUB1" << std::endl;
    if (!towersOH3) std::cerr << " - TOWERINFO_CALIB_HCALOUT_SUB1" << std::endl;
    if (!geomEM) std::cerr << " - TOWERGEOM_HCALIN" << std::endl;
    if (!geomOH) std::cerr << " - TOWERGEOM_HCALOUT" << std::endl;
    if (!bg) std::cerr << " - TowerInfoBackground_Sub2" << std::endl;
    if (!cent_node) std::cerr << " - CentralityInfo" << std::endl;
    
    return Fun4AllReturnCodes::ABORTRUN;
  }  
  m_centrality = cent_node->get_centile(CentralityInfo::PROP::bimp);
  m_impactparam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
  float v2 = bg->get_v2();
  float psi2 = bg->get_Psi2();

  // Perform matching first:
  MatchJets1to1(jets, jetsMC, 0.2);  // You can tune dRMax here

  std::cout << "Matched pairs: recoToTruth.size() = " << recoToTruth.size()
          << ", truthToReco.size() = " << truthToReco.size() << std::endl;

  // Then analyze based on the established match maps:
  AnalyzeMatchedJets(jets, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2);
  //AnalyzeTruthJets(jetsMC, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2);

  PHG4TruthInfoContainer *truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthInfo)
    {
      std::cerr << "Missing G4TruthInfo node!" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  AnalyzeTruthJets(jetsMC, truthInfo);  

  m_T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::Init(PHCompositeNode*) {
  
  std::cout << "DEBUG: Opening file with name = " << m_outputFileName << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  m_T = new TTree("T", "Jet Tree");
  m_T->Branch("event", &m_event, "event/I");
  m_T->Branch("cent", &m_centrality, "cent/F");
  m_T->Branch("b", &m_impactparam, "b/F");
  m_T->Branch("pt", &m_pt);
  m_T->Branch("eta", &m_eta);
  m_T->Branch("phi", &m_phi);
  m_T->Branch("pt_truth", &m_pt_truth);
  m_T->Branch("eta_truth", &m_eta_truth);
  m_T->Branch("phi_truth", &m_phi_truth);

  _h_R04_z_sj_10 = new TH1F("z_sj", "z_{sj};z;Entries", 20, 0, 0.5);
  _h_R04_theta_sj_10 = new TH1F("theta_sj", "#theta_{sj};#theta;Entries", 20, 0, 0.5);
  _h_R04_z_g_10_01 = new TH1F("z_g", "z_{g};z;Entries", 20, 0, 0.5);
  _h_R04_theta_g_10_01 = new TH1F("theta_g", "#theta_{g};#theta;Entries", 20, 0, 0.5);

  _h_R04_truth_z_sj_10 = new TH1F("truth_z_sj", "Truth z_{sj};z;Entries", 20, 0, 0.5);
  _h_R04_truth_theta_sj_10 = new TH1F("truth_theta_sj", "Truth #theta_{sj};#theta;Entries", 20, 0, 0.5);
  _h_R04_truth_z_g_10_01 = new TH1F("truth_z_g", "Truth z_{g};z;Entries", 20, 0, 0.5);
  _h_R04_truth_theta_g_10_01 = new TH1F("truth_theta_g", "Truth #theta_{g};#theta;Entries", 20, 0, 0.5);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetMatchingSubjets::End(PHCompositeNode*) {
  PHTFileServer::get().cd(m_outputFileName);
  
  m_T->Write();
  _h_R04_z_sj_10->Write();
  _h_R04_theta_sj_10->Write();
  _h_R04_z_g_10_01->Write();
  _h_R04_theta_g_10_01->Write();

  _h_R04_truth_z_sj_10->Write();
  _h_R04_truth_theta_sj_10->Write();
  _h_R04_truth_z_g_10_01->Write();
  _h_R04_truth_theta_g_10_01->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
  int JetMatchingSubjets::Reset(PHCompositeNode *topNode)
  {
    std::cout << "JetMatchingSubjets::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //____________________________________________________________________________..
  void JetMatchingSubjets::Print(const std::string &what) const
  {
     std::cout << "JetMatchingSubjets::Print(const std::string &what) const Printing info for " << what << std::endl;
  }
    

