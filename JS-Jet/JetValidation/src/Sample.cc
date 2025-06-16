#include "JetUnfoldingSubjets.h"

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
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>

using namespace std;
using namespace fastjet;

// -- GLOBAL pt bins --
std::vector<float> pt_bins = {30, 35, 40, 45, 50, 55, 60};

// -- FindPtBin --
int FindPtBin(float pt) {
  for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
    if (pt >= pt_bins[i] && pt < pt_bins[i + 1]) return i;
  }
  return -1;
}

// -- Helper deltaR --
float deltaR(Jet* a, Jet* b) {
  float deta = a->get_eta() - b->get_eta();
  float dphi = TVector2::Phi_mpi_pi(a->get_phi() - b->get_phi());
  return sqrt(deta * deta + dphi * dphi);
}

// -- Constructor / Destructor --
JetUnfoldingSubjets::JetUnfoldingSubjets(const std::string& recojetname,
                                         const std::string& truthjetname,
                                         const std::string& outputfilename)
  : SubsysReco("JetUnfoldingSubjets"),
    m_recoJetName(recojetname),
    m_truthJetName(truthjetname),
    m_outputFileName(outputfilename),
    m_etaRange(-0.7, 0.7),
    m_ptRange(30, 60),
    m_doTruthJets(1),
    m_T(nullptr),
    m_event(-1)
{}

JetUnfoldingSubjets::~JetUnfoldingSubjets() {}

// -- Jet Matching --
void JetUnfoldingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax) {
  recoToTruth.clear();
  truthToReco.clear();
  std::vector<std::tuple<float, Jet*, Jet*>> pairs;

  for (auto reco : *recoJets) {
    if (reco->get_pt() < 5.0) continue;
    if (fabs(reco->get_eta()) > 0.7) continue;

    for (auto truth : *truthJets) {
      if (truth->get_pt() < 30.0) continue;
      if (fabs(truth->get_eta()) > 0.7) continue;

      float dR = deltaR(reco, truth);
      if (dR < dRMax) {
        pairs.emplace_back(dR, reco, truth);
      }
    }
  }

  std::sort(pairs.begin(), pairs.end());
  std::set<Jet*> matchedReco, matchedTruth;

  for (const auto& tup : pairs) {
    Jet* reco = std::get<1>(tup);
    Jet* truth = std::get<2>(tup);

    if (matchedReco.count(reco) == 0 && matchedTruth.count(truth) == 0) {
      recoToTruth[reco] = truth;
      truthToReco[truth] = reco;
      matchedReco.insert(reco);
      matchedTruth.insert(truth);
    }
  }
}

// -- Build PseudoJets (reco) --
std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildPseudoJets(Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                                                                    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
                                                                    TowerBackground* bg, float v2, float psi2, bool doUnsub) {
  std::vector<fastjet::PseudoJet> particles;

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
    float energy = tower->get_energy();
    if (doUnsub) energy -= UE;

    float pt = energy / cosh(eta);
    float px = pt * cos(phi);
    float py = pt * sin(phi);
    float pz = pt * sinh(eta);

    particles.emplace_back(px, py, pz, energy);
  }

  return particles;
}

// -- Build Truth PseudoJets --
std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildTruthPseudoJets(Jet* truthJet, PHG4TruthInfoContainer* truthInfo) {
  std::vector<fastjet::PseudoJet> particles;

  for (auto comp : truthJet->get_comp_vec()) {
    unsigned int truth_id = comp.second;
    PHG4Particle *particle = truthInfo ? truthInfo->GetParticle(truth_id) : nullptr;
    if (!particle) continue;

    float px = particle->get_px();
    float py = particle->get_py();
    float pz = particle->get_pz();
    float e  = particle->get_e();

    particles.emplace_back(px, py, pz, e);
  }

  return particles;
}

// -- Analyze Matched Jets --
void JetUnfoldingSubjets::AnalyzeMatchedJets(JetContainer* recoJets,
                                             TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
                                             RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
                                             TowerBackground* bg, float v2, float psi2,
                                             PHG4TruthInfoContainer* truthInfo) {
  JetDefinition jetDefAKT_R04(antikt_algorithm, 0.4);
  JetDefinition jetDefAKT_R01(antikt_algorithm, 0.1);

  for (auto& pair : recoToTruth) {
    Jet* recoJet = pair.first;
    Jet* truthJet = pair.second;

    if (!recoJet || !truthJet) continue;
    if (recoJet->get_pt() < 5.0) continue;
    if (truthJet->get_pt() < 30.0) continue;
    if (fabs(recoJet->get_eta()) > 0.7 || fabs(truthJet->get_eta()) > 0.7) continue;

    m_pt.push_back(recoJet->get_pt());
    m_eta.push_back(recoJet->get_eta());
    m_phi.push_back(recoJet->get_phi());
    m_pt_truth.push_back(truthJet->get_pt());
    m_eta_truth.push_back(truthJet->get_eta());
    m_phi_truth.push_back(truthJet->get_phi());

    hRecoJetPtMatched->Fill(recoJet->get_pt());
    hTruthJetPtMatched->Fill(truthJet->get_pt());
    m_response1D->Fill(recoJet->get_pt(), truthJet->get_pt());

    // Subjets z_sj
    auto particles = BuildPseudoJets(recoJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false);
    ClusterSequence clustSeq(particles, jetDefAKT_R04);
    auto jets = sorted_by_pt(clustSeq.inclusive_jets());
    if (jets.empty()) continue;

    PseudoJet leading = jets[0];
    ClusterSequence subClust(leading.constituents(), jetDefAKT_R01);
    auto subjets = sorted_by_pt(subClust.inclusive_jets());
    if (subjets.size() >= 2 && subjets[0].pt() >= 3 && subjets[1].pt() >= 3) {
      double reco_z_sj = subjets[1].pt() / (subjets[0].pt() + subjets[1].pt());
      int bin = FindPtBin(truthJet->get_pt());
      if (bin >= 0) {
        m_hRecoZsjMatched[bin]->Fill(reco_z_sj);
        m_responseZsj[bin]->Fill(reco_z_sj, reco_z_sj); // Diagonal fill for reco/truth
      }
    }
  }
}

// -- Analyze Truth Jets --
void JetUnfoldingSubjets::AnalyzeTruthJets(JetContainer* truthJets,
                                           PHG4TruthInfoContainer* truthInfo) {
  JetDefinition jetDefAKT_R04(antikt_algorithm, 0.4);
  JetDefinition jetDefAKT_R01(antikt_algorithm, 0.1);

  for (auto truthJet : *truthJets) {
    if (truthJet->get_pt() < 30.0) continue;
    if (fabs(truthJet->get_eta()) > 0.7) continue;

    if (truthToReco.find(truthJet) != truthToReco.end()) continue; // Skip matched

    hTruthJetPtMatched->Fill(truthJet->get_pt());
    m_response1D->Miss(truthJet->get_pt());

    auto truth_particles = BuildTruthPseudoJets(truthJet, truthInfo);
    ClusterSequence truthClustSeq(truth_particles, jetDefAKT_R04);
    auto truthjets = sorted_by_pt(truthClustSeq.inclusive_jets());
    if (truthjets.empty()) continue;

    PseudoJet truth_leading = truthjets[0];
    ClusterSequence truthSubClust(truth_leading.constituents(), jetDefAKT_R01);
    auto truthsubjets = sorted_by_pt(truthSubClust.inclusive_jets());
    if (truthsubjets.size() >= 2 && truthsubjets[0].pt() >= 3 && truthsubjets[1].pt() >= 3) {
      double truth_z_sj = truthsubjets[1].pt() / (truthsubjets[0].pt() + truthsubjets[1].pt());
      int bin = FindPtBin(truthJet->get_pt());
      if (bin >= 0) {
        m_hTruthZsjMatched[bin]->Fill(truth_z_sj);
        m_responseZsj[bin]->Miss(truth_z_sj);
      }
    }
  }
}

int JetUnfoldingSubjets::Init(PHCompositeNode*) {
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  // Setup tree
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

  // 1D pT response
  float reco_ptmin = 0;
  float reco_ptmax = 50;
  float truth_ptmin = 30;
  float truth_ptmax = 60;
  int nbins_reco = 20;
  int nbins_truth = 10;

  m_response1D = new RooUnfoldResponse(nbins_reco, reco_ptmin, reco_ptmax,
                                   nbins_truth, truth_ptmin, truth_ptmax);

  hRecoJetPtMatched  = new TH1F("hRecoJetPtMatched", "Reco Jet pT (Matched);p_{T} [GeV];Jets",
                                nbins_reco, reco_ptmin, reco_ptmax);

  hTruthJetPtMatched = new TH1F("hTruthJetPtMatched", "Truth Jet pT (Matched);p_{T} [GeV];Jets",
                                nbins_truth, truth_ptmin, truth_ptmax);

  // z_sj histograms and responses per pt bin
  for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
    float ptlow = pt_bins[i];
    float pthigh = pt_bins[i + 1];
    std::string label = Form("ptbin_%d_%d", (int)ptlow, (int)pthigh);

    TH1F* hReco = new TH1F(("hRecoZsj_" + label).c_str(),
                           Form("Reco z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                           20, 0, 0.5);

    TH1F* hTruth = new TH1F(("hTruthZsj_" + label).c_str(),
                            Form("Truth z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                            20, 0, 0.5);

    RooUnfoldResponse* resp = new RooUnfoldResponse(hReco, hTruth);

    m_hRecoZsjMatched.push_back(hReco);
    m_hTruthZsjMatched.push_back(hTruth);
    m_responseZsj.push_back(resp);
    m_hRecoZsjUnfolded.push_back(nullptr);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::process_event(PHCompositeNode* topNode) {
  ++m_event;
  m_pt.clear(); m_eta.clear(); m_phi.clear();
  m_pt_truth.clear(); m_eta_truth.clear(); m_phi_truth.clear();
  recoToTruth.clear();
  truthToReco.clear();

  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  JetContainer* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  TowerBackground* bg = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  PHG4TruthInfoContainer* truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!jets || !jetsMC || !towersEM3 || !towersIH3 || !towersOH3 || !geomEM || !geomOH || !bg || !cent_node || !truthInfo) {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_centrality = cent_node->get_centile(CentralityInfo::PROP::bimp);
  m_impactparam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
  float v2 = bg->get_v2();
  float psi2 = bg->get_Psi2();

  // Matching
  MatchJets1to1(jets, jetsMC, 0.2);

  // Analyze matched jets
  AnalyzeMatchedJets(jets, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2, truthInfo);

  // Analyze truth jets
  AnalyzeTruthJets(jetsMC, truthInfo);

  // Fill tree
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::End(PHCompositeNode*) {
  PHTFileServer::get().cd(m_outputFileName);

  // Write tree
  m_T->Write();

  // Unfold pT
  RooUnfoldBayes unfold(m_response1D, hRecoJetPtMatched, 4);
  unfold.SetNToys(0);
  hRecoJetPtUnfolded = (TH1F*) unfold.Hunfold(RooUnfolding::kErrors);

  hRecoJetPtMatched->Write();
  hTruthJetPtMatched->Write();
  hRecoJetPtUnfolded->Write();

  TH2D* hResponseMatrix = (TH2D*) m_response1D->Hresponse();
  hResponseMatrix->Write("responsePt");

  // Unfold z_sj
  for (size_t i = 0; i < m_responseZsj.size(); ++i) {
    RooUnfoldBayes unfold(m_responseZsj[i], m_hRecoZsjMatched[i], 4);
    unfold.SetNToys(0);
    m_hRecoZsjUnfolded[i] = (TH1F*) unfold.Hunfold(RooUnfolding::kErrors);

    if (m_hRecoZsjUnfolded[i]) {
      m_hRecoZsjUnfolded[i]->Write(Form("hRecoZsjUnfolded_ptbin_%zu", i));
    }

    m_hRecoZsjMatched[i]->Write();
    m_hTruthZsjMatched[i]->Write();
    m_responseZsj[i]->Hresponse()->Write(Form("m_responseZsj_%zu", i));
  }

  std::cout << "Unfolded z_sj for " << m_responseZsj.size() << " pt bins written." << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
int JetUnfoldingSubjets::Reset(PHCompositeNode* /*topNode*/) {
  std::cout << "JetUnfoldingSubjets::Reset(PHCompositeNode* topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void JetUnfoldingSubjets::Print(const std::string& what) const {
  std::cout << "JetUnfoldingSubjets::Print(const std::string& what) Printing info for " << what << std::endl;
}
