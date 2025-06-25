#include "JetUnfoldingSubjets.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>
#include <centrality/CentralityInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <jetbackground/TowerBackground.h>
#include <fastjet/ClusterSequence.hh>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include <TVector2.h>
#include <TH1F.h>
#include <TTree.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <tuple>
#include <iostream>

// -- Constants --
namespace {
constexpr float kTruthJetPtMin = 10.0;   // 10 GeV truth jet cut
constexpr float kRecoJetPtMin  = 5.0;    // reco jet cut
constexpr float kJetEtaMax     = 0.7;    // |eta| acceptance
constexpr float kSubjetPtMin   = 3.0;    // subjet pt cut
const std::vector<float> pt_bins = {10, 15, 20, 25, 30, 35};

// Helper: pt bin finder
inline int FindPtBin(float pt) {
    for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
        if (pt >= pt_bins[i] && pt < pt_bins[i + 1]) return static_cast<int>(i);
    }
    return -1;
}

// Helper: deltaR between two jets
inline float deltaR(const Jet* a, const Jet* b) {
    float deta = a->get_eta() - b->get_eta();
    float dphi = TVector2::Phi_mpi_pi(a->get_phi() - b->get_phi());
    return sqrt(deta * deta + dphi * dphi);
}
} // namespace

JetUnfoldingSubjets::JetUnfoldingSubjets(const std::string& recojetname,
                                         const std::string& truthjetname,
                                         const std::string& outputfilename)
    : SubsysReco("JetUnfoldingSubjets"),
      m_recoJetName(recojetname),
      m_truthJetName(truthjetname),
      m_outputFileName(outputfilename),
      m_T(nullptr),
      m_event(-1),
      m_centrality(-1),
      m_impactparam(-1),
      m_response1D(nullptr),
      hRecoJetPtMatched(nullptr),
      hTruthJetPtMatched(nullptr),
      hRecoJetPtUnfolded(nullptr)
{}
/*
JetUnfoldingSubjets::~JetUnfoldingSubjets() {
  // RooUnfoldResponse clean-up
  if (m_response1D) {
    delete m_response1D;
    m_response1D = nullptr;
  }

  for (auto& r : m_responseZsj) {
    delete r;
    r = nullptr;
  }

  // Delete matched z_sj histograms (we own them)
  for (auto* h : m_hRecoZsjMatched) {
    delete h;
  }
  m_hRecoZsjMatched.clear();

  for (auto* h : m_hTruthZsjMatched) {
    delete h;
  }
  m_hTruthZsjMatched.clear();

  // These were cloned from Hunfold and written, do not delete
  m_hRecoZsjUnfolded.clear();

  // Jet pT matched histograms (owned by us)
  delete hRecoJetPtMatched;
  hRecoJetPtMatched = nullptr;

  delete hTruthJetPtMatched;
  hTruthJetPtMatched = nullptr;

  // Do not delete hRecoJetPtUnfolded: written clone
  hRecoJetPtUnfolded = nullptr;
}
*/
JetUnfoldingSubjets::~JetUnfoldingSubjets() {
  /*  if (m_response1D) { delete m_response1D; m_response1D = nullptr; }
  std::cout << "Reached line " << __LINE__ << std::endl;
  for (auto& r : m_responseZsj) { delete r; r = nullptr; }
  std::cout << "Reached line " << __LINE__ << std::endl;
  delete hRecoJetPtMatched; hRecoJetPtMatched = nullptr;
  std::cout << "Reached line " << __LINE__ << std::endl;
  delete hTruthJetPtMatched; hTruthJetPtMatched = nullptr;
  std::cout << "Reached line " << __LINE__ << std::endl;
  delete hRecoJetPtUnfolded; hRecoJetPtUnfolded = nullptr;
  std::cout << "Reached line " << __LINE__ << std::endl;
  // DO NOT delete m_hRecoZsjUnfolded[i]; let ROOT manage it or let them leak safely
  m_hRecoZsjUnfolded.clear();
  std::cout << "Reached line " << __LINE__ << std::endl;*/
}

void JetUnfoldingSubjets::MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax) {
    recoToTruth.clear();
    truthToReco.clear();
    std::vector<std::tuple<float, Jet*, Jet*>> pairs;
    for (auto reco : *recoJets) {
        if (reco->get_pt() < kRecoJetPtMin || fabs(reco->get_eta()) > kJetEtaMax) continue;
        for (auto truth : *truthJets) {
            if (truth->get_pt() < kTruthJetPtMin || fabs(truth->get_eta()) > kJetEtaMax) continue;
            float dR = deltaR(reco, truth);
            if (dR < dRMax) pairs.emplace_back(dR, reco, truth);
        }
    }
    std::sort(pairs.begin(), pairs.end());
    std::set<Jet*> matchedReco, matchedTruth;
    for (const auto& tup : pairs) {
        Jet* reco = std::get<1>(tup);
        Jet* truth = std::get<2>(tup);
        if (!matchedReco.count(reco) && !matchedTruth.count(truth)) {
            recoToTruth[reco] = truth;
            truthToReco[truth] = reco;
            matchedReco.insert(reco);
            matchedTruth.insert(truth);
        }
    }
}

std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildPseudoJets(
    Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
    TowerBackground* bg, float v2, float psi2, bool doUnsub) const
{
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

std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildTruthPseudoJets(
    Jet* truthJet, PHG4TruthInfoContainer* truthInfo) const
{
    std::vector<fastjet::PseudoJet> particles;
    for (auto comp : truthJet->get_comp_vec()) {
        unsigned int truth_id = comp.second;
        PHG4Particle *particle = truthInfo ? truthInfo->GetParticle(truth_id) : nullptr;
        if (!particle) continue;
        particles.emplace_back(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());
    }
    return particles;
}

void JetUnfoldingSubjets::AnalyzeMatchedJets(JetContainer* recoJets,
    TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
    TowerBackground* bg, float v2, float psi2,
    PHG4TruthInfoContainer* truthInfo)
{
    fastjet::JetDefinition jetDefAKT_R04(fastjet::antikt_algorithm, 0.4);
    fastjet::JetDefinition jetDefAKT_R01(fastjet::antikt_algorithm, 0.1);
    for (const auto& pair : recoToTruth) {
        Jet* recoJet = pair.first;
        Jet* truthJet = pair.second;
        if (!recoJet || !truthJet) continue;
        if (recoJet->get_pt() < kRecoJetPtMin || truthJet->get_pt() < kTruthJetPtMin) continue;
        if (fabs(recoJet->get_eta()) > kJetEtaMax || fabs(truthJet->get_eta()) > kJetEtaMax) continue;
        m_pt.push_back(recoJet->get_pt());
        m_eta.push_back(recoJet->get_eta());
        m_phi.push_back(recoJet->get_phi());
        m_pt_truth.push_back(truthJet->get_pt());
        m_eta_truth.push_back(truthJet->get_eta());
        m_phi_truth.push_back(truthJet->get_phi());
        hRecoJetPtMatched->Fill(recoJet->get_pt());
        hTruthJetPtMatched->Fill(truthJet->get_pt());
        m_response1D->Fill(recoJet->get_pt(), truthJet->get_pt());
        auto particles = BuildPseudoJets(recoJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false);
        fastjet::ClusterSequence clustSeq(particles, jetDefAKT_R04);
        auto jets = sorted_by_pt(clustSeq.inclusive_jets());
        if (jets.empty()) continue;
        fastjet::PseudoJet leading = jets[0];
        fastjet::ClusterSequence subClust(leading.constituents(), jetDefAKT_R01);
        auto subjets = sorted_by_pt(subClust.inclusive_jets());
        if (subjets.size() >= 2 && subjets[0].pt() >= kSubjetPtMin && subjets[1].pt() >= kSubjetPtMin) {
            double reco_z_sj = subjets[1].pt() / (subjets[0].pt() + subjets[1].pt());
            int bin = FindPtBin(truthJet->get_pt());
            if (bin >= 0) {
                m_hRecoZsjMatched[bin]->Fill(reco_z_sj);
                m_responseZsj[bin]->Fill(reco_z_sj, reco_z_sj); // Diagonal fill for reco/truth
            }
        }
    }
}

void JetUnfoldingSubjets::AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo) {
    fastjet::JetDefinition jetDefAKT_R04(fastjet::antikt_algorithm, 0.4);
    fastjet::JetDefinition jetDefAKT_R01(fastjet::antikt_algorithm, 0.1);
    for (auto truthJet : *truthJets) {
        if (truthJet->get_pt() < kTruthJetPtMin || fabs(truthJet->get_eta()) > kJetEtaMax) continue;
        if (truthToReco.count(truthJet)) continue; // Skip matched
        hTruthJetPtMatched->Fill(truthJet->get_pt());
        m_response1D->Miss(truthJet->get_pt());
        auto truth_particles = BuildTruthPseudoJets(truthJet, truthInfo);
        fastjet::ClusterSequence truthClustSeq(truth_particles, jetDefAKT_R04);
        auto truthjets = sorted_by_pt(truthClustSeq.inclusive_jets());
        if (truthjets.empty()) continue;
        fastjet::PseudoJet truth_leading = truthjets[0];
        fastjet::ClusterSequence truthSubClust(truth_leading.constituents(), jetDefAKT_R01);
        auto truthsubjets = sorted_by_pt(truthSubClust.inclusive_jets());
        if (truthsubjets.size() >= 2 && truthsubjets[0].pt() >= kSubjetPtMin && truthsubjets[1].pt() >= kSubjetPtMin) {
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
    constexpr float reco_ptmin = kRecoJetPtMin, reco_ptmax = 60;
    constexpr float truth_ptmin = kTruthJetPtMin, truth_ptmax = 35;
    constexpr int nbins_reco = 20, nbins_truth = 10;

    m_response1D = new RooUnfoldResponse(nbins_reco, reco_ptmin, reco_ptmax, nbins_truth, truth_ptmin, truth_ptmax);
    m_response1D->Hresponse()->SetName("responsePt");

    hRecoJetPtMatched  = new TH1F("hRecoJetPtMatched", "Reco Jet pT (Matched);p_{T} [GeV];Jets", nbins_reco, reco_ptmin, reco_ptmax);
    hTruthJetPtMatched = new TH1F("hTruthJetPtMatched", "Truth Jet pT (Matched);p_{T} [GeV];Jets", nbins_truth, truth_ptmin, truth_ptmax);
   
    hRecoJetPtMatched->SetDirectory(nullptr);
    hTruthJetPtMatched->SetDirectory(nullptr);

    m_hRecoZsjMatched.clear();
    m_hTruthZsjMatched.clear();
    m_responseZsj.clear();
    m_hRecoZsjUnfolded.clear();

    // Loop over z_sj pt bins
    for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
        float ptlow = pt_bins[i];
        float pthigh = pt_bins[i + 1];
        std::string label = Form("ptbin_%d_%d", static_cast<int>(ptlow), static_cast<int>(pthigh));

        TH1F* hReco = new TH1F(("hRecoZsj_" + label).c_str(),
                               Form("Reco z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                               20, 0, 0.5);
        TH1F* hTruth = new TH1F(("hTruthZsj_" + label).c_str(),
                                Form("Truth z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                                20, 0, 0.5);

        hReco->SetDirectory(nullptr);
        hTruth->SetDirectory(nullptr);

        RooUnfoldResponse* resp = new RooUnfoldResponse(hReco, hTruth);
        resp->Hresponse()->SetName(Form("m_responseZsj_%zu", i));

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
    recoToTruth.clear(); truthToReco.clear();

    auto* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
    auto* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
    auto* towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
    auto* towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
    auto* towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
    auto* geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    auto* geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    auto* bg = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
    auto* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
    auto* truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (!jets || !jetsMC || !towersEM3 || !towersIH3 || !towersOH3 || !geomEM || !geomOH || !bg || !cent_node || !truthInfo)
        return Fun4AllReturnCodes::ABORTRUN;

    m_centrality = cent_node->get_centile(CentralityInfo::PROP::bimp);
    m_impactparam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
    float v2 = bg->get_v2(), psi2 = bg->get_Psi2();

    MatchJets1to1(jets, jetsMC, 0.2);
    AnalyzeMatchedJets(jets, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2, truthInfo);
    AnalyzeTruthJets(jetsMC, truthInfo);
    m_T->Fill();
    return Fun4AllReturnCodes::EVENT_OK;
}
/*
int JetUnfoldingSubjets::End(PHCompositeNode*) {
  PHTFileServer::get().cd(m_outputFileName);
  std::cout << " about to write" << std::endl;
  if (m_T) m_T->Write();
  std::cout << " should have written" << std::endl;

  // ---- 1D pT Unfolding ----
  if (m_response1D && hRecoJetPtMatched && hTruthJetPtMatched &&
      hRecoJetPtMatched->GetEntries() > 0 && hTruthJetPtMatched->GetEntries() > 0) {
    
    RooUnfoldBayes unfold(m_response1D, hRecoJetPtMatched, 4);
    unfold.SetNToys(0);

    hRecoJetPtUnfolded = dynamic_cast<TH1F*>(unfold.Hunfold(RooUnfolding::kErrors));
    if (!hRecoJetPtUnfolded) {
      std::cout << "WARNING: Hunfold() returned null for 1D pT" << std::endl;
    } else {
      hRecoJetPtUnfolded->SetDirectory(nullptr);
      hRecoJetPtUnfolded->Write("hRecoJetPtUnfolded");
    }
  }

  if (hRecoJetPtMatched) hRecoJetPtMatched->Write();
  if (hTruthJetPtMatched) hTruthJetPtMatched->Write();
  if (m_response1D && m_response1D->Hresponse()) {
    m_response1D->Hresponse()->Write("m_response1D");
  }

  // ---- z_sj Unfolding in pt bins ----
  size_t nbins = std::min({m_hRecoZsjMatched.size(), m_hTruthZsjMatched.size(),
                           m_responseZsj.size(), m_hRecoZsjUnfolded.size()});
  
  for (size_t i = 0; i < nbins; ++i) {
    auto* reco = m_hRecoZsjMatched[i];
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    auto* truth = m_hTruthZsjMatched[i];
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    auto* resp = m_responseZsj[i];
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    if (!reco || !truth || !resp) {
      std::cout << "WARNING: Missing input for pt bin " << i << ", skipping." << std::endl;
      continue;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    if (reco->GetEntries() == 0 || truth->GetEntries() == 0) {
      std::cout << "WARNING: No entries in pt bin " << i << ", skipping unfolding." << std::endl;
      continue;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    if (!resp->Hresponse() || resp->Hresponse()->GetEntries() == 0) {
      std::cout << "WARNING: Empty response matrix in pt bin " << i << ", skipping." << std::endl;
      continue;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    RooUnfoldBayes unfold(resp, reco, 4);
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    unfold.SetNToys(0);
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    
    TH1F* unfolded = dynamic_cast<TH1F*>(unfold.Hunfold(RooUnfolding::kErrors));
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    if (!unfolded) {
      std::cout << "WARNING: Hunfold() returned null for pt bin " << i << std::endl;
      continue;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    unfolded->SetDirectory(nullptr);
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    unfolded->Write(Form("hRecoZsjUnfolded_ptbin_%zu", i));
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;

    // Replace old hist safely
    if (m_hRecoZsjUnfolded[i]) {
      delete m_hRecoZsjUnfolded[i];
      std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    m_hRecoZsjUnfolded[i] = unfolded;
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    
    // Always write matched histograms and response matrix
    reco->Write(Form("hRecoZsjMatched_ptbin_%zu", i));
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    truth->Write(Form("hTruthZsjMatched_ptbin_%zu", i));
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    if (resp->Hresponse()) {
      resp->Hresponse()->Write(Form("m_responseZsj_%zu", i));
      std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    }
    std::cout << "DEBUG: About to unfold z_sj bin " << i << " at " << __FILE__ << ":" << __LINE__ << std::endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}
*/

int JetUnfoldingSubjets::End(PHCompositeNode*) {
  PHTFileServer::get().cd(m_outputFileName);
  std::cout << " about to write" << std::endl;
  if (m_T) m_T->Write();
  std::cout << " should have written" << std::endl;

  // Unfold 1D pT
  if (m_response1D && hRecoJetPtMatched && hTruthJetPtMatched) {
    if (hRecoJetPtMatched->GetEntries() > 0 && hTruthJetPtMatched->GetEntries() > 0) {
      RooUnfoldBayes unfold(m_response1D, hRecoJetPtMatched, 4);
      unfold.SetNToys(0);
      TH1* hUnfoldedRaw = unfold.Hunfold(RooUnfolding::kErrors);
      
      if (hUnfoldedRaw) {
	std::cout << "Unfolded histogram type: " << hUnfoldedRaw->ClassName() << std::endl;
	hRecoJetPtUnfolded = (TH1F*) hUnfoldedRaw->Clone("hRecoJetPtUnfolded");  // Try cast
	
	if (!hRecoJetPtUnfolded) {
	  std::cerr << "ERROR: Failed to clone hRecoJetPtUnfolded as TH1F*. Writing original as generic TH1." << std::endl;
	  hUnfoldedRaw->SetName("hRecoJetPtUnfolded");
	  hUnfoldedRaw->SetDirectory(nullptr);
	  hUnfoldedRaw->Write();
	} else {
	  hRecoJetPtUnfolded->SetDirectory(nullptr);
	  hRecoJetPtUnfolded->Write();
	  delete hUnfoldedRaw;
	}
      }
    } else {
      std::cerr << "WARNING: 1D unfolding skipped due to empty reco/truth." << std::endl;
    }
  }
  
  if (hRecoJetPtMatched) {
    hRecoJetPtMatched->SetDirectory(nullptr);
    hRecoJetPtMatched->Write();
  }
  if (hTruthJetPtMatched) {
    hTruthJetPtMatched->SetDirectory(nullptr);
    hTruthJetPtMatched->Write();
  }
  
  if (m_response1D && m_response1D->Hresponse()) {
    auto* hResp = m_response1D->Hresponse();
    hResp->SetDirectory(nullptr);
    hResp->Write("m_response1D");
  }

  for (size_t i = 0; i < m_responseZsj.size(); ++i) {
    auto reco = m_hRecoZsjMatched[i];
    auto truth = m_hTruthZsjMatched[i];
    auto resp = m_responseZsj[i];

    if (!reco || !truth || !resp) continue;
    if (reco->GetEntries() == 0 || truth->GetEntries() == 0) continue;

    auto* hResp = resp->Hresponse();
    if (!hResp || hResp->GetEntries() == 0) continue;

    RooUnfoldBayes unfold(resp, reco, 4);
    unfold.SetNToys(0);

    auto* hUnfolded = unfold.Hunfold(RooUnfolding::kErrors);
    if (hUnfolded) {
      auto* clone = (TH1F*) hUnfolded->Clone(Form("hRecoZsjUnfolded_ptbin_%zu", i));
      clone->SetDirectory(nullptr);
      clone->Write();
      m_hRecoZsjUnfolded[i] = clone;
      std::cout << "Wrote z_sj bin " << i << " successfully" << std::endl;

    }

    reco->SetDirectory(nullptr);
    truth->SetDirectory(nullptr);
    hResp->SetDirectory(nullptr);

    reco->Write(Form("hRecoZsjMatched_ptbin_%zu", i));
    truth->Write(Form("hTruthZsjMatched_ptbin_%zu", i));
    hResp->Write(Form("m_responseZsj_%zu", i));
  }
  std::cout << "JetUnfoldingSubjets::End() finished without crash." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::Reset(PHCompositeNode*) {
    std::cout << "JetUnfoldingSubjets::Reset(PHCompositeNode*) being Reset" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
}

void JetUnfoldingSubjets::Print(const std::string& what) const {
    std::cout << "JetUnfoldingSubjets::Print(const std::string& what) Printing info for " << what << std::endl;
}
