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

inline int FindPtBin(float pt) {
    for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
        if (pt >= pt_bins[i] && pt < pt_bins[i + 1]) return static_cast<int>(i);
    }
    return -1;
}

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
      m_pt(),
      m_eta(),
      m_phi(),
      m_pt_truth(),
      m_eta_truth(),
      m_phi_truth(),
      m_response1D(),
      hRecoJetPtMatched(),
      hTruthJetPtMatched(),
      hRecoJetPtUnfolded()
{
}

JetUnfoldingSubjets::~JetUnfoldingSubjets() {
    // Do NOT delete any histogram or tree here. Let ROOT handle all cleanup.
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

    m_T->SetDirectory(nullptr);
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
    hRecoJetPtMatched->Sumw2();
    hTruthJetPtMatched = new TH1F("hTruthJetPtMatched", "Truth Jet pT (Matched);p_{T} [GeV];Jets", nbins_truth, truth_ptmin, truth_ptmax);

    m_hRecoZsjMatched.clear();
    m_hTruthZsjMatched.clear();
    m_responseZsj.clear();
    m_hRecoZsjUnfolded.clear();
    for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
      float ptlow = pt_bins[i], pthigh = pt_bins[i + 1];
      std::string label = Form("ptbin_%d_%d", (int)ptlow, (int)pthigh);
      
      // Create unique histogram names
      std::string hRecoName = Form("hRecoZsj_%s_%zu", label.c_str(), i);
      std::string hTruthName = Form("hTruthZsj_%s_%zu", label.c_str(), i);
      std::string responseName = Form("responseZsj_%s_%zu", label.c_str(), i);
      
      // Allocate and set names to ensure uniqueness
      auto hReco = new TH1F(hRecoName.c_str(),
			    Form("Reco z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
			    20, 0, 0.5);
      hReco->SetName(hRecoName.c_str());
      
      auto hTruth = new TH1F(hTruthName.c_str(),
                           Form("Truth z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
			     20, 0, 0.5);
      hTruth->SetName(hTruthName.c_str());
      
      // Construct RooUnfoldResponse with a unique name
      RooUnfoldResponse* resp = new RooUnfoldResponse(hReco, hTruth, responseName.c_str());
      resp->SetName(Form("m_responseZsj_%zu", i));
      
      // Store for later use
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
    assert(m_T != nullptr);
    m_T->Fill();
    return Fun4AllReturnCodes::EVENT_OK;
}
int JetUnfoldingSubjets::End(PHCompositeNode*) {
  PHTFileServer::get().cd(m_outputFileName);
  std::cout << "JetUnfoldingSubjets::End -- start writing output..." << std::endl;
  
  if (m_T) {
    std::cout << "Writing tree..." << std::endl;
        m_T->Write();
        std::cout << "Tree written." << std::endl;
  }
  
 
  // Unfold 1D pT
  if (m_response1D && hRecoJetPtMatched) {
    std::cout << "Unfolding 1D pT..." << std::endl;
    RooUnfoldBayes unfold(m_response1D, hRecoJetPtMatched, 4);
    unfold.SetNToys(0);
    hRecoJetPtUnfolded = static_cast<TH1F*>(unfold.Hunfold(RooUnfolding::kErrors));

    if (hRecoJetPtUnfolded) {
      hRecoJetPtUnfolded->SetName("hRecoJetPtUnfolded");
      hRecoJetPtUnfolded->SetDirectory(nullptr);
      std::cout << "Writing hRecoJetPtUnfolded..." << std::endl;
      hRecoJetPtUnfolded->Write();
      std::cout << "hRecoJetPtUnfolded written." << std::endl;
    }
  }
  
  if (hRecoJetPtMatched) {
    hRecoJetPtMatched->SetDirectory(nullptr);
    std::cout << "Writing hRecoJetPtMatched..." << std::endl;
    hRecoJetPtMatched->Write();
    std::cout << "hRecoJetPtMatched written." << std::endl;
  }
  
  if (hTruthJetPtMatched) {
    hTruthJetPtMatched->SetDirectory(nullptr);
    std::cout << "Writing hTruthJetPtMatched..." << std::endl;
    hTruthJetPtMatched->Write();
    std::cout << "hTruthJetPtMatched written." << std::endl;
  }
  
  if (m_response1D && m_response1D->Hresponse()) {
    auto* hist = m_response1D->Hresponse();
    hist->SetDirectory(nullptr);
    std::cout << "Writing m_response1D->Hresponse()..." << std::endl;
    hist->Write("m_response1D");
    std::cout << "m_response1D written." << std::endl;
  }
  
  // Unfold z_sj in bins
  for (size_t i = 0; i < m_responseZsj.size(); ++i) {
    auto& reco = m_hRecoZsjMatched[i];
    auto& truth = m_hTruthZsjMatched[i];
    auto& resp = m_responseZsj[i];
    
    if (!reco || !truth || !resp) {
      std::cout << "WARNING: Missing input for pt bin " << i << ", skipping." << std::endl;
      continue;
    }
    
    if (reco->GetEntries() == 0 || truth->GetEntries() == 0) {
      std::cout << "WARNING: No entries in pt bin " << i << ", skipping unfolding." << std::endl;
      continue;
    }
    
    if (!resp->Hresponse() || resp->Hresponse()->GetEntries() == 0) {
      std::cout << "WARNING: Empty response matrix in pt bin " << i << ", skipping." << std::endl;
      continue;
    }
    
    std::cout << "Unfolding z_sj in pt bin " << i << "..." << std::endl;
    RooUnfoldBayes unfoldZsj(resp, reco, 4);
    unfoldZsj.SetNToys(0);
    m_hRecoZsjUnfolded[i] = static_cast<TH1F*>(unfoldZsj.Hunfold(RooUnfolding::kErrors));
    m_hRecoZsjMatched[i]->Sumw2();
    std::cout << "Unfolded histogram stored in m_hRecoZsjUnfolded[" << i << "] = " << m_hRecoZsjUnfolded[i] << std::endl;  
    if (m_hRecoZsjUnfolded[i]) {
      m_hRecoZsjUnfolded[i]->SetDirectory(nullptr);
      std::cout << "Writing hRecoZsjUnfolded_ptbin_" << i << "..." << std::endl;
      m_hRecoZsjUnfolded[i]->Write(Form("hRecoZsjUnfolded_ptbin_%zu", i));
      std::cout << "hRecoZsjUnfolded_ptbin_" << i << " written." << std::endl;
      if (m_hRecoZsjUnfolded[i]) m_hRecoZsjUnfolded[i]->Write();
      else std::cerr << "WARNING: m_hRecoZsjUnfolded[" << i << "] is null during Write()\n";
      
    }
    
    reco->SetDirectory(nullptr);
    std::cout << "Writing hRecoZsjMatched_ptbin_" << i << "..." << std::endl;
    reco->Write(Form("hRecoZsjMatched_ptbin_%zu", i));
    
    truth->SetDirectory(nullptr);
    std::cout << "Writing hTruthZsjMatched_ptbin_" << i << "..." << std::endl;
    truth->Write(Form("hTruthZsjMatched_ptbin_%zu", i));
    
    auto* hist = resp->Hresponse();
    if (hist) {
      hist->SetDirectory(nullptr);
      std::cout << "Writing m_responseZsj_" << i << "..." << std::endl;
      hist->Write(Form("m_responseZsj_%zu", i));
      std::cout << "m_responseZsj_" << i << " written." << std::endl;
    }
  }
  
  std::cout << "JetUnfoldingSubjets::End -- finished writing output." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::Reset(PHCompositeNode*) {
    std::cout << "JetUnfoldingSubjets::Reset(PHCompositeNode*) being Reset" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
}

void JetUnfoldingSubjets::Print(const std::string& what) const {
    std::cout << "JetUnfoldingSubjets::Print(const std::string& what) Printing info for " << what << std::endl;
}
