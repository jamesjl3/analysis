#include "JetUnfoldingSubjets.h"

#include <fun4all/SubsysReco.h>
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
#include <memory>
#include <vector>
#include <map>

namespace {
constexpr float kTruthJetPtMin = 10.0;
constexpr float kRecoJetPtMin  = 5.0;
constexpr float kJetEtaMax     = 0.7;
constexpr float kSubjetPtMin   = 3.0;
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
      m_event(-1),
      m_centrality(-1),
      m_impactparam(-1)
{}

// Modern, robust destructor: unique_ptr for all histograms/objects, so no manual delete
JetUnfoldingSubjets::~JetUnfoldingSubjets() = default;

int JetUnfoldingSubjets::Init(PHCompositeNode*) {
    PHTFileServer::get().open(m_outputFileName, "RECREATE");
    m_T = std::make_unique<TTree>("T", "Jet Tree");
    m_T->Branch("event", &m_event, "event/I");
    m_T->Branch("cent", &m_centrality, "cent/F");
    m_T->Branch("b", &m_impactparam, "b/F");
    m_T->Branch("pt", &m_pt);
    m_T->Branch("eta", &m_eta);
    m_T->Branch("phi", &m_phi);
    m_T->Branch("pt_truth", &m_pt_truth);
    m_T->Branch("eta_truth", &m_eta_truth);
    m_T->Branch("phi_truth", &m_phi_truth);

    constexpr float reco_ptmin = kRecoJetPtMin, reco_ptmax = 60;
    constexpr float truth_ptmin = kTruthJetPtMin, truth_ptmax = 35;
    constexpr int nbins_reco = 20, nbins_truth = 10;

    m_response1D = std::make_unique<RooUnfoldResponse>(
        nbins_reco, reco_ptmin, reco_ptmax, nbins_truth, truth_ptmin, truth_ptmax);
    m_response1D->Hresponse()->SetName("responsePt");

    hRecoJetPtMatched = std::make_unique<TH1F>("hRecoJetPtMatched", "Reco Jet pT (Matched);p_{T} [GeV];Jets", nbins_reco, reco_ptmin, reco_ptmax);
    hTruthJetPtMatched = std::make_unique<TH1F>("hTruthJetPtMatched", "Truth Jet pT (Matched);p_{T} [GeV];Jets", nbins_truth, truth_ptmin, truth_ptmax);
    hRecoJetPtMatched->SetDirectory(nullptr);
    hTruthJetPtMatched->SetDirectory(nullptr);

    m_hRecoZsjMatched.clear();
    m_hTruthZsjMatched.clear();
    m_responseZsj.clear();
    m_hRecoZsjUnfolded.clear();

    for (size_t i = 0; i < pt_bins.size() - 1; ++i) {
        float ptlow = pt_bins[i];
        float pthigh = pt_bins[i + 1];
        std::string label = Form("ptbin_%d_%d", static_cast<int>(ptlow), static_cast<int>(pthigh));
        auto hReco = std::make_unique<TH1F>(("hRecoZsj_" + label).c_str(),
                           Form("Reco z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                           20, 0, 0.5);
        auto hTruth = std::make_unique<TH1F>(("hTruthZsj_" + label).c_str(),
                            Form("Truth z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                            20, 0, 0.5);
        hReco->SetDirectory(nullptr);
        hTruth->SetDirectory(nullptr);
        auto resp = std::make_unique<RooUnfoldResponse>(hReco.get(), hTruth.get());
        resp->Hresponse()->SetName(Form("m_responseZsj_%zu", i));
        m_hRecoZsjMatched.push_back(std::move(hReco));
        m_hTruthZsjMatched.push_back(std::move(hTruth));
        m_responseZsj.push_back(std::move(resp));
        m_hRecoZsjUnfolded.push_back(nullptr);
    }

	std::cout << "[INFO] m_response1D: nbins_truth=" << nbins_truth
		  << ", truth_ptmin=" << truth_ptmin
		  << ", truth_ptmax=" << truth_ptmax << std::endl;
   
    return Fun4AllReturnCodes::EVENT_OK;
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
            auto geom = geomEM->get_tower_geometry(key);
            if (!geom) continue;
            eta = geom->get_eta();
            phi = geom->get_phi();
            UE = bg->get_UE(0).at(ieta);
        } else if (comp.first == 15 || comp.first == 30) {
            tower = ih->get_tower_at_channel(ch);
            if (!tower || !geomEM) continue;
            auto calokey = ih->encode_key(ch);
            int ieta = ih->getTowerEtaBin(calokey);
            int iphi = ih->getTowerPhiBin(calokey);
            auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALIN, ieta, iphi);
            auto geom = geomEM->get_tower_geometry(key);
            if (!geom) continue;
            eta = geom->get_eta();
            phi = geom->get_phi();
            UE = bg->get_UE(1).at(ieta);
        } else if (comp.first == 16 || comp.first == 31) {
            tower = oh->get_tower_at_channel(ch);
            if (!tower || !geomOH) continue;
            auto calokey = oh->encode_key(ch);
            int ieta = oh->getTowerEtaBin(calokey);
            int iphi = oh->getTowerPhiBin(calokey);
            auto key = RawTowerDefs::encode_towerid(RawTowerDefs::HCALOUT, ieta, iphi);
            auto geom = geomOH->get_tower_geometry(key);
            if (!geom) continue;
            eta = geom->get_eta();
            phi = geom->get_phi();
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
void JetUnfoldingSubjets::AnalyzeMatchedJets(
    JetContainer* recoJets,
    TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
    TowerBackground* bg, float v2, float psi2,
    PHG4TruthInfoContainer* truthInfo)
{
    constexpr float reco_ptmin = kRecoJetPtMin;
    constexpr float reco_ptmax = 60.0;
    constexpr float truth_ptmin = kTruthJetPtMin;
    constexpr float truth_ptmax = 35.0;

    fastjet::JetDefinition jetDefAKT_R04(fastjet::antikt_algorithm, 0.4);
    fastjet::JetDefinition jetDefAKT_R01(fastjet::antikt_algorithm, 0.1);

    for (const auto& pair : recoToTruth) {
        Jet* recoJet = pair.first;
        Jet* truthJet = pair.second;
        if (!recoJet || !truthJet) continue;
        float reco_pt = recoJet->get_pt();
        float truth_pt = truthJet->get_pt();
        if (reco_pt < reco_ptmin || truth_pt < truth_ptmin) continue;
        if (fabs(recoJet->get_eta()) > kJetEtaMax || fabs(truthJet->get_eta()) > kJetEtaMax) continue;

        m_pt.push_back(reco_pt);
        m_eta.push_back(recoJet->get_eta());
        m_phi.push_back(recoJet->get_phi());
        m_pt_truth.push_back(truth_pt);
        m_eta_truth.push_back(truthJet->get_eta());
        m_phi_truth.push_back(truthJet->get_phi());

        if (hRecoJetPtMatched) hRecoJetPtMatched->Fill(reco_pt);
        if (hTruthJetPtMatched) hTruthJetPtMatched->Fill(truth_pt);

        // --- Bounds check for response1D Fill
        if (m_response1D) {
            if (reco_pt >= reco_ptmin && reco_pt < reco_ptmax &&
                truth_pt >= truth_ptmin && truth_pt < truth_ptmax)
            {
                m_response1D->Fill(reco_pt, truth_pt);
            } else {
                std::cout << "[WARNING] Skipped Fill(reco_pt=" << reco_pt << ", truth_pt=" << truth_pt
                          << "), out of bounds: reco[" << reco_ptmin << "," << reco_ptmax
                          << "), truth[" << truth_ptmin << "," << truth_ptmax << ")\n";
            }
        }

        // Subjet z_sj
        auto particles = BuildPseudoJets(recoJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false);
        fastjet::ClusterSequence clustSeq(particles, jetDefAKT_R04);
        auto jets = sorted_by_pt(clustSeq.inclusive_jets());
        if (jets.empty()) continue;
        fastjet::PseudoJet leading = jets[0];
        fastjet::ClusterSequence subClust(leading.constituents(), jetDefAKT_R01);
        auto subjets = sorted_by_pt(subClust.inclusive_jets());
        if (subjets.size() >= 2 && subjets[0].pt() >= kSubjetPtMin && subjets[1].pt() >= kSubjetPtMin) {
            double reco_z_sj = subjets[1].pt() / (subjets[0].pt() + subjets[1].pt());
            int bin = FindPtBin(truth_pt);
            if (bin >= 0 && bin < static_cast<int>(m_hRecoZsjMatched.size()) && bin < static_cast<int>(m_responseZsj.size())) {
                if (m_hRecoZsjMatched[bin]) m_hRecoZsjMatched[bin]->Fill(reco_z_sj);
                // For z_sj, only fill in-range
                constexpr float zsj_min = 0.0;
                constexpr float zsj_max = 0.5;
                if (m_responseZsj[bin] && reco_z_sj >= zsj_min && reco_z_sj < zsj_max) {
                    m_responseZsj[bin]->Fill(reco_z_sj, reco_z_sj); // Diagonal fill for reco/truth
                } else if (m_responseZsj[bin]) {
                    std::cout << "[WARNING] Skipped z_sj Fill(" << reco_z_sj << "), out of [" << zsj_min << "," << zsj_max << ")\n";
                }
            }
        }
    }
}

void JetUnfoldingSubjets::AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo) {
    constexpr float truth_ptmin = kTruthJetPtMin;
    constexpr float truth_ptmax = 35.0;
    constexpr float zsj_min = 0.0;
    constexpr float zsj_max = 0.5;

    fastjet::JetDefinition jetDefAKT_R04(fastjet::antikt_algorithm, 0.4);
    fastjet::JetDefinition jetDefAKT_R01(fastjet::antikt_algorithm, 0.1);

    for (auto truthJet : *truthJets) {
        float pt = truthJet->get_pt();
        float eta = truthJet->get_eta();

        if (pt < truth_ptmin || std::abs(eta) > kJetEtaMax) continue;
        if (truthToReco.count(truthJet)) continue; // Skip matched

        if (hTruthJetPtMatched) hTruthJetPtMatched->Fill(pt);

        // --- Safe Miss call for response1D
        if (m_response1D) {
            if (pt >= truth_ptmin && pt < truth_ptmax) {
                m_response1D->Miss(pt);
            } else {
                std::cout << "[WARNING] Skipped Miss(pt=" << pt << "), out of [" << truth_ptmin << "," << truth_ptmax << ")\n";
            }
        }

        // Subjet z_sj
        auto truth_particles = BuildTruthPseudoJets(truthJet, truthInfo);
        fastjet::ClusterSequence truthClustSeq(truth_particles, jetDefAKT_R04);
        auto truthjets = sorted_by_pt(truthClustSeq.inclusive_jets());
        if (truthjets.empty()) continue;

        fastjet::PseudoJet truth_leading = truthjets[0];
        fastjet::ClusterSequence truthSubClust(truth_leading.constituents(), jetDefAKT_R01);
        auto truthsubjets = sorted_by_pt(truthSubClust.inclusive_jets());

        if (truthsubjets.size() < 2) continue;
        if (truthsubjets[0].pt() < kSubjetPtMin || truthsubjets[1].pt() < kSubjetPtMin) continue;

        double truth_z_sj = truthsubjets[1].pt() / (truthsubjets[0].pt() + truthsubjets[1].pt());
        int bin = FindPtBin(pt);

        if (bin >= 0 && bin < static_cast<int>(m_hTruthZsjMatched.size()) && bin < static_cast<int>(m_responseZsj.size())) {
            if (m_hTruthZsjMatched[bin]) m_hTruthZsjMatched[bin]->Fill(truth_z_sj);

            // --- Safe Miss call for z_sj
            if (m_responseZsj[bin]) {
                if (truth_z_sj >= zsj_min && truth_z_sj < zsj_max) {
                    m_responseZsj[bin]->Miss(truth_z_sj);
                } else {
                    std::cout << "[WARNING] Skipped Miss(z_sj=" << truth_z_sj << "), out of [" << zsj_min << "," << zsj_max << ")\n";
                }
            }
        }
    }
}

/*
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
        if (hRecoJetPtMatched) hRecoJetPtMatched->Fill(recoJet->get_pt());
        if (hTruthJetPtMatched) hTruthJetPtMatched->Fill(truthJet->get_pt());
        if (m_response1D) m_response1D->Fill(recoJet->get_pt(), truthJet->get_pt());
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
            if (bin >= 0 && bin < (int)m_hRecoZsjMatched.size() && bin < (int)m_responseZsj.size()) {
                if (m_hRecoZsjMatched[bin]) m_hRecoZsjMatched[bin]->Fill(reco_z_sj);
                if (m_responseZsj[bin]) m_responseZsj[bin]->Fill(reco_z_sj, reco_z_sj); // Diagonal fill for reco/truth
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
        if (hTruthJetPtMatched) hTruthJetPtMatched->Fill(truthJet->get_pt());
        if (m_response1D) m_response1D->Miss(truthJet->get_pt());
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
            if (bin >= 0 && bin < (int)m_hTruthZsjMatched.size() && bin < (int)m_responseZsj.size()) {
                if (m_hTruthZsjMatched[bin]) m_hTruthZsjMatched[bin]->Fill(truth_z_sj);
                if (m_responseZsj[bin]) m_responseZsj[bin]->Miss(truth_z_sj);
            }
        }
    }
}
*/
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
    {
        std::cerr << "ERROR: Missing input node(s), aborting event." << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    m_centrality = cent_node->get_centile(CentralityInfo::PROP::bimp);
    m_impactparam = cent_node->get_quantity(CentralityInfo::PROP::bimp);
    float v2 = bg->get_v2(), psi2 = bg->get_Psi2();

    MatchJets1to1(jets, jetsMC, 0.2);
    AnalyzeMatchedJets(jets, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2, truthInfo);
    AnalyzeTruthJets(jetsMC, truthInfo);
    if (m_T) m_T->Fill();
    return Fun4AllReturnCodes::EVENT_OK;
}
int JetUnfoldingSubjets::End(PHCompositeNode*) {
    PHTFileServer::get().cd(m_outputFileName);
    if (m_T) m_T->Write();
    /*
    // Unfold 1D pT
    if (m_response1D && hRecoJetPtMatched && hTruthJetPtMatched) {
        if (hRecoJetPtMatched->GetEntries() > 0 && hTruthJetPtMatched->GetEntries() > 0) {
            RooUnfoldBayes unfold(m_response1D.get(), hRecoJetPtMatched.get(), 4);
            unfold.SetNToys(0);
            std::unique_ptr<TH1> hUnfoldedRaw(unfold.Hunfold(RooUnfolding::kErrors));
            if (hUnfoldedRaw) {
                hUnfoldedRaw->SetDirectory(nullptr);
                hUnfoldedRaw->SetName("hRecoJetPtUnfolded");
                hUnfoldedRaw->Write();
                hRecoJetPtUnfolded.reset(dynamic_cast<TH1F*>(hUnfoldedRaw.release()));
            }
        }
    }
    if (hRecoJetPtMatched) hRecoJetPtMatched->Write();
    if (hTruthJetPtMatched) hTruthJetPtMatched->Write();
    if (m_response1D && m_response1D->Hresponse()) m_response1D->Hresponse()->Write("m_response1D");

    // Zsj response unfolding
    for (size_t i = 0; i < m_responseZsj.size(); ++i) {
        auto& reco = m_hRecoZsjMatched[i];
        auto& truth = m_hTruthZsjMatched[i];
        auto& resp = m_responseZsj[i];
        if (!reco || !truth || !resp) continue;
        if (reco->GetEntries() == 0 || truth->GetEntries() == 0) continue;
        auto* hResp = resp->Hresponse();
        if (!hResp || hResp->GetEntries() == 0) continue;

        RooUnfoldBayes unfold(resp.get(), reco.get(), 4);
        unfold.SetNToys(0);
        std::unique_ptr<TH1> hUnfoldedRaw(unfold.Hunfold(RooUnfolding::kErrors));
        if (hUnfoldedRaw) {
            hUnfoldedRaw->SetDirectory(nullptr);
            hUnfoldedRaw->SetName(Form("hRecoZsjUnfolded_ptbin_%zu", i));
            hUnfoldedRaw->Write();
            // Use unique_ptr for robust memory management
            m_hRecoZsjUnfolded[i].reset(dynamic_cast<TH1F*>(hUnfoldedRaw.release()));
        }
        reco->Write(Form("hRecoZsjMatched_ptbin_%zu", i));
        truth->Write(Form("hTruthZsjMatched_ptbin_%zu", i));
        if (hResp) hResp->Write(Form("m_responseZsj_%zu", i));
	}   
    */
    
    return Fun4AllReturnCodes::EVENT_OK;
}
/*
int JetUnfoldingSubjets::End(PHCompositeNode*) {
    PHTFileServer::get().cd(m_outputFileName);
    if (m_T) m_T->Write();

    // Unfold 1D pT
    if (m_response1D && hRecoJetPtMatched && hTruthJetPtMatched) {
        if (hRecoJetPtMatched->GetEntries() > 0 && hTruthJetPtMatched->GetEntries() > 0) {
            RooUnfoldBayes unfold(m_response1D.get(), hRecoJetPtMatched.get(), 4);
            unfold.SetNToys(0);
            std::unique_ptr<TH1> hUnfoldedRaw(unfold.Hunfold(RooUnfolding::kErrors));
            TH1F* hUnfoldedFloat = dynamic_cast<TH1F*>(hUnfoldedRaw.get());
            if (hUnfoldedFloat) {
                hUnfoldedFloat->SetDirectory(nullptr);
                hUnfoldedFloat->Write();
                hRecoJetPtUnfolded.reset(dynamic_cast<TH1F*>(hUnfoldedRaw.release()));  // Take ownership
            } else if (hUnfoldedRaw) {
                hUnfoldedRaw->SetName("hRecoJetPtUnfolded");
                hUnfoldedRaw->SetDirectory(nullptr);
                hUnfoldedRaw->Write();
            }
        }
    }
    if (hRecoJetPtMatched) hRecoJetPtMatched->Write();
    if (hTruthJetPtMatched) hTruthJetPtMatched->Write();
    if (m_response1D && m_response1D->Hresponse()) m_response1D->Hresponse()->Write("m_response1D");

    for (size_t i = 0; i < m_responseZsj.size(); ++i) {
        auto& reco = m_hRecoZsjMatched[i];
        auto& truth = m_hTruthZsjMatched[i];
        auto& resp = m_responseZsj[i];
        if (!reco || !truth || !resp) continue;
        if (reco->GetEntries() == 0 || truth->GetEntries() == 0) continue;
        auto* hResp = resp->Hresponse();
        if (!hResp || hResp->GetEntries() == 0) continue;
        RooUnfoldBayes unfold(resp.get(), reco.get(), 4);
        unfold.SetNToys(0);
        std::unique_ptr<TH1> hUnfoldedRaw(unfold.Hunfold(RooUnfolding::kErrors));
        TH1F* clone = dynamic_cast<TH1F*>(hUnfoldedRaw ? hUnfoldedRaw->Clone(Form("hRecoZsjUnfolded_ptbin_%zu", i)) : nullptr);
        if (clone) {
            clone->SetDirectory(nullptr);
            clone->Write();
            if (m_hRecoZsjUnfolded[i]) delete m_hRecoZsjUnfolded[i];
            m_hRecoZsjUnfolded[i] = clone;
        }
        reco->Write(Form("hRecoZsjMatched_ptbin_%zu", i));
        truth->Write(Form("hTruthZsjMatched_ptbin_%zu", i));
        if (hResp) hResp->Write(Form("m_responseZsj_%zu", i));
    }
    return Fun4AllReturnCodes::EVENT_OK;
}
*/
// Reset and Print unchanged
int JetUnfoldingSubjets::Reset(PHCompositeNode*) {
    std::cout << "JetUnfoldingSubjets::Reset(PHCompositeNode*) being Reset" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
}
void JetUnfoldingSubjets::Print(const std::string& what) const {
    std::cout << "JetUnfoldingSubjets::Print(const std::string& what) Printing info for " << what << std::endl;
}
