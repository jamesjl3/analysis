#include "JetUnfoldingSubjets.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <centrality/CentralityInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerDefs.h>
#include <jetbackground/TowerBackground.h>

#include <fastjet/ClusterSequence.hh>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include <TTree.h>
#include <TH1F.h>
#include <TVector2.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <tuple>

namespace
{
constexpr float kTruthJetPtMin = 10.f;
constexpr float kRecoJetPtMin  = 5.f;
constexpr float kJetEtaMax     = 0.7f;
constexpr float kSubjetPtMin   = 3.f;
const std::vector<float> pt_bins = {10, 15, 20, 25, 30, 35};

inline int FindPtBin( float pt )
{
  for( size_t i = 0; i < pt_bins.size() - 1; ++i )
  {
    if( pt >= pt_bins[i] && pt < pt_bins[i+1] ) return static_cast<int>(i);
  }
  return -1;
}

inline float deltaR( const Jet* a, const Jet* b )
{
  const float deta = a->get_eta() - b->get_eta();
  const float dphi = TVector2::Phi_mpi_pi( a->get_phi() - b->get_phi() );
  return std::sqrt( deta*deta + dphi*dphi );
}
} // namespace

JetUnfoldingSubjets::JetUnfoldingSubjets( const std::string& recojetname,
                                          const std::string& truthjetname,
                                          const std::string& outputfilename )
  : SubsysReco("JetUnfoldingSubjets")
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_event(-1)
  , m_centrality(-1)
  , m_impactparam(-1)
{}

JetUnfoldingSubjets::~JetUnfoldingSubjets() = default;

int JetUnfoldingSubjets::Init( PHCompositeNode* )
{
  PHTFileServer::get().open( m_outputFileName, "RECREATE" );

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

  constexpr float reco_ptmin = kRecoJetPtMin, reco_ptmax = 60.f;
  constexpr float truth_ptmin = kTruthJetPtMin, truth_ptmax = 35.f;
  constexpr int nbins_reco = 20, nbins_truth = 10;

  m_response1D = std::make_unique<RooUnfoldResponse>( nbins_reco, reco_ptmin, reco_ptmax,
                                                      nbins_truth, truth_ptmin, truth_ptmax );
  m_response1D->Hresponse()->SetName("responsePt");

  hRecoJetPtMatched = std::make_unique<TH1F>("hRecoJetPtMatched",
                                             "Reco Jet pT (Matched);p_{T} [GeV];Jets",
                                             nbins_reco, reco_ptmin, reco_ptmax );
  hTruthJetPtMatched = std::make_unique<TH1F>("hTruthJetPtMatched",
                                              "Truth Jet pT (Matched);p_{T} [GeV];Jets",
                                              nbins_truth, truth_ptmin, truth_ptmax );
  hRecoJetPtMatched->SetDirectory(nullptr);
  hTruthJetPtMatched->SetDirectory(nullptr);

  for( size_t i = 0; i < pt_bins.size() - 1; ++i )
  {
    const float ptlow = pt_bins[i];
    const float pthigh = pt_bins[i+1];
    const std::string label = Form("ptbin_%d_%d", static_cast<int>(ptlow), static_cast<int>(pthigh));

    auto hReco = std::make_unique<TH1F>( ("hRecoZsj_" + label).c_str(),
                                         Form("Reco z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                                         20, 0, 0.5 );
    auto hTruth = std::make_unique<TH1F>( ("hTruthZsj_" + label).c_str(),
                                          Form("Truth z_{sj} [%d-%d GeV];z_{sj};Entries", (int)ptlow, (int)pthigh),
                                          20, 0, 0.5 );
    hReco->SetDirectory(nullptr);
    hTruth->SetDirectory(nullptr);

    auto resp = std::make_unique<RooUnfoldResponse>( hReco.get(), hTruth.get() );
    resp->Hresponse()->SetName( Form("m_responseZsj_%zu", i) );

    m_hRecoZsjMatched.push_back( std::move(hReco) );
    m_hTruthZsjMatched.push_back( std::move(hTruth) );
    m_responseZsj.push_back( std::move(resp) );
    m_hRecoZsjUnfolded.emplace_back( nullptr );
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetUnfoldingSubjets::MatchJets1to1( JetContainer* recoJets, JetContainer* truthJets, float dRMax )
{
  recoToTruth.clear();
  truthToReco.clear();

  std::vector<std::tuple<float, Jet*, Jet*>> pairs;
  for( auto reco : *recoJets )
  {
    if( reco->get_pt() < kRecoJetPtMin || std::fabs(reco->get_eta()) > kJetEtaMax ) continue;
    for( auto truth : *truthJets )
    {
      if( truth->get_pt() < kTruthJetPtMin || std::fabs(truth->get_eta()) > kJetEtaMax ) continue;
      const float dR = deltaR( reco, truth );
      if( dR < dRMax ) pairs.emplace_back( dR, reco, truth );
    }
  }

  std::sort( pairs.begin(), pairs.end() );
  std::set<Jet*> matchedReco, matchedTruth;
  for( const auto& tup : pairs )
  {
    Jet* reco = std::get<1>(tup);
    Jet* truth = std::get<2>(tup);
    if( matchedReco.count(reco) || matchedTruth.count(truth) ) continue;
    recoToTruth[reco] = truth;
    truthToReco[truth] = reco;
    matchedReco.insert(reco);
    matchedTruth.insert(truth);
  }
}

std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildPseudoJets( Jet* jet,
    TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
    TowerBackground* bg, float v2, float psi2, bool doUnsub ) const
{
  std::vector<fastjet::PseudoJet> particles;
  for( auto comp : jet->get_comp_vec() )
  {
    TowerInfo* tower = nullptr;
    const unsigned int ch = comp.second;
    float eta = 0, phi = 0, UE = 0;
    if( comp.first == 14 || comp.first == 29 )
    {
      tower = em->get_tower_at_channel( ch );
      if( !tower || !geomEM ) continue;
      const auto calokey = em->encode_key( ch );
      const int ieta = em->getTowerEtaBin( calokey );
      const int iphi = em->getTowerPhiBin( calokey );
      const auto key = RawTowerDefs::encode_towerid( RawTowerDefs::HCALIN, ieta, iphi );
      auto geom = geomEM->get_tower_geometry( key );
      if( !geom ) continue;
      eta = geom->get_eta();
      phi = geom->get_phi();
      UE = bg->get_UE(0).at( ieta );
    }
    else if( comp.first == 15 || comp.first == 30 )
    {
      tower = ih->get_tower_at_channel( ch );
      if( !tower || !geomEM ) continue;
      const auto calokey = ih->encode_key( ch );
      const int ieta = ih->getTowerEtaBin( calokey );
      const int iphi = ih->getTowerPhiBin( calokey );
      const auto key = RawTowerDefs::encode_towerid( RawTowerDefs::HCALIN, ieta, iphi );
      auto geom = geomEM->get_tower_geometry( key );
      if( !geom ) continue;
      eta = geom->get_eta();
      phi = geom->get_phi();
      UE = bg->get_UE(1).at( ieta );
    }
    else if( comp.first == 16 || comp.first == 31 )
    {
      tower = oh->get_tower_at_channel( ch );
      if( !tower || !geomOH ) continue;
      const auto calokey = oh->encode_key( ch );
      const int ieta = oh->getTowerEtaBin( calokey );
      const int iphi = oh->getTowerPhiBin( calokey );
      const auto key = RawTowerDefs::encode_towerid( RawTowerDefs::HCALOUT, ieta, iphi );
      auto geom = geomOH->get_tower_geometry( key );
      if( !geom ) continue;
      eta = geom->get_eta();
      phi = geom->get_phi();
      UE = bg->get_UE(2).at( ieta );
    }
    else continue;

    UE *= ( 1 + 2*v2*std::cos( 2*(phi - psi2) ) );
    float energy = tower->get_energy();
    if( doUnsub ) energy -= UE;
    const float pt = energy / std::cosh( eta );
    const float px = pt * std::cos( phi );
    const float py = pt * std::sin( phi );
    const float pz = pt * std::sinh( eta );
    particles.emplace_back( px, py, pz, energy );
  }
  return particles;
}

std::vector<fastjet::PseudoJet> JetUnfoldingSubjets::BuildTruthPseudoJets( Jet* truthJet,
    PHG4TruthInfoContainer* truthInfo ) const
{
  std::vector<fastjet::PseudoJet> particles;
  for( auto comp : truthJet->get_comp_vec() )
  {
    const unsigned int tid = comp.second;
    PHG4Particle* particle = truthInfo ? truthInfo->GetParticle( tid ) : nullptr;
    if( !particle ) continue;
    particles.emplace_back( particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e() );
  }
  return particles;
}

void JetUnfoldingSubjets::AnalyzeMatchedJets( JetContainer* /*recoJets*/,
    TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
    RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
    TowerBackground* bg, float v2, float psi2,
    PHG4TruthInfoContainer* truthInfo )
{
  fastjet::JetDefinition jetDefAKT_R04( fastjet::antikt_algorithm, 0.4 );
  fastjet::JetDefinition jetDefAKT_R01( fastjet::antikt_algorithm, 0.1 );

  for( const auto& pair : recoToTruth )
  {
    Jet* recoJet = pair.first;
    Jet* truthJet = pair.second;
    if( !recoJet || !truthJet ) continue;
    const float reco_pt = recoJet->get_pt();
    const float truth_pt = truthJet->get_pt();
    if( reco_pt < kRecoJetPtMin || truth_pt < kTruthJetPtMin ) continue;
    if( std::fabs(recoJet->get_eta()) > kJetEtaMax || std::fabs(truthJet->get_eta()) > kJetEtaMax ) continue;

    m_pt.push_back( reco_pt );
    m_eta.push_back( recoJet->get_eta() );
    m_phi.push_back( recoJet->get_phi() );
    m_pt_truth.push_back( truth_pt );
    m_eta_truth.push_back( truthJet->get_eta() );
    m_phi_truth.push_back( truthJet->get_phi() );

    if( hRecoJetPtMatched ) hRecoJetPtMatched->Fill( reco_pt );
    if( hTruthJetPtMatched ) hTruthJetPtMatched->Fill( truth_pt );
    if( m_response1D ) m_response1D->Fill( reco_pt, truth_pt );

    auto particles = BuildPseudoJets( recoJet, em, ih, oh, geomEM, geomOH, bg, v2, psi2, false );
    fastjet::ClusterSequence clustSeq( particles, jetDefAKT_R04 );
    auto jets = fastjet::sorted_by_pt( clustSeq.inclusive_jets() );
    if( jets.empty() ) continue;
    fastjet::PseudoJet leading = jets[0];
    fastjet::ClusterSequence subClust( leading.constituents(), jetDefAKT_R01 );
    auto subjets = fastjet::sorted_by_pt( subClust.inclusive_jets() );
    if( subjets.size() < 2 ) continue;
    if( subjets[0].pt() < kSubjetPtMin || subjets[1].pt() < kSubjetPtMin ) continue;

    const double reco_z_sj = subjets[1].pt() / ( subjets[0].pt() + subjets[1].pt() );
    const int bin = FindPtBin( truth_pt );
    if( bin >= 0 && bin < static_cast<int>(m_hRecoZsjMatched.size()) )
    {
      if( m_hRecoZsjMatched[bin] ) m_hRecoZsjMatched[bin]->Fill( reco_z_sj );
      if( m_responseZsj[bin] ) m_responseZsj[bin]->Fill( reco_z_sj, reco_z_sj );
    }
  }
}

void JetUnfoldingSubjets::AnalyzeTruthJets( JetContainer* truthJets,
    PHG4TruthInfoContainer* truthInfo )
{
  fastjet::JetDefinition jetDefAKT_R04( fastjet::antikt_algorithm, 0.4 );
  fastjet::JetDefinition jetDefAKT_R01( fastjet::antikt_algorithm, 0.1 );

  for( auto truthJet : *truthJets )
  {
    const float pt = truthJet->get_pt();
    if( pt < kTruthJetPtMin || std::fabs(truthJet->get_eta()) > kJetEtaMax ) continue;
    if( truthToReco.count( truthJet ) ) continue; // already matched

    if( hTruthJetPtMatched ) hTruthJetPtMatched->Fill( pt );
    if( m_response1D ) m_response1D->Miss( pt );

    auto truth_particles = BuildTruthPseudoJets( truthJet, truthInfo );
    fastjet::ClusterSequence truthClustSeq( truth_particles, jetDefAKT_R04 );
    auto truthjets = fastjet::sorted_by_pt( truthClustSeq.inclusive_jets() );
    if( truthjets.empty() ) continue;
    fastjet::PseudoJet truth_leading = truthjets[0];
    fastjet::ClusterSequence truthSubClust( truth_leading.constituents(), jetDefAKT_R01 );
    auto truthsubjets = fastjet::sorted_by_pt( truthSubClust.inclusive_jets() );
    if( truthsubjets.size() < 2 ) continue;
    if( truthsubjets[0].pt() < kSubjetPtMin || truthsubjets[1].pt() < kSubjetPtMin ) continue;

    const double truth_z_sj = truthsubjets[1].pt() / ( truthsubjets[0].pt() + truthsubjets[1].pt() );
    const int bin = FindPtBin( pt );
    if( bin >= 0 && bin < static_cast<int>(m_hTruthZsjMatched.size()) )
    {
      if( m_hTruthZsjMatched[bin] ) m_hTruthZsjMatched[bin]->Fill( truth_z_sj );
      if( m_responseZsj[bin] ) m_responseZsj[bin]->Miss( truth_z_sj );
    }
  }
}

int JetUnfoldingSubjets::process_event( PHCompositeNode* topNode )
{
  ++m_event;
  m_pt.clear();
  m_eta.clear();
  m_phi.clear();
  m_pt_truth.clear();
  m_eta_truth.clear();
  m_phi_truth.clear();
  recoToTruth.clear();
  truthToReco.clear();

  auto* jets = findNode::getClass<JetContainer>( topNode, m_recoJetName );
  auto* jetsMC = findNode::getClass<JetContainer>( topNode, m_truthJetName );
  auto* towersEM3 = findNode::getClass<TowerInfoContainer>( topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1" );
  auto* towersIH3 = findNode::getClass<TowerInfoContainer>( topNode, "TOWERINFO_CALIB_HCALIN_SUB1" );
  auto* towersOH3 = findNode::getClass<TowerInfoContainer>( topNode, "TOWERINFO_CALIB_HCALOUT_SUB1" );
  auto* geomEM = findNode::getClass<RawTowerGeomContainer>( topNode, "TOWERGEOM_HCALIN" );
  auto* geomOH = findNode::getClass<RawTowerGeomContainer>( topNode, "TOWERGEOM_HCALOUT" );
  auto* bg = findNode::getClass<TowerBackground>( topNode, "TowerInfoBackground_Sub2" );
  auto* cent_node = findNode::getClass<CentralityInfo>( topNode, "CentralityInfo" );
  auto* truthInfo = findNode::getClass<PHG4TruthInfoContainer>( topNode, "G4TruthInfo" );

  if( !jets || !jetsMC || !towersEM3 || !towersIH3 || !towersOH3 || !geomEM || !geomOH || !bg || !cent_node || !truthInfo )
  {
    std::cerr << "ERROR: Missing input node(s), aborting event." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_centrality = cent_node->get_centile( CentralityInfo::PROP::bimp );
  m_impactparam = cent_node->get_quantity( CentralityInfo::PROP::bimp );
  const float v2 = bg->get_v2();
  const float psi2 = bg->get_Psi2();

  MatchJets1to1( jets, jetsMC, 0.2 );
  AnalyzeMatchedJets( jets, towersEM3, towersIH3, towersOH3, geomEM, geomOH, bg, v2, psi2, truthInfo );
  AnalyzeTruthJets( jetsMC, truthInfo );
  if( m_T ) m_T->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::End( PHCompositeNode* )
{
  PHTFileServer::get().cd( m_outputFileName );
  if( m_T ) m_T->Write();

  if( m_response1D && hRecoJetPtMatched )
  {
    RooUnfoldBayes unfold( m_response1D.get(), hRecoJetPtMatched.get(), 4 );
    unfold.SetNToys(0);
    std::unique_ptr<TH1> hRaw( unfold.Hunfold( RooUnfolding::kErrors ) );
    if( hRaw )
    {
      hRaw->SetDirectory(nullptr);
      hRaw->SetName("hRecoJetPtUnfolded");
      hRaw->Write();
      hRecoJetPtUnfolded.reset( dynamic_cast<TH1F*>( hRaw.release() ) );
    }
  }
  if( hRecoJetPtMatched ) hRecoJetPtMatched->Write();
  if( hTruthJetPtMatched ) hTruthJetPtMatched->Write();
  if( m_response1D && m_response1D->Hresponse() ) m_response1D->Hresponse()->Write("m_response1D");

  for( size_t i = 0; i < m_responseZsj.size(); ++i )
  {
    auto& reco = m_hRecoZsjMatched[i];
    auto& truth = m_hTruthZsjMatched[i];
    auto& resp = m_responseZsj[i];
    if( !reco || !truth || !resp ) continue;
    if( reco->GetEntries() == 0 || truth->GetEntries() == 0 ) continue;
    auto* hResp = resp->Hresponse();
    if( !hResp || hResp->GetEntries() == 0 ) continue;
    RooUnfoldBayes unfold( resp.get(), reco.get(), 4 );
    unfold.SetNToys(0);
    std::unique_ptr<TH1> hRaw( unfold.Hunfold( RooUnfolding::kErrors ) );
    if( hRaw )
    {
      hRaw->SetDirectory(nullptr);
      hRaw->SetName( Form("hRecoZsjUnfolded_ptbin_%zu", i) );
      hRaw->Write();
      m_hRecoZsjUnfolded[i].reset( dynamic_cast<TH1F*>( hRaw.release() ) );
    }
    reco->Write( Form("hRecoZsjMatched_ptbin_%zu", i) );
    truth->Write( Form("hTruthZsjMatched_ptbin_%zu", i) );
    if( hResp ) hResp->Write( Form("m_responseZsj_%zu", i) );
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetUnfoldingSubjets::Reset( PHCompositeNode* )
{
  std::cout << "JetUnfoldingSubjets::Reset(PHCompositeNode*) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void JetUnfoldingSubjets::Print( const std::string& what ) const
{
  std::cout << "JetUnfoldingSubjets::Print(const std::string& what) Printing info for "
            << what << std::endl;
}

