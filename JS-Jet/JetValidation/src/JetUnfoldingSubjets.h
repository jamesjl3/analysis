#ifndef JETUNFOLDINGSUBJETS_H
#define JETUNFOLDINGSUBJETS_H

#include <fun4all/SubsysReco.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <string>
#include <vector>
#include <utility>
#include <TH1F.h>
#include <TH2D.h>
#include <set>
#include <map>

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include <fastjet/PseudoJet.hh>

class PHCompositeNode;
class TFile;
class TTree;
class TH1F;
class Jet;
class JetContainer;
class JetMap;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
//class PHG4TruthInfoContainer

class JetUnfoldingSubjets : public SubsysReco
{
public:
  JetUnfoldingSubjets(const std::string& recojetname,
		     const std::string& truthjetname,
		     const std::string& outputfilename);
  
  ~JetUnfoldingSubjets() override;
  
  void setEtaRange(double low, double high) {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
  void setPtRange(double low, double high) {
    m_ptRange.first = low;
    m_ptRange.second = high;
  }
  
  void doTruth(int flag) {
    m_doTruthJets = flag;
  }
  
  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* topNode) override;
  void Print(const std::string& what) const override;

  void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax = 0.2);

  void AnalyzeMatchedJets(JetContainer* recoJets,
			  TowerInfoContainer* towersEM3,
			  TowerInfoContainer* towersIH3,
			  TowerInfoContainer* towersOH3,
			  RawTowerGeomContainer* tower_geom,
			  RawTowerGeomContainer* tower_geomOH,
			  TowerBackground* background,
			  float background_v2,
			  float background_Psi2);

  void AnalyzeTruthJets(JetContainer* truthJets,
                      PHG4TruthInfoContainer* truthInfo);
  
private:
  std::string m_recoJetName;
  std::string m_truthJetName;
  std::string m_outputFileName;

  // Output file and TTree
  
  TTree* m_T = nullptr;

  // TTree variables
  int m_event = -1;
  float m_centrality = 0;
  float m_impactparam = 0;
  std::vector<float> m_pt;
  std::vector<float> m_eta;
  std::vector<float> m_phi;

  std::vector<float> m_pt_truth;
  std::vector<float> m_eta_truth;
  std::vector<float> m_phi_truth;

  // Histograms
  TH1F* _h_R04_z_sj_30 = nullptr;
  TH1F* _h_R04_theta_sj_30 = nullptr;
  TH1F* _h_R04_z_g_30_01 = nullptr;
  TH1F* _h_R04_theta_g_30_01 = nullptr;
  TH1F* _h_R04_truth_z_sj_30 = nullptr;
  TH1F* _h_R04_truth_theta_sj_30 = nullptr;
  TH1F* _h_R04_truth_z_g_30_01 = nullptr;
  TH1F* _h_R04_truth_theta_g_30_01 = nullptr;

  TH1F* hRecoJetPtMatched = nullptr;
  TH1F* hTruthJetPtMatched = nullptr;
  RooUnfoldResponse* response = nullptr;
  TH1F* hRecoJetPtUnfolded = nullptr;
  
  std::pair<double, double> m_etaRange { -1.1, 1.1 };
  std::pair<double, double> m_ptRange  { 5.0, 60.0 };
  bool m_doTruthJets = false;
  
  std::map<Jet*, Jet*> recoToTruth;
  std::map<Jet*, Jet*> truthToReco;
};

#endif
