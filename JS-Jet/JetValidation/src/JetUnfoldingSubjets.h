#ifndef JETUNFOLDINGSUBJETS_H
#define JETUNFOLDINGSUBJETS_H

#include <fun4all/SubsysReco.h>
#include <fastjet/PseudoJet.hh>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward declarations for ROOT classes
class TTree;
class TH1F;
class RooUnfoldResponse;
class JetContainer;
class Jet;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;
class PHCompositeNode;

class JetUnfoldingSubjets : public SubsysReco {
public:
    JetUnfoldingSubjets(const std::string& recojetname,
                        const std::string& truthjetname,
                        const std::string& outputfilename);

  ~JetUnfoldingSubjets() override;
  
    int Init(PHCompositeNode*) override;
    int process_event(PHCompositeNode*) override;
    int End(PHCompositeNode*) override;
    int Reset(PHCompositeNode*) override;
    void Print(const std::string& what = "ALL") const override;
protected:
    // Matching and subjet analysis
    void MatchJets1to1(JetContainer* recoJets, JetContainer* truthJets, float dRMax);
    std::vector<fastjet::PseudoJet> BuildPseudoJets(
        Jet* jet, TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
        RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
        TowerBackground* bg, float v2, float psi2, bool doUnsub) const;
    std::vector<fastjet::PseudoJet> BuildTruthPseudoJets(
        Jet* truthJet, PHG4TruthInfoContainer* truthInfo) const;
    void AnalyzeMatchedJets(JetContainer* recoJets,
        TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
        RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomOH,
        TowerBackground* bg, float v2, float psi2,
        PHG4TruthInfoContainer* truthInfo);
    void AnalyzeTruthJets(JetContainer* truthJets, PHG4TruthInfoContainer* truthInfo);

    // I/O and configuration
    std::string m_recoJetName;
    std::string m_truthJetName;
    std::string m_outputFileName;

    // Output tree and event-level info
    std::unique_ptr<TTree> m_T;
    int m_event = -1;
    float m_centrality = -1.0f;
    float m_impactparam = -1.0f;

    // Jet tree variables
    std::vector<float> m_pt;
    std::vector<float> m_eta;
    std::vector<float> m_phi;
    std::vector<float> m_pt_truth;
    std::vector<float> m_eta_truth;
    std::vector<float> m_phi_truth;

    // Jet matching maps
    std::map<Jet*, Jet*> recoToTruth;
    std::map<Jet*, Jet*> truthToReco;

    // 1D response and spectra
    std::unique_ptr<RooUnfoldResponse> m_response1D;
    std::unique_ptr<TH1F> hRecoJetPtMatched;
    std::unique_ptr<TH1F> hTruthJetPtMatched;
    std::unique_ptr<TH1F> hRecoJetPtUnfolded;

    // z_sj per-pt-bin histograms and responses
    std::vector<std::unique_ptr<TH1F>> m_hRecoZsjMatched;
    std::vector<std::unique_ptr<TH1F>> m_hTruthZsjMatched;
    std::vector<std::unique_ptr<TH1F>> m_hRecoZsjUnfolded;
    std::vector<std::unique_ptr<RooUnfoldResponse>> m_responseZsj;
};

#endif // JETUNFOLDINGSUBJETS_H
