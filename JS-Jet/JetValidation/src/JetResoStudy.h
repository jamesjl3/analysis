// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETRESOSTUDY_H
#define JETRESOSTUDY_H

#include <fun4all/SubsysReco.h>
#include <jetbase/Jetv1.h>
#include <jetbase/Jetv2.h>
#include <fastjet/PseudoJet.hh>
#include <calobase/RawClusterv1.h>

#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <deque>
#include <CLHEP/Vector/ThreeVector.h>

#include <g4jets/TruthJetInput.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Hit.h>

class PHCompositeNode;
class TTree;
class GlobalVertexMap;
class Fun4AllHistoManager;
class RawClusterv1;
class PHG4TruthInfoContainer;

struct EventCluster {
  std::vector<RawClusterv1> clusters;
  CLHEP::Hep3Vector vertex;
};

class JetResoStudy : public SubsysReco
{
 public:

  JetResoStudy(const std::string &recojetname = "AntiKt_Tower_r04",
		const std::string &truthjetname = "AntiKt_Truth_r04",
		const std::string &outputfilename = "myjetanalysis.root");

  ~JetResoStudy() override;

  void
    setEtaRange(double low, double high)
  {
    m_etaRange.first = low;
    m_etaRange.second = high;
  }
 void
   setPtRange(double low, double high)
 {
   m_ptRange.first = low;
   m_ptRange.second = high;
 }
 void
   doTruth(int flag)
 {
   m_doTruthJets = flag;
 }

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
 float calculateDeltaR(float eta1, float phi1, float eta2, float phi2);
 
 int Init(PHCompositeNode *topNode) override;
 
 /** Called for first event when run number is known.
     Typically this is where you may want to fetch data from
     database, because you know the run number. A place
     to book histograms which have to know the run number.
   */
 int InitRun(PHCompositeNode *topNode) override;
 
  /** Called for each event.
      This is where you do the real work.
  */
 int process_event(PHCompositeNode *topNode) override;
 
 /// Clean up internals after each event.
 int ResetEvent(PHCompositeNode *topNode) override;
 
 /// Called at the end of each run.
 int EndRun(const int runnumber) override;
 
 /// Called at the end of all processing.
 int End(PHCompositeNode *topNode) override;
 
 /// Reset
 int Reset(PHCompositeNode * /*topNode*/) override;
 
 void Print(const std::string &what = "ALL") const override;
 
 void set_doSim(bool doSim) { m_doSim = doSim; };
 
 private:
 Fun4AllHistoManager *hm;
 TFile *outfile{nullptr};
 // Histograms                                                                                                                                                                                                                                                                  
 TH1F *hDeltaR;
 TH1F *hRecoJetPt;
 TH1F *hTruthJetPt;
 TH1F *hFakeRecoJetPt;
 TH1F *hFakeTruthJetPt;
 TH1F *hEfficiency;
 TH1F *hPurity;
 TH1F *hFakeRate;
 TH1F *hMomentumResolution;
 TH1F *hScaleFactor;

 std::string m_recoJetName;
 std::string m_truthJetName;
 std::string m_outputFileName;
 std::string m_histoFileName;
 std::pair<double, double> m_etaRange;
 std::pair<double, double> m_ptRange;
 
 bool m_doSim{true};
 
 int m_doTruthJets;
 float m_maxZvtx;

 //! Output Tree variables                                                                                                                                                         
 TTree *m_T;

 //! eventwise quantities                                                                       
 float m_jetPt;
 float m_jetEta;
 float m_jetPhi;
 int m_event;
 int m_nJet;
 int m_nTruthJet;
 int m_centrality;
 float m_impactparam;
 float m_zvtx;

 float m_deltaR;

 //spectra
 TH1F *spectrum;
 TH1F *energy_spectrum;

 TH1F *spectrum_truth;
 TH1F *energy_spectrum_truth;
 
 TH2D *eta_pt;
 TH2D *phi_pt;
 TH2D *M_pt;
 
 TH2D *eta_pt_truth;
 TH2D *phi_pt_truth;
 TH2D *M_pt_truth;

 TH1F *jet_AN;
 
 //////////////////// Substructure hists                                                                                                                                                                   // data anti-kT hists
 TH1F *_h_R04_z_sj_30_35_data;
 TH1F *_h_R04_theta_sj_30_35_data;
 TH1F *_h_R04_z_sj_35_40_data;
 TH1F *_h_R04_theta_sj_35_40_data;
 TH1F *_h_R04_z_sj_40_45_data;
 TH1F *_h_R04_theta_sj_40_45_data;      
 TH1F *_h_R04_z_sj_45_data;
 TH1F * _h_R04_theta_sj_45_data;
// truth anti-kT hists
 TH1F *_h_R04_z_sj_30_35_truth;
 TH1F *_h_R04_theta_sj_30_35_truth;
 TH1F *_h_R04_z_sj_35_40_truth;
 TH1F *_h_R04_theta_sj_35_40_truth;
 TH1F *_h_R04_z_sj_40_45_truth;
 TH1F *_h_R04_theta_sj_40_45_truth;      
 TH1F *_h_R04_z_sj_45_truth;
 TH1F * _h_R04_theta_sj_45_truth;
 // softdrop hists zcut = 0.1 ,data          
 TH1F *_h_R04_z_g_30_35_data;
 TH1F *_h_R04_theta_g_30_35_data;
 TH1F *_h_R04_z_g_35_40_data;
 TH1F *_h_R04_theta_g_35_40_data;
 TH1F *_h_R04_z_g_40_45_data;             
 TH1F *_h_R04_theta_g_40_45_data;
 TH1F *_h_R04_z_g_45_data;
 TH1F *_h_R04_theta_g_45_data;                                       
 // softdrop hists zcut = 0.1 ,truth                                                                                                                                                            
 TH1F *_h_R04_z_g_30_35_truth;
 TH1F *_h_R04_theta_g_30_35_truth;
 TH1F *_h_R04_z_g_35_40_truth;
 TH1F *_h_R04_theta_g_35_40_truth;
 TH1F *_h_R04_z_g_40_45_truth;
 TH1F *_h_R04_theta_g_40_45_truth;
 TH1F *_h_R04_z_g_45_truth;
 TH1F *_h_R04_theta_g_45_truth;
 // correlation hists
 TH2D* correlation_theta_30_35_data;
 TH2D* correlation_theta_35_40_data;
 TH2D* correlation_theta_40_45_data;
 TH2D* correlation_theta_45_data;
 TH2D* correlation_theta_30_35_truth;
 TH2D* correlation_theta_35_40_truth;
 TH2D* correlation_theta_40_45_truth;
 TH2D* correlation_theta_45_truth;
 TH2D* correlation_z_30_35_data;
 TH2D* correlation_z_35_40_data;
 TH2D* correlation_z_40_45_data;
 TH2D* correlation_z_45_data;
 TH2D* correlation_z_30_35_truth;
 TH2D* correlation_z_35_40_truth;
 TH2D* correlation_z_40_45_truth;
 TH2D* correlation_z_45_truth;
 
 //!trigger info
 std::vector<bool> m_triggerVector;
 
 //! reconstructed jets
 std::vector<int> m_id;
 std::vector<int> m_nComponent;
 std::vector<float> m_eta;
 std::vector<float> m_phi;
 std::vector<float> m_e;
 std::vector<float> m_pt;

 //! truth jets                                      
 std::vector<int> m_id_truth;
 std::vector<int> m_nComponent_truth;
 std::vector<float> m_eta_truth;
 std::vector<float> m_phi_truth;
 std::vector<float> m_e_truth;
 std::vector<float> m_pt_truth;
};

#endif // JETRESOSTUDY_H

