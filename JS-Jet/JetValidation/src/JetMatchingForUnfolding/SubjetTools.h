#pragma once
#include <jetbase/Jet.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <jetbackground/TowerBackground.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <vector>

// Forward declarations from sPHENIX
class Jet;
class TowerInfoContainer;
class RawTowerGeomContainer;
class TowerBackground;
class PHG4TruthInfoContainer;

namespace fastjet {
  class PseudoJet;
}

/// Compute subjet observables (z_sj, theta_sj) for a *reco jet*
/// built from calorimeter towers (all calos).
///
/// \param rraw             The reco jet to analyze (usually the matched jet).
/// \param em, ih, oh       TowerInfoContainers for CEMC, HCALIN, HCALOUT.
/// \param geomEM, geomIH, geomOH  Geometry containers for each calo.
/// \param bg               TowerBackground for UE subtraction (can be nullptr).
/// \param v2, psi2         Flow background parameters (ignored if !doUnsub).
/// \param doUnsub          Whether to subtract UE.
/// \param etaCalMax        Fiducial |eta| cut for parent jet.
/// \param R_jet            Jet radius (e.g. 0.4).
/// \param z_cut, beta      Not used here (for possible SoftDrop extension).
/// \param pt_min_subjet    Minimum subjet pT threshold.
/// \param[out] z_sj, theta_sj  Computed subjet momentum fraction and opening angle.
///
/// \return true if successfully computed; false otherwise.
bool ComputeSubjetForJet(
		      Jet* rraw,
		      TowerInfoContainer* em,
		      TowerInfoContainer* ih,
		      TowerInfoContainer* oh,
		      RawTowerGeomContainer* geomEM,
		      RawTowerGeomContainer* geomIH,
		      RawTowerGeomContainer* geomOH,
		      TowerBackground* bg,
		      float v2, float psi2, bool doUnsub,
		      double etaCalMax, double R_jet,
		      double z_cut, double beta, double pt_min_subjet,
		      double& z_sj, double& theta_sj
		      );

/// Compute subjet observables (z_sj, theta_sj) for a *truth jet*
/// built from PHG4 truth particles.
///
/// \param truthJet         The truth jet to analyze (usually the matched jet).
/// \param truthInfo        Truth info container to retrieve PHG4Particles.
/// \param R_jet            Jet radius (e.g. 0.4).
/// \param etaCalMax        Fiducial |eta| cut for parent jet.
/// \param pt_min_subjet    Minimum subjet pT threshold.
/// \param[out] z_sj, theta_sj  Computed subjet momentum fraction and opening angle.
///
/// \return true if successfully computed; false otherwise.

bool ComputeSubjetForTruthJet(
			   Jet* truthJet,
			   PHG4TruthInfoContainer* truthInfo,
			   double R_jet, double etaCalMax,
			   double pt_min_subjet,
			   double& z_sj, double& theta_sj
);

// Get leading/subleading subjet pT (R=0.1) for a reco jet.
// Returns true if >= 2 subjets and both pass ptMinSubjet.
bool LeadingSubjetPts_RecoJet(
  Jet* rraw,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  TowerBackground* bg,
  float v2, float psi2, bool doUnsub,
  double etaCalMax, double R_jet, double ptMinSubjet,
  double& pt1, double& pt2);

// Same for a truth jet.
bool LeadingSubjetPts_TruthJet(
  Jet* truthJet, PHG4TruthInfoContainer* truthInfo,
  double R_jet, double etaCalMax, double ptMinSubjet,
  double& pt1, double& pt2);


// Build z_sj, theta_sj for a "truth" jet by taking ALL towers from
// CEMC/HCALIN/HCALOUT within ΔR<R_jet of the truth-jet axis.
bool ComputeSubjetForTruthJet_FromTowersCone(
  Jet* truthJet,
  TowerInfoContainer* emSim, TowerInfoContainer* ihSim, TowerInfoContainer* ohSim,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  double R_jet, double etaCalMax, double pt_min_subjet,
  double& z_sj, double& theta_sj);

// Same, but just return leading/subleading subjet pT (R=0.1)
bool LeadingSubjetPts_TruthJet_FromTowersCone(
  Jet* truthJet,
  TowerInfoContainer* emSim, TowerInfoContainer* ihSim, TowerInfoContainer* ohSim,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  double R_jet, double etaCalMax, double ptMinSubjet,
  double& pt1, double& pt2);
