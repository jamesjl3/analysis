#pragma once
#include <jetbase/Jet.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <jetbackground/TowerBackground.h>
#include <g4main/PHG4TruthInfoContainer.h>

// Long, fully-configurable form (matches the symbol your other code expects)
bool ComputeZsjForJet(
  Jet* jet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  TowerBackground* bg,
  float v2, float psi2, bool doUnsub,
  double etaCalMax, double R_jet, double z_cut, double beta, double pt_min_subjet,
  double& theta_sj, double& z_sj);

// Short convenience form used by your current JetMatchingSubjets.cc calls
bool ComputeZsjForJet(
  Jet* jet,
  TowerInfoContainer* em, TowerInfoContainer* ih, TowerInfoContainer* oh,
  RawTowerGeomContainer* geomEM, RawTowerGeomContainer* geomIH, RawTowerGeomContainer* geomOH,
  TowerBackground* bg,
  double R_jet, double etaCalMax,
  double& z_sj, double& theta_sj);

// Truth helper (what you’re already calling)
bool ComputeZsjForTruthJet(
  Jet* truthJet,
  PHG4TruthInfoContainer* truthInfo,
  double R_jet, double etaCalMax,
  double& z_sj, double& theta_sj);
