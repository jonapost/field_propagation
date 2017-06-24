#ifndef GEANT_CMSApp
#define GEANT_CMSApp

#include <map>

#include "globals.hh"

class G4Step;
class TH1F;

/** @brief CMSApp class */
class CMSApp {
  static const G4int kNvolumes     = 4900;
  static const G4int kNECALModules = 36;
  static const G4int kNHCALModules = 100;

private:
  G4bool  fInitialized;                           /** Initialized flag */
  G4bool  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  G4double fEdepECAL[kNECALModules];              /** Energy deposition in ECAL */
  G4double fEdepHCAL[kNHCALModules];              /** Energy deposition in HCAL */
  G4int   fECALid[kNECALModules];                 /** ECAL volume id's */
  G4int   fHCALid[kNHCALModules];                 /** HCAL volume id's */
  G4double fCubicVolumes[kNvolumes];              /** Precomputed volume capacity */ 
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */
  TH1F   *fFluxElec;                              /** Flux histogram for electrons */
  TH1F   *fFluxGamma;                             /** Flux histogram for gammas */
  TH1F   *fFluxP;                                 /** Flux histogram for protons */
  TH1F   *fFluxPi;                                /** Flux histogram for pions */
  TH1F   *fFluxK;                                 /** Flux histogram for kaons */
  TH1F   *fEdepElec;                              /** Edep histogram for electrons */
  TH1F   *fEdepGamma;                             /** Edep histogram for gammas */
  TH1F   *fEdepP;                                 /** Edep histogram for protons */
  TH1F   *fEdepPi;                                /** Edep histogram for pions */
  TH1F   *fEdepK;                                 /** Edep histogram for kaons */ 

  /**
   * @brief Copy constructor CMSApp
   * * @todo Still not implemented
   */
  CMSApp(const CMSApp &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  CMSApp &operator=(const CMSApp &);
public:

  /** @brief Constructor CMSApp */
  CMSApp();

  /** @brief Destructor CMSApp */
  ~CMSApp() {}

  /**
   * @brief Function of initialization
   */
  G4bool Initialize();

  /**
   * @brief Step by step to do
   */
  void SteppingAction(const G4Step *);

  /**
   * @brief End of run to do
   */
  void EndOfRunAction(G4int);

public:
  static G4bool fgIsScoreActive; 

};
#endif

