// Purpose:
//    Single process to replace all (G4) physics processes using tabulations 
//      of x-section and sampled reactions.
// 
// This must be the *only* physics process registered for a particle 
//    (it uses the total cross-section, summed over the G4/real physics processes)

#ifndef TotalPhysicsProcess_HH
#define TotalPhysicsProcess_HH 1

#include "G4VRestContinuousDiscreteProcess.hh"
#include "G4ProcessType.hh"


class TabulatedDataManager;
class TTabPhysMgr; // Alternative to using TabulatedDataManager

class G4VParticleChange;

class TotalPhysicsProcess : public G4VRestContinuousDiscreteProcess {
  
public:

  TotalPhysicsProcess(G4String processName);
  virtual ~TotalPhysicsProcess();	

  // CONTINUOUS PART OF THE TOTAL PROCESS:
  // continuous step limit: will be set to a high value in order to exculde from
  // the step selection (physics step will be determined by the Discrete part) 
  // PreStepPoint momentumDirection will be changed according to the MSC model
  G4double GetContinuousStepLimit(const G4Track& track, G4double /*previousStepSize*/,
                        G4double /*currentMinimumStep*/, G4double& /*currentSafety*/);
  // continuous process action: compute energy loss for charged particles from
  // dedx
  G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step& step);

  // DISCRETE PART OF THE TOTAL PROCESS
  // mean free path for the total Discrete process part 
  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                           G4ForceCondition* condition);
  // samples, transforms and return with a final state for the total discrete 
  // process part
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& /*step*/);

  // At REST PART OF THE TOTAL PROCESS
  // mean life time for at rest processes: will be zero to force it since we have
  // only one at rest process (NuclearCaptureAtRest)
  G4double GetMeanLifeTime(const G4Track& /*track*/, G4ForceCondition* condition){
    *condition = NotForced; // will be forced other way i.e. setting to 0.0
    return 0.0;
  }
  // at rest process action: samples and return final state in the particle 
  // change for NuclearCaptureAtRest
  G4VParticleChange* AtRestDoIt(const G4Track& track, const G4Step& /*step*/);


  int SetupForMaterial(const G4Track& track);
   // For material from the track, find the corresponding Root material
   // Return (and store) its index.

  void Print(const G4Step& astep);
  
  const G4ParticleDefinition* ParticleDefinition(G4int ipdg);

  // for setting the cut in energy; particles will be stopped below this 
  static void SetEnergyLimit(G4double energyLimit){fgEnergyLimit = energyLimit;}

private:
  TotalPhysicsProcess(const TotalPhysicsProcess&); // Not implemented
  TotalPhysicsProcess& operator=(const TotalPhysicsProcess&); // Not implemented

  G4int                 fMaterialIndex;
  G4int                 fParticleId; 
  G4ParticleChange*     fParticleChange;
  // G4ParticleDefinition* fSecDefinition;
  TabulatedDataManager* theDataManager;
  // TTabPhysMgr *pTTabPhysMgr;
  static G4double       fgEnergyLimit;  //particles will be stopped below this
};
#endif

