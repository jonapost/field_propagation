//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <time.h>
// #include <sys/resource.h>  //  For use of 'rusage'

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class CMSApp;
class G4Step;

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void FillPerEvent(unsigned long nphyssteps, unsigned long nprims, 
                    unsigned long nsecs, unsigned long ntotal, unsigned long nall);

  void FillHist(const G4Step *);

  static G4bool isTabPhys;     // running with TABPHYS ?

  // 0 -> just run-time is reported 
  // 1 -> run-time and step statistics
  // 2 -> run-time, step statistics and CMS ECAL histograms
  static G4int fgScoreTypeFlag; 

private:
  clock_t  fRunTime; // the total runtime

  // rusage   fStartRunUsage; // -- could use to show more detailed stats

  unsigned long    fNumPysLimStepsRun;   // number of steps limited by physics
  unsigned long    fNumPrimsRun;         // number of secondaries  
  unsigned long    fNumSecsRun;          // number of primaries
  unsigned long    fNumTotalStepsRun;    // number of total steps 
  unsigned long    fNumAllStepsRun;      // number of ALL steps

  // G4double  fSumTime; // Unused

  // The CMS application 
  CMSApp   *fCMSApp; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

