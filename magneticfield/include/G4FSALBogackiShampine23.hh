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
//  Bogacki-Shampine 4-stage 3(2) non-FSAL implementation
//  Author:  Somnath Banerjee
//             sponsored by Google in Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
// 
// First version: 20 May 2015
//
//  History
// -----------------------------
//  Created by Somnath on 20 May 2015
//
//  TODO: there are only small differences with non-FSAL version => look to unify
/////////////////////////////////////////////////////////////////////////////////


#ifndef G4FSAL_BOGACKI_SHAMPINE23_HH
#define G4FSAL_BOGACKI_SHAMPINE23_HH

#include "G4VFSALIntegrationStepper.hh"

class G4FSALBogackiShampine23 : public G4VFSALIntegrationStepper
{
 public:
 	//constructor
 	G4FSALBogackiShampine23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~G4FSALBogackiShampine23() ;

 	//Stepper
 	 void Stepper(const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[],
                      	G4double nextDydx[]) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }
     G4bool isFSAL() const{ return true; }
	G4double *getLastDydx();
    
	G4FSALBogackiShampine23(const G4FSALBogackiShampine23&);
   G4FSALBogackiShampine23& operator=(const G4FSALBogackiShampine23&);


 private:

	   G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    G4double *pseudoDydx_for_DistChord;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    G4FSALBogackiShampine23* fAuxStepper;

//	G4int No_of_vars;
//	G4double hinit, tinit, tmax, *yinit;
//	double hmax, hmin, safe_const, err0, Step_factor;
//	void (*derivs)(double, double *, double *);


};

#endif /* G4FSAL_BOGACKI_SHAMPINE23_HH */





























