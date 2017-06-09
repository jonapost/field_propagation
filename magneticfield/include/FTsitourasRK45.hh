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
//  FTsitourasRK45.hh
//  Geant4
//
//  An implementation of an embedded RK method from the paper 
//     C. Tsitouras, “Runge–Kutta pairs of order 5(4) satisfying only the first
//    column simplifying assumption,”
//  Comput. Math. with Appl., vol. 62, no. 2, pp. 770–775, Jul. 2011.
//
//  This version utilises the FSAL property of the method, which
//  allows reuse, in the next step, of the derivative of the endpoint
//  of a successful step.
//
//  Implemented by Somnath Banerjee (hackabot) on 1 August 2015.
//      Work supported by the Google as part of Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
//
///////////////////////////////////////////////////////////////////////////////


#ifndef FTsitouras_RK45
#define FTsitouras_RK45

#include "G4VFSALIntegrationStepper.hh"

class FTsitourasRK45 : public G4VFSALIntegrationStepper
{
public:
    FTsitourasRK45(G4EquationOfMotion *EqRhs,
                  G4int numberOfVariables = 6,
                  G4bool primary =  true);
    ~FTsitourasRK45();
    
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[],
                 G4double nextDydx[] ) ;
    
    void interpolate( const G4double yInput[],
                     const G4double dydx[],
                     G4double yOut[],
                     G4double Step,
                     G4double tau);
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    // G4bool isFSAL() const{ return false; }
    
    private :
    
    FTsitourasRK45(const FTsitourasRK45&);
    FTsitourasRK45& operator=(const FTsitourasRK45&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double *fake_Next_Dydx;
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    FTsitourasRK45* fAuxStepper;
};



#endif /* defined(__Geant4__FTsitourasRK45__) */
