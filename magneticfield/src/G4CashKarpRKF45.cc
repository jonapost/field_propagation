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
// $Id: G4CashKarpRKF45.cc 97572 2016-06-03 21:52:00Z japost $
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.
//
//  [ref. Numerical Recipes in C, 2nd Edition]
//
// -------------------------------------------------------------------

#include <cassert>

#include "G4CashKarpRKF45.hh"
#include "G4LineSection.hh"

#define ncomp G4FieldTrack::ncompSVEC

/////////////////////////////////////////////////////////////////////
//
// Constructor

G4CashKarpRKF45::G4CashKarpRKF45(G4EquationOfMotion *EqRhs,
                                 G4int noIntegrationVariables,
                                 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, std::min(noIntegrationVariables,G4int(ncomp))),
    fLastStepLength(0.), fAuxStepper(nullptr)
{
  // Ensure that the temporary arrays are large enough 
  assert( noIntegrationVariables <= G4FieldTrack::ncompSVEC );

  // Must ensure space extra 'state' variables exists - i.e. yIn[7]
  const G4int numStateMax  = std::max(GetNumberOfStateVariables(), 8);  
  const G4int numStateVars = std::max(noIntegrationVariables,
                                      numStateMax );
  assert( numStateVars <= G4FieldTrack::ncompSVEC );
  
#if 0  
<<<<<<< HEAD
  const G4int numberOfVariables =
      std::max( noIntegrationVariables,
               ( ( (noIntegrationVariables-1)/4 + 1 ) * 4 ) );
  // For better alignment with cache-line

  ak2 = new G4double[numberOfVariables] ;
  ak3 = new G4double[numberOfVariables] ;
  ak4 = new G4double[numberOfVariables] ;
  ak5 = new G4double[numberOfVariables] ;
  ak6 = new G4double[numberOfVariables] ;



  // yTemp = new G4double[numStateVars] ;
  // yIn = new G4double[numStateVars] ;

  // fLastInitialVector = new G4double[numStateVars] ;
  // fLastFinalVector = new G4double[numStateVars] ;
  // fLastDyDx = new G4double[numberOfVariables];

  fMidVector = new G4double[numStateVars];
  fMidError =  new G4double[numStateVars];
=======
#endif
  if( primary )
  {
    fAuxStepper = new G4CashKarpRKF45(EqRhs,
                                      std::min(noIntegrationVariables, G4int(ncomp)),
                                      !primary);
  }
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4CashKarpRKF45::~G4CashKarpRKF45()
{
  delete fAuxStepper;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for n = 6 variables yIn[0,...,n-1]
// known  at x, use the fifth-order Cash-Karp Runge-
// Kutta-Fehlberg-4-5 method to advance the solution over an interval
// Step and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 4th-order method. The user supplies routine
// RightHandSide(y,dydx), which returns derivatives dydx for y .

void
G4CashKarpRKF45::Stepper(const G4double yInput[],
                         const G4double dydx[],
                               G4double Step,
                               G4double yOutput[],
                               G4double yErr[])
{
 // const G4int nvar = 6 ;
 // const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
 G4int i;

 static const G4double
         b21 = 1./5.,
         b31 = 3./40., b32 = 9./40.,
         b41 = 3./10., b42 = -9./10. , b43 = 6./5.,
         b51 = -11./54., b52 = 5./2., b53 = -70./27., b54 = 35./27.,

         b61 = 1631./55296., b62 = 175./512., b63 =  575./13824.,
         b64 = 44275./110592., b65 = 253./4096.,

         c1 = 37./378., c3 = 250./621., c4 = 125./594., c6 = 512./1771.,

         dc1 = c1 - 2825./27648.,  dc3 = c3 - 18575./48384.,
         dc4 = c4 - 13525./55296., dc5 = -277./14336., dc6 = c6 - 1./4.;

 // Initialise time to t0, needed when it is not updated by the integration.
 //        [ Note: Only for time dependent fields (usually electric)
 //                  is it neccessary to integrate the time.]
 yTemp[7] = yIn[7];

 // The number of variables to be integrated over
 const G4int nvar = GetNumberOfVariables();

 //Save yInput and dydx for chord calculation
 memcpy(yIn,yInput,sizeof(G4double)*ncomp);
 memcpy(dydxIn,dydx,sizeof(G4double)*ncomp);


 // RightHandSide(yIn, dydx);              // 1st Step

 for(i=0; i < nvar; ++i)
 {
   yTemp[i] = yIn[i] + b21*Step*dydx[i];
 }
 RightHandSide(yTemp, ak2);              // 2nd Step

 for(i=0; i < nvar; ++i)
 {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]);
 }
 RightHandSide(yTemp, ak3);              // 3rd Step

 for(i=0; i < nvar; ++i)
 {
    yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]);
 }
 RightHandSide(yTemp, ak4);              // 4th Step

 for(i=0; i < nvar; ++i)
 {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                      b54*ak4[i]);
 }
 RightHandSide(yTemp, ak5);              // 5th Step

 for(i=0; i < nvar; ++i)
 {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                      b64*ak4[i] + b65*ak5[i]);
 }
 RightHandSide(yTemp, ak6);              // 6th Step

 for(i=0; i < nvar; ++i)
 {
    // Accumulate increments with proper weights
    yOutput[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);

    // Estimate error as difference between 4th and 5th order methods
    yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
              dc5*ak5[i] + dc6*ak6[i]);
 }
 // copy non-invegrated vars to yOutput
 memcpy(yOutput+nvar, yIn+nvar, sizeof(G4double)*(ncomp-nvar));

 // Store Final values, for chord calculation
 memcpy(yOut,yOutput,sizeof(G4double)*ncomp);


 // NormaliseTangentVector( yOut ); // Not wanted

 fLastStepLength = Step;

}

///////////////////////////////////////////////////////////////////////////////

void
G4CashKarpRKF45::StepWithEst( const G4double*,
                              const G4double*,
                                    G4double,
                                    G4double*,
                                    G4double&,
                                    G4double&,
                              const G4double*,
                                    G4double*  )
{
  G4Exception("G4CashKarpRKF45::StepWithEst()", "GeomField0001",
              FatalException, "Method no longer used.");
}

/////////////////////////////////////////////////////////////////

G4double  G4CashKarpRKF45::DistChord() const
{
  // Do half a step using fAuxStepper
  fAuxStepper->Stepper(yIn, dydxIn, 0.5 * fLastStepLength, yMid, yError);

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord
  G4ThreeVector initialPoint = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector finalPoint   = G4ThreeVector(yOut[0], yOut[1], yOut[2]);
  G4ThreeVector midPoint = G4ThreeVector(yMid[0], yMid[1], yMid[2]);


  G4double distChord;
  if (initialPoint != finalPoint)
  {
     distChord  = G4LineSection::Distline(midPoint, initialPoint, finalPoint);
  }
  else
  {
     distChord = (midPoint-initialPoint).mag();
  }
  return distChord;
}


//inline
//void G4MagIntegratorStepper::RightHandSide( const  double y[], double dydx[] )
//{
//    fEquation_Rhs-> RightHandSide(y, dydx);
//    increasefNORHSCalls();
//}
