//  Bogacki-Shampine - 8 - 5(4) implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 20 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DormandPrince745.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 20/05/15.
//

/*

This contains the stepper function of the BogackiShampine23 class

The Bogacki shampine method has the following Butcher's tableau

0  |
1/2|1/2
3/4|0	3/4
1  |2/9	1/3	4/9
-------------------
   |2/9	1/3	4/9	0
   |7/24 1/4 1/3 1/8

*/

#include "BogackiShampine23.hh"
#include "G4LineSection.hh"

using namespace std;

//Constructor
BogackiShampine23::BogackiShampine23(G4EquationOfMotion *EqRhs,
				 G4int noIntegrationVariables,
				 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    fLastStepLength(0.), fAuxStepper(0)
{
  const G4int numberOfVariables = noIntegrationVariables;

//New Chunk of memory being created for use by the stepper
  ak2 = new G4double[numberOfVariables] ;
  ak3 = new G4double[numberOfVariables] ;
  ak4 = new G4double[numberOfVariables] ;

  yTemp = new G4double[numberOfVariables] ;
  yIn = new G4double[numberOfVariables] ;

  fLastInitialVector = new G4double[numberOfVariables] ;
  fLastFinalVector = new G4double[numberOfVariables] ;
  fLastDyDx = new G4double[numberOfVariables];

  fMidVector = new G4double[numberOfVariables];
  fMidError =  new G4double[numberOfVariables];
  if( primary )
  {
    fAuxStepper = new BogackiShampine23(EqRhs, numberOfVariables, !primary);
  }
}


//Destructor
BogackiShampine23::~BogackiShampine23()
{


//clear all previously allocated memory for stepper and DistChord

  delete[] ak2;
  delete[] ak3;
  delete[] ak4;

  delete[] yTemp;
  delete[] yIn;

  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError;

  delete fAuxStepper;
}

void BogackiShampine23::Stepper(const G4double yInput[],
                         const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[])
{
    
 G4int i;
    

 const G4double  b21 = 0.5 ,
                 b31 = 0 , b32 = 3.0/4.0 ,
                 b41 = 2.0/9.0, b42 = 1.0/3.0 , b43 = 4.0/9.0 ,
				 c1 = 2.0/9.0, c2 = 1.0/3.0 , c3 = 4.0/9.0 , c4 = 0 ;


 const G4double  dc1 = c1 - 7.0/24.0 ,  dc2 = c2 - 1.0/4.0 ,
  				 dc3 = c3 - 1.0/3.0 , dc4 = c4 - 0.125 ;

 // Initialise time to t0, needed when it is not updated by the integration.
 //        [ Note: Only for time dependent fields (usually electric)
 //                  is it neccessary to integrate the time.]
 yOut[7] = yTemp[7]   = yIn[7];

 const G4int numberOfVariables= this->GetNumberOfVariables();
    
   // The number of variables to be integrated over

   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<numberOfVariables;i++)
   {
     yIn[i]=yInput[i];
   }
 // RightHandSide(yIn, dydx) ;              // 1st Step --Not doing, getting passed

	 for(i=0;i<numberOfVariables;i++)
	 {
	   yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
	 }
	 RightHandSide(yTemp, ak2) ;              // 2nd Step

	 for(i=0;i<numberOfVariables;i++)
	 {
		yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
	 }
	 RightHandSide(yTemp, ak3) ;              // 3rd Step

	 for(i=0;i<numberOfVariables;i++)
	 {
		yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
	 }
	 RightHandSide(yTemp, ak4) ;              // 4th Step

	 for(i=0;i<numberOfVariables;i++)
	 {
         yOut[i] = yIn[i] + Step*(c1*dydx[i]+ c2*ak2[i] + c3*ak3[i] + c4*ak4[i]);
         
         yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] +
		          dc4*ak4[i] ) ;

		// Store Input and Final values, for possible use in calculating chord
		fLastInitialVector[i] = yIn[i] ;
		fLastFinalVector[i]   = yOut[i];
		fLastDyDx[i]          = dydx[i];
	 }
 // NormaliseTangentVector( yOut ); // Not wanted

 fLastStepLength =Step;

 return ;
}


//The following has not been tested
//The DistChord() function - must define it here
G4double  BogackiShampine23::DistChord() const
{
  G4double distLine, distChord;
  G4ThreeVector initialPoint, finalPoint, midPoint;

  // Store last initial and final points (they will be overwritten in self-Stepper call!)
  initialPoint = G4ThreeVector( fLastInitialVector[0],
                                fLastInitialVector[1], fLastInitialVector[2]);
  finalPoint   = G4ThreeVector( fLastFinalVector[0],
                                fLastFinalVector[1],  fLastFinalVector[2]);

  // Do half a step using StepNoErr

  fAuxStepper->Stepper( fLastInitialVector, fLastDyDx, 0.5 * fLastStepLength,
           fMidVector,   fMidError );

  midPoint = G4ThreeVector( fMidVector[0], fMidVector[1], fMidVector[2]);

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord


  if (initialPoint != finalPoint)
  {
     distLine  = G4LineSection::Distline( midPoint, initialPoint, finalPoint );
     distChord = distLine;
  }
  else
  {
     distChord = (midPoint-initialPoint).mag();
  }
  return distChord;
}

























