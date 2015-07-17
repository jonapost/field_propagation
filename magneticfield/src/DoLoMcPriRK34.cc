// Dormand-Lockyer-McGorrigan-Prince-6-3-4 non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 7 July 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DoLoMcPri34.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 7 July 2015





//
//
//    This is the source file of DoLoMcPri34 class containing the
//    definition of the stageper() method that evaluates one stage in
//    field propagation.
//    The Butcher table of the Dormand-Lockyer-McGorrigan-Prince-6-3-4 method is as follows :
//



#include "DoLoMcPri34.hh"
#include "G4LineSection.hh"
// #include <cmath>


//Constructor
DoLoMcPri34::DoLoMcPri34(G4EquationOfMotion *EqRhs,
                                  G4int noIntegrationVariables,
                                  G4bool primary)
: FSALMagIntegratorStepper(EqRhs, noIntegrationVariables){
   
   const G4int numberOfVariables = noIntegrationVariables;
   
   //New Chunk of memory being created for use by the stageper
   
   //aki - for storing intermediate RHS
   ak2 = new G4double[numberOfVariables];
   ak3 = new G4double[numberOfVariables];
   ak4 = new G4double[numberOfVariables];
   ak5 = new G4double[numberOfVariables];
   ak6 = new G4double[numberOfVariables];
   
   yTemp = new G4double[numberOfVariables] ;
   yIn = new G4double[numberOfVariables] ;
   
   pseudoDydx_for_DistChord = new G4double[numberOfVariables];
   
   fLastInitialVector = new G4double[numberOfVariables] ;
   fLastFinalVector = new G4double[numberOfVariables] ;
   fLastDyDx = new G4double[numberOfVariables];
   
   fMidVector = new G4double[numberOfVariables];
   fMidError =  new G4double[numberOfVariables];
   if( primary )
   {
       fAuxStepper = new DoLoMcPri34(EqRhs, numberOfVariables,
                                          !primary);
   }
}


//Destructor
DoLoMcPri34::~DoLoMcPri34(){
   //clear all previously allocated memory for stageper and DistChord
   delete[] ak2;
   delete[] ak3;
   delete[] ak4;
   delete[] ak5;
   delete[] ak6;
   
   delete[] yTemp;
   delete[] yIn;
   
   delete[] fLastInitialVector;
   delete[] fLastFinalVector;
   delete[] fLastDyDx;
   delete[] fMidVector;
   delete[] fMidError;
   
   delete fAuxStepper;
   
   delete[] pseudoDydx_for_DistChord;
   
}


//stageper :

// Passing in the value of yInput[],the first time dydx[] and stage length
// Giving back yOut and yErr arrays for output and error respectively

void DoLoMcPri34::Stepper(const G4double yInput[],
                          const G4double dydx[],
                          	    G4double stage,
                            	G4double yOut[],
                              	G4double yErr[] )                              )
{
   G4int i;
   
   //The various constants defined on the basis of butcher tableu
   const G4double  //G4double - only once
  
   b21 = 0.0 ,
   
   b21 = 7.0/27.0 , 
   

   b31 = 7.0/72.0 , 
   b32 = 7.0/24.0 ,
   
   b41 = 3043.0/3528.0 , 
   b42 = -3757.0/1176.0 ,
   b43 = 1445.0/441.0,
   
   b51 =  17617.0/11662.0 ,
   b52 = -4023.0/686.0 , 
   b53 =  9372.0/1715.0 ,
   b54 = -66.0/595.0 ,
   
   b61 =  29.0/238.0 ,
   b62 =  0.0 , 
   b63 =  216.0/385.0 ,
   b64 =  54.0/85.0 ,
   b65 =  -7.0/22.0 ;
   

   
//    c1 = 35.0/384.0, c2 = .0,
//    c3 = 500.0/1113.0, c4 = 125.0/192.0,
//    c5 = -2187.0/6784.0, c6 = 11.0/84.0,
//    c7 = 0,
   
   dc1 = b61 - 363.0/2975.0 ,
	dc2 = b62 - 0.0 ,
   dc3 = b63 - 981.0/1750.0 ,
   dc4 = b64 - 2709.0/4250.0 ,
   dc5 = b65 + 3.0/10.0 ,
   dc6 = -1.0/50.0 ; //end of declaration
   
   
   const G4int numberOfVariables= this->GetNumberOfVariables();
   G4double *DyDx = new G4double[numberOfVariables];
   
   // The number of variables to be integrated over
   yOut[7] = yTemp[7]  = yIn[7];
   //  Saving yInput because yInput and yOut can be aliases for same array
   
   for(i=0;i<numberOfVariables;i++)
   {
       yIn[i]=yInput[i];
       DyDx[i] = dydx[i];
   }
   
   
   
   // RightHandSide(yIn, DyDx) ;
   // 1st stage - Not doing, getting passed
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + b21*stage*DyDx[i] ;
   }
   RightHandSide(yTemp, ak2) ;              // 2nd stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + stage*(b31*DyDx[i] + b32*ak2[i]) ;
   }
   RightHandSide(yTemp, ak3) ;              // 3rd stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + stage*(b41*DyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
   }
   RightHandSide(yTemp, ak4) ;              // 4th stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yTemp[i] = yIn[i] + stage*(b51*DyDx[i] + b52*ak2[i] + b53*ak3[i] +
                                 b54*ak4[i]) ;
   }
   RightHandSide(yTemp, ak5) ;              // 5th stage
   
   for(i=0;i<numberOfVariables;i++)
   {
       yOut[i] = yIn[i] + stage*(b61*DyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                 b64*ak4[i] + b65*ak5[i]) ;
   }
   RightHandSide(yOut, ak6) ;              // 6th and Final stage
   

   
   for(i=0;i<numberOfVariables;i++)
   {
       
       yErr[i] = stage*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                       dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] ) ;
       

       // Store Input and Final values, for possible use in calculating chord
       fLastInitialVector[i] = yIn[i] ;
       fLastFinalVector[i]   = yOut[i];
       fLastDyDx[i]          = DyDx[i];
       
       
   }
   
   fLaststageLength = stage;
   
   return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  DoLoMcPri34::DistChord() const
{
   G4double distLine, distChord;
   G4ThreeVector initialPoint, finalPoint, midPoint;
   
   // Store last initial and final points (they will be overwritten in self-stageper call!)
   initialPoint = G4ThreeVector( fLastInitialVector[0],
                                fLastInitialVector[1], fLastInitialVector[2]);
   finalPoint   = G4ThreeVector( fLastFinalVector[0],
                                fLastFinalVector[1],  fLastFinalVector[2]);
   
   // Do half a stage using stageNoErr
   
   fAuxstageper->stageper( fLastInitialVector, fLastDyDx, 0.5 * fLaststageLength,
                        fMidVector,   fMidError, pseudoDydx_for_DistChord );
   
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


void DoLoMcPri34::interpolate(  const G4double yInput[],
                                    const G4double dydx[],
                                    G4double yOut[],
                                    G4double stage,
                                	  G4double tau){
   G4double
   bf1, bf2, bf3, bf4, bf5, bf6;    

   
   const G4int numberOfVariables= this->GetNumberOfVariables();
   
   for(int i=0;i<numberOfVariables;i++)
   {
       yIn[i]=yInput[i];
   }
   
   G4double
   tau_2 = tau*tau ,
   tau_3 = tau*tau_2,
   tau_4 = tau_2*tau_2;
   
   bf1 = -(162.0*tau_3 - 504.0*tau_2 + 551.0*tau - 238.0)/238.0 ,
   bf2 =  0.0 ,
   bf3 =  27.0*tau*(27.0*tau_2 - 70.0*tau + 51.0 )/385.0 ,
   bf4 = -27*tau*(27.0*tau_2 - 50.0*tau + 21.0)/85.0 ,
   bf5 =  7.0*tau*(2232.0*tau_2 - 4166.0*tau + 1785.0 )/3278.0 
   bf6 = tau*(tau - 1)*(387.0*tau - 238.0)/149.0 ;
   
   for( int i=0; i<numberOfVariables; i++){
       yOut[i] = yIn[i] + stage*tau*(bf1*dydx[i] + bf2*ak2[i] + bf3*ak3[i] + 
								       	bf4*ak4[i] + bf5*ak5[i] + bf6*ak6[i] ) ;
   }
   

   
}




