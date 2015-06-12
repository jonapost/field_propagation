//  Tsitouras- 7 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
//  Sponsored by Google in Google Summer of Code 2015.
// 
//  First version:  11 June 2015
//
//  This code is made available subject to the Geant4 license, a copy of
//  which is available at
//  http://geant4.org/license
//  
//  History
// -----------------------------
//  Created by Somnath on 11 June 2015.
//
//////////////////////////////////////////////////////////////////////////////


#include "TsitourasRK45.hh"
#include "G4LineSection.hh"

TsitourasRK45::TsitourasRK45(G4EquationOfMotion *EqRhs, 
				 G4int noIntegrationVariables, 
				 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    fLastStepLength(0.), fAuxStepper(0)
{
  const G4int numberOfVariables = noIntegrationVariables;

  ak2 = new G4double[numberOfVariables] ;  
  ak3 = new G4double[numberOfVariables] ; 
  ak4 = new G4double[numberOfVariables] ; 
  ak5 = new G4double[numberOfVariables] ; 
  ak6 = new G4double[numberOfVariables] ; 
  ak7 = new G4double[numberOfVariables] ;
  ak8 = new G4double[numberOfVariables] ;
  yTemp = new G4double[numberOfVariables] ; 
  yIn = new G4double[numberOfVariables] ;

  fLastInitialVector = new G4double[numberOfVariables] ;
  fLastFinalVector = new G4double[numberOfVariables] ;
  fLastDyDx = new G4double[numberOfVariables];

  fMidVector = new G4double[numberOfVariables];
  fMidError =  new G4double[numberOfVariables];
  if( primary )
  { 
    fAuxStepper = new TsitourasRK45(EqRhs, numberOfVariables, !primary);
  }
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

TsitourasRK45::~TsitourasRK45()
{
  delete[] ak2;
  delete[] ak3;
  delete[] ak4;
  delete[] ak5;
  delete[] ak6;
  delete[] ak7;
    delete[] ak8;

  delete[] yTemp;
  delete[] yIn;

  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError; 

  delete fAuxStepper;
}

void
TsitourasRK45::Stepper(const G4double yInput[],
                         const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[])
{

 G4int i;

    const G4double
//    b21 = 0.2 ,
//    
//    b31 = 3.0/40.0, b32 = 9.0/40.0 ,
//    
//    b41 = 44.0/45.0, b42 = -56.0/15.0, b43 = 32.0/9.0,
//    
//    b51 = 19372.0/6561.0, b52 = -25360.0/2187.0, b53 = 64448.0/6561.0,
//    b54 = -212.0/729.0 ,
//    
//    b61 = 9017.0/3168.0 , b62 =   -355.0/33.0,
//    b63 =  46732.0/5247.0    , b64 = 49.0/176.0 ,
//    b65 = -5103.0/18656.0 ,
//    
//    b71 = 35.0/384.0, b72 = .0,
//    b73 = 500.0/1113.0, b74 = 125.0/192.0,
//    b75 = -2187.0/6784.0, b76 = 11.0/84.0,

    
    b21 = 0.161 ,
    
    b31 = 0.327 - 0.3354806554923570 , 
    b32 = 0.3354806554923570 ,
    
    b41 = 0.9 - 4.362295432869581 + 6.359448489975075 , 
    b42 = -6.359448489975075 , 
    b43 = 4.362295432869581 ,

    b51 = 0.9800255409045097 + 11.74888356406283 - 7.495539342889836 + 0.09249506636175525, 
    b52 = -11.74888356406283 , 
    b53 = 7.495539342889836 ,
    b54 = -0.09249506636175525,

    b61 = 1 + 12.92096931784711 - 8.159367898576159 + 0.07158497328140100 + 0.02826905039406838, 
    b62 = -12.92096931784711 ,
    b63 = 8.159367898576159 , 
    b64 = -0.07158497328140100,
    b65 = -0.02826905039406838,

    
    b71 = 0.09646076681806523 ,
    b72 = 0.01, 
    b73 = 0.4798896504144996, 
    b74 = 1.379008574103742, 
    b75 = -3.290069515436081, 
    b76 = 2.324710524099774,
//
//    dc1 = b71 - 5179.0/57600.0,
//    dc2 = b72 - 0.,
//    dc3 = b73 - 7571.0/16695.0,
//    dc4 = b74 - 393.0/640.0,
//    dc5 = b75 + 92097.0/339200.0,
//    dc6 = b76 - 187.0/2100.0,
//    dc7 = - 1.0/40.0 ; //end of declaration
//
    dc1 = b71 - 0.001780011052226,
    dc2 = b72 - 0.000816434459657,
    dc3 = b73 + 0.007880878010262,
    dc4 = b74 - 0.144711007173263,
    dc5 = b75 + 0.582357165452555,
    dc6 = b76 - 0.458082105929187,
    dc7 = - 1.0/66.0 ;

    
//    dc1 = -3.0/1280.0,
//    dc2 = 0.0,
//    dc3 = 6561.0/632320.0,
//    dc4 = -343.0/20800.0,
//    dc5 = 243.0/12800.0,
//    dc6 = -1.0/95.0,
//    dc7 = 0.0   ;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    //  Saving yInput because yInput and yOut can be aliases for same array
    
    for(i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    
    
    // RightHandSide(yIn, dydx) ;
    // 1st Step - Not doing, getting passed
    
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
        yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b71*dydx[i] + b72*ak2[i] + b73*ak3[i] +
                                 b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yOut, ak7);				//7th Step
    
    
    //Calculate the error in the step:
    for(i=0;i<numberOfVariables;i++)
    {
        
        yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] ) ;
        
        
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}

///////////////////////////////////////////////////////////////////////////////


G4double  TsitourasRK45::DistChord() const
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


