//  DormandPrince745.cc
//  Geant4
//
//  Created by Somnath on 25/05/15.
//
//    This is the source file of DormandPrince745 class containing the
//    definition of the stepper() method that evaluates one step in
//    field propagation.
//    The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
//
//    0   |
//    1/5 | 1/5
//    3/10| 3/40        9/40
//    4/5 | 44/45      −56/15      32/9
//    8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
//    1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
//    1   | 35/384      0         500/1113    125/192 −2187/6784    11/84
//    ------------------------------------------------------------------------
//          35/384       0        500/1113    125/192  −2187/6784    11/84   0
//          5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40


#include "DormandPrince745.hh"
#include "G4LineSection.hh"
#include <cmath>


//Constructor
DormandPrince745::DormandPrince745(G4EquationOfMotion *EqRhs,
                                   G4int noIntegrationVariables,
                                   G4bool primary)
: FSALMagIntegratorStepper(EqRhs, noIntegrationVariables){
    
    const G4int numberOfVariables = noIntegrationVariables;
    
    //New Chunk of memory being created for use by the stepper
    
    //aki - for storing intermediate RHS
    ak2 = new G4double[numberOfVariables];
    ak3 = new G4double[numberOfVariables];
    ak4 = new G4double[numberOfVariables];
    ak5 = new G4double[numberOfVariables];
    ak6 = new G4double[numberOfVariables];
    ak7 = new G4double[numberOfVariables];
    ak8 = 0;
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
        fAuxStepper = new DormandPrince745(EqRhs, numberOfVariables,
                                           !primary);
    }
}


//Destructor
DormandPrince745::~DormandPrince745(){
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    
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


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void DormandPrince745::Stepper(const G4double yInput[],
                               const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[],
                               G4double nextDydx[]
                               )
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    b21 = 0.2 ,
    
    b31 = 3.0/40.0, b32 = 9.0/40.0 ,
    
    b41 = 44.0/45.0, b42 = -56.0/15.0, b43 = 32.0/9.0,
    
    b51 = 19372.0/6561.0, b52 = -25360.0/2187.0, b53 = 64448.0/6561.0,
    b54 = -212.0/729.0 ,
    
    b61 = 9017.0/3168.0 , b62 =   -355.0/33.0,
    b63 =  46732.0/5247.0    , b64 = 49.0/176.0 ,
    b65 = -5103.0/18656.0 ,
    
    b71 = 35.0/384.0, b72 = .0,
    b73 = 500.0/1113.0, b74 = 125.0/192.0,
    b75 = -2187.0/6784.0, b76 = 11.0/84.0,
    
//    c1 = 35.0/384.0, c2 = .0,
//    c3 = 500.0/1113.0, c4 = 125.0/192.0,
//    c5 = -2187.0/6784.0, c6 = 11.0/84.0,
//    c7 = 0,
    
    dc1 = b71 - 5179.0/57600.0,
    dc2 = b72 - .0,
    dc3 = b73 - 7571.0/16695.0,
    dc4 = b74 - 393.0/640.0,
    dc5 = b75 + 92097.0/339200.0,
    dc6 = b76 - 187.0/2100.0,
    dc7 = - 1.0/40.0 ; //end of declaration
    
    
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
    // 1st Step - Not doing, getting passed
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + b21*Step*DyDx[i] ;
    }
    RightHandSide(yTemp, ak2) ;              // 2nd Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*DyDx[i] + b32*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3) ;              // 3rd Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b41*DyDx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4) ;              // 4th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*DyDx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*DyDx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yOut[i] = yIn[i] + Step*(b71*DyDx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yOut, ak7);				//7th and Final step
    
    for(i=0;i<numberOfVariables;i++)
    {
        
        yErr[i] = Step*(dc1*DyDx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] ) ;
        

        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = DyDx[i];
        nextDydx[i] = ak7[i];
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  DormandPrince745::DistChord() const
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


void DormandPrince745::interpolate(  const G4double yInput[],
                                     const G4double dydx[],
                                     G4double yOut[],
                                     G4double Step,
                                 	  G4double tau){
    G4double
    a81 = -24018683.0/8152320000.0 ,
    a82 = 25144.0/43425.0 ,
    a83 = -76360723.0/337557000.0 ,
    a84 = 349808429.0/2445696000.0 ,
    a85 =  -13643731773.0/144024320000.0 ,
    a86 = 1.0/20.0 ,
    a87 = -12268567.0/254760000.0 ,
    
    bf1 = 2104901.0/9204000.0 ,
    bf2 = 0.0 ,
    bf3 = -76360723.0/337557000.0 ,
    bf4 = 349808429.0/2445696000.0 ,
    bf5 = -13268529.0/162604000.0 ,
    bf6 = 13486.0/402675.0 ,
    bf7 = -3162.0/95875.0 ,
    bf8 = -1737.0/3068.0 ;
    
//    G4double *ak8;
    const G4int numberOfVariables= this->GetNumberOfVariables();
    ak8 = new G4double[numberOfVariables];
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(a81*dydx[i] + a82*ak2[i] + a83*ak3[i] +
                                  a84*ak4[i] + a85*ak5[i] + a86*ak6[i] +
                                  a87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);
    
    for( int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*tau*(bf1*dydx[i] + bf2*ak2[i]*pow(tau, 1) +
                                     bf3*ak3[i]*pow(tau, 2) + bf4*ak3[i]*pow(tau, 3)
                                     + bf5*ak5[i]*pow(tau, 4) + bf6*ak6[i]*pow(tau, 5)
                                     + bf7*ak7[i]*pow(tau, 6) + bf8*ak8[i]*pow(tau, 7));
    }
    
}







