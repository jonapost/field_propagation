//  BogackiShampine45.cc
//  Geant4
//
//  Created by Somnath on 25/05/15.
//
//This is the source file of BogackiShampine45 class containing the
//definition of the stepper() method that evaluates one step in
//field propagation.
//The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
//
//0   |
//1/6 | 1/6
//2/9 | 2/27	   	  4/27
//3/7 | 183/1372  -162/343   1053/1372
//2/3 | 68/297      -4/11      42/143 1960/3861
//3/4 | 597/22528     81/352     63099/585728     58653/366080   4617/20480
//1	  | 174197/959244 -30942/79937 8152137/19744439 666106/1039181 -29421/29068 482048/414219
//1   | 587/8064           0       4440339/15491840 24353/124800   387/44800     2152/5985   7267/94080
//-------------------------------------------------------------------------------------------------------------------
//      587/8064           0       4440339/15491840 24353/124800    387/44800     2152/5985   7267/94080       0
//      2479/34992         0          123/416       612941/3411720  43/1440       2272/6561  79937/1113912  3293/556956


#include "BogackiShampine45.hh"
#include "G4LineSection.hh"


//Constructor
BogackiShampine45::BogackiShampine45(G4EquationOfMotion *EqRhs,
                                     G4int noIntegrationVariables,
                                     G4bool primary)
: G4MagIntegratorStepper(EqRhs, noIntegrationVariables){
    
    const G4int numberOfVariables = noIntegrationVariables;
    
    //New Chunk of memory being created for use by the stepper
    
    //aki - for storing intermediate RHS
    ak2 = new G4double[numberOfVariables];
    ak3 = new G4double[numberOfVariables];
    ak4 = new G4double[numberOfVariables];
    ak5 = new G4double[numberOfVariables];
    ak6 = new G4double[numberOfVariables];
    ak7 = new G4double[numberOfVariables];
    ak8 = new G4double[numberOfVariables];
    yTemp = new G4double[numberOfVariables];
    yIn = new G4double[numberOfVariables] ;
    
    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fLastDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    if( primary )
    {
        fAuxStepper = new BogackiShampine45(EqRhs, numberOfVariables,
                                            !primary);
    }
}


//Destructor
BogackiShampine45::~BogackiShampine45(){
    //clear all previously allocated memory for stepper and DistChord
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


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively

void BogackiShampine45::Stepper(const G4double yInput[],
                                const G4double dydx[],
                                G4double Step,
                                G4double yOut[],
                                G4double yErr[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
    b21 = 1.0/6.0 ,
    b31 = 2.0/27.0 , b32 = 4.0/27.0,
    
    b41 = 183.0/1372.0 , b42 = -162.0/343.0, b43 = 1053.0/1372.0,
    
    b51 = 68.0/297.0, b52 = -4.0/11.0,
    b53 = 42.0/143.0, b54 = 1960.0/3861.0,
    
    b61 = 597.0/22528.0, b62 = 81.0/352.0,
    b63 = 63099.0/585728.0, b64 = 58653.0/366080.0,
    b65 = 4617.0/20480.0,
    
    b71 = 174197.0/959244.0, b72 = -30942.0/79937.0,
    b73 = 8152137.0/19744439.0, b74 = 666106.0/1039181.0,
    b75 = -29421.0/29068.0,  b76 = 482048.0/414219.0,
    
    b81 = 587.0/8064.0,  b82 = 0.0,
    b83 = 4440339.0/15491840.0,  b84 = 24353.0/124800.0,
    b85 = 387.0/44800.0, b86 = 2152.0/5985.0,
    b87 = 7267.0/94080.0,
    
    
    c1 = 2479.0/34992.0,
    c2 = 0.0,
    c3 = 123.0/416.0,
    c4 = 612941.0/3411720.0,
    c5 = 43.0/1440.0,
    c6 = 2272.0/6561.0,
    c7 = 79937.0/1113912.0,
    c8 = 3293.0/556956.0,
    
    //For the embedded higher order method only the difference of values
    // taken and is used directly later instead of defining the last row
    // of butcher table in a separate set of variables and taking the
    // difference there
    
    dc1 = -3.0/1280.0,
    dc2 = 0.0,
    dc3 = 6561.0/632320.0,
    dc4 = -343.0/20800.0,
    dc5 = 243.0/12800.0,
    dc6 = -1.0/95.0,
    dc7 = 0.0   ;   //end of declaration
    
    
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
        yTemp[i] = yIn[i] + Step*(b71*dydx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yTemp, ak7);				//7th Step
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] + b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);				//7th Step

    
    for(i=0;i<numberOfVariables;i++)
    {
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak4[i] + c3*ak3[i] +
                                 c4*ak4[i] + c5*ak5[i] + c6*ak6[i] + c7*ak7[i] + c8*ak8[i]) ;
        
        // Estimate error as difference between 4th and
        // 5th order methods
        
        yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] ) ;
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  BogackiShampine45::DistChord() const
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










