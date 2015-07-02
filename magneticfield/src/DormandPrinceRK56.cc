#include "DormandPrinceRK56.hh"
#include "G4LineSection.hh"


//Constructor
DormandPrinceRK56::DormandPrinceRK56(G4EquationOfMotion *EqRhs,
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
//    ak9 = new G4double[numberOfVariables];
    
    yTemp = new G4double[numberOfVariables];
    yIn = new G4double[numberOfVariables] ;
    
    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fLastDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    if( primary )
    {
        fAuxStepper = new DormandPrinceRK56(EqRhs, numberOfVariables,
                                     !primary);
    }
}


//Destructor
DormandPrinceRK56::~DormandPrinceRK56(){
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    delete[] ak8;
//    delete[] ak9;
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

void DormandPrinceRK56::Stepper(const G4double yInput[],
                         const G4double dydx[],
                         G4double Step,
                         G4double yOut[],
                         G4double yErr[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
    
    b21 =  1.0/10.0 ,
    
    b31 =  -2.0/81.0 ,
    b32 =  20.0/81.0 ,
    
    b41 =  615.0/1372.0 ,
    b42 =  -270.0/343.0 ,
    b43 =  1053.0/1372.0 ,
    
    b51 =  3243.0/5500.0 ,
    b52 =  -54.0/55.0 ,
    b53 = 50949.0/71500.0 ,
    b54 =  4998.0/17875.0 ,
    
    b61 = -26492.0/37125.0 ,
    b62 =  72.0/55.0 ,
    b63 =  2808.0/23375.0 ,
    b64 = -24206.0/37125.0 ,
    b65 =  338.0/459.0 ,
    
    b71 = 5561.0/2376.0 ,
    b72 =  -35.0/11.0 ,
    b73 =  -24117.0/31603.0 ,
    b74 = 899983.0/200772.0 ,
    b75 =  -5225.0/1836.0 ,
    b76 = 3925.0/4056.0 ,
    
    b81 = 465467.0/266112.0 ,
    b82 = -2945.0/1232.0 ,
    b83 = -5610201.0/14158144.0 ,
    b84 =  10513573.0/3212352.0 ,
    b85 = -424325.0/205632.0 ,
    b86 = 376225.0/454272.0 ,
    b87 = 0.0 ,
    
    c1 =  61.0/864.0 ,
    c2 =  0.0 ,
    c3 =  98415.0/321776.0 ,
    c4 =  16807.0/146016.0 ,
    c5 =  1375.0/7344.0 ,
    c6 =  1375.0/5408.0 ,
    c7 = -37.0/1120.0 ,
    c8 =  1.0/10.0 ,
    
    dc1 =  c1  - 821.0/10800.0 ,
    dc2 =  c2 - 0.0 ,
    dc3 =  c3 - 19683.0/71825,
    dc4 =  c4 - 175273.0/912600.0 ,
    dc5 =  c5 - 395.0/3672.0 ,
    dc6 =  c6 - 785.0/2704.0 ,
    dc7 =  c7 - 3.0/50.0 ,
    dc8 =  c8 - 0.0 ;
        //end of declaration
    
    
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
    RightHandSide(yTemp, ak2) ;              // 2nd Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3) ;              // 3rd Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4) ;              // 4th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                                  b54*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5) ;              // 5th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6) ;              // 6th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b71*dydx[i] + b72*ak2[i] + b73*ak3[i] +
                                  b74*ak4[i] + b75*ak5[i] + b76*ak6[i]);
    }
    RightHandSide(yTemp, ak7);				//7th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                  b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);				//8th Stage
    
//    for(i=0;i<numberOfVariables;i++)
//    {
//        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
//                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
//                                  b97*ak7[i] + b98*ak8[i] );
//    }
//    RightHandSide(yTemp, ak9);          //9th Stage
    
    
    for(i=0;i<numberOfVariables;i++)
    {
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak2[i] + c3*ak3[i] +
                                 c4*ak4[i] + c5*ak5[i] + c6*ak6[i] +
                                 c7*ak7[i] + c8*ak8[i]  ) ;
        
        // Estimate error as difference between 5th and
        // 6th order methods
        
        yErr[i] = Step*( dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                         dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i]  ) ;
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  DormandPrinceRK56::DistChord() const
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
