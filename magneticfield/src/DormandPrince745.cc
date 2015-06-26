
//  DormandPrince7 - 5(4) implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 25 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DormandPrince745.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 25/05/15.


//This is the source file of DormandPrince745 class containing the 
//definition of the stepper() method that evaluates one step in 
//field propagation.
//The Butcher table of the Dormand-Prince-7-4-5 method is as follows :
// 
// 0   |
// 1/5 | 1/5
// 3/10| 3/40        9/40
// 4/5 | 44/45      −56/15      32/9
// 8/9 | 19372/6561 −25360/2187 64448/6561 −212/729
// 1   | 9017/3168  −355/33    46732/5247  49/176  −5103/18656
// 1   | 35/384      0         500/1113    125/192 −2187/6784    11/84   
// ------------------------------------------------------------------------     
//      35/384       0        500/1113    125/192  −2187/6784    11/84   0
//      5179/57600   0       7571/16695  393/640  −92097/339200 187/2100 1/40



#include "DormandPrince745.hh"
#include "G4LineSection.hh"


//Constructor
DormandPrince745::DormandPrince745(G4EquationOfMotion *EqRhs,
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
    ak8 = 0;
    yTemp = new G4double[numberOfVariables] ;
    yIn = new G4double[numberOfVariables] ;
    
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
    
}


//Stepper :

// Passing in the value of yInput[],the first time dydx[] and Step length
// Giving back yOut and yErr arrays for output and error respectively
void DormandPrince745::Stepper(const G4double yInput[],
                         const G4double dydx[],
                         G4double Step,
                         G4double yOut[],
                         G4double yErr[])
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

    c1 = 35.0/384.0, c2 = .0,
    c3 = 500.0/1113.0, c4 = 125.0/192.0,
    c5 = -2187.0/6784.0, c6 = 11.0/84.0,
    c7 = 0,

    //For the embedded higher order method only the difference of values
    // taken and is used directly later instead of defining the last row
    // of butcher table in a separate set of variables and taking the 
    // difference there
    dc1 = c1 - 5179.0/57600.0,
    dc2 = c2 - .0,
    dc3 = c3 - 7571.0/16695.0,
    dc4 = c4 - 393.0/640.0,
    dc5 = c5 + 92097.0/339200.0,
    dc6 = c6 - 187.0/2100.0,
    dc7 = c7 - 1.0/40.0 ; //end of declaration
    
    
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
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak4[i] + c3*ak3[i] + 
                                 c4*ak4[i] + c5*ak5[i] + c6*ak6[i] + c7*ak7[i] ) ;
        
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

/*
    Code for interpolation. Formula from :
    (Ref.)
    Dormand, J. R., & Prince, P. J. (1986). Runge-Kutta triples. 
    Computers & Mathematics with Applications, 12(9), 1007–1017. 
    doi:10.1016/0898-1221(86)90025-8
*/


void DormandPrince745::interpolate(  const G4double yInput[],
                                     const G4double dydx[],
                                     G4double yOut[],       
                                     G4double Step,
                                      G4double tau){
    

    
    
    const G4int numberOfVariables= this->GetNumberOfVariables();

    G4double tau0 = tau;        //Just saving the value of tau - for changes later
    
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    G4double
    tau_2 = tau0*tau0 ,
    tau_3 = tau0*tau_2,
    tau_4 = tau_2*tau_2,
    
    //Calculating the stage coefficient polynomials bi(tau)
    bi1 = (157015080.0*tau_4 - 13107642775.0*tau_3+ 34969693132.0*tau_2- 32272833064.0*tau
           + 11282082432.0)/11282082432.0,
    bi2 = 0.0 ,
    bi3 = - 100.0*tau*(15701508.0*tau_3 - 914128567.0*tau_2 + 2074956840.0*tau
                 - 1323431896.0)/32700410799.0,
    bi4 = 25.0*tau*(94209048.0*tau_3- 1518414297.0*tau_2+ 2460397220.0*tau - 889289856.0)/5641041216.0 ,
    bi5 = -2187.0*tau*(52338360.0*tau_3 - 451824525.0*tau_2 + 687873124.0*tau - 259006536.0)/199316789632.0 ,
    bi6 =  11.0*tau*(106151040.0*tau_3- 661884105.0*tau_2 + 946554244.0*tau - 361440756.0)/2467955532.0 ,
    bi7 = tau*(1.0 - tau)*(8293050.0*tau_2 - 82437520.0*tau + 44764047.0)/ 29380423.0 ;
    
    //Getting the value of yOut using the formula for y(n+tau)
    for( int i=0; i<numberOfVariables; i++){
        yOut[i] = yIn[i] + Step*tau*(bi1*dydx[i] + bi2*ak2[i] + bi3*ak3[i] + bi4*ak4[i]
                                     + bi5*ak5[i] + bi6*ak6[i] + bi7*ak7[i]  ) ;
    }
        
}






