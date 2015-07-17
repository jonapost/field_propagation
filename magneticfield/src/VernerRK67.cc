/*
 ********************************************************************
 // Acknowledgement
 //
 // The following code uses the work of Jim Verner, obtained from
 // http://people.math.sfu.ca/~jverner/
 //
 // Sets of all coefficients provided in attachments are copyrighted
 // as such by the author. They many not be published for general
 // distribution. They may be used for any research, industrial
 // application or development of software provided that any
 // product arising using any set of coefficients acknowledges this
 // source and includes the URL for this site within the produced
 // item.
 
 ********************************************************************
 */


#include "VernerRK67.hh"
#include "G4LineSection.hh"


//Constructor
VernerRK67::VernerRK67(G4EquationOfMotion *EqRhs,
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
    ak9 = new G4double[numberOfVariables];
    ak10 = new G4double[numberOfVariables];
    yTemp = new G4double[numberOfVariables];
    yIn = new G4double[numberOfVariables] ;
    
    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fLastDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    if( primary )
    {
        fAuxStepper = new VernerRK67(EqRhs, numberOfVariables,
                                            !primary);
    }
}


//Destructor
VernerRK67::~VernerRK67(){
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    delete[] ak8;
    delete[] ak9;
    delete[] ak10;
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

void VernerRK67::Stepper(const G4double yInput[],
                         const G4double dydx[],
                                G4double Step,
                                G4double yOut[],
                                G4double yErr[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
  
    b21 =  .5e-2 ,

    b31 = -1.076790123456790123456790123456790123457 ,
    b32 =  1.185679012345679012345679012345679012346 ,

    b41 =  .4083333333333333333333333333333333333333e-1 ,
    b42 =  0. ,
    b43 =  .1225 ,

    b51 =  .6389139236255726780508121615993336109954 ,
    b52 =  0. ,
    b53 = -2.455672638223656809662640566430653894211 ,
    b54 =  2.272258714598084131611828404831320283215 ,

    b61 = -2.661577375018757131119259297861818119279 ,
    b62 =  0. ,
    b63 =  10.80451388645613769565396655365532838482 ,
    b64 = -8.353914657396199411968048547819291691541 ,
    b65 =  .8204875949566569791420417341743839209619 ,

    b71 =  6.067741434696770992718360183877276714679 ,
    b72 =  0. ,
    b73 = -24.71127363591108579734203485290746001803 ,
    b74 =  20.42751793078889394045773111748346612697 ,
    b75 = -1.906157978816647150624096784352757010879 ,
    b76 =  1.006172249242068014790040335899474187268 ,

    b81 =  12.05467007625320299509109452892778311648 ,
    b82 =  0. ,
    b83 = -49.75478495046898932807257615331444758322 ,
    b84 =  41.14288863860467663259698416710157354209 ,
    b85 = -4.461760149974004185641911603484815375051 ,
    b86 =  2.042334822239174959821717077708608543738 ,
    b87 = -0.9834843665406107379530801693870224403537e-1 ,

    b91 =  10.13814652288180787641845141981689030769 ,
    b92 =  0. ,
    b93 = -42.64113603171750214622846006736635730625 ,
    b94 =  35.76384003992257007135021178023160054034 ,
    b95 = -4.348022840392907653340370296908245943710 ,
    b96 =  2.009862268377035895441943593011827554771 ,
    b97 =  .3487490460338272405953822853053145879140 ,
    b98 = -.2714390051048312842371587140910297407572 ,

    b101 = -45.03007203429867712435322405073769635151 ,
    b102 =  0. ,
    b103 =  187.3272437654588840752418206154201997384 ,
    b104 = -154.0288236935018690596728621034510402582 ,
    b105 =  18.56465306347536233859492332958439136765 ,
    b106 = -7.141809679295078854925420496823551192821 ,
    b107 =  1.308808578161378625114762706007696696508 ,
    b108 =  0. ,
    b109 =  0. ,

    c1 =  .4715561848627222170431765108838175679569e-1 ,
    c2 =  0. ,
    c3 =  0. ,
    c4 =  .2575056429843415189596436101037687580986 ,
    c5 =  .2621665397741262047713863095764527711129 ,
    c6 =  .1521609265673855740323133199165117535523 ,
    c7 =  .4939969170032484246907175893227876844296 ,
    c8 = -.2943031171403250441557244744092703429139 ,
    c9 =  .8131747232495109999734599440136761892478e-1 ,
    c10 =  0. ,

    
//Remove the stupidity from below ci - XYZ
    dc1 = .4715561848627222170431765108838175679569e-1 - .4460860660634117628731817597479197781432e-1 ,
    dc2 =  0. ,
    dc3 =  0. ,
    dc4 = .2575056429843415189596436101037687580986 - .2671640378571372680509102260943837899738 ,
    dc5 = .2621665397741262047713863095764527711129 - .2201018300177293019979715776650753096323 ,
    dc6 = .1521609265673855740323133199165117535523 - .2188431703143156830983120833512893824578 ,
    dc7 = .4939969170032484246907175893227876844296 - .2289871705411202883378173889763552365362 ,
    dc8 = -.2943031171403250441557244744092703429139 ,
    dc9 = .8131747232495109999734599440136761892478e-1 ,
    dc10 = - .2029518466335628222767054793810430358554e-1 ;

//end of declaration !
    
    
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
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] + 
                                  b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);				//8th Step

     for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] + 
                                  b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yTemp, ak9);          //9th Step 

     for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b101*dydx[i] + b102*ak2[i] + b103*ak3[i] +
                                  b104*ak4[i] + b105*ak5[i] + b106*ak6[i] + 
                                  b107*ak7[i] + b108*ak8[i] + b109*ak9[i]);
    }
    RightHandSide(yTemp, ak10);          //10th Step 


    
    for(i=0;i<numberOfVariables;i++)
    {
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak2[i] + c3*ak3[i] +
                                 c4*ak4[i] + c5*ak5[i] + c6*ak6[i] + 
                                 c7*ak7[i] + c8*ak8[i] +c9*ak9[i] + c10*ak10[i]) ;
        
        // Estimate error as difference between 4th and
        // 5th order methods
        
        yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i] + 
                        dc9*ak9[i] + dc10*ak10[i]) ;
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  VernerRK67::DistChord() const
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










