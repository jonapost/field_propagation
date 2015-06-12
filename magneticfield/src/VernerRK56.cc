/*
 ********************************************************************
// Acknowledgement
// 
// The following code uses the work of J.H.Verner, obtained from
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



#include "VernerRK56.hh"
#include "G4LineSection.hh"


//Constructor
VernerRK56::VernerRK56(G4EquationOfMotion *EqRhs,
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

    yTemp = new G4double[numberOfVariables];
    yIn = new G4double[numberOfVariables] ;
    
    fLastInitialVector = new G4double[numberOfVariables] ;
    fLastFinalVector = new G4double[numberOfVariables] ;
    fLastDyDx = new G4double[numberOfVariables];
    
    fMidVector = new G4double[numberOfVariables];
    fMidError =  new G4double[numberOfVariables];
    if( primary )
    {
        fAuxStepper = new VernerRK56(EqRhs, numberOfVariables,
                                     !primary);
    }
}


//Destructor
VernerRK56::~VernerRK56(){
    //clear all previously allocated memory for stepper and DistChord
    delete[] ak2;
    delete[] ak3;
    delete[] ak4;
    delete[] ak5;
    delete[] ak6;
    delete[] ak7;
    delete[] ak8;
    delete[] ak9;
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

void VernerRK56::Stepper(const G4double yInput[],
                         const G4double dydx[],
                         G4double Step,
                         G4double yOut[],
                         G4double yErr[])
{
    G4int i;
    
    //The various constants defined on the basis of butcher tableu
    const G4double  //G4double - only once
    
    
    b21 =  .6e-1 ,
    
    b31 =  .1923996296296296296296296296296296296296e-1 ,
    b32 =  .7669337037037037037037037037037037037037e-1 ,
    
    b41 =  .35975e-1 ,
    b42 =  0. ,
    b43 =  .107925 ,
    
    b51 =  1.318683415233148260919747276431735612861 ,
    b52 =  0. ,
    b53 = -5.042058063628562225427761634715637693344 ,
    b54 =  4.220674648395413964508014358283902080483 ,
    
    b61 = -41.87259166432751461803757780644346812905 ,
    b62 =  0. ,
    b63 =  159.4325621631374917700365669070346830453 ,
    b64 = -122.1192135650100309202516203389242140663 ,
    b65 =  5.531743066200053768252631238332999150076 ,
    
    b71 = -54.43015693531650433250642051294142461271 ,
    b72 =  0. ,
    b73 =  207.0672513650184644273657173866509835987 ,
    b74 = -158.6108137845899991828742424365058599469 ,
    b75 =  6.991816585950242321992597280791793907096 ,
    b76 = -.1859723106220323397765171799549294623692e-1 ,
    
    b81 = -54.66374178728197680241215648050386959351 ,
    b82 =  0. ,
    b83 =  207.9528062553893734515824816699834244238 ,
    b84 = -159.2889574744995071508959805871426654216 ,
    b85 =  7.018743740796944434698170760964252490817 ,
    b86 = -.1833878590504572306472782005141738268361e-1 ,
    b87 = -.5119484997882099077875432497245168395840e-3 ,
    
    b91 =  .3438957868357036009278820124728322386520e-1 ,
    b92 =  0. ,
    b93 =  0. ,
    b94 =  .2582624555633503404659558098586120858767 ,
    b95 =  .4209371189673537150642551514069801967032 ,
    b96 =  4.405396469669310170148836816197095664891 ,
    b97 = -176.4831190242986576151740942499002125029 ,
    b98 =  172.3641334014150730294022582711902413315 ,
    
    c1 =  .3438957868357036009278820124728322386520e-1 ,
    c2 =  0. ,
    c3 =  0. ,
    c4 =  .2582624555633503404659558098586120858767 ,
    c5 =  .4209371189673537150642551514069801967032 ,
    c6 =  4.405396469669310170148836816197095664891 ,
    c7 = -176.4831190242986576151740942499002125029 ,
    c8 =  172.3641334014150730294022582711902413315 ,
    c9 =  0. ,
    
    dc1 =  .3438957868357036009278820124728322386520e-1  - .4909967648382489730906854927971225836479e-1 ,
    dc2 =  0. ,
    dc3 =  0. ,
    dc4 =  .2582624555633503404659558098586120858767 - .2251112229516524153401395320539875329485 ,
    dc5 =  .4209371189673537150642551514069801967032 - .4694682253029562039431948525047387412553 ,
    dc6 =  4.405396469669310170148836816197095664891 - .8065792249988867707634161808995217981443 ,
    dc7 = -176.4831190242986576151740942499002125029 ,
    dc8 =  172.3641334014150730294022582711902413315 + .6071194891777959797672951465256217122488 ,
    dc9 =  - .5686113944047569241147603178766138153594e-1 ;
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
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                  b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yTemp, ak9);          //9th Stage
    
    
    
    for(i=0;i<numberOfVariables;i++)
    {
        // Accumulate increments with proper weights
        
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak2[i] + c3*ak3[i] +
                                 c4*ak4[i] + c5*ak5[i] + c6*ak6[i] +
                                 c7*ak7[i] + c8*ak8[i] +c9*ak9[i] ) ;
        
        // Estimate error as difference between 5th and
        // 6th order methods
        
        yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i] +
                        dc5*ak5[i] + dc6*ak6[i] + dc7*ak7[i] + dc8*ak8[i] +
                        dc9*ak9[i] ) ;
        
        
    }
    
    fLastStepLength = Step;
    
    return ;
}


//The following has not been tested

//The DistChord() function fot the class - must define it here.
G4double  VernerRK56::DistChord() const
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










