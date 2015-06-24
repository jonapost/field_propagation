/**************************************************************************
                        |   Acknowledgement  |
 
     The following code uses the work of J.H.Verner, obtained from
     http://people.math.sfu.ca/~jverner/
     
     Sets of all coefficients provided in attachments are copyrighted
     as such by the author. They many not be published for general 
     distribution. They may be used for any research, industrial 
     application or development of software provided that any
     product arising using any set of coefficients acknowledges this 
     source and includes the URL for this site within the produced 
     item.
 
**************************************************************************/

//  Verner - 9 - 6(5) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
//  Sponsored by Google in Google Summer of Code 2015.
// 
//  First version:  9 June 2015
//
//  This code is made available subject to the Geant4 license, a copy of
//  which is available at
//  http://geant4.org/license
//  
//  History
// ---------------------------------------
//  Created                 09 June 2015    Somnath
//  Added interpolate()     23 June 2015    Somnath
//  Added interpolate6()    - do -          - do -



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
    
    
    //Redundancy here :
    c1 =  .3438957868357036009278820124728322386520e-1 ,
    c2 =  0. ,
    c3 =  0. ,
    c4 =  .2582624555633503404659558098586120858767 ,
    c5 =  .4209371189673537150642551514069801967032 ,
    c6 =  4.405396469669310170148836816197095664891 ,
    c7 = -176.4831190242986576151740942499002125029 ,
    c8 =  172.3641334014150730294022582711902413315 ,
    c9 =  0. ,
    
    
    //Redundancy here :
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
    RightHandSide(yTemp, ak7);              //7th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b81*dydx[i] + b82*ak2[i] + b83*ak3[i] +
                                  b84*ak4[i] + b85*ak5[i] + b86*ak6[i] +
                                  b87*ak7[i]);
    }
    RightHandSide(yTemp, ak8);              //8th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b91*dydx[i] + b92*ak2[i] + b93*ak3[i] +
                                  b94*ak4[i] + b95*ak5[i] + b96*ak6[i] +
                                  b97*ak7[i] + b98*ak8[i] );
    }
    RightHandSide(yTemp, ak9);          //9th Stage
    
//  -------  FSAL NOT IMPLEMENTED  ----------  //
    
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






void VernerRK56::interpolate( const G4double yInput[],
                             const G4double dydx[],
                             G4double yOut[],
                             G4double Step,
                             G4double tau
                             ){
    
    
    G4double *ak10;
    
    G4double
    a101 =  35289331988986254405692535758830683.0/2135620454874580332949729350544993288.0 ,
    a102 = 0.0 ,
    a103 = 0.0 ,
    a104 =  313937014583068512255490687992212890625.0/1028247080705354654473994781524199691557.0 ,
    a105 =  1309307687253621245836726130885318359375.0/6321490412177191231557635904400612215708.0 ,
    a106 = -35295844079877524186147726060781875.0/27279088881521314684841470427640876.0 ,
    a107 =  794353492803973228770716697389421875.0/13906777037439977359946774228636361.0 ,
    a108 = -15228408956329265381787438679500067.0/272520859345009876882656783678732.0 ,
    a109 =  28587810357600962662801.0/1151340224617184234295192.0 ;
    
    
    //
    //  --------------------------------------------------------
    //  COEFFICIENTS FOR INTERPOLANT  bi5  WITH  10  STAGES
    //  --------------------------------------------------------
    //
    
    G4double bi5[11][7], b[11];
    //  COEFFICIENTS OF bi5[1]
    bi5[1][1] =  1.0 ,
    bi5[1][2] = -2834058897718490495086218793721472473.0/533905113718645083237432337636248322.0 ,
    bi5[1][3] =  2718025628974094767211106485595747024.0/266952556859322541618716168818124161.0 ,
    bi5[1][4] = -2007493102587435133656511668645819580.0/266952556859322541618716168818124161.0 ,
    bi5[1][5] =  249346645146318025711596899739877112.0/266952556859322541618716168818124161.0 ,
    bi5[1][6] =  199378106425839009650374224000000000.0/266952556859322541618716168818124161.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[2]
    bi5[2][1] =  0.0 ,
    bi5[2][2] =  0.0 ,
    bi5[2][3] =  0.0 ,
    bi5[2][4] =  0.0 ,
    bi5[2][5] =  0.0 ,
    bi5[2][6] =  0.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[3]
    bi5[3][1] =  0.0 ,
    bi5[3][2] =  0.0 ,
    bi5[3][3] =  0.0 ,
    bi5[3][4] =  0.0 ,
    bi5[3][5] =  0.0 ,
    bi5[3][6] =  0.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[4]
    bi5[4][1] =  0.0 ,
    bi5[4][2] =  2149739120967678287896284375471359375000.0/342749026901784884824664927174733230519.0 ,
    bi5[4][3] = -5492958105397111152037592078122406250000.0/342749026901784884824664927174733230519.0 ,
    bi5[4][4] =  4402390631408885178088267799125656250000.0/342749026901784884824664927174733230519.0 ,
    bi5[4][5] = -56249742645646759461666034659375000000.0/48964146700254983546380703882104747217.0 ,
    bi5[4][6] = -1730712729390963625437276250000000000000.0/1028247080705354654473994781524199691557.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[5]
    bi5[5][1] =  0.0 ,
    bi5[5][2] =  3622473030746576800982284292813464843750.0/526790867681432602629802992033384351309.0 ,
    bi5[5][3] = -38933691634017210049674664360163992187500.0/1580372603044297807889408976100153053927.0 ,
    bi5[5][4] =  17495139028182305773126471135900867187500.0/526790867681432602629802992033384351309.0 ,
    bi5[5][5] = -9216003564492900591852706378813281250000.0/526790867681432602629802992033384351309.0 ,
    bi5[5][6] =  432678182347740906359319062500000000000.0/175596955893810867543267664011128117103.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[6]
    bi5[6][1] =  0.0 ,
    bi5[6][2] = -80801688121532406876813280779008750.0/2273257406793442890403455868970073.0 ,
    bi5[6][3] =  1130047284618441598167544907799477500.0/6819772220380328671210367606910219.0 ,
    bi5[6][4] = -876257846328841227135521923123077500.0/2273257406793442890403455868970073.0 ,
    bi5[6][5] =  1005762761452595148951569429250170000.0/2273257406793442890403455868970073.0 ,
    bi5[6][6] = -138457018351277090034982100000000000.0/757752468931147630134485289656691.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[7]
    bi5[7][1] =  0.0 ,
    bi5[7][2] =  8894101767966865321149886325974625000.0/4635592345813325786648924742878787.0 ,
    bi5[7][3] = -128889699381092513087660977440685250000.0/13906777037439977359946774228636361.0 ,
    bi5[7][4] =  96690747476972701449592103439602750000.0/4635592345813325786648924742878787.0 ,
    bi5[7][5] = -104976825419006157083194997793935000000.0/4635592345813325786648924742878787.0 ,
    bi5[7][6] =  13845701835127709003498210000000000000.0/1545197448604441928882974914292929.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[8]
    bi5[8][1] =  0.0 ,
    bi5[8][2] = -85529300113974351208051144641213185.0/45420143224168312813776130613122.0 ,
    bi5[8][3] =  620054801234124026686518242620266725.0/68130214836252469220664195919683.0 ,
    bi5[8][4] = -464947578142702593618050980843076970.0/22710071612084156406888065306561.0 ,
    bi5[8][5] =  72055052308090805849478606896950124.0/3244295944583450915269723615223.0 ,
    bi5[8][6] = -598331009666258752869007208000000000.0/68130214836252469220664195919683.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[9]
    bi5[9][1] =  0.0 ,
    bi5[9][2] =  5709918156918632012901.0/47972509359049343095633.0 ,
    bi5[9][3] = -17993572040875704216709.0/143917528077148029286899.0 ,
    bi5[9][4] =  85388999974381230343470.0/47972509359049343095633.0 ,
    bi5[9][5] = -223596609894610627617468.0/47972509359049343095633.0 ,
    bi5[9][6] =  415486647330808000000000.0/143917528077148029286899.0 ,
    //  --------------------------------------------------------
    //
    //  COEFFICIENTSOF bi5[10]
    bi5[10][1] =  0.0 ,
    bi5[10][2] = -8.0 ,
    bi5[10][3] =  32.0 ,
    bi5[10][4] = -40.0 ,
    bi5[10][5] =  16.0 ,
    bi5[10][6] =  0.0 ;
    

    
    const G4int numberOfVariables= this->GetNumberOfVariables();

//  Saving yInput because yInput and yOut can be aliases for same array
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    ak10 = new G4double[numberOfVariables];
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    

    
    //    calculating extra stage functions
    for(int i=0; i<6; i++){
        yTemp[i] = yIn[i] + Step*(a101*dydx[i] + a102*ak2[i] + a103*ak3[i] +
                                  a104*ak4[i] + a105*ak5[i] + a106*ak6[i] +
                                  a107*ak7[i] + a108*ak8[i] + a109*ak9[i] );
    }
    
    RightHandSide(yTemp, ak10);
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    for(int i=1; i<=10; i++){   //Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = tau0;
        for(int j=1; j<=6; j++){
            b[i] += bi5[i][j]*tau;
            tau*=tau0;
        }
    }
    
    for(int i=0; i<6; i++){
        yOut[i] = yIn[i] + Step*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i]);
    }
    
}


void VernerRK56::interpolate6( const G4double yInput[],
                             const G4double dydx[],
                             G4double yOut[],
                             G4double Step,
                             G4double tau
                             ){
    
    
    G4double *ak10, *ak11, *ak12;
    
//  ********************************************************
//
//  ADDITIONAL STAGES FOR INTERPOLANT OF ORDER  6
//  ********************************************************

    G4double
    a101 =  35289331988986254405692535758830683.0/2135620454874580332949729350544993288.0 ,
    a102 = 0.0 ,
    a103 = 0.0 ,
    a104 =  313937014583068512255490687992212890625.0/1028247080705354654473994781524199691557.0 ,
    a105 =  1309307687253621245836726130885318359375.0/6321490412177191231557635904400612215708.0 ,
    a106 = -35295844079877524186147726060781875.0/27279088881521314684841470427640876.0 ,
    a107 =  794353492803973228770716697389421875.0/13906777037439977359946774228636361.0 ,
    a108 = -15228408956329265381787438679500067.0/272520859345009876882656783678732.0 ,
    a109 =  28587810357600962662801.0/1151340224617184234295192.0 ,
    

    //  ********************************************************

    //  Coupling coefficients for   c11 =  207/250
    //  --------------------------------------------------------
    
    a111 =  2486392061981208591025761263164027224438868971.0/65173964076983042387381877152862343994140625000.0 ,
    a112 =  0.0 ,
    a113 =  0.0 ,
    a114 =  2330654500023704838558579323179918419669.0/9313832252765893609365894760182968220625.0 ,
    a115 =  5283259505481013273874688940942473187741.0/16258977397575080328080339260289640472500.0 ,
    a116 =  9989685106081485386057729811605187743723.0/5481427003263510055949691042076757812500.0 ,
    a117 = -65815640423883764662985178413751186161.0/971969007022721623945108012714453125.0 ,
    a118 =  183066350554023250298437927498791289370414247.0/2772225538584491748887703284492309570312500.0 ,
    a119 = -426178927623072052719640507155669.0/11712038417736656029207275390625000.0 ,
    a1110 =  3248339841.0/30517578125.0 ,
    
    //  ********************************************************
    //
    //  Coupling coefficients for   c12 =  7.0/25.0 ,
    //  --------------------------------------------------------
    
    a121 =  4676747786898097735038451956075910033997933945857.0/41838231186922043164464169766109251031526972656250.0 ,
    a122 =  0.0 ,
    a123 =  0.0 ,
    a124 =  1320032412954312695441306548681592444623240.0/51248457773784347881352490499724836575577977.0 ,
    a125 =  2087002134582726310861746540254017903014374710.0/551367099344274428347227263044005314054687829.0 ,
    a126 =  3432932836484348829479408524345545011748570706.0/37176735450871998946806722732624135633015625.0 ,
    a127 = -2316434358511265475362584844804601519943610264.0/606481922490173339581866127622363581143375.0 ,
    a128 =  82514605285282414051716141603447021470923168793.0/22107104196177512751528507591142367597656250.0 ,
    a129 = -7560161019374651900153317984708038834.0/7028170531590816328729091157353515625.0 ,
    a1210 = -21655450552377696842870155771710589332.0/6701278878958685336695179940732421875.0 ,
    a1211 = -3194830887993202085244614477336220.0/678662636676110315314332975245759.0 ;


    
    
    
    //  --------------------------------------------------------
    //  COEFFICIENTS FOR INTERPOLANT   bi6   WITH   12  STAGES
    //  --------------------------------------------------------
    //
    
    G4double bi6[13][7], b[13];
    

    //  COEFFICIENTS OF bi6[1]
    bi6[1][1] =  1.0,
    bi6[1][2] = -940811006205413129.0/120948724610397495.0 ,
    bi6[1][3] =  88342864458754360181.0/3265615564480732365.0 ,
    bi6[1][4] = -99667000922033025307.0/2177077042987154910.0 ,
    bi6[1][5] =  7995049273203130972.0/217707704298715491.0 ,
    bi6[1][6] = -7303903485456272500.0/653123112896146473.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[2]
    bi6[2][1] =  0.0 ,
    bi6[2][2] =  0.0 ,
    bi6[2][3] =  0.0 ,
    bi6[2][4] =  0.0 ,
    bi6[2][5] =  0.0 ,
    bi6[2][6] =  0.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[3]
    bi6[3][1] =  0.0 ,
    bi6[3][2] =  0.0 ,
    bi6[3][3] =  0.0 ,
    bi6[3][4] =  0.0 ,
    bi6[3][5] =  0.0 ,
    bi6[3][6] =  0.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[4]
    bi6[4][1] =  0.0 ,
    bi6[4][2] =  2214248281250000.0/133130993475189.0 ,
    bi6[4][3] = -49918013252500000000.0/578720428636646583.0 ,
    bi6[4][4] =  1440368506953125000.0/8387252588936907.0 ,
    bi6[4][5] = -28873797587500000000.0/192906809545548861.0 ,
    bi6[4][6] =  27678103515625000000.0/578720428636646583.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[5]
    bi6[5][1] =  0.0 ,
    bi6[5][2] =  893038428789062500.0/32943296570459319.0 ,
    bi6[5][3] = -125047567320625000000.0/889469007402401613.0 ,
    bi6[5][4] =  82988785418183593750.0/296489669134133871.0 ,
    bi6[5][5] = -72330565909375000000.0/296489669134133871.0 ,
    bi6[5][6] =  69335281738281250000.0/889469007402401613.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[6]
    bi6[6][1] =  0.0 ,
    bi6[6][2] =  40331864555500.0/142160006043.0 ,
    bi6[6][3] = -5647463071672000.0/3838320163161.0 ,
    bi6[6][4] =  3747982556193250.0/1279440054387.0 ,
    bi6[6][5] = -3266630520520000.0/1279440054387.0 ,
    bi6[6][6] =  3131355943750000.0/3838320163161.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[7]
    bi6[7][1] =  0.0 ,
    bi6[7][2] = -143250206750000.0/12603936879.0 ,
    bi6[7][3] =  461347522996000000.0/7827044801859.0 ,
    bi6[7][4] = -13312037070125000.0/113435431911.0 ,
    bi6[7][5] =  266854670860000000.0/2609014933953.0 ,
    bi6[7][6] = -255803940625000000.0/7827044801859.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[8]
    bi6[8][1] =  0.0 ,
    bi6[8][2] =  3753451420391.0/338141155.0 ,
    bi6[8][3] = -3679035166143248.0/63908678295.0 ,
    bi6[8][4] =  4883240297928691.0/42605785530.0 ,
    bi6[8][5] = -425608752364336.0/4260578553.0 ,
    bi6[8][6] =  407983850042500.0/12781735659.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[9]
    bi6[9][1] =  0.0 ,
    bi6[9][2] = -69713.0/23220.0 ,
    bi6[9][3] =  4685161.0/313470.0 ,
    bi6[9][4] = -135239.0/4860.0 ,
    bi6[9][5] =  228046.0/10449.0 ,
    bi6[9][6] = -186250.0/31347.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[10]
    bi6[10][1] =  0.0 ,
    bi6[10][2] = -132664.0/6765.0 ,
    bi6[10][3] =  17011336.0/182655.0 ,
    bi6[10][4] = -10067296.0/60885.0 ,
    bi6[10][5] =  1579832.0/12177.0 ,
    bi6[10][6] = -1385000.0/36531.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[11]
    bi6[11][1] =  0.0 ,
    bi6[11][2] = -2734375000.0/149990751.0 ,
    bi6[11][3] =  391796875000.0/4049750277.0 ,
    bi6[11][4] = -6250000000.0/31393413.0 ,
    bi6[11][5] =  244140625000.0/1349916759.0 ,
    bi6[11][6] = -244140625000.0/4049750277.0 ,
    //  -------------------------------------------------------
    //
    //  COEFFICIENTS OF bi6[12]
    bi6[12][1] =  0.0 ,
    bi6[12][2] = -15453125.0/1139292.0 ,
    bi6[12][3] =  1393796875.0/15380442.0 ,
    bi6[12][4] = -2092203125.0/10253628.0 ,
    bi6[12][5] =  488281250.0/2563407.0 ,
    bi6[12][6] = -488281250.0/7690221.0 ;
    
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    //  Saving yInput because yInput and yOut can be aliases for same array
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    ak10 = new G4double[numberOfVariables];
    ak11 = new G4double[numberOfVariables];
    ak12 = new G4double[numberOfVariables];
    
    // The number of variables to be integrated over
    yOut[7] = yTemp[7]  = yIn[7];
    
    
    
    //    calculating extra stage functions
    for(int i=0; i<6; i++){
        yTemp[i] = yIn[i] + Step*(a101*dydx[i] + a102*ak2[i] + a103*ak3[i] +
                                  a104*ak4[i] + a105*ak5[i] + a106*ak6[i] +
                                  a107*ak7[i] + a108*ak8[i] + a109*ak9[i] );
    }
    
    RightHandSide(yTemp, ak10);
    
    for(int i=0; i<6; i++){
        yTemp[i] = yIn[i] + Step*(a111*dydx[i] + a112*ak2[i] + a113*ak3[i] +
                                  a114*ak4[i] + a115*ak5[i] + a116*ak6[i] +
                                  a117*ak7[i] + a118*ak8[i] + a119*ak9[i] +
                                  a1110*ak10[i]);
    }
    
    RightHandSide(yTemp, ak11);
    
    for(int i=0; i<6; i++){
        yTemp[i] = yIn[i] + Step*(a121*dydx[i] + a122*ak2[i] + a123*ak3[i] +
                                  a124*ak4[i] + a125*ak5[i] + a126*ak6[i] +
                                  a127*ak7[i] + a128*ak8[i] + a129*ak9[i] +
                                  a1210*ak10[i] + a1211*ak11[i]);
    }
    RightHandSide(yTemp, ak12);
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    for(int i=1; i<=12; i++){   //Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = tau0;
        for(int j=1; j<=6; j++){
            b[i] += bi6[i][j]*tau;
            tau*=tau0;
        }
    }
    
    for(int i=0; i<6; i++){
        yOut[i] = yIn[i] + Step*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i] + b[11]*ak11[i] + b[12]*ak12[i]);
    }
    
}


