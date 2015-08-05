//  Verner's RK 6(5) non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 9 June 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  VernerRK56.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 9 June 2015
//
//
///////////////////////////////////////////////////////////////////////////////


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


void VernerRK67::SetupInterpolate_low( const G4double yInput[],
                                  const G4double dydx[],
                                  const G4double Step ){
    
    G4double
    b111 =  .4715561848627222170431765108838175679569e-1 ,
    b112 =  0. ,
    b113 =  0. ,
    b114 =  .2575056429843415189596436101037687580986 ,
    b115 =  .2621665397741262047713863095764527711129 ,
    b116 =  .1521609265673855740323133199165117535523 ,
    b117 =  .4939969170032484246907175893227876844296 ,
    b118 = -.2943031171403250441557244744092703429139 ,
    b119 =  .8131747232495109999734599440136761892478e-1 ,
    b1110 =  0. ,
    
    b121 =  .5232227691599689815470932256735029887614e-1 ,
    b122 =  0. ,
    b123 =  0. ,
    b124 =  .2249586182670571550244187743667190903405 ,
    b125 =  .1744370924877637539031751304611402542578e-1 ,
    b126 = -.7669379876829393188009028209348812321417e-2 ,
    b127 =  .3435896044073284645684381456417912794447e-1 ,
    b128 = -.4102097230093949839125144540100346681769e-1 ,
    b129 =  .2565113300520561655297104906598973655221e-1 ,
    b1210 =  0. ,
    b1211 = -.1604434570000000000000000000000000000000e-1 ,
    
    b131 =  .5305334125785908638834747243817578898946e-1 ,
    b132 =  0. ,
    b133 =  0. ,
    b134 =  .1219530101140188607092225622195251463666 ,
    b135 =  .1774684073760249704011573985936092552347e-1 ,
    b136 = -.5928372667681494328907467430302313286925e-3 ,
    b137 =  .8381833970853750873624781948796072714855e-2 ,
    b138 = -.1293369259698611956700998079778496462996e-1 ,
    b139 =  .9412056815253860804791356641605087829772e-2 ,
    b1310 =  0. ,
    b1311 = -.5353253107275676032399320754008272222345e-2 ,
    b1312 = -.6666729992455811078380186481263955324311e-1 ;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    G4int i;
    
    //  Saving yInput because yInput and yOut can be aliases for same array
    for( i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    yTemp[7]  = yIn[7];
    
    ak11 = new G4double[numberOfVariables];
    ak12 = new G4double[numberOfVariables];
    ak13 = new G4double[numberOfVariables];
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b111*dydx[i] + b112*ak2[i] + b113*ak3[i] +
                                  b114*ak4[i] + b115*ak5[i] + b116*ak6[i] +
                                  b117*ak7[i] + b118*ak8[i] + b119*ak9[i] +
                                  b1110*ak10[i]);
    }
    RightHandSide(yTemp, ak11);			//11th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b121*dydx[i] + b122*ak2[i] + b123*ak3[i] +
                                  b124*ak4[i] + b125*ak5[i] + b126*ak6[i] +
                                  b127*ak7[i] + b128*ak8[i] + b129*ak9[i] +
                                  b1210*ak10[i] + b1211*ak11[i]);
    }
    RightHandSide(yTemp, ak12);			//12th Stage
    
    for(i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b131*dydx[i] + b132*ak2[i] + b133*ak3[i] +
                                  b134*ak4[i] + b135*ak5[i] + b136*ak6[i] +
                                  b137*ak7[i] + b138*ak8[i] + b139*ak9[i] +
                                  b1310*ak10[i] + b1311*ak11[i] + b1312*ak12[i]);
    }
    RightHandSide(yTemp, ak13);			//13th Stage
    
    
}


void VernerRK67::Interpolate_low( const G4double yInput[],
                     const G4double dydx[],
                     const G4double Step,
                     G4double yOut[],
                     G4double tau ){
    G4double bi6[14][7], b[14];
    G4int numberOfVariables = this->GetNumberOfVariables();
    
    bi6[1][1] =  1. ,
    bi6[1][2] = -7.579486522562013856370387489358761917136 ,
    bi6[1][3] =  24.84859042701758998305254381572114767883 ,
    bi6[1][4] = -38.85067748922540724491934058373805286459 ,
    bi6[1][5] =  28.75646349856373329638850558427800090199 ,
    bi6[1][6] = -8.127734295307629956447003675813952042307 ,
    
    bi6[2][1] =  0. ,
    bi6[2][2] =  0. ,
    bi6[2][3] =  0. ,
    bi6[2][4] =  0. ,
    bi6[2][5] =  0. ,
    bi6[2][6] =  0. ,
    
    bi6[3][1] =  0. ,
    bi6[3][2] =  0. ,
    bi6[3][3] =  0. ,
    bi6[3][4] =  0. ,
    bi6[3][5] =  0. ,
    bi6[3][6] =  0. ,
    
    bi6[4][1] =  0. ,
    bi6[4][2] =  4.551232240400547830289289496617659221741 ,
    bi6[4][3] = -41.74306197989451542904654379577378740985 ,
    bi6[4][4] =  125.9208614672431114278470271124781332255 ,
    bi6[4][5] = -143.2724320984988187759537191634829285599 ,
    bi6[4][6] =  54.80090601373401646582358996026469228056 ,
    
    bi6[5][1] =  0. ,
    bi6[5][2] =  2.347388766837311148481149416205748182220 ,
    bi6[5][3] = -21.42965001437327943727021057869980940196 ,
    bi6[5][4] =  64.01102018753913070410229000503296635821 ,
    bi6[5][5] = -71.54964616066291046169022808133078061275 ,
    bi6[5][6] =  26.88305376043387425114838554836832824539 ,
    
    bi6[6][1] =  0. ,
    bi6[6][2] =  .6629628602943922170984684021172226061984 ,
    bi6[6][3] = -5.991791209485099663334522931412659054138 ,
    bi6[6][4] =  17.51358806151491209369363907599542681959 ,
    bi6[6][5] = -18.79068837634778062158370304672269638025 ,
    bi6[6][6] =  6.758089590590961548158431819939217762152 ,
    
    bi6[7][1] =  0. ,
    bi6[7][2] = -1.894931483197030289667134961741273285056 ,
    bi6[7][3] =  17.84551891193619310038245879582238106213 ,
    bi6[7][4] = -56.77440614878830815775145830865327266897 ,
    bi6[7][5] =  70.55596299657564872116838561274122139834 ,
    bi6[7][6] = -29.23814735952325494944153354884626882202 ,
    
    bi6[8][1] =  0. ,
    bi6[8][2] =  1.489077233496668566079506992903531187388 ,
    bi6[8][3] = -13.95068086791358924071990604089950620891 ,
    bi6[8][4] =  43.93573738643896559384998172064294274652 ,
    bi6[8][5] = -53.74155980596578799479262013665711367333 ,
    bi6[8][6] =  21.97312293680341803142731298960087560542 ,
    
    bi6[9][1] =  0. ,
    bi6[9][2] = -.5934749977615343231780548830714052210176 ,
    bi6[9][3] =  5.532214575131257019966932786240153623467 ,
    bi6[9][4] = -17.25057884887540359301582894972755244557 ,
    bi6[9][5] =  20.76631879735288001882715503942847061836 ,
    bi6[9][6] = -8.373162053522248022602857998468298956311 ,
    
    bi6[10][1] =  0. ,
    bi6[10][2] =  0. ,
    bi6[10][3] =  0. ,
    bi6[10][4] =  0. ,
    bi6[10][5] =  0. ,
    bi6[10][6] =  0. ,
    
    bi6[11][1] =  0. ,
    bi6[11][2] =  .5229705513661462700547131555553062786922 ,
    bi6[11][3] = -4.935085895387111016017254136967884797119 ,
    bi6[11][4] =  15.77686659404020360951757145807650520227 ,
    bi6[11][5] = -19.84035770738365925120223312747058112796 ,
    bi6[11][6] =  8.475606457364420387647202650806654444113 ,
    
    bi6[12][1] =  0. ,
    bi6[12][2] = -8.546914399459162196852655523910775929846 ,
    bi6[12][3] =  71.42491845557421929764191001202432467195 ,
    bi6[12][4] = -182.3743002138117003976306073059834582214 ,
    bi6[12][5] =  184.6615026587373916897461066715370461463 ,
    bi6[12][6] = -65.16520650104074839290475385366713666699 ,
    
    bi6[13][1] =  0. ,
    bi6[13][2] =  9.041175750584674634065105394682748876817 ,
    bi6[13][3] = -31.60097240260566461465540792605436016441 ,
    bi6[13][4] =  28.09188900392449596430672577587636184837 ,
    bi6[13][5] =  2.454436197629303379092350647679361289224 ,
    bi6[13][6] = -7.986528549532809362808773892184111850000 ;

    for(G4int i = 0; i< numberOfVariables; i++)
        yIn[i] = yInput[i];
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    
    for(int i=1; i<=13; i++){	//Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = tau0;
        for(int j=1; j<=6; j++){
            b[i] += bi6[i][j]*tau;
            tau*=tau0;
        }
    }
    for(int i=0; i<numberOfVariables; i++){		//Here is IS the cooridnate no.
        yOut[i] = yIn[i] + Step*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i] + b[11]*ak11[i] +
                                 b[12]*ak12[i] + b[13]*ak13[i] );
    }

}

void VernerRK67::SetupInterpolate_high( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step ){
    SetupInterpolate_low(yInput, dydx, Step);
    		//Sets the values of *ak11,*ak12,*ak13
    //Defining the coefficients for the new stages
    G4double
    b141 =  .3887903257436303686399931060834951327899e-1 ,
    b142 =  0. ,
    b143 =  0. ,
    b144 = -.2440320330830131517910045090190069290791e-2 ,
    b145 = -.1392891721467262281273220992320214734208e-2 ,
    b146 = -.4744629155868013465038358934145339168472e-3 ,
    b147 =  .3920793241315951369383517310870803393356e-3 ,
    b148 = -.4055473328512800136385880031750264996936e-3 ,
    b149 =  .1989709314771672628794304728258886009267e-3 ,
    b1410 =  0. ,
    b1411 = -.1027819879317916884712606136811051029682e-3 ,
    b1412 =  .3385661513870266715302548402957613704604e-1 ,
    b1413 =  .1814893063199928004309543737509423302792 ,
    
    //c15 =  53/100 ,
    
    b151 =  .5723681204690012909606837582140921695189e-1 ,
    b152 =  0. ,
    b153 =  0. ,
    b154 =  .2226594806676118099285816235023183680020 ,
    b155 =  .1234486420018689904911221497830317287757 ,
    b156 =  .4006332526666490875113688731927762275433e-1 ,
    b157 = -.5269894848581452066926326838943832327366e-1 ,
    b158 =  .4765971214244522856887315416093212596338e-1 ,
    b159 = -.2138895885042213036387863538386958914368e-1 ,
    b1510 =  0. ,
    b1511 =  .1519389106403640165459624646184297766866e-1 ,
    b1512 =  .1206054671628965554251364472502413614358 ,
    b1513 = -.2277942301618737288237298052574548913451e-1 ,
    b1514 =  0. ,
    
    //c16 =  79/100 ,
    
    b161 =  .5137203880275681426595607279552927584506e-1 ,
    b162 =  0. ,
    b163 =  0. ,
    b164 =  .5414214473439405582401399378307410450482 ,
    b165 =  .3503998066921840081154745647747846804810 ,
    b166 =  .1419311226969218216861835872156617148040 ,
    b167 =  .1052737747842942254816302629823570359198 ,
    b168 = -.3108184780587401700842726199589213259835e-1 ,
    b169 = -.7401883149519145061791854716430279714483e-2 ,
    b1610 =  0. ,
    b1611 = -.6377932504865363437569726480040013149706e-2 ,
    b1612 = -.1732549590836186403386348310205265959935 ,
    b1613 = -.1822815677762202619429607513861847306420 ,
    b1614 =  0. ,
    b1615 =  0. ;
    
    const G4int numberOfVariables= this->GetNumberOfVariables();
    
    //  Saving yInput because yInput and yOut can be aliases for same array
    for(int i=0;i<numberOfVariables;i++)
    {
        yIn[i]=yInput[i];
    }
    
    yTemp[7]  = yIn[7];
    
    
    ak14 = new G4double[numberOfVariables];
    ak15 = new G4double[numberOfVariables];
    ak16 = new G4double[numberOfVariables];
    
    
    
    //    calculating extra stage functions
    
    //Evaluate the stages :
    
    for(G4int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b141*dydx[i] + b142*ak2[i] + b143*ak3[i] +
                                  b144*ak4[i] + b145*ak5[i] + b146*ak6[i] +
                                  b147*ak7[i] + b148*ak8[i] + b149*ak9[i] +
                                  b1410*ak10[i] + b1411*ak11[i] + b1412*ak12[i] +
                                  b1413*ak13[i] );
    }
    RightHandSide(yTemp, ak14);			//14th Stage
    
    for(G4int i=0; i<numberOfVariables; i++)
    {
        yTemp[i] = yIn[i] + Step*(b151*dydx[i] + b152*ak2[i] + b153*ak3[i] +
                                  b154*ak4[i] + b155*ak5[i] + b156*ak6[i] +
                                  b157*ak7[i] + b158*ak8[i] + b159*ak9[i] +
                                  b1510*ak10[i] + b1511*ak11[i] + b1512*ak12[i] +
                                  b1513*ak13[i] + b1514*ak14[i]);
    }
    RightHandSide(yTemp, ak15);			//15th Stage
    
    
    for(G4int i=0;i<numberOfVariables;i++)
    {
        yTemp[i] = yIn[i] + Step*(b161*dydx[i] + b162*ak2[i] + b163*ak3[i] +
                                  b164*ak4[i] + b165*ak5[i] + b166*ak6[i] +
                                  b167*ak7[i] + b168*ak8[i] + b169*ak9[i] +
                                  b1610*ak10[i] + b1611*ak11[i] + b1612*ak12[i] +
                                  b1613*ak13[i] + b1614*ak14[i] +b1615*ak15[i]
                                  );
    }
    RightHandSide(yTemp, ak16);			//16th Stage
    
}



//For calculating the output at the tau fraction of Step
void VernerRK67::Interpolate_high( const G4double yInput[],
                                   const G4double dydx[],
                                   const G4double Step,
                      					 G4double yOut[],
                    					 G4double tau ){
    
    G4double bi7[17][8], b[17];
    G4int numberOfVariables = this->GetNumberOfVariables();

//    Coefficients for the interpolant bi7 with 16 stages
    
    bi7[1][1] =  1. ,
    bi7[1][2] = -8.413387198332767469319987751201351965810 ,
    bi7[1][3] =  33.67550888449089654479469983556967202215 ,
    bi7[1][4] = -70.80159089484886164618905961010838757357 ,
    bi7[1][5] =  80.64695108301297872968868805293298389704 ,
    bi7[1][6] = -47.19413969837521580145883430419406103536 ,
    bi7[1][7] =  11.13381344253924186418881142808952641234 ,
    
    bi7[2][1] =  0. ,
    bi7[2][2] =  0. ,
    bi7[2][3] =  0. ,
    bi7[2][4] =  0. ,
    bi7[2][5] =  0. ,
    bi7[2][6] =  0. ,
    bi7[2][7] =  0. ,
    
    bi7[3][1] =  0. ,
    bi7[3][2] =  0. ,
    bi7[3][3] =  0. ,
    bi7[3][4] =  0. ,
    bi7[3][5] =  0. ,
    bi7[3][6] =  0. ,
    bi7[3][7] =  0. ,
    
    bi7[4][1] =  0. ,
    bi7[4][2] =  8.754921980674397160629587282876763437696 ,
    bi7[4][3] = -88.45968286997709426134300934922618655402 ,
    bi7[4][4] =  346.9017638429916309499891288356321692825 ,
    bi7[4][5] = -629.2580030059837046812187141184986252218 ,
    bi7[4][6] =  529.6773755604192983874116479833480529304 ,
    bi7[4][7] = -167.3588698651401860365089970240284051167 ,
    
    bi7[5][1] =  0. ,
    bi7[5][2] =  8.913387586637921662996190126913331844214 ,
    bi7[5][3] = -90.06081846893217794712014609702916991513 ,
    bi7[5][4] =  353.1807459217057824951538014683541349020 ,
    bi7[5][5] = -640.6476819744374433668701027882567716886 ,
    bi7[5][6] =  539.2646279047155261551781390920363285084 ,
    bi7[5][7] = -170.3880944299154827945664954924414008798 ,
    
    bi7[6][1] =  0. ,
    bi7[6][2] =  5.173312029847800338889849068990984974299 ,
    bi7[6][3] = -52.27111590005538823385270070373176751689 ,
    bi7[6][4] =  204.9853867374073094711024260808085419491 ,
    bi7[6][5] = -371.8306118563602890875634623992262437796 ,
    bi7[6][6] =  312.9880934374529000210073972654145891826 ,
    bi7[6][7] = -98.89290352172494693555119599233959305606 ,
    
    bi7[7][1] =  0. ,
    bi7[7][2] =  16.79537744079695986364946329034055578253 ,
    bi7[7][3] = -169.7004000005972744435739149730966805754 ,
    bi7[7][4] =  665.4937727009246303131700313781960584913 ,
    bi7[7][5] = -1207.163889233600728395392916633015853882 ,
    bi7[7][6] =  1016.129151581854603280159105697386989470 ,
    bi7[7][7] = -321.0600155723749421933210511704882816019 ,
    
    bi7[8][1] =  0. ,
    bi7[8][2] = -10.00599753609866476866352971232058330270 ,
    bi7[8][3] =  101.1005433052275068199636113246449312792 ,
    bi7[8][4] = -396.4739151237843754958939772727577263768 ,
    bi7[8][5] =  719.1787707014182914108130834128646525498 ,
    bi7[8][6] = -605.3681033918824350795711030652978269725 ,
    bi7[8][7] =  191.2743989279793520691961908384572824802 ,
    
    bi7[9][1] =  0. ,
    bi7[9][2] =  2.764708833638599139713222853969606774131 ,
    bi7[9][3] = -27.93460263739046178114640484830267988046 ,
    bi7[9][4] =  109.5477918613789217803046856340175757800 ,
    bi7[9][5] = -198.7128113064482116421691972646370773711 ,
    bi7[9][6] =  167.2663357164031670694252647113936863857 ,
    bi7[9][7] = -52.85010499525706346613022509203974406942 ,
    
    bi7[10][1] =  0. ,
    bi7[10][2] =  0. ,
    bi7[10][3] =  0. ,
    bi7[10][4] =  0. ,
    bi7[10][5] =  0. ,
    bi7[10][6] =  0. ,
    bi7[10][7] =  0. ,
    
    bi7[11][1] =  0. ,
    bi7[11][2] = -2.169632028016350481156919876642428429100 ,
    bi7[11][3] =  22.01669603756987625585768587320929912766 ,
    bi7[11][4] = -86.90152427798948350846176288615482496306 ,
    bi7[11][5] =  159.2238897386147443720253338471077193471 ,
    bi7[11][6] = -135.9618306534587908363115231453760181702 ,
    bi7[11][7] =  43.79240118328000419804718618785625308759 ,
    
    bi7[12][1] =  0. ,
    bi7[12][2] = -4.890070188793803933769786966428026149549 ,
    bi7[12][3] =  22.75407737425176120799532459991506803585 ,
    bi7[12][4] = -30.78034218537730965082079824005797506535 ,
    bi7[12][5] = -2.797194317207249021142015125037024035537 ,
    bi7[12][6] =  31.36945663750840183161406140272783187147 ,
    bi7[12][7] = -15.65592732038180043387678567111987465689 ,
    
    bi7[13][1] =  0. ,
    bi7[13][2] =  10.86217092955196715517224349929627754387 ,
    bi7[13][3] = -50.54297141782710697188187875653305700081 ,
    bi7[13][4] =  68.37148040407511827604242008548181691494 ,
    bi7[13][5] =  6.213326521632409162585500428935637861213 ,
    bi7[13][6] = -69.68006323194158104163196358466588618336 ,
    bi7[13][7] =  34.77605679450919341971367832748521086414 ,
    
    bi7[14][1] =  0. ,
    bi7[14][2] = -11.37286691922922915922346687401389055763 ,
    bi7[14][3] =  130.7905807824671644130452602841032046030 ,
    bi7[14][4] = -488.6511367778560207543260583489312609826 ,
    bi7[14][5] =  832.2148793276440873476229585070779183432 ,
    bi7[14][6] = -664.7743368554426242883314487337054193624 ,
    bi7[14][7] =  201.7928804424166224412127551654694479565 ,
    
    bi7[15][1] =  0. ,
    bi7[15][2] = -5.919778732715006698693070786679427540601 ,
    bi7[15][3] =  63.27679965889218829298274978013773800731 ,
    bi7[15][4] = -265.4326820887379575820873554556433306580 ,
    bi7[15][5] =  520.1009254140610824835871087519714692468 ,
    bi7[15][6] = -467.4121095339020118993777963241667608460 ,
    bi7[15][7] =  155.3868452824017054035883640343803117904 ,
    
    bi7[16][1] =  0. ,
    bi7[16][2] = -10.49214619796182281022379415510181241136 ,
    bi7[16][3] =  105.3553852518801101042787230303396283676 ,
    bi7[16][4] = -409.4397501198893846479834816688367917005 ,
    bi7[16][5] =  732.8314489076540326880337353277812147333 ,
    bi7[16][6] = -606.3044574733512377981129469949015057785 ,
    bi7[16][7] =  188.0495196316683024640077644607192667895 ;
    
    for(G4int i = 0; i< numberOfVariables; i++)
        yIn[i] = yInput[i];
    
    G4double tau0 = tau;
    //    Calculating the polynomials :
    
    for(int i=1; i<=16; i++){	//Here i is NOT the coordinate no. , it's stage no.
        b[i] = 0;
        tau = tau0;
        for(int j=1; j<=7; j++){
            b[i] += bi7[i][j]*tau;
            tau*=tau0;
        }
    }
    
    for(int i=0; i<numberOfVariables; i++){		//Here is IS the cooridnate no.
        yOut[i] = yIn[i] + Step*(b[1]*dydx[i] + b[2]*ak2[i] + b[3]*ak3[i] +
                                 b[4]*ak4[i] + b[5]*ak5[i] + b[6]*ak6[i] +
                                 b[7]*ak7[i] + b[8]*ak8[i] + b[9]*ak9[i] +
                                 b[10]*ak10[i] + b[11]*ak11[i] + b[12]*ak12[i] +
                                 b[13]*ak13[i] + b[14]*ak14[i] + b[15]*ak15[i] +
                                 b[16]*ak16[i] );
    }
    
}

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










