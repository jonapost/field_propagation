// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
#include "MuruaRKN5459.hh"

#include "G4Mag_EqRhs.hh"
#include "MagEqRhsbyTimestoreB.hh"

#include "G4CachedMagneticField.hh"

// Needed for DistChord()
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"



MuruaRKN5459::MuruaRKN5459(G4EquationOfMotion *EqRhs,
                           G4int numberOfVariables,
                           G4bool primary_status)
   : G4MagIntegratorStepper(EqRhs, numberOfVariables),
     fLastStepLength(0.), fAuxStepper(0) {

   primary = primary_status;

   myField_as_storeBfield = dynamic_cast<MagEqRhs_byTime_storeB*>(GetEquationOfMotion()); // So we don't have to keep recasting it over and over again.

   position_interpolant = new Interpolant();

   c2 = 771./3847., c3 = 771./3847., c4 = 3051./6788., c5 = 4331./6516., c6 = 4331./6516., c7 = 10463./11400., c8 = 10463./11400., c9 = 1./1.;

   a21 = 771./3847.,
   a31 = 771./7694., a32 = 771./7694.,
   a41 = -264272222./4845505509., a42 = -9458865980./12714902623., a43 = 17133165419./13729279360.,
   a51 = 1943604853./18116134489., a52 = -2470367896./7636570485., a53 = 1733951147./3918733571., a54 = 4613437932./10523350595.,
   a61 = 369952551./2046485744., a62 = 281630106828./143708239525., a63 = -9868262031./5606899429., a64 = 208606720./5597531799., a65 = 792516107./3209667255.,
   a71 = -2089737154./15083636501., a72 = -39924138556./8175090533., a73 = 72922890855./14010113917., a74 = 9484193351./15493195043., a75 = -17895265139./12412283353., a76 = 278232./177835.,
   a81 = -1762013041./13188190032., a82 = -22636373880./4795132451., a83 = 30527401913./6048941340., a84 = 11564353310./19632283007., a85 = -50677425731./36595197965., a86 = 12408./8167., a87 = 10722067./5782709432.,
   a91 = 8034174097./12261534992., a92 = 72032427203./6782716235., a93 = -90566218637./8185393121., a94 = 18770105843./41171085325., a95 = 28010344030./6199889941., a96 = -21917292279./4540377286., a97 = -236637914115./8183370127., a98 = 71217630373./2409299224.;



   alpha21 = 594441./29598818.,
   alpha31 = 594441./29598818., alpha32 = 0./1.,
   alpha41 = -311625081./28869248936., alpha42 = 128./8219., alpha43 = 1015645524./10554116159.,
   alpha51 = 1852480471./26299626569., alpha52 = -247./14069., alpha53 = 648800762./5897141541., alpha54 = 519849979./8963946221.,
   alpha61 = 1852480471./26299626569., alpha62 = -247./14069., alpha63 = 648800762./5897141541., alpha64 = 519849979./8963946221., alpha65 = 0./1.,
   alpha71 = 229929851./7158517178., alpha72 = 113395809./8665398238., alpha73 = 4865737279./19748497543., alpha74 = 340133672./10137556453., alpha75 = 738./11587., alpha76 = 509108839./15737542787.,
   alpha81 = 229929851./7158517178., alpha82 = 113395809./8665398238., alpha83 = 4865737279./19748497543., alpha84 = 340133672./10137556453., alpha85 = 738./11587., alpha86 = 509108839./15737542787., alpha87 = 0./1.,
   alpha91 = 164505448./2653157365., alpha92 = 0./1., alpha93 = 9357192./40412735., alpha94 = 736403089./7677655029., alpha95 = 960089./17896194., alpha96 = 482653907./11393392643., alpha97 = -47281957./150822000., alpha98 = 6715245221./20471724521.;


   b1 = 164505448./2653157365., b2 = 0./1., b3 = 3042./10505., b4 = 1586146904./9104113535., b5 = 4394./27465., b6 = 2081836558./16479128289., b7 = -50461./13230., b8 = 13928550541./3490062596., b9 = 91464477./8242174145.;

   //beta1 = 164505448./2653157365., beta2 = 0./1., beta3 = 3042./10505., beta4 = 1586146904./9104113535., beta5 = 4394./27465., beta6 = 2081836558./16479128289., beta7 = -50461./13230., beta8 = 13928550541./3490062596., beta9 = 91464477./8242174145.;

   beta1 = 164505448./2653157365., beta2 = 0./1., beta3 = 9357192./40412735., beta4 = 736403089./7677655029., beta5 = 960089./17896194., beta6 = 482653907./11393392643., beta7 = -47281957./150822000., beta8 = 6715245221./20471724521., beta9 = 0./1.;

   berr1 = b1 - 53757362./127184461., berr2 = b2 - 0./1., berr3 = b3 - -138687950./204047369., berr4 = b4 - 161961633./188152853., berr5 = b5 - 36242723./103243418., berr6 = b6 - 1./2., berr7 = b7 - 1147554103./9981952., berr8 = b8 - -2395015001./20532034., berr9 = b9 - 1./1., berr10 = - 23./100.;

   //bhat1 = 53757362./127184461., bhat2 = 0./1., bhat3 = -138687950./204047369., bhat4 = 161961633./188152853., bhat5 = 36242723./103243418., bhat6 = 1./2., bhat7 = 1147554103./9981952., bhat8 = -2395015001./20532034., bhat9 = 1./1., bhat10 = 23./100.;
   //beta_hat1 = 53757362./127184461., beta_hat2 = 0./1., beta_hat3 = -426604134200./784970228543., beta_hat4 = 605250622521./1277181566164., beta_hat5 = 79190349755./672734111688., beta_hat6 = 2185./13032., beta_hat7 = 1075258194511./113794252800., beta_hat8 = -2244129055937./234065187600., beta_hat9 = 0./1., beta_hat10 = 0./1.;

   beta_err1 = beta1 - 53757362./127184461., beta_err2 = beta2 - 0./1., beta_err3 = beta3 - -426604134200./784970228543., beta_err4 = beta4 - 605250622521./1277181566164., beta_err5 = beta5 - 79190349755./672734111688., beta_err6 = beta6 - 2185./13032., beta_err7 = beta7 - 1075258194511./113794252800., beta_err8 = beta8 - -2244129055937./234065187600., beta_err9 = beta9 - 0./1., beta_err10 = - 0./1.;



   fLastInitialVector = new G4double[numberOfVariables] ;
   fLastFinalVector = new G4double[numberOfVariables] ;

   fMidVector = new G4double[numberOfVariables];
   fMidError =  new G4double[numberOfVariables];


   if( primary_status )
  {
    fAuxStepper = new MuruaRKN5459(EqRhs, numberOfVariables, !primary_status);
  }
}


void MuruaRKN5459::SetEquationOfMotion(G4EquationOfMotion* newEquation) {

   if( newEquation != 0 ){

      G4MagIntegratorStepper::SetEquationOfMotion( newEquation );
      if (primary)
         fAuxStepper -> G4MagIntegratorStepper::SetEquationOfMotion( newEquation );
   }
}



void MuruaRKN5459::set_MagEqRhs_storedBfield( G4EquationOfMotion *EqRhs ) {
   myField_as_storeBfield = dynamic_cast<MagEqRhs_byTime_storeB*>(GetEquationOfMotion());
   if ( primary )
      fAuxStepper -> set_MagEqRhs_storedBfield( EqRhs );
}


MuruaRKN5459::~MuruaRKN5459() {
   delete position_interpolant;

   delete fLastInitialVector; delete fLastFinalVector;
   delete fMidVector; delete fMidError;
}

void MuruaRKN5459::ComputeRhsWithStoredB( const G4double y[],
                                          //const G4double B[3],
                                          G4double dydx[] )
{
   G4double cof = myField_as_storeBfield -> FCof();
   cof /= myField_as_storeBfield -> FMass();   // Correction because we are using dynamic variables
                                // position and velocity (not momentum).
                                // (This is just a copy from Mag_UsualEqRhs_IntegrateByTime.)

   B = myField_as_storeBfield -> get_last_B_eval();

   // Ignore 0 through 3 components. We don't use them for a Nystrom stepper.
   /*
   dydx[0] = y[3];
   dydx[1] = y[4];
   dydx[2] = y[5];
   */

   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   return ;
}


G4double  MuruaRKN5459::DistChord()   const {


   //G4int no_function_calls_before_aux_stepper =
   //      (( G4CachedMagneticField* )( mTracker -> getStepper() -> GetEquationOfMotion() -> GetFieldObj() ))
   //                                            -> GetCountCalls();

   // Implementation borrowed from G4CashKarpRK45:

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


     //mTracker -> no_function_calls_used_by_DistChord +=
     //      (( G4CachedMagneticField* )( mTracker -> getStepper() -> GetEquationOfMotion() -> GetFieldObj() ))
     //      -> GetCountCalls() - no_function_calls_before_aux_stepper;

     return distChord;

   /* Ideal implementation, using interpolant. But we don't have the interpolant working yet.
   G4double distLine, distChord;
   // Store last initial and final points (they will be overwritten in self-Stepper call!)

   if (! position_interpolant -> IsInitialized_Position() ) {
      position_interpolant -> Initialize( fLastInitialVector,
                                          fLastFinalVector, fLastDyDx,
                                          fNextDyDx, fLastStepLength );
   }

   position_interpolant -> InterpolatePosition( 0.5, fMidVector );


   // Should probably make these vectors class variables.
   G4ThreeVector midPoint( fMidVector[0], fMidVector[1], fMidVector[2] );
   G4ThreeVector initialPoint( fLastInitialVector[0], fLastInitialVector[1], fLastInitialVector[2] );
   G4ThreeVector finalPoint( fLastFinalVector[0], fLastFinalVector[1], fLastFinalVector[2] );

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

   */

}


void MuruaRKN5459::Stepper(  const G4double y[],
                  const G4double dydx[],
                  G4double h,
                  G4double yout[],
                  G4double yerr[] ) {

   int k;

   for (k = 0; k < 3; k ++) {
      q0[k] = y[k];
      v0[k] = y[k + 3];
      Vdot1[k + 3] = dydx[k + 3];
   }


   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c2*v0[k] + h*h*(alpha21*Vdot1[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a21*Vdot1[k+3]);
   }
   RightHandSide(ytemp, Vdot2);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a31*Vdot1[k+3] + a32*Vdot2[k+3]);
   }
   // Use last B field Evaluation

   ComputeRhsWithStoredB(ytemp, Vdot3);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c4*v0[k] + h*h*(alpha41*Vdot1[k+3]+alpha42*Vdot2[k+3]+alpha43*Vdot3[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a41*Vdot1[k+3] + a42*Vdot2[k+3] + a43*Vdot3[k+3]);
   }
   RightHandSide(ytemp, Vdot4);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c5*v0[k] + h*h*(alpha51*Vdot1[k+3]+alpha52*Vdot2[k+3]+alpha53*Vdot3[k+3]+alpha54*Vdot4[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a51*Vdot1[k+3] + a52*Vdot2[k+3] + a53*Vdot3[k+3] + a54*Vdot4[k+3]);
   }
   RightHandSide(ytemp, Vdot5);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a61*Vdot1[k+3] + a62*Vdot2[k+3] + a63*Vdot3[k+3] + a64*Vdot4[k+3] + a65*Vdot5[k+3]);
   }
   // Use last B field Evaluation

   ComputeRhsWithStoredB(ytemp, Vdot6);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c7*v0[k] + h*h*(alpha71*Vdot1[k+3]+alpha72*Vdot2[k+3]+alpha73*Vdot3[k+3]+alpha74*Vdot4[k+3]+alpha75*Vdot5[k+3]+alpha76*Vdot6[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a71*Vdot1[k+3] + a72*Vdot2[k+3] + a73*Vdot3[k+3] + a74*Vdot4[k+3] + a75*Vdot5[k+3] + a76*Vdot6[k+3]);
   }
   RightHandSide(ytemp, Vdot7);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a81*Vdot1[k+3] + a82*Vdot2[k+3] + a83*Vdot3[k+3] + a84*Vdot4[k+3] + a85*Vdot5[k+3] + a86*Vdot6[k+3] + a87*Vdot7[k+3]);
   }
   // Use last B field Evaluation

   ComputeRhsWithStoredB(ytemp, Vdot8);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c9*v0[k] + h*h*(alpha91*Vdot1[k+3]+alpha92*Vdot2[k+3]+alpha93*Vdot3[k+3]+alpha94*Vdot4[k+3]+alpha95*Vdot5[k+3]+alpha96*Vdot6[k+3]+alpha97*Vdot7[k+3]+alpha98*Vdot8[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a91*Vdot1[k+3] + a92*Vdot2[k+3] + a93*Vdot3[k+3] + a94*Vdot4[k+3] + a95*Vdot5[k+3] + a96*Vdot6[k+3] + a97*Vdot7[k+3] + a98*Vdot8[k+3]);
   }
   RightHandSide(ytemp, Vdot9);

   G4double temp_var;

   // Just copy position values into yout[0..2]:
   for (k = 0; k < 3; k ++) {
      yout[k] = ytemp[k];
      //temp_var = beta1*Vdot1[k] + beta2*Vdot2[k] + beta3*Vdot3[k] + beta4*Vdot4[k] + beta5*Vdot5[k] + beta6*Vdot6[k] + beta7*Vdot7[k] + beta8*Vdot8[k] + beta9*Vdot9[k];
      //yout[k - 3] = q0[k - 3] + h*v0[k - 3] + h*h*temp_var;
   }
   for (k = 3; k < 6; k ++) {
      temp_var = b1*Vdot1[k] + b2*Vdot2[k] + b3*Vdot3[k] + b4*Vdot4[k] + b5*Vdot5[k] + b6*Vdot6[k] + b7*Vdot7[k] + b8*Vdot8[k] + b9*Vdot9[k];
      yout[k] = v0[k - 3] + h*temp_var;
   }

   // Error:
   ComputeRhsWithStoredB(yout, Vdot10);

   for (k = 3; k < 6; k ++) {
         temp_var = beta_err1*Vdot1[k] + beta_err2*Vdot2[k] + beta_err3*Vdot3[k] + beta_err4*Vdot4[k] + beta_err5*Vdot5[k] + beta_err6*Vdot6[k] + beta_err7*Vdot7[k] + beta_err8*Vdot8[k] + beta_err9*Vdot9[k] + beta_err10*Vdot10[k];
         yerr[k - 3] = h*h*temp_var;
      }
   for (k = 3; k < 6; k ++) {
      temp_var = berr1*Vdot1[k] + berr2*Vdot2[k] + berr3*Vdot3[k] + berr4*Vdot4[k] + berr5*Vdot5[k] + berr6*Vdot6[k] + berr7*Vdot7[k] + berr8*Vdot8[k] + berr9*Vdot9[k] + berr10*Vdot10[k];
      yerr[k] = h*temp_var;
   }

   //last_step_len = h;
   for(int i = 0; i < 6; i ++ )
   {
      fLastInitialVector[i] = y[i] ;
      fLastFinalVector[i]   = yout[i];
   }
   for (int i = 0; i < 3; i ++ ) {
      fLastDyDx[i + 3]  = dydx[i + 3];        // Switch back to i from "i + 3" when going back to using interpolation.
      fNextDyDx[i + 3]  = Vdot10[i + 3];      // Same comment as above.
   }
   // NormaliseTangentVector( yOut ); // Not wanted

   fLastStepLength = h;



}
