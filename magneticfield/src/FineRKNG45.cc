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

#include "FineRKNG45.hh"

#include "G4ThreeVector.hh"
#include "G4LineSection.hh"


//#define JTS_DEBUG

#ifdef JTS_DEBUG
#include <iostream>
using namespace std;
#endif


FineRKNG45::~FineRKNG45() {

   for (int i = 0; i < 7; i ++) {
      delete[] f[i];
   }
   delete[] f;
   for (int i = 0; i < 7; i ++) {
      delete[] a[i];
      delete[] aprime[i];
   }
   delete[] a;
   delete[] aprime;
   delete position_interpolant;

   delete fLastInitialVector; delete fLastFinalVector; delete fLastDyDx;
   delete fNextDyDx; delete fMidVector; delete fMidError;

}


FineRKNG45::FineRKNG45(G4EquationOfMotion *EqRhs,
      G4int numberOfVariables)
: G4MagIntegratorStepper(EqRhs, numberOfVariables){

   position_interpolant = new Interpolant();

   f = new G4double*[7];
   for (int i = 0; i < 7; i ++) {
      f[i] = new G4double[6]; // 6 components for RHS evaluation (we only care about the last 3 components).
   }
   a = new G4double*[7];
   aprime = new G4double*[7];
   for (int i = 0; i < 7; i ++) {
      a[i] = new G4double[6];
      aprime[i] = new G4double[6];
   }
   for (int i = 0; i < 7; i ++) {
      for (int j = 0; j < 6; j ++) {
         a[i][j] = 0.;
         aprime[i][j] = 0.;
      }
   }

   a[1][0] = 32. / 1521.;
   a[2][0] = 4. / 169.; a[2][1] = 4. / 169.;
   a[3][0] = 175. / 5184.; a[3][1] = 0. / 1.; a[3][2] = 1625. / 5184.;
   a[4][0] = -342497279. / 5618900760.; a[4][1] = 6827067. / 46824173.; a[4][2] = 35048741. / 102161832.; a[4][3] = -2201514. / 234120865.;
   a[5][0] = -7079. / 52152.; a[5][1] = 767. / 2173.; a[5][2] = 14027. / 52152.; a[5][3] = 30. / 2173.; a[5][4] = 0. / 1.;
   a[6][0] = 4817. / 51600.; a[6][1] = 0. / 1.; a[6][2] = 388869. / 1216880.; a[6][3] = 3276. / 23575.; a[6][4] = -1142053. / 22015140.; a[6][5] = 0. / 1.;


   aprime[1][0] = 8. / 39.;
   aprime[2][0] = 1. / 13.; aprime[2][1] = 3. / 13.;
   aprime[3][0] = 7385. / 6912.; aprime[3][1] = -9425. / 2304.; aprime[3][2] = 13325. / 3456.;
   aprime[4][0] = 223324757. / 91364240.; aprime[4][1] = -174255393. / 18272848.; aprime[4][2] = 382840094. / 46824173.; aprime[4][3] = -39627252. / 234120865.;
   aprime[5][0] = 108475. / 36464.; aprime[5][1] = -9633. / 848.; aprime[5][2] = 7624604. / 806183.; aprime[5][3] = 8100. / 49979.; aprime[5][4] = -4568212. / 19446707.;
   aprime[6][0] = 4817. / 51600.; aprime[6][1] = 0. / 1.; aprime[6][2] = 1685099. / 3650640.; aprime[6][3] = 19656. / 23575.; aprime[6][4] = -53676491. / 88060560.; aprime[6][5] = 53. / 240.;

   for (int i = 1; i < 7; i ++) {
      c[i] = 0.;
      for (int j = 0; j < 6; j ++) {
         c[i] += aprime[i][j];
      }
   }

   b[0] = 4817. / 51600.; b[1] = 0. / 1.; b[2] = 388869. / 1216880.; b[3] = 3276. / 23575.; b[4] = -1142053. / 22015140.; b[5] = 0. / 1.; b[6] = 0. / 1.;
   bprime[0] = 4817. / 51600.; bprime[1] = 0. / 1.; bprime[2] = 1685099. / 3650640.; bprime[3] = 19656. / 23575.; bprime[4] = -53676491. / 88060560.; bprime[5] = 53. / 240.; bprime[6] = 0. / 1.;

   b_error[0] = 8151. / 2633750.; b_error[1] = 0. / 1.; b_error[2] = -1377519. / 186334750.; b_error[3] = 586872. / 28879375.; b_error[4] = -36011118. / 2247378875.; b_error[5] = 0. / 1.; b_error[6] = 0. / 1.;
   bprime_error[0] = 8151. / 2633750.; bprime_error[1] = 0. / 1.; bprime_error[2] = -5969249. / 559004250.; bprime_error[3] = 3521232. / 28879375.; bprime_error[4] = -846261273. / 4494757750.; bprime_error[5] = 4187. / 36750.; bprime_error[6] = -1. / 25.;


   fLastInitialVector = new G4double[numberOfVariables] ;
   fLastFinalVector = new G4double[numberOfVariables] ;

   fMidVector = new G4double[numberOfVariables];
   fMidError =  new G4double[numberOfVariables];

}

////////////////////////////////////////////////////////////////////////////////////////////////

void FineRKNG45::Stepper( const G4double y[],
                    const G4double dydx[],
                    G4double h,
                    G4double yout[],
                    G4double yerr[] ) {

   // const G4int numberOfVariables= this->GetNumberOfVariables();

   position_interpolant->DeInitialize(); // DeInitialize since we are on a new interval.
                                         // Will have to take care of what to do when there
                                         // is a rejected step later.

   G4double ytemp[10]; // temp value of numberOfVariables

   G4int i, j, k;

   for (k = 0; k < 3; k ++) {
      f[0][k + 3] = fInitial[k + 3] = dydx[k + 3];
   }

   for (i = 1; i < 7; i++) {
      for (k = 0; k < 3; k ++) {

         ytemp[k] = 0.;
         for (j = 0; j < i; j ++) {
            ytemp[k] += a[i][j] * f[j][k + 3];
         }
         ytemp[k] *= h*h;
         ytemp[k] += y[k] + c[i]*h*y[k + 3];

         ytemp[k + 3] = 0.;
         for (j = 0; j < i; j ++) {
            ytemp[k + 3] += aprime[i][j] * f[j][k + 3];
         }
         ytemp[k + 3] *= h;
         ytemp[k + 3] +=  y[k + 3];
      }
      RightHandSide(ytemp, f[i]);

   }

   // yout:
   for (k = 0; k < 3; k ++) {
      yout[k] = 0.;
      for (i = 0; i < 7; i ++) {
         yout[k] += b[i]*f[i][k + 3];
      }
      yout[k] *= h*h;
      yout[k] += y[k] + h*y[k + 3];
   }

   for (k = 3; k < 6; k ++) {
      yout[k] = 0.;
      for (i = 0; i < 7; i ++) {
         yout[k] += bprime[i]*f[i][k];
      }
      yout[k] *= h;
      yout[k] += y[k]; // y'_n
   }

   // yerr:
   for (k = 0; k < 3; k ++) {
      yerr[k] = 0.;
      for (i = 0; i < 7; i ++) {
         yerr[k] += b_error[i]*f[i][k + 3];
      }
      yerr[k] *= h*h;
   }
   for (k = 3; k < 6; k ++) {
      yerr[k] = 0.;
      for (i = 0; i < 7; i ++) {
         yerr[k] += bprime_error[i]*f[i][k];
      }
      yerr[k] *= h;
   }

   for (k = 0; k < 6; k ++) {
      yInitial[k] = y[k];
      yNext[k] = yout[k];
   }
   for (k = 0; k < 3; k ++) {
      fNext[k] = f[4][k];
   }

   //last_step_len = h;
   for( i = 0; i < 6; i ++ )
   {
      fLastInitialVector[i] = y[i] ;
      fLastFinalVector[i]   = yout[i];
   }
   for ( i = 0; i < 3; i ++ ) {
      fLastDyDx[i]          = dydx[i + 3];
      fNextDyDx[i]          = f[6][i + 3];
   }
   // NormaliseTangentVector( yOut ); // Not wanted

   fLastStepLength = h;

}



G4double  FineRKNG45::DistChord()   const {
   G4double distLine, distChord;
   // Store last initial and final points (they will be overwritten in self-Stepper call!)

#ifdef JTS_DEBUG
   cout << "DistChord" << endl;
#endif

   if (! position_interpolant -> IsInitialized() ) {
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

}


