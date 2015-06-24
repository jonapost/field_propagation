/*
 * FineRKNG34.cc
 *
 *  Created on: Jun 21, 2015
 *      Author: Jason Suagee
 */

#include "FineRKNG34.hh"
#include <iostream>
using namespace std;


FineRKNG34::~FineRKNG34() {
   for (int i = 0; i < 5; i ++) {
      delete[] f[i];
   }
   delete[] f;
   for (int i = 0; i < 5; i ++) {
      delete[] a[i];
      delete[] aprime[i];
   }
   delete[] a;
   delete[] aprime;
   delete position_interpolant;
}

FineRKNG34::FineRKNG34(G4EquationOfMotion *EqRhs,
      G4int numberOfVariables,
      G4int numberOfStateVariables)
: G4MagIntegratorStepper(EqRhs, numberOfVariables){

   position_interpolant = new Interpolant();

   f = new G4double*[5];
   for (int i = 0; i < 5; i ++) {
      f[i] = new G4double[6];
   }
   a = new G4double*[5];
   aprime = new G4double*[5];
   for (int i = 0; i < 5; i ++) {
      a[i] = new G4double[4];
      aprime[i] = new G4double[4];
   }
   for (int i = 0; i < 5; i ++) {
      for (int j = 0; j < 5; j ++) {
         a[i][j] = 0.;
         aprime[i][j] = 0.;
      }
   }

   a[1][0] = 2. / 81.;
   a[2][0] = 1. / 36.; a[2][1] = 1. / 36.;
   a[3][0] = 9. / 128; a[3][1] = 0. / 1.; a[3][2] = 27. / 128.;
   a[4][0] = 11. / 60.; a[4][1] = -3. / 20.; a[4][2] = 9. / 25.; a[4][3] = 8. / 75.;

   aprime[1][0] = 2. / 9.;
   aprime[2][0] = 1. / 12.; aprime[2][1] = 1. / 4.;
   aprime[3][0] = 69. / 128.; aprime[3][1] = -243. / 128; aprime[3][2] = 135. / 64.;
   aprime[4][0] = -17. / 12.; aprime[4][1] = 27. / 4.; aprime[4][2] = -27. / 5.; aprime[4][3] = 16. / 15.;

   for (int i = 1; i < 5; i ++) {
      c[i] = 0.;
      for (int j = 0; j < 5; j ++) {
         c[i] += aprime[i][j];
      }
   }

   b[0] = 19. / 180.; b[1] = 0. / 1.; b[2] = 63. / 200.; b[3] = 16. / 225.; b[4] = 1. / 120.;
   bprime[0] = 1. / 9.; bprime[1] = 0. / 1.; bprime[2] = 9. / 20.; bprime[3] = 16. / 45.; bprime[4] = 1. / 12.;

   b_error[0] = 25. / 1116.; b_error[1] = 0. / 1.; b_error[2] = -63. / 1240.; b_error[3] = 64. / 1395.; b_error[4] = -13. / 744.;
   bprime_error[0] = 2. / 125.; bprime_error[1] = 0. / 1.; bprime_error[2] = -27. / 625.; bprime_error[3] = 32. / 625.; bprime_error[4] = -3. / 125.;
}


void FineRKNG34::Stepper( const G4double y[],
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

   for (i = 1; i < 5; i++) {
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
      ComputeRightHandSide(ytemp, f[i]);
   }

   // yout:
   for (k = 0; k < 3; k ++) {
      yout[k] = 0.;
      for (i = 0; i < 5; i ++) {
         yout[k] += b[i]*f[i][k + 3];
      }
      yout[k] *= h*h;
      yout[k] += y[k] + h*y[k + 3];
   }

   for (k = 3; k < 6; k ++) {
      yout[k] = 0.;
      for (i = 0; i < 5; i ++) {
         yout[k] += bprime[i]*f[i][k];
      }
      yout[k] *= h;
      yout[k] += y[k]; // y'_n
   }

   // yerr:
   for (k = 0; k < 3; k ++) {
      yerr[k] = 0.;
      for (i = 0; i < 5; i ++) {
         yerr[k] += b_error[i]*f[i][k + 3];
      }
      yerr[k] *= h*h;
   }
   for (k = 3; k < 6; k ++) {
      yerr[k] = 0.;
      for (i = 0; i < 5; i ++) {
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

   last_step_len = h;

}

void FineRKNG34::InterpolatePosition(G4double xi, G4double yout[]) {
   if (! position_interpolant->IsInitialized() ) {
      position_interpolant->Initialize(yInitial, yNext, fInitial, fNext, last_step_len);
   }
   position_interpolant->InterpolatePosition(xi, yout);
}


G4double  FineRKNG34::DistChord()   const {
   // temporary function stub
   return 1.;
}

