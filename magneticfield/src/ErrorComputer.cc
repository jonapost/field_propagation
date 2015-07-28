/*
 * ErrorComputer.cc
 *
 *  Created on: Jul 5, 2015
 *      Author: jason
 */

#include "ErrorComputer.hh"
#include "Interpolant.hh"

#define TIME_SLOT 0
#define ARCLENGTH_SLOT 1

#define POSITION_SLOT 2
#define VELOCITY_SLOT 5
#define RHS_SLOT 8
#define BUFFER_COLUMN_LEN 11

#define NUMBER_INTERPOLATION_VARIABLES 6

ErrorComputer::ErrorComputer(G4double **in_bufferA,
                             G4int lenBufferA,
                             G4double **in_bufferB,
                             G4int lenBufferB) {
   counter = 1;
   minterpolant = new Interpolant();

   bufferA = in_bufferA;
   bufferB = in_bufferB;

   len_bufferA = lenBufferA;
   len_bufferB = lenBufferB;
}

ErrorComputer::~ErrorComputer() {
   delete minterpolant;
}


G4bool ErrorComputer::PathA_Interpolant(G4int pathB_index,
                                        G4double *pathA_interpolant) {

   G4double h;
   G4double t = bufferB[pathB_index][TIME_SLOT];
   while (t > bufferA[counter][TIME_SLOT]) {

      counter ++;
      minterpolant -> DeInitialize();

      if ( counter >= len_bufferA ) // We have passed the possible
                                    // interpolation zone
         return false;
   }

   if (t == bufferA[counter][TIME_SLOT]) {
      for (int i = POSITION_SLOT; i < POSITION_SLOT + 3; i ++) {
         pathA_interpolant[i] = bufferA[counter][i];
      }
      for (int i = VELOCITY_SLOT; i < VELOCITY_SLOT + 3; i ++) {
         pathA_interpolant[i] = bufferA[counter][i];
      }
      counter ++;
   }
   else {
      if ( ! minterpolant -> IsInitialized_Position() ) {
         for (int i = 0; i < 3; i ++) {
            F0[i] = bufferA[counter - 1][i + RHS_SLOT];
            F1[i] = bufferA[counter][i + RHS_SLOT];
         }
         h = bufferA[counter][TIME_SLOT] - bufferA[counter - 1][TIME_SLOT];

         minterpolant -> Initialize( &( bufferA[counter - 1][POSITION_SLOT] ),
                                     &( bufferA[counter][POSITION_SLOT] ),
                                     F0, F1, h );
      }
      G4double xi = ( t - bufferA[counter - 1][TIME_SLOT] ) / h;
      minterpolant -> InterpolatePosition( xi, pathA_interpolant );
      minterpolant -> InterpolateVelocity( xi, &(pathA_interpolant[3]) );
      // InterpolateVelocity() will initialize the velocity interpolating polynomials
      // if they are not initialized.
   }
   return true;
}

G4int ErrorComputer::ErrorArray(G4double **err) {

   G4double pathA_interpolant[NUMBER_INTERPOLATION_VARIABLES]; // Interpolate has only been written
                                  // to handle position (not velocity yet)

   for ( int pathB_index = 0; pathB_index < len_bufferB; pathB_index ++ ) {
      if ( PathA_Interpolant(pathB_index, pathA_interpolant) == false )
         return pathB_index; // We were only able to interpolate up to pathB_index
                             // number of rows in the B array.

      err[pathB_index][0] = bufferB[pathB_index][0]; // Record time.
      err[pathB_index][1] = bufferB[pathB_index][1]; // Record arclength.

      for (int i = 0; i < NUMBER_INTERPOLATION_VARIABLES; i ++) { // 3 because we are only computing error for position.

         err[pathB_index][i + POSITION_SLOT] =
               bufferB[pathB_index][i + POSITION_SLOT] - pathA_interpolant[i];
      }
   }
   // Were able to interpolate error for all values in B array
   return len_bufferB;
}
