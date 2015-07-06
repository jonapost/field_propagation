/*
 * ErrorComputer.cc
 *
 *  Created on: Jul 5, 2015
 *      Author: jason
 */

#include <ErrorComputer.hh>
#include "Interpolant.hh"
#include "assert.h"

#define TIME_SLOT 9

ErrorComputer::ErrorComputer(G4double **bufferA, G4int lenBufferA, G4double **bufferB,
                              G4int lenBufferB) {
   counter = 1;
   minterpolant = new Interpolant();
   len_bufferA = lenBufferA;
   len_bufferB = lenBufferB;
}

ErrorComputer::~ErrorComputer() {
   delete minterpolant;
}


void ErrorComputer::PathA_Interpolant(G4int pathB_index,  G4double *pathA_interpolant) {

   G4double t = bufferB[pathB_index][TIME_SLOT];
   while (t < bufferA[counter][TIME_SLOT]) {
      counter ++;
      minterpolant -> DeInitialize();
      assert( counter < len_bufferA );
   }
   if (t == bufferA[counter][TIME_SLOT]) {
      for (int i = 0; i < 6; i ++) {
         pathA_interpolant[i] = bufferA[counter][i];
      }
   }
   else {
      if ( ! minterpolant -> IsInitialized() ) {
         for (int i = 0; i < 3; i ++) {
            F0[i] = bufferA[counter - 1][i + 6];
            F1[i] = bufferA[counter][i + 6];
         }
         G4double h = bufferA[counter][9] - bufferA[counter - 1][TIME_SLOT];
         minterpolant -> Initialize(bufferA[counter - 1], bufferA[counter], F0, F1, h );
      }
      G4double xi = bufferB[pathB_index][TIME_SLOT] - bufferA[counter - 1][TIME_SLOT];
      minterpolant -> InterpolatePosition( xi, pathA_interpolant );
   }
}

void ErrorComputer::ErrorArray(G4double **err) {

   G4double pathA_interpolant[3]; // Interpolate has only been written to handle position (not velocity yet)

   for ( int pathB_index = 0; pathB_index < len_bufferB; pathB_index ++ ) {
      PathA_Interpolant(pathB_index, pathA_interpolant);
      for (int i = 0; i < 3; i ++) { // 3 because only for position yet.

         err[pathB_index][i] = bufferB[pathB_index][i] - pathA_interpolant[i];
      }
   }



}
