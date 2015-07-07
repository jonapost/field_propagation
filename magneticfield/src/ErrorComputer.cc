/*
 * ErrorComputer.cc
 *
 *  Created on: Jul 5, 2015
 *      Author: jason
 */

#include "ErrorComputer.hh"
#include "Interpolant.hh"
#include "assert.h"

#include <iostream>
using namespace std;


#define TIME_SLOT 9

ErrorComputer::ErrorComputer(G4double **in_bufferA, G4int lenBufferA, G4double **in_bufferB,
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


void ErrorComputer::PathA_Interpolant(G4int pathB_index,  G4double *pathA_interpolant) {

   G4double h;
   G4double t = bufferB[pathB_index][TIME_SLOT];
   while (t > bufferA[counter][TIME_SLOT]) {
      /*
      for (int i = 0; i < 3; i ++) {
         cout << bufferA[counter][i] << ",";
      }
      cout << bufferA[counter][TIME_SLOT] << " (Original) " << endl;
      */

      counter ++;
      minterpolant -> DeInitialize();
      assert( counter < len_bufferA );
   }
   if (t == bufferA[counter][TIME_SLOT]) {
      for (int i = 0; i < 3; i ++) {
         pathA_interpolant[i] = bufferA[counter][i];
      }
      counter ++;
   }
   else {
      if ( ! minterpolant -> IsInitialized() ) {
         for (int i = 0; i < 3; i ++) {
            F0[i] = bufferA[counter - 1][i + 6];
            F1[i] = bufferA[counter][i + 6];
         }
         h = bufferA[counter][TIME_SLOT] - bufferA[counter - 1][TIME_SLOT];
         minterpolant -> Initialize(bufferA[counter - 1], bufferA[counter], F0, F1, h );
      }
      G4double xi = ( t - bufferA[counter - 1][TIME_SLOT] ) / h;
      //cout << "xi,h " << xi << ", " << h << endl;
      minterpolant -> InterpolatePosition( xi, pathA_interpolant );
   }
}

void ErrorComputer::ErrorArray(G4double **err) {

   G4double pathA_interpolant[3]; // Interpolate has only been written to handle position (not velocity yet)

   for ( int pathB_index = 0; pathB_index < len_bufferB; pathB_index ++ ) {
      PathA_Interpolant(pathB_index, pathA_interpolant);

      /*
      for (int i = 0; i < 3; i ++) {
         cout << pathA_interpolant[i] << ",";
      }
      cout << bufferB[pathB_index][TIME_SLOT] << endl;
      */

      //cout << bufferA[counter - 1][TIME_SLOT] << " < " << bufferB[pathB_index][TIME_SLOT] << " < " << bufferA[counter][TIME_SLOT] << endl;

      for (int i = 0; i < 3; i ++) { // 3 because only for position yet.

         err[pathB_index][i] = bufferB[pathB_index][i] - pathA_interpolant[i];
      }
   }



}
