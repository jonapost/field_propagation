/*
 * ErrorComputer.hh
 *
 *  Created on: Jul 5, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_ERRORCOMPUTER_HH_
#define MAGNETICFIELD_INCLUDE_ERRORCOMPUTER_HH_

#include "Interpolant.hh"

class ErrorComputer {
public:
   ErrorComputer(G4double **in_bufferA, G4int lenBufferA, G4double **in_bufferB,
                             G4int lenBufferB);
   virtual ~ErrorComputer();


   G4bool PathA_Interpolant(G4int pathB_index,  G4double pathA_interpolant[]);

   //void ErrorFromInterpolant(G4double t, G4double err[]);

   G4int ErrorArray(G4double **err);

private:
   Interpolant *minterpolant;
   G4double **bufferA, **bufferB;
   G4int counter;
   G4double F0[3], F1[3];

   G4int len_bufferB, len_bufferA;

};

#endif /* MAGNETICFIELD_INCLUDE_ERRORCOMPUTER_HH_ */
