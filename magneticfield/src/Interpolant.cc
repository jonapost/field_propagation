/*
 * Interpolant.cc
 *
 *  Created on: Jun 22, 2015
 *      Author: jason
 */

#include "Interpolant.hh"

#include <iostream>
using namespace std;


Interpolant::Interpolant(const G4double y0in[],
                         const G4double y1in[],
                         const G4double F0[],
                         const G4double F1[],
                         G4double step) {
   Initialize(y0in, y1in, F0, F1, step);
}

void Interpolant::Initialize( const G4double y0in[],
                              const G4double y1in[],
                              const G4double F0[],
                              const G4double F1[],
                              G4double step) {

   for (int k = 0; k < 3; k ++) {
            y1[k] = y0in[k];
            y1prime[k] = y0in[k + 3];
            y1prime2[k] = F0[k]; // was F0[k + 3]
            y2[k] = y1in[k];
            y2prime[k] = y1in[k + 3];
            y2prime2[k] = F1[k]; // was F1[k + 3]
      }
      h = step;
      construct_position_polynomials();
      is_initialized = true;
}

void Interpolant::construct_position_polynomials() {

   G4double h2 = h*h;

   for (int k = 0; k < 3; k ++) {
      p1[k] = h * y1prime[k];
      p2[k] = (h2 / 2.) * y1prime2[k];
      p3[k] = -10.*y1[k] - 6*h*y1prime[k] - 3*(h2/2.)*y1prime2[k] + (h2/2.)*y2prime2[k] - 4.*h*y2prime[k] + 10.*y2[k];
      p4[k] = 15.*y1[k] + 8.*h*y1prime[k] + 3.*(h2/2.)*y1prime2[k] - 2.*(h2/2.)*y2prime2[k] + 7.*h*y2prime[k] - 15.*y2[k];
      p5[k] = -6.*y1[k] - 3.*h*y1prime[k] - (h2/2.)*y1prime2[k] + (h2/2.)*y2prime2[k] - 3.*h*y2prime[k] + 6.*y2[k];
   }
}


void Interpolant::InterpolatePosition(G4double xi, G4double yout[]) {



   G4double xi1 = xi, xi2 = xi*xi;
   G4double xi3 = xi*xi2, xi4 = xi2*xi2;
   G4double xi5 = xi3*xi2;
   for (int k = 0; k < 3; k ++) {
      yout[k] = y1[k] + p1[k]*xi1 + p2[k]*xi2 + p3[k]*xi3 + p4[k]*xi4 + p5[k]*xi5;
   }
}

void Interpolant::InterpolateMomentum(G4double xi, G4double yout[]) {
   // Dummy implementation
   yout[0] = 0.; yout[1] = 0.; yout[2] = 0.;
}



