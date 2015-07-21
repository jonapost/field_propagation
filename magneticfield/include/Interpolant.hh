/*
 * Interpolant.hh
 *
 *  Created on: Jun 22, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_INTERPOLANT_HH_
#define MAGNETICFIELD_INCLUDE_INTERPOLANT_HH_


#include "G4Types.hh"


class Interpolant {
public:
   Interpolant() { is_initialized = false; } // For beginning of run.

   Interpolant(const G4double y0in[],
               const G4double y1in[],
               const G4double F0[],
               const G4double F1[],
               G4double h);


   virtual ~Interpolant() {};

   void Initialize(  const G4double y0in[],
                     const G4double y1in[],
                     const G4double F0[],
                     const G4double F1[],
                     G4double step);

   inline void DeInitialize();

   void InterpolatePosition(G4double xi, G4double yout[]);
   void InterpolateMomentum(G4double xi, G4double yout[]);

   inline bool IsInitialized();

private:

   void construct_position_polynomials();

   G4double h;
   G4double y1[3], y2[3], y1prime[3], y2prime[3], y1prime2[3], y2prime2[3];
   G4double p1[3], p2[3], p3[3], p4[3], p5[3];

   bool is_initialized;

};

inline
bool Interpolant::IsInitialized() { return is_initialized; }

inline void Interpolant::DeInitialize() { is_initialized = false; }

//inline Interpolant::Interpolant() { is_initialized = false; }

#endif /* MAGNETICFIELD_INCLUDE_INTERPOLANT_HH_ */
