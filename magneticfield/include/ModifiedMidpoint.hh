// modified midpoint method implementation by Dmitry Sorokin
// Implementation is based on modified_midpoint.hpp from boost
// Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2016
//
//
// This code is made available subject to the Geant4 license, a copy of
//  which is available at http://www.geant4.org/geant4/license/
//
//  History
// -----------------------------
//  Created by Dmitry Sorokin 2016
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ModifiedMidpoint_HH
#define ModifiedMidpoint_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4EquationOfMotion.hh"

class ModifiedMidpoint
{
public :
    ModifiedMidpoint(G4EquationOfMotion* Equation,
                     G4int nvar = 6,
                     G4int steps = 2);

    ~ModifiedMidpoint();

    void do_step(const G4double xIn[], const G4double dxdtIn[] ,G4double xOut[], G4double dt);

    inline void set_steps(G4int steps);

    inline G4int steps() const;

private:

    G4EquationOfMotion* Equation_Rhs;
    G4int fnvar;
    G4int fsteps;

    G4double x0[G4FieldTrack::ncompSVEC];
    G4double x1[G4FieldTrack::ncompSVEC];
    G4double dxdt[G4FieldTrack::ncompSVEC];
    G4double tmp[G4FieldTrack::ncompSVEC];

};


#include "ModifiedMidpoint.icc"

#endif
