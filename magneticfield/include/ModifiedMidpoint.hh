//
// classes ModifiedMidpoint and ModifiedMidpointDenseOut
//
// Class description:
//
// This is a copy of modified_midpoint.hpp
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

#ifndef MODIFIEDMIDPOINT_HH
#define MODIFIEDMIDPOINT_HH

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
