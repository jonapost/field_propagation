//
// classes ModifiedMidpoint and ModifiedMidpointDenseOut
//
// Class description:
//      Modified midpoint which stores derivatives and state at dt/2
//      in some external storage for later usage in dense output calculation
//      This Stepper is for use in Bulirsch Stoer only.
//
//      This is a copy of modified_midpoint.hpp
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------


#ifndef ModifiedMidpointDenseOut_HH
#define ModifiedMidpointDenseOut_HH



#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4EquationOfMotion.hh"

class ModifiedMidpointDenseOut
{
public :
    ModifiedMidpointDenseOut(G4EquationOfMotion* Equation,
                             unsigned int nvar = 6,
                             unsigned int steps = 2);

    ~ModifiedMidpointDenseOut();

    /*
     * performs a modified midpoint step with fsteps intermediate points
     * stores approximation for x(t+dt/2) in x_mp and all evaluated function results in derivs
     */
    void do_step(const G4double xIn[], const G4double dxdtIn[],
                 G4double xOut[], G4double dt,
                 G4double x_mp[], G4double* derivs[]);


    void set_steps(G4int steps);


    G4int steps() const;

private:

    G4EquationOfMotion* Equation_Rhs;
    G4int fnvar;
    G4int fsteps;

    G4double x0[G4FieldTrack::ncompSVEC];
    G4double x1[G4FieldTrack::ncompSVEC];
    G4double dxdt[G4FieldTrack::ncompSVEC];
    G4double tmp[G4FieldTrack::ncompSVEC];

};

#include "ModifiedMidpointDenseOut.icc"

#endif

