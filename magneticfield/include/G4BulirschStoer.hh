// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: $
//
// Helper namespace 'magneticfield', 'internal'
//
// Description:
// G4BulirschStoer class implementation by Dmitry Sorokin
// Implementation is based on bulirsch_stoer.hpp from boost
//
// The Bulirsch-Stoer is a controlled driver that adjusts both step size
// and order of the method. The algorithm uses the modified midpoint and
// a polynomial extrapolation compute the solution.
//
//    Implementation by Dmitry Sorokin - GSoC 2016
//       Work supported by Google as part of Google Summer of Code 2016.
//    Supervision / code review: John Apostolakis
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4BULIRSCH_STOER_HH
#define G4BULIRSCH_STOER_HH

#include "G4ModifiedMidpoint.hh"

#include "G4FieldTrack.hh"

enum class step_result{
    success,
    fail
};

namespace internal {

class G4BulirschStoer {
public:
    const static G4int m_k_max = 8;

    G4BulirschStoer(
        G4EquationOfMotion* equation,
        G4int nvar,
        G4double eps_rel,
        G4double max_dt);

    inline void set_max_dt(G4double max_dt);
    inline void set_max_relative_error(G4double eps_rel);

    //stepper method
    step_result try_step(const G4double in[],
        const G4double dxdt[],
        G4double& t,
        G4double out[],
        G4double& dt);

    //Reset the internal state of the stepper
    void reset();


private:
    void extrapolate(size_t k, G4double xest[]);
    G4double calc_h_opt(G4double h, G4double error, size_t k) const;

    G4bool set_k_opt(size_t k, G4double& dt);
    G4bool in_convergence_window(G4int k) const;
    G4bool should_reject(G4double error, G4int k) const;

private:
    //pointer to the system class
    // G4EquationOfMotion* fequation;

    //number of vars to be integrated
    G4int fnvar;

    //relative tolerance
    G4double m_eps_rel;

    //modified midpoint algorithm
    G4ModifiedMidpoint m_midpoint;

    G4bool m_last_step_rejected;
    G4bool m_first;

    G4double m_dt_last;
    // G4double m_t_last;

    //max allowed time step
    G4double m_max_dt;

    G4int m_current_k_opt;

    // G4double m_xnew[G4FieldTrack::ncompSVEC];
    G4double m_err[G4FieldTrack::ncompSVEC];
    // G4double m_dxdt[G4FieldTrack::ncompSVEC];

    // stores the successive interval counts
    G4int m_interval_sequence[m_k_max+1];

    //Extrapolation coeffs (Neville’s algorithm)
    G4double m_coeff[m_k_max+1][m_k_max];

    // costs for interval count
    G4int m_cost[m_k_max+1];

    // sequence of states for extrapolation
    G4double m_table[m_k_max][G4FieldTrack::ncompSVEC];

    //oprimal step size
    G4double h_opt[m_k_max+1];

    //Work per unit step
    G4double work[m_k_max+1];
};

} // internal

#include "G4BulirschStoer.icc"

#endif

