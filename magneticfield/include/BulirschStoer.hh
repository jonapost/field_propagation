// BulirschStoer class implementation by Dmitry Sorokin
// Implementation is based on bulirsch_stoer.hpp from boost
//
// The Bulirsch-Stoer is a controlled driver that adjusts both step size
// and order of the method. The algorithm uses the modified midpoint and
// a polynomial extrapolation compute the solution.
//
//
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

#ifndef BulirschStoer_HH
#define BulirschStoer_HH

#include "G4EquationOfMotion.hh"
#include "G4Types.hh"
#include "ModifiedMidpoint.hh"

inline G4double check_error(const G4double yOut[], const G4double yError[],
                            const G4double hstep, const G4double eps){

    G4double errPos2 = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    G4double errMom2 = sqr(yError[3]) + sqr(yError[4]) + sqr(yError[5]);
    G4double Mom2 = sqr(yOut[3]) + sqr(yOut[4]) + sqr(yOut[5]);
    errMom2 /= Mom2;
    errPos2 /= (hstep*hstep);

    return std::sqrt(std::max(errPos2, errMom2)/(eps*eps));
}

enum step_result{
    success,
    fail
};


class BulirschStoer {

public:
    const static G4int m_k_max = 8;

    BulirschStoer(G4EquationOfMotion* equation,G4int nvar,G4double eps_rel,G4double max_dt);

    inline void set_max_dt(G4double max_dt);
    inline void set_max_relative_error(G4double eps_rel);

    //stepper method
    step_result try_step(const G4double in[] , const G4double dxdt[] , G4double& t , G4double out[] , G4double& dt);

    //Reset the internal state of the stepper
    void reset();


private:
    void extrapolate(size_t k , G4double xest[]);
    G4double calc_h_opt(G4double h , G4double error , size_t k ) const;

    G4bool set_k_opt(size_t k , G4double &dt );
    G4bool in_convergence_window(G4int k ) const;
    G4bool should_reject(G4double error , G4int k) const;

private:

    //pointer to the system class
    G4EquationOfMotion* fequation;

    //number of vars to be integrated
    G4int fnvar;

    //relative tolerance
    G4double m_eps_rel;

    //modified midpoint algorithm
    ModifiedMidpoint m_midpoint;

    G4bool m_last_step_rejected;
    G4bool m_first;

    G4double m_dt_last;
    G4double m_t_last;

    //max allowed time step
    G4double m_max_dt;

    G4int m_current_k_opt;

    G4double m_xnew[G4FieldTrack::ncompSVEC];
    G4double m_err[G4FieldTrack::ncompSVEC];
    G4double m_dxdt[G4FieldTrack::ncompSVEC];

    // stores the successive interval counts
    G4int m_interval_sequence[m_k_max+1];

    //Extrapolation coeffs (Neville’s algorithm)
    G4double m_coeff[m_k_max+1][m_k_max];

    // costs for interval count
    G4int m_cost[m_k_max+1];

    // sequence of states for extrapolation
    G4double m_table[m_k_max][G4FieldTrack::ncompSVEC];

    G4double STEPFAC1 , STEPFAC2 , STEPFAC3 , STEPFAC4 , KFAC1 , KFAC2;

    //oprimal step size
    G4double h_opt[m_k_max+1];

    //Work per unit step
    G4double work[m_k_max+1];
};

#include "BulirschStoer.icc"

#endif

