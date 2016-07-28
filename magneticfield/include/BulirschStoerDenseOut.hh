
#ifndef BulirschStoerDenseOut_HH
#define BulirschStoerDenseOut_HH

#include "G4EquationOfMotion.hh"
#include "G4Types.hh"
#include "ModifiedMidpointDenseOut.hh"
#include "BulirschStoer.hh"

inline G4double factorial(G4int val)
{
    G4double res = 1;
    while(val > 1)
    {
        res *= val;
        --val;
    }
    return res;
}

inline G4double binomial_coefficient(G4int kappa , G4int c)
{//kappa!/(c!*(kappa-c)!)
    return factorial(kappa)/factorial(c)/factorial(kappa-c);
}

class BulirschStoerDenseOut {

public:

    const static G4int m_k_max = 8;

    BulirschStoerDenseOut(G4EquationOfMotion* equation,G4int nvar,G4double eps,
                          G4bool control_interpolation = false);

    //stepper method
    step_result try_step(const G4double in[], const G4double dxdt[], G4double& t, G4double out[], G4double &dt);

    // Resets the internal state of the stepper.
    void reset();

    void set_max_relative_error(G4double eps)
    { eps_rel = eps;}

    void SetEquationOfMotion(G4EquationOfMotion* newEquation)
    {fequation = newEquation; }

    G4EquationOfMotion* GetEquationOfMotion()
    {return fequation;}

    G4int GetNumberOfVariables()
    {return fnvar;}

    void SetNumberOfVariables(G4int nvar)
    {fnvar = nvar;}

    G4int GetCurrentOrder()
    {return pow(2,m_current_k_opt);}

    void do_interpolation(G4double t, G4double out[]) const;


private:

    void extrapolate(G4int k, G4double xest[]);


    void extrapolate_dense_out(G4int k, G4double table[][G4FieldTrack::ncompSVEC],
    const G4double coeff[m_k_max+1][m_k_max], G4int order_start_index = 0);

    G4double calc_h_opt(G4double h, G4double error, G4int k) const;

    G4bool in_convergence_window(G4int k) const;

    G4bool should_reject(G4double error, G4int k) const;

    G4double prepare_dense_output(G4int k, const G4double dxdt_start[] ,
                                  const G4double x_end[],G4double dt );

    void calculate_finite_difference(G4int j, G4int kappa, G4double fac, const G4double dxdt[]);




    G4EquationOfMotion* fequation;
    G4int fnvar;
    G4double eps_rel;

    ModifiedMidpointDenseOut m_midpoint;


    G4bool m_control_interpolation;

    G4bool m_last_step_rejected;
    G4bool m_first;

    G4double m_t;
    G4double m_dt;
    G4double m_dt_last;
    G4double m_t_last;

    G4int m_current_k_opt;
    G4int m_k_final;

    //oprimal step size
    G4double h_opt[m_k_max+1];

    //Work per unit step
    G4double work[m_k_max+1];

    // errors of repeated midpoint steps and extrapolations
    G4double m_err[G4FieldTrack::ncompSVEC];

    // stores the successive interval counts
    G4int m_interval_sequence[m_k_max+1];

    //Extrapolation coeffs (Neville’s algorithm)
    G4double m_coeff[m_k_max+1][m_k_max];

    // costs for interval count
    G4int m_cost[m_k_max+1];

    // sequence of states for extrapolation
    G4double m_table[m_k_max][G4FieldTrack::ncompSVEC]; // sequence of states for extrapolation

    //for dense output:
    G4double m_mp_states[m_k_max+1][G4FieldTrack::ncompSVEC]; // sequence of approximations of x at distance center
    G4double m_derivs[m_k_max+1][2+4*m_k_max][G4FieldTrack::ncompSVEC]; // table of function values
    G4double m_diffs[2*m_k_max+2][m_k_max+1][G4FieldTrack::ncompSVEC]; // table of function values

    G4double STEPFAC1 , STEPFAC2 , STEPFAC3 , STEPFAC4 , KFAC1 , KFAC2;
};

#endif
