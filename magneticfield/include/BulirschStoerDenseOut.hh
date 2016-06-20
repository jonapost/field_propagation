

#ifndef BULIRSCH_STOER_DENSE_OUT_HH
#define BULIRSCH_STOER_DENSE_OUT_HH

#include "G4EquationOfMotion.hh"
#include "ModifiedMidpointDenseOut.hh"

#include "BulirschStoer.hh"

class bulirsch_stoer_dense_out {

public:

    const static size_t m_k_max = 8;

    bulirsch_stoer_dense_out(
            G4EquationOfMotion* equation,
            G4int nvar,
            G4double eps_rel,
            G4double max_dt,
            G4bool control_interpolation = false ):
        fequation(equation),
        fnvar(nvar),
        m_midpoint(equation,nvar),
        m_max_dt(max_dt) ,
        m_control_interpolation( control_interpolation) ,
        m_last_step_rejected( false ) , m_first( true ) ,
        m_current_state_x1( true ) ,
        STEPFAC1( 0.65 ) , STEPFAC2( 0.94 ) , STEPFAC3( 0.02 ) , STEPFAC4( 4.0 ) , KFAC1( 0.8 ) , KFAC2( 0.9 )
    {
        for(G4int i = 0; i < m_k_max+1; ++i)
        {
            /* only this specific sequence allows for dense output */
            m_interval_sequence[i] = 2 + 4*i;  // 2 6 10 14 ...
            if( i == 0 )
                m_cost[i] = m_interval_sequence[i];
            else
                m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
            for(G4int k = 0 ; k < i ; ++k  )
            {
                const G4double r = static_cast<G4double>( m_interval_sequence[i] ) / static_cast<G4double>( m_interval_sequence[k] );
                m_coeff[i][k] = 1.0 / (r*r - 1.0); // coefficients for extrapolation
            }
            // crude estimate of optimal order

            m_current_k_opt = 4;
            /* no calculation because log10 might not exist for value_type!
            const value_type logfact( -log10( max BOOST_PREVENT_MACRO_SUBSTITUTION( eps_rel , static_cast< value_type >( 1.0E-12 ) ) ) * 0.6 + 0.5 );
            m_current_k_opt = max BOOST_PREVENT_MACRO_SUBSTITUTION( 1 , min BOOST_PREVENT_MACRO_SUBSTITUTION( static_cast<int>( m_k_max-1 ) , static_cast<int>( logfact ) ));
            */
        }
    }

    G4bool try_step(const G4double in[], const G4double dxdt[],
                    G4double& t, G4double out[],
                    G4double dxdt_new[], G4double& dt)
    {
        if(m_max_dt < dt)
        {
            // given step size is bigger then max_dt
            // set limit and return fail
            dt = m_max_dt;
            return false;
        }

        static const value_type val1( 1.0 );

        bool reject = true;

        m_k_final = 0;
        G4double new_h = dt;

        //std::cout << "t=" << t <<", dt=" << dt << ", k_opt=" << m_current_k_opt << ", first: " << m_first << std::endl;

        for(G4int k = 0 ; k <= m_current_k_opt+1 ; ++k)
        {
            m_midpoint.set_steps(m_interval_sequence[k]);
            if( k == 0 )
            {
                m_midpoint.do_step(in, dxdt, out, dt, m_mp_states[k], m_derivs[k]);
            }
            else
            {
                m_midpoint.do_step(in, dxdt, m_table[k-1], dt , m_mp_states[k], m_derivs[k]);
                extrapolate(k , m_table , m_coeff , out);
                // get error estimate
                for (G4int i = 0; i < fnvar; ++i){
                    m_err[i] = out[i] - m_table[0][i];
                }
                const G4double error = check_error(out, m_err, dt, feps_rel);
                h_opt[k] = calc_h_opt( dt , error , k );
                work[k] = static_cast<G4double>( m_cost[k] ) / h_opt[k];

                m_k_final = k;

                if( (k == m_current_k_opt-1) || m_first )
                { // convergence before k_opt ?
                    if( error < 1.0 )
                    {
                        //convergence
                        reject = false;
                        if( (work[k] < KFAC2*work[k-1]) || (m_current_k_opt <= 2) )
                        {
                            // leave order as is (except we were in first round)
                            m_current_k_opt = std::min(m_k_max-1, std::max(2, k+1));
                            new_h = h_opt[k] * static_cast<G4double>(m_cost[k+1]) / static_cast<G4double>( m_cost[k] );
                        }
                        else
                        {
                            m_current_k_opt = std::min(m_k_max-1 , std::max(2, k));
                            new_h = h_opt[k];
                        }
                        break;
                    }
                    else if( should_reject( error , k ) && !m_first )
                    {
                        reject = true;
                        new_h = h_opt[k];
                        break;
                    }
                }
                if( k == m_current_k_opt )
                { // convergence at k_opt ?
                    if( error < 1.0 )
                    {
                        //convergence
                        reject = false;
                        if( (work[k-1] < KFAC2*work[k]) )
                        {
                            m_current_k_opt = std::max(2, m_current_k_opt-1);
                            new_h = h_opt[m_current_k_opt];
                        }
                        else if( (work[k] < KFAC2*work[k-1]) && !m_last_step_rejected )
                        {
                            m_current_k_opt = std::min(m_k_max-1, m_current_k_opt+1 );
                            new_h = h_opt[k]*static_cast<G4double>( m_cost[m_current_k_opt] ) / static_cast<G4double>( m_cost[k] );
                        }
                        else
                            new_h = h_opt[m_current_k_opt];
                        break;
                    }
                    else if( should_reject( error , k ) )
                    {
                        reject = true;
                        new_h = h_opt[m_current_k_opt];
                        break;
                    }
                }
                if( k == m_current_k_opt+1 )
                { // convergence at k_opt+1 ?
                    if( error < 1.0 )
                    {   //convergence
                        reject = false;
                        if( work[k-2] < KFAC2*work[k-1] )
                            m_current_k_opt = std::max(2, m_current_k_opt-1);
                        if( (work[k] < KFAC2*work[m_current_k_opt]) && !m_last_step_rejected )
                            m_current_k_opt = std::min(m_k_max-1, k);
                        new_h = h_opt[m_current_k_opt];
                    }
                    else
                    {
                        reject = true;
                        new_h = h_opt[m_current_k_opt];
                    }
                    break;
                }
            }
        }

        if( !reject )
        {

            //calculate dxdt for next step and dense output
            fequation->RightHandSide(out,dxdt_new);

            //prepare dense output
            G4double error = prepare_dense_output( m_k_final , in , dxdt , out , dxdt_new , dt );

            if( error > 10. ) // we are not as accurate for interpolation as for the steps
            {
                reject = true;
                new_h = dt * std::pow( error , -1./(2*m_k_final+2) );
            }
            else
            {
                t += dt;
            }
        }
        //set next stepsize
        if( !m_last_step_rejected || (new_h < dt) )
        {
            // limit step size
            new_h = std::min(m_max_dt,new_h);
            dt = new_h;
        }

        m_last_step_rejected = reject;
        if( reject )
            return false;
        else
            return true;
    }

    void initialize(const G4double x0[], const G4double& t0, const G4double& dt0)
    {
        memcpy(get_current_state(),x0,sizeof(G4double)*fnvar);
        m_t = t0;
        m_dt = dt0;
        reset();
    }


    /*  =======================================================
     *  the actual step method that should be called from outside (maybe make try_step private?)
     */
    std::pair<G4double, G4double> do_step()
    {
        if(m_first)
        {
            fequation->RightHandSide(get_current_state(), get_current_deriv());
        }

        G4bool res = false;
        m_t_last = m_t;
        while( res == false )
        {
            res = try_step(get_current_state(), get_current_deriv(), m_t,
                           get_old_state(), get_old_deriv(), m_dt );
            m_first = false;
        }
        toggle_current_state();
        return std::make_pair( m_t_last , m_t );
    }

    /* performs the interpolation from a calculated step */
    void calc_state(G4double t , G4double x[]) const
    {
        do_interpolation( t , x );
    }

    const G4double* current_state() const
    {
        return get_current_state();
    }

    G4double current_time() const
    {
        return m_t;
    }

    const G4double* previous_state() const
    {
        return get_old_state();
    }

    G4double previous_time() const
    {
        return m_t_last;
    }

    G4double current_time_step() const
    {
        return m_dt;
    }

    /** \brief Resets the internal state of the stepper. */
    void reset()
    {
        m_first = true;
        m_last_step_rejected = false;
    }


private:

    void extrapolate(G4int k, StateInOut &xest, G4int order_start_index = 0 )
    //polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
    {
        static const value_type val1( 1.0 );
        for( int j=k-1 ; j>0 ; --j )
        {
            for (G4int i = 0; i < fnvar; ++i){
                m_table[j-1][i] = m_table[j][i]*(1 + coeff[k + order_start_index][j + order_start_index] )
                        - table[j-1][i]*coeff[k + order_start_index][j + order_start_index];
            }
        }
        for (G4int i = 0; i < fnvar; ++i){
            xest[i] = table[0][i]*val1 + coeff[k + order_start_index][0 + order_start_index]
                    - xest[i]*coeff[k + order_start_index][0 + order_start_index];
        }
    }


    void extrapolate_dense_out(G4int k , size_t order_start_index = 0 )
    //polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
    {
        // result is written into table[0]
        for(G4int j=k ; j>1 ; --j )
        {
            for (G4int i = 0; i < fnvar; ++i){
                table[j-1][i] = table[j][i]*(1. + coeff[k + order_start_index][j + order_start_index - 1])
                        -table[j-1][i]*coeff[k + order_start_index][j + order_start_index - 1];
            }
        }
        for (G4int i = 0; i < fnvar; ++i){
            table[0][i] = table[1][i]*(1. + coeff[k + order_start_index][order_start_index])
                    -table[0][i]*coeff[k + order_start_index][order_start_index];
        }
    }

    G4double calc_h_opt(G4double h, G4double error, G4int k) const
    {

        G4double expo = 1./(m_interval_sequence[k-1]);
        G4double facmin = std::pow( STEPFAC3 , expo );
        G4double fac;
        if (error == 0.0)
            fac = 1./facmin;
        else
        {
            fac = STEPFAC2 / std::pow( error / STEPFAC1 , expo );
            fac = std::max(facmin/STEPFAC4, std::min(1./facmin , fac) );
        }
        return h*fac;
    }

    G4bool in_convergence_window(G4int k ) const
    {
        if( (k == m_current_k_opt-1) && !m_last_step_rejected )
            return true; // decrease order only if last step was not rejected
        return ( (k == m_current_k_opt) || (k == m_current_k_opt+1) );
    }

    G4bool should_reject(G4double error, G4int k ) const
    {
        if( k == m_current_k_opt-1 )
        {
            const G4double d = m_interval_sequence[m_current_k_opt] * m_interval_sequence[m_current_k_opt+1] /
                (m_interval_sequence[0]*m_interval_sequence[0]);
            //step will fail, criterion 17.3.17 in NR
            return ( error > d*d );
        }
        else if( k == m_current_k_opt )
        {
            const G4double d = m_interval_sequence[m_current_k_opt+1] / m_interval_sequence[0];
            return ( error > d*d );
        } else
            return error > 1.0;
    }

    G4double prepare_dense_output(G4nt k , const G4double x_start[], const G4double dxdt_start[],
                                  const G4double /* x_end */[] , const G4double /*dxdt_end */[] ,
                                  G4double dt )
    /* k is the order to which the result was approximated */
    {

        /* compute the coefficients of the interpolation polynomial
         * we parametrize the interval t .. t+dt by theta = -1 .. 1
         * we use 2k+3 values at the interval center theta=0 to obtain the interpolation coefficients
         * the values are x(t+dt/2) and the derivatives dx/dt , ... d^(2k+2) x / dt^(2k+2) at the midpoints
         * the derivatives are approximated via finite differences
         * all values are obtained from interpolation of the results from the increasing orders of the midpoint calls
         */

        // calculate finite difference approximations to derivatives at the midpoint
        for( int j = 0 ; j<=k ; j++ )
        {
            /* not working with boost units... */
            const G4double d = m_interval_sequence[j] / ( 2. * dt );
            G4double f = 1.0; //factor 1/2 here because our interpolation interval has length 2 !!!
            for( int kappa = 0 ; kappa <= 2*j+1 ; ++kappa )
            {
                calculate_finite_difference( j , kappa , f , dxdt_start );
                f *= d;
            }

            if( j > 0 )
                extrapolate_dense_out( j , m_mp_states , m_coeff );
        }

        time_type d = dt/2;

        // extrapolate finite differences
        for( int kappa = 0 ; kappa<=2*k+1 ; kappa++ )
        {
            for( int j=1 ; j<=(k-kappa/2) ; ++j )
                extrapolate_dense_out( j , m_diffs[kappa] , m_coeff , kappa/2 );

            // extrapolation results are now stored in m_diffs[kappa][0]

            // divide kappa-th derivative by kappa because we need these terms for dense output interpolation
            m_algebra.for_each1( m_diffs[kappa][0].m_v , typename operations_type::template scale< time_type >( static_cast<time_type>(d) ) );

            d *= dt/(2*(kappa+2));
        }

        // dense output coefficients a_0 is stored in m_mp_states[0], a_i for i = 1...2k are stored in m_diffs[i-1][0]

        // the error is just the highest order coefficient of the interpolation polynomial
        // this is because we use only the midpoint theta=0 as support for the interpolation (remember that theta = -1 .. 1)

        value_type error = 0.0;
        if( m_control_interpolation )
        {
            boost::numeric::odeint::copy( m_diffs[2*k+1][0].m_v , m_err.m_v );
            error = m_error_checker.error( m_algebra , x_start , dxdt_start , m_err.m_v , dt );
        }

        return error;
    }

    template< class DerivIn >
    void calculate_finite_difference( size_t j , size_t kappa , value_type fac , const DerivIn &dxdt )
    {
        const int m = m_interval_sequence[j]/2-1;
        if( kappa == 0) // no calculation required for 0th derivative of f
        {
            m_algebra.for_each2( m_diffs[0][j].m_v , m_derivs[j][m].m_v ,
                                 typename operations_type::template scale_sum1< value_type >( fac ) );
        }
        else
        {
            // calculate the index of m_diffs for this kappa-j-combination
            const int j_diffs = j - kappa/2;

            m_algebra.for_each2( m_diffs[kappa][j_diffs].m_v , m_derivs[j][m+kappa].m_v ,
                                 typename operations_type::template scale_sum1< value_type >( fac ) );
            value_type sign = -1.0;
            int c = 1;
            //computes the j-th order finite difference for the kappa-th derivative of f at t+dt/2 using function evaluations stored in m_derivs
            for( int i = m+static_cast<int>(kappa)-2 ; i >= m-static_cast<int>(kappa) ; i -= 2 )
            {
                if( i >= 0 )
                {
                    m_algebra.for_each3( m_diffs[kappa][j_diffs].m_v , m_diffs[kappa][j_diffs].m_v , m_derivs[j][i].m_v ,
                                         typename operations_type::template scale_sum2< value_type , value_type >( 1.0 ,
                                                                                                                   sign * fac * boost::math::binomial_coefficient< value_type >( kappa , c ) ) );
                }
                else
                {
                    m_algebra.for_each3( m_diffs[kappa][j_diffs].m_v , m_diffs[kappa][j_diffs].m_v , dxdt ,
                                         typename operations_type::template scale_sum2< value_type , value_type >( 1.0 , sign * fac ) );
                }
                sign *= -1;
                ++c;
            }
        }
    }

    template< class StateOut >
    void do_interpolation( time_type t , StateOut &out ) const
    {
        // interpolation polynomial is defined for theta = -1 ... 1
        // m_k_final is the number of order-iterations done for the last step - it governs the order of the interpolation polynomial
        const value_type theta = 2 * get_unit_value( (t - m_t_last) / (m_t - m_t_last) ) - 1;
        // we use only values at interval center, that is theta=0, for interpolation
        // our interpolation polynomial is thus of order 2k+2, hence we have 2k+3 terms

        boost::numeric::odeint::copy( m_mp_states[0].m_v , out );
        // add remaining terms: x += a_1 theta + a2 theta^2 + ... + a_{2k} theta^{2k}
        value_type theta_pow( theta );
        for( size_t i=0 ; i<=2*m_k_final+1 ; ++i )
        {
            m_algebra.for_each3( out , out , m_diffs[i][0].m_v ,
                                 typename operations_type::template scale_sum2< value_type >( static_cast<value_type>(1) , theta_pow ) );
            theta_pow *= theta;
        }
    }



    state_type& get_current_state( void )
    {
        return m_current_state_x1 ? m_x1.m_v : m_x2.m_v ;
    }

    const state_type& get_current_state( void ) const
    {
        return m_current_state_x1 ? m_x1.m_v : m_x2.m_v ;
    }

    state_type& get_old_state( void )
    {
        return m_current_state_x1 ? m_x2.m_v : m_x1.m_v ;
    }

    const state_type& get_old_state( void ) const
    {
        return m_current_state_x1 ? m_x2.m_v : m_x1.m_v ;
    }

    deriv_type& get_current_deriv( void )
    {
        return m_current_state_x1 ? m_dxdt1.m_v : m_dxdt2.m_v ;
    }

    const deriv_type& get_current_deriv( void ) const
    {
        return m_current_state_x1 ? m_dxdt1.m_v : m_dxdt2.m_v ;
    }

    deriv_type& get_old_deriv( void )
    {
        return m_current_state_x1 ? m_dxdt2.m_v : m_dxdt1.m_v ;
    }

    const deriv_type& get_old_deriv( void ) const
    {
        return m_current_state_x1 ? m_dxdt2.m_v : m_dxdt1.m_v ;
    }


    void toggle_current_state( void )
    {
        m_current_state_x1 = ! m_current_state_x1;
    }


private:
    G4EquationOfMotion* fequation;
    G4int fnvar;
    ModifiedMidpointDenseOut m_midpoint;
    G4double feps_rel;


    G4double m_max_dt;

    G4bool m_control_interpolation;

    G4bool m_last_step_rejected;
    G4bool m_first;

    G4double m_t;
    G4double m_dt;
    G4double m_dt_last;
    G4double m_t_last;

    G4int m_current_k_opt;
    G4int m_k_final;


    G4double m_x1[G4FieldTrack::ncompSVEC] , m_x2[G4FieldTrack::ncompSVEC];
    G4double m_dxdt1[G4FieldTrack::ncompSVEC] , m_dxdt2[G4FieldTrack::ncompSVEC];
    G4double m_err[G4FieldTrack::ncompSVEC];
    G4bool m_current_state_x1;



    G4double m_error[m_k_max]; // errors of repeated midpoint steps and extrapolations
    G4int m_interval_sequence[m_k_max+1]; // stores the successive interval counts
    G4double m_coeff[m_k_max+1][m_k_max];
    G4int m_cost[m_k_max+1]; // costs for interval count

    G4double m_table[m_k_max][G4FieldTrack::ncompSVEC]; // sequence of states for extrapolation

    //for dense output:
    G4double m_mp_states[m_k_max+1][G4FieldTrack::ncompSVEC]; // sequence of approximations of x at distance center
    G4double m_derivs[m_k_max+1][2+4*(m_k_max+1)]; // table of function values
    G4double m_diffs[2*m_k_max+2][2*m_k_max+2][G4FieldTrack::ncompSVEC]; // table of function values

    //wrapped_state_type m_a1 , m_a2 , m_a3 , m_a4;

    G4double STEPFAC1 , STEPFAC2 , STEPFAC3 , STEPFAC4 , KFAC1 , KFAC2;


    G4double h_opt[m_k_max+1];
    G4double work[m_k_max+1];
};



/********** DOXYGEN **********/

/**
 * \class bulirsch_stoer_dense_out
 * \brief The Bulirsch-Stoer algorithm.
 *
 * The Bulirsch-Stoer is a controlled stepper that adjusts both step size
 * and order of the method. The algorithm uses the modified midpoint and
 * a polynomial extrapolation compute the solution. This class also provides
 * dense output facility.
 *
 * \tparam State The state type.
 * \tparam Value The value type.
 * \tparam Deriv The type representing the time derivative of the state.
 * \tparam Time The time representing the independent variable - the time.
 * \tparam Algebra The algebra type.
 * \tparam Operations The operations type.
 * \tparam Resizer The resizer policy type.
 */

    /**
     * \fn bulirsch_stoer_dense_out::bulirsch_stoer_dense_out( value_type eps_abs , value_type eps_rel , value_type factor_x , value_type factor_dxdt , bool control_interpolation )
     * \brief Constructs the bulirsch_stoer class, including initialization of
     * the error bounds.
     *
     * \param eps_abs Absolute tolerance level.
     * \param eps_rel Relative tolerance level.
     * \param factor_x Factor for the weight of the state.
     * \param factor_dxdt Factor for the weight of the derivative.
     * \param control_interpolation Set true to additionally control the error of
     * the interpolation.
     */

    /**
     * \fn bulirsch_stoer_dense_out::try_step( System system , const StateIn &in , const DerivIn &dxdt , time_type &t , StateOut &out , DerivOut &dxdt_new , time_type &dt )
     * \brief Tries to perform one step.
     *
     * This method tries to do one step with step size dt. If the error estimate
     * is to large, the step is rejected and the method returns fail and the
     * step size dt is reduced. If the error estimate is acceptably small, the
     * step is performed, success is returned and dt might be increased to make
     * the steps as large as possible. This method also updates t if a step is
     * performed. Also, the internal order of the stepper is adjusted if required.
     *
     * \param system The system function to solve, hence the r.h.s. of the ODE.
     * It must fulfill the Simple System concept.
     * \param in The state of the ODE which should be solved.
     * \param dxdt The derivative of state.
     * \param t The value of the time. Updated if the step is successful.
     * \param out Used to store the result of the step.
     * \param dt The step size. Updated.
     * \return success if the step was accepted, fail otherwise.
     */

    /**
     * \fn bulirsch_stoer_dense_out::initialize( const StateType &x0 , const time_type &t0 , const time_type &dt0 )
     * \brief Initializes the dense output stepper.
     *
     * \param x0 The initial state.
     * \param t0 The initial time.
     * \param dt0 The initial time step.
     */

    /**
     * \fn bulirsch_stoer_dense_out::do_step( System system )
     * \brief Does one time step. This is the main method that should be used to
     * integrate an ODE with this stepper.
     * \note initialize has to be called before using this method to set the
     * initial conditions x,t and the stepsize.
     * \param system The system function to solve, hence the r.h.s. of the
     * ordinary differential equation. It must fulfill the Simple System concept.
     * \return Pair with start and end time of the integration step.
     */

    /**
     * \fn bulirsch_stoer_dense_out::calc_state( time_type t , StateOut &x ) const
     * \brief Calculates the solution at an intermediate point within the last step
     * \param t The time at which the solution should be calculated, has to be
     * in the current time interval.
     * \param x The output variable where the result is written into.
     */

    /**
     * \fn bulirsch_stoer_dense_out::current_state( void ) const
     * \brief Returns the current state of the solution.
     * \return The current state of the solution x(t).
     */

    /**
     * \fn bulirsch_stoer_dense_out::current_time( void ) const
     * \brief Returns the current time of the solution.
     * \return The current time of the solution t.
     */

    /**
     * \fn bulirsch_stoer_dense_out::previous_state( void ) const
     * \brief Returns the last state of the solution.
     * \return The last state of the solution x(t-dt).
     */

    /**
     * \fn bulirsch_stoer_dense_out::previous_time( void ) const
     * \brief Returns the last time of the solution.
     * \return The last time of the solution t-dt.
     */

    /**
     * \fn bulirsch_stoer_dense_out::current_time_step( void ) const
     * \brief Returns the current step size.
     * \return The current step size.
     */

    /**
     * \fn bulirsch_stoer_dense_out::adjust_size( const StateIn &x )
     * \brief Adjust the size of all temporaries in the stepper manually.
     * \param x A state from which the size of the temporaries to be resized is deduced.
     */



#endif
