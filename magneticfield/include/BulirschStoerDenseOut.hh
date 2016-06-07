#ifndef BulirschStoerDenseOut_HH
#define BulirschStoerDenseOut_HH

#include <iostream>
#include <algorithm>
#include <vector>

#include "ModifiedMidpoint.hh"

typedef std::vector<int> int_vector;
typedef std::vector<state_vector_type> state_table_type;
typedef std::vector<double> double_vector;
typedef std::vector<std::vector<double>> double_matrix;

class DefaultErrorChecker{
public:
    DefaultErrorChecker(double eps):m_eps_rel(eps){
    }
    double operator() (const state_type& in, const state_type& err, double hstep){
        double errh2 = 0,errp2 = 0,p2 = 0;
        for (int i = 0; i < 3; ++i){
            errh2 += err[i]*err[i];
            errp2 += err[i+3]*err[i+3];
            p2 += in[i+3]*in[i+3];
        }
        errh2 /= hstep;
        errp2 /= p2;
        return std::max(errh2,errp2)/(m_eps_rel*m_eps_rel);
    }

private:
    double m_eps_rel;
};

enum controlled_step_result{
    fail,
    success
};

class BulirschStoerDenseOut{
public:

    const static int m_k_max = 8;

    BulirschStoerDenseOut(
        int nvar,
        double eps_rel,
        double max_dt = 0,
        bool control_interpolation = false )
        : fnvar(nvar),
          m_error_checker(eps_rel),
          m_max_dt(max_dt) ,
          m_control_interpolation(control_interpolation) ,
          m_last_step_rejected( false ) , m_first( true ) ,
          m_current_state_x1( true ) ,
          m_error( m_k_max ) ,
          m_interval_sequence( m_k_max+1 ) ,
          m_coeff( m_k_max+1 ) ,
          m_cost( m_k_max+1 ) ,
          m_table( m_k_max ) ,
          m_mp_states( m_k_max+1 ) ,
          m_derivs( m_k_max+1 ) ,
          m_diffs( 2*m_k_max+2 ) ,
          STEPFAC1( 0.65 ) , STEPFAC2( 0.94 ) , STEPFAC3( 0.02 ) , STEPFAC4( 4.0 ) , KFAC1( 0.8 ) , KFAC2( 0.9 )
    {
        for(unsigned short i = 0; i < m_k_max+1; ++i)
        {
            /* only this specific sequence allows for dense output */
            m_interval_sequence[i] = 2 + 4*i;  // 2 6 10 14 ...
            m_derivs[i].resize(m_interval_sequence[i]);
            if(i == 0)
                m_cost[i] = m_interval_sequence[i];
            else
                m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
            m_coeff[i].resize(i);
            for( int k = 0 ; k < i ; ++k  )
            {
                const double r = m_interval_sequence[i] / m_interval_sequence[k];
                m_coeff[i][k] = 1.0 / (r*r - 1.0); // coefficients for extrapolation
            }
            // crude estimate of optimal order

            m_current_k_opt = 4;
            /* no calculation because log10 might not exist for double!
            const double logfact( -log10( max BOOST_PREVENT_MACRO_SUBSTITUTION( eps_rel , static_cast< double >( 1.0E-12 ) ) ) * 0.6 + 0.5 );
            m_current_k_opt = max BOOST_PREVENT_MACRO_SUBSTITUTION( 1 , min BOOST_PREVENT_MACRO_SUBSTITUTION( static_cast<int>( m_k_max-1 ) , static_cast<int>( logfact ) ));
            */
        }
        int num = 1;
        for( int i = 2*(m_k_max)+1 ; i >=0  ; i-- )
        {
            m_diffs[i].resize( num );
            num += (i+1)%2;
        }
    }

    template< class System>
    controlled_step_result try_step( System system , const state_type& in , const state_type& dxdt,
                                     double &t , state_type& out , state_type& dxdt_new , double& dt )
    {
        if(fabs(dt) > m_max_dt)
        {
            // given step size is bigger then max_dt
            // set limit and return fail
            dt = m_max_dt;
            return fail;
        }

        bool reject = true;

        double_vector h_opt( m_k_max+1 );
        double_vector work( m_k_max+1 );

        m_k_final = 0;
        double new_h = dt;

        //std::cout << "t=" << t <<", dt=" << dt << ", k_opt=" << m_current_k_opt << ", first: " << m_first << std::endl;

        for( int k = 0 ; k <= m_current_k_opt+1 ; k++ )
        {
            m_midpoint.set_steps(m_interval_sequence[k]);
            if( k == 0 )
            {
                m_midpoint.do_step( system , in , dxdt, out , dt , m_mp_states[k] , m_derivs[k]);
            }
            else
            {
                m_midpoint.do_step( system , in , dxdt, m_table[k-1], dt , m_mp_states[k], m_derivs[k] );
                extrapolate( k , m_table , m_coeff , out );
                // get error estimate
                for (int i = 0; i < fnvar; ++i){
                    m_err[i] = out[i] - m_table[0][i];
                }
                const double error = m_error_checker(in, m_err, dt);
                h_opt[k] = calc_h_opt( dt , error , k );
                work[k] = m_cost[k]  / h_opt[k];

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
                            m_current_k_opt = std::min(m_k_max-1 , std::max(2 ,k+1));
                            new_h = h_opt[k] * m_cost[k+1]/m_cost[k];
                        } else {
                            m_current_k_opt = std::min(m_k_max-1 , std::max (2 ,k));
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
                            m_current_k_opt = std::max(2 ,m_current_k_opt-1);
                            new_h = h_opt[m_current_k_opt];
                        }
                        else if( (work[k] < KFAC2*work[k-1]) && !m_last_step_rejected )
                        {
                            m_current_k_opt = std::min(m_k_max-1 , m_current_k_opt+1 );
                            new_h = h_opt[k]*m_cost[m_current_k_opt] / m_cost[k];
                        } else
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
                            m_current_k_opt = std::max(2 ,m_current_k_opt-1);
                        if( (work[k] < KFAC2*work[m_current_k_opt]) && !m_last_step_rejected )
                            m_current_k_opt = std::min(m_k_max-1 , k);
                        new_h = h_opt[m_current_k_opt];
                    } else
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
            system( out , dxdt_new);

            //prepare dense output
            double error = prepare_dense_output( m_k_final , in , dxdt , out , dxdt_new , dt );

            if( error > 10) // we are not as accurate for interpolation as for the steps
            {
                reject = true;
                new_h = dt * std::pow( error , -1/(2*m_k_final+2) );
            } else {
                t += dt;
            }
        }
        //set next stepsize
        if( !m_last_step_rejected || (new_h < dt) )
        {
            // limit step size
            if( m_max_dt != 0)
            {
                new_h = std::min(fabs(m_max_dt),fabs(new_h));
            }
            dt = new_h;
        }

        m_last_step_rejected = reject;
        if( reject )
            return fail;
        else
            return success;
    }

    void initialize( const state_type &x0 , const double &t0 , const double &dt0 )
    {
        //m_resizer.adjust_size( x0 , detail::bind( &controlled_error_bs_type::template resize_impl< StateType > , detail::ref( *this ) , detail::_1 ) );
        get_current_state() = x0;
        m_t = t0;
        m_dt = dt0;
        reset();
    }


    /*  =======================================================
     *  the actual step method that should be called from outside (maybe make try_step private?)
     */
    template< class System >
    std::pair< double , double > do_step(System system)
    {
        if( m_first )
        {
            system( get_current_state() , get_current_deriv());
        }

        controlled_step_result res = fail;
        m_t_last = m_t;
        while( res == fail )
        {
            res = try_step( system , get_current_state() , get_current_deriv() , m_t , get_old_state() , get_old_deriv() , m_dt );
            m_first = false;
        }
        toggle_current_state();
        return std::make_pair( m_t_last , m_t );
    }

    /* performs the interpolation from a calculated step */
    void calc_state(double t , state_type &x) const
    {
        do_interpolation( t , x );
    }

    const state_type& current_state() const
    {
        return get_current_state();
    }

    double current_time() const
    {
        return m_t;
    }

    const state_type& previous_state() const
    {
        return get_old_state();
    }

    double previous_time() const
    {
        return m_t_last;
    }

    double current_time_step() const
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

    void extrapolate( int k , state_vector_type &table , const double_matrix &coeff , state_type &xest , int order_start_index = 0 )
    //polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
    {
        for( int j=k-1 ; j>0 ; --j )
        {
            for (int i = 0; i < fnvar; ++i){
                table[j-1][i] = table[j][i] * (1 + coeff[k + order_start_index][j + order_start_index]) +
                                table[j-1][i] * (-coeff[k + order_start_index][j + order_start_index]);
            }

        }
        for (int i = 0; i < fnvar; ++i){
            xest[i] = table[0][i] * (1 + coeff[k + order_start_index][0 + order_start_index]) +
                      xest[i] * (-coeff[k + order_start_index][0 + order_start_index]);
        }
    }


    void extrapolate_dense_out( int k , state_vector_type &table , const double_matrix &coeff , int order_start_index = 0 )
    //polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
    {
        // result is written into table[0]
        for( int j=k ; j>1 ; --j )
        {
            for (int i = 0; i < fnvar; ++i){
                table[j-1][i] = table[j][i] * (1 + coeff[k + order_start_index][j + order_start_index - 1]) +
                                table[j-1][i] * (-coeff[k + order_start_index][j + order_start_index - 1]);
            }
        }
        for (int i = 0; i < fnvar; ++i){
            table[0][i] = table[1][i] * (1 + coeff[k + order_start_index][order_start_index]) +
                          table[0][i] * (-coeff[k + order_start_index][order_start_index]);
        }

    }

    double calc_h_opt( double h , double error , int k ) const
    {
        double expo = 1./(m_interval_sequence[k-1]);
        double facmin = std::pow( STEPFAC3 , expo );
        double fac;
        if (error == 0.0)
            fac = 1./facmin;
        else
        {
            fac = STEPFAC2 / std::pow( error / STEPFAC1 , expo );
            fac = std::max(facmin/STEPFAC4 , std::min(1.0/facmin , fac) );
        }
        return h*fac;
    }

    bool in_convergence_window( int k ) const
    {
        if( (k == m_current_k_opt-1) && !m_last_step_rejected )
            return true; // decrease order only if last step was not rejected
        return ( (k == m_current_k_opt) || (k == m_current_k_opt+1) );
    }

    bool should_reject( double error , int k ) const
    {
        if( k == m_current_k_opt-1 )
        {
            const double d = m_interval_sequence[m_current_k_opt] * m_interval_sequence[m_current_k_opt+1] /
                (m_interval_sequence[0]*m_interval_sequence[0]);
            //step will fail, criterion 17.3.17 in NR
            return ( error > d*d );
        }
        else if( k == m_current_k_opt )
        {
            const double d = m_interval_sequence[m_current_k_opt+1] / m_interval_sequence[0];
            return ( error > d*d );
        } else
            return error > 1.0;
    }

    double prepare_dense_output( int k , const state_type &x_start , const state_type &dxdt_start ,
                                 const state_type & /* x_end */ , const state_type & /*dxdt_end */ ,
                                 double dt )
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
            const double d = m_interval_sequence[j] / ( 2 * dt );
            double f = 1.0; //factor 1/2 here because our interpolation interval has length 2 !!!
            for( int kappa = 0 ; kappa <= 2*j+1 ; ++kappa )
            {
                calculate_finite_difference( j , kappa , f , dxdt_start );
                f *= d;
            }

            if( j > 0 )
                extrapolate_dense_out( j , m_mp_states , m_coeff );
        }

        double d = dt/2;

        // extrapolate finite differences
        for( int kappa = 0 ; kappa<=2*k+1 ; kappa++ )
        {
            for( int j=1 ; j<=(k-kappa/2) ; ++j )
                extrapolate_dense_out( j , m_diffs[kappa] , m_coeff , kappa/2 );

            // extrapolation results are now stored in m_diffs[kappa][0]

            // divide kappa-th derivative by kappa because we need these terms for dense output interpolation
            for (int i = 0; i < fnvar; ++i){
                m_diffs[kappa][0][i] *= d;
            }
            d *= dt/(2*(kappa+2));
        }

        // dense output coefficients a_0 is stored in m_mp_states[0], a_i for i = 1...2k are stored in m_diffs[i-1][0]

        // the error is just the highest order coefficient of the interpolation polynomial
        // this is because we use only the midpoint theta=0 as support for the interpolation (remember that theta = -1 .. 1)

        double error = 0.0;
        if( m_control_interpolation )
        {
            m_err = m_diffs[2*k+1][0];
            error = m_error_checker(x_start, m_err , dt );
        }

        return error;
    }

    void calculate_finite_difference( int j , int kappa , double fac , const state_type &dxdt )
    {
        const int m = m_interval_sequence[j]/2-1;
        if( kappa == 0) // no calculation required for 0th derivative of f
        {
            for (int i = 0; i < fnvar; ++i){
                m_diffs[0][j][i] = m_derivs[j][m][i] * fac;
            }
        }
        else
        {
            // calculate the index of m_diffs for this kappa-j-combination
            const int j_diffs = j - kappa/2;

            for (int i = 0; i < fnvar; ++i){
                m_diffs[kappa][j_diffs][i] = m_derivs[j][m+kappa][i] * fac;
            }

            double sign = -1.0;
            int c = 1;
            //computes the j-th order finite difference for the kappa-th derivative of f at t+dt/2 using function evaluations stored in m_derivs
            for( int i = m + kappa - 2 ; i >= m - kappa ; i -= 2 )
            {
                if( i >= 0 )
                {
                    for (int idx = 0; idx < fnvar; ++idx){
                        m_diffs[kappa][j_diffs][idx] =  m_diffs[kappa][j_diffs][idx] +
                                                   m_derivs[j][i][idx] * sign * fac * boost::math::binomial_coefficient< double >( kappa , c );
                    }

                }
                else
                {
                    for (int idx = 0; idx < fnvar; ++idx){
                        m_diffs[kappa][j_diffs][idx] = m_diffs[kappa][j_diffs][idx] +
                                                  dxdt[idx] * sign * fac;
                    }
                }
                sign *= -1;
                ++c;
            }
        }
    }

    void do_interpolation( double t , state_type &out ) const
    {
        // interpolation polynomial is defined for theta = -1 ... 1
        // m_k_final is the number of order-iterations done for the last step - it governs the order of the interpolation polynomial
        const double theta = 2 *  (t - m_t_last) / (m_t - m_t_last)  - 1;
        // we use only values at interval center, that is theta=0, for interpolation
        // our interpolation polynomial is thus of order 2k+2, hence we have 2k+3 terms

        out = m_mp_states[0];
        // add remaining terms: x += a_1 theta + a2 theta^2 + ... + a_{2k} theta^{2k}
        double theta_pow( theta );
        for( int i=0 ; i<=2*m_k_final+1 ; ++i )
        {
            for (int idx = 0; idx < fnvar; ++idx){
                out[idx] = out[idx] + m_diffs[i][0][idx] * theta_pow;
            }
            theta_pow *= theta;
        }
    }


    state_type& get_current_state()
    {
        return m_current_state_x1 ? m_x1 : m_x2;
    }

    const state_type& get_current_state() const
    {
        return m_current_state_x1 ? m_x1 : m_x2;
    }

    state_type& get_old_state()
    {
        return m_current_state_x1 ? m_x2 : m_x1;
    }

    const state_type& get_old_state() const
    {
        return m_current_state_x1 ? m_x2 : m_x1;
    }

    state_type& get_current_deriv()
    {
        return m_current_state_x1 ? m_dxdt1 : m_dxdt2;
    }

    const state_type& get_current_deriv() const
    {
        return m_current_state_x1 ? m_dxdt1 : m_dxdt2;
    }

    state_type& get_old_deriv( void )
    {
        return m_current_state_x1 ? m_dxdt2 : m_dxdt1;
    }

    const state_type& get_old_deriv() const
    {
        return m_current_state_x1 ? m_dxdt2 : m_dxdt1;
    }


    void toggle_current_state()
    {
        m_current_state_x1 = ! m_current_state_x1;
    }

    const int fnvar;

    DefaultErrorChecker m_error_checker;
    ModifiedMidpointDenseOut m_midpoint;

    double m_max_dt;

    bool m_control_interpolation;

    bool m_last_step_rejected;
    bool m_first;

    double m_t;
    double m_dt;
    double m_dt_last;
    double m_t_last;

    int m_current_k_opt;
    int m_k_final;

    state_type m_x1 , m_x2;
    state_type m_dxdt1 , m_dxdt2;
    state_type m_err;
    bool m_current_state_x1;



    double_vector m_error; // errors of repeated midpoint steps and extrapolations
    int_vector m_interval_sequence; // stores the successive interval counts
    double_matrix m_coeff;
    int_vector m_cost; // costs for interval count

    state_vector_type m_table; // sequence of states for extrapolation

    //for dense output:
    state_vector_type m_mp_states; // sequence of approximations of x at distance center
    state_table_type m_derivs; // table of function values
    state_table_type m_diffs; // table of function values

    //wrapped_state_type m_a1 , m_a2 , m_a3 , m_a4;

    double STEPFAC1 , STEPFAC2 , STEPFAC3 , STEPFAC4 , KFAC1 , KFAC2;
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
     * \fn bulirsch_stoer_dense_out::bulirsch_stoer_dense_out( double eps_abs , double eps_rel , double factor_x , double factor_dxdt , bool control_interpolation )
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
     * \fn bulirsch_stoer_dense_out::try_step( System system , const StateIn &in , const DerivIn &dxdt , double &t , StateOut &out , DerivOut &dxdt_new , double &dt )
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
     * \fn bulirsch_stoer_dense_out::initialize( const StateType &x0 , const double &t0 , const double &dt0 )
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
     * \fn bulirsch_stoer_dense_out::calc_state( double t , StateOut &x ) const
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



#endif // BOOST_NUMERIC_ODEINT_STEPPER_BULIRSCH_STOER_HPP_INCLUDED
