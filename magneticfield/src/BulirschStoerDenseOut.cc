#include "BulirschStoerDenseOut.hh"

BulirschStoerDenseOut::BulirschStoerDenseOut(G4EquationOfMotion* equation,
                                             G4int nvar,G4double
                                             eps,G4bool
                                             control_interpolation):
    // fequation(equation),
    fnvar(nvar),
    eps_rel(eps),
    m_midpoint(equation,nvar),
    m_control_interpolation(control_interpolation),
    m_last_step_rejected(false), m_first(true),
    STEPFAC1(0.65), STEPFAC2(0.94), STEPFAC3(0.02), STEPFAC4(4.0), /*KFAC1(0.8),*/  KFAC2(0.9)
{

    for(G4int i = 0; i < m_k_max+1; ++i)
    {
        /* only this specific sequence allows for dense output */
        m_interval_sequence[i] = 2 + 4*i;  // 2 6 10 14 ...
        if(i == 0)
            m_cost[i] = m_interval_sequence[i];
        else
            m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
        for(G4int k = 0; k < i; ++k)
        {
            const G4double r = G4double(m_interval_sequence[i]) / G4double(m_interval_sequence[k]);
            m_coeff[i][k] = 1.0 / (r*r - 1.0 ); // coefficients for extrapolation
        }
        // crude estimate of optimal order

        m_current_k_opt = 4;
        /* no calculation because log10 might not exist for G4double!
        const G4double logfact( -log10( max BOOST_PREVENT_MACRO_SUBSTITUTION( eps_rel , static_cast< G4double >( 1.0E-12 ) ) ) * 0.6 + 0.5 );
        m_current_k_opt = max BOOST_PREVENT_MACRO_SUBSTITUTION( 1 , min BOOST_PREVENT_MACRO_SUBSTITUTION( static_cast<int>( m_k_max-1 ) , static_cast<int>( logfact ) ));
        */
    }
}


step_result BulirschStoerDenseOut::try_step(const G4double in[], const G4double dxdt[], G4double& t,
                           G4double out[], G4double &dt)
{

    if(dt != m_dt_last)
    {
        reset(); // step size changed from outside -> reset
    }

    G4bool reject = true;

    m_k_final = 0;
    G4double new_h = dt;

    //std::cout << "t=" << t <<", dt=" << dt << ", k_opt=" << m_current_k_opt << ", first: " << m_first << std::endl;

    for(G4int k = 0; k <= m_current_k_opt+1; ++k)
    {
        m_midpoint.set_steps(m_interval_sequence[k]);
        if( k == 0 )
        {
            m_midpoint.do_step(in , dxdt , out , dt , m_mp_states[k] , m_derivs[k]);
        }
        else
        {
            m_midpoint.do_step(in, dxdt, m_table[k-1], dt, m_mp_states[k], m_derivs[k]);
            extrapolate(k , out);
            // get error estimate
            for (G4int ivar = 0; ivar < fnvar; ++ivar)
            {
                m_err[ivar] = out[ivar] - m_table[0][ivar];
            }
            const G4double error = check_error(out, m_err, dt, eps_rel);
            h_opt[k] = calc_h_opt(dt, error, k);
            work[k] = G4double(m_cost[k]) / h_opt[k];

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
                        new_h = h_opt[k] * G4double(m_cost[k+1]) / G4double(m_cost[k]);
                    }
                    else
                    {
                        m_current_k_opt = std::min(m_k_max-1, std::max(2, k));
                        new_h = h_opt[k];
                    }
                    break;
                }
                else if(should_reject( error , k ) && !m_first)
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
                    else if((work[k] < KFAC2*work[k-1]) && !m_last_step_rejected)
                    {
                        m_current_k_opt = std::min(m_k_max-1, m_current_k_opt+1);
                        new_h = h_opt[k];
                        new_h *= G4double(m_cost[m_current_k_opt]) / G4double(m_cost[k]);
                    }
                    else
                        new_h = h_opt[m_current_k_opt];
                    break;
                }
                else if(should_reject(error, k))
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
                    {
                        m_current_k_opt = std::max(2, m_current_k_opt-1);
                    }
                    if( (work[k] < KFAC2*work[m_current_k_opt]) && !m_last_step_rejected )
                    {
                        m_current_k_opt = std::min(m_k_max-1, k);
                    }
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


    if(!reject)
    {
        G4double error = prepare_dense_output(m_k_final, dxdt, out, dt);
        if (error > eps_rel)
        {
            reject = true;
            m_t_last = m_t = 0;
            new_h = dt * std::pow(error , -1./(2*m_k_final+2));
        }
        else
        {
            m_t_last = t;
            m_t = t + dt;
            t += dt;
        }
    }
    else
    {
        m_t_last = m_t = 0;
    }
    //set next stepsize
    if(!m_last_step_rejected || (new_h < dt))
    {
        m_dt_last = new_h;
        dt = new_h;
    }

    m_last_step_rejected = reject;
    m_first = false;

    if(reject)
        return fail;
    else
        return success;
}

void BulirschStoerDenseOut::reset()
{
    m_first = true;
    m_last_step_rejected = false;
}

void BulirschStoerDenseOut::extrapolate(G4int k, G4double xest[])
//polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
{
    for(G4int j=k-1 ; j>0 ; --j)
    {
        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            m_table[j-1][ivar] = m_table[j][ivar]*(1. + m_coeff[k][j])
                    - m_table[j-1][ivar]*(m_coeff[k][j]);
        }
    }
    for (G4int ivar = 0; ivar < fnvar; ++ivar)
    {
        xest[ivar] = m_table[0][ivar]*(1. + m_coeff[k][0]) - xest[ivar]*m_coeff[k][0];
    }
}

G4double BulirschStoerDenseOut::calc_h_opt(G4double h, G4double error, G4int k) const
{
    G4double expo = 1./(m_interval_sequence[k-1]);
    G4double facmin = std::pow(STEPFAC3, expo);
    G4double fac;
    if (error == 0.0)
        fac = 1./facmin;
    else
    {
        fac = STEPFAC2 / std::pow(error / STEPFAC1, expo);
        fac = std::max(facmin/STEPFAC4, std::min(1./facmin, fac));
    }
    return h*fac;
}

G4bool BulirschStoerDenseOut::in_convergence_window(G4int k) const
{
    if( (k == m_current_k_opt-1) && !m_last_step_rejected )
        return true; // decrease order only if last step was not rejected
    return ( (k == m_current_k_opt) || (k == m_current_k_opt+1) );
}

G4bool BulirschStoerDenseOut::should_reject(G4double error, G4int k) const
{
    if(k == m_current_k_opt-1)
    {
        const G4double d = m_interval_sequence[m_current_k_opt] * m_interval_sequence[m_current_k_opt+1] /
            (m_interval_sequence[0]*m_interval_sequence[0]);
        //step will fail, criterion 17.3.17 in NR
        return ( error > d*d );
    }
    else if(k == m_current_k_opt)
    {
        const G4double d = m_interval_sequence[m_current_k_opt+1] / m_interval_sequence[0];
        return (error > d*d);
    }
    else
        return error > 1.0;
}




void BulirschStoerDenseOut::extrapolate_dense_out(G4int k, G4double table[][G4FieldTrack::ncompSVEC],
const G4double coeff[m_k_max+1][m_k_max], G4int order_start_index)
//polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
{
    // result is written into table[0]
    for(G4int j=k; j>1; --j)
    {
        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            table[j-1][ivar] = table[j][ivar]*(1. + coeff[k + order_start_index][j + order_start_index - 1])
                    - table[j-1][ivar]*(coeff[k + order_start_index][j + order_start_index - 1]);
        }
    }
    for (G4int ivar = 0; ivar < fnvar; ++ivar)
    {
        table[0][ivar] = table[1][ivar]*(1. + coeff[k + order_start_index][order_start_index])
                - table[0][ivar]*(coeff[k + order_start_index][order_start_index]);
    }
}




G4double BulirschStoerDenseOut::prepare_dense_output(G4int k, const G4double dxdt_start[] ,
                              const G4double x_end[],G4double dt )
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
    for(G4int j = 0; j<=k; ++j)
    {
        /* not working with boost units... */
        const G4double d = G4double(m_interval_sequence[j]) / (2.* dt);
        G4double f = 1.0; //factor 1/2 here because our interpolation interval has length 2 !!!
        for(G4int kappa = 0 ; kappa <= 2*j+1 ; ++kappa)
        {
            calculate_finite_difference(j, kappa, f , dxdt_start);
            f *= d;
        }

        if(j > 0)
            extrapolate_dense_out(j, m_mp_states, m_coeff);
    }

    G4double d = dt/2;

    // extrapolate finite differences
    for(G4int kappa = 0 ; kappa<=2*k+1 ; ++kappa)
    {
        for(G4int j=1; j <= (k-kappa/2); ++j)
            extrapolate_dense_out(j, m_diffs[kappa], m_coeff, kappa/2);

        // extrapolation results are now stored in m_diffs[kappa][0]

        // divide kappa-th derivative by kappa because we need these terms for dense output interpolation
        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            m_diffs[kappa][0][ivar] *= d;
        }
        d *= dt/(2*(kappa+2));
    }

    // dense output coefficients a_0 is stored in m_mp_states[0], a_i for i = 1...2k are stored in m_diffs[i-1][0]

    // the error is just the highest order coefficient of the interpolation polynomial
    // this is because we use only the midpoint theta=0 as support for the interpolation (remember that theta = -1 .. 1)

    G4double error = 0.0;
    if(m_control_interpolation)
    {
        memcpy(m_err, m_diffs[2*k+1][0], sizeof(G4double)*fnvar);
        error = check_error(x_end, m_err, dt, eps_rel);
    }
    return error;
}



void BulirschStoerDenseOut::calculate_finite_difference(G4int j, G4int kappa, G4double fac, const G4double dxdt[])
{
    const G4int m = m_interval_sequence[j]/2-1;
    if(kappa == 0) // no calculation required for 0th derivative of f
    {
        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            m_diffs[0][j][ivar] = m_derivs[j][m][ivar]*fac;
        }
    }
    else
    {
        // calculate the index of m_diffs for this kappa-j-combination
        const G4int j_diffs = j - kappa/2;

        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            m_diffs[kappa][j_diffs][ivar] = m_derivs[j][m+kappa][ivar]*fac;
        }

        G4double sign = -1.0;
        G4int c = 1;
        //computes the j-th order finite difference for the kappa-th derivative of f at t+dt/2 using function evaluations stored in m_derivs
        for(G4int i = m+kappa-2; i >= m-kappa; i -= 2)
        {
            if( i >= 0 )
            {
                for (G4int ivar = 0; ivar < fnvar; ++ivar)
                {
                    m_diffs[kappa][j_diffs][ivar] +=
                            m_derivs[j][i][ivar] * sign * fac * binomial_coefficient( kappa , c );
                }

            }
            else
            {
                for (G4int ivar = 0; ivar < fnvar; ++ivar)
                {
                    m_diffs[kappa][j_diffs][ivar] += dxdt[ivar]*sign * fac;
                }
            }
            sign *= -1;
            ++c;
        }
    }
}


void BulirschStoerDenseOut::do_interpolation(G4double t, G4double out[]) const
{
    // interpolation polynomial is defined for theta = -1 ... 1
    // m_k_final is the number of order-iterations done for the last step - it governs the order of the interpolation polynomial
    const G4double theta = 2 * (t - m_t_last) / (m_t - m_t_last) - 1;
    // we use only values at interval center, that is theta=0, for interpolation
    // our interpolation polynomial is thus of order 2k+2, hence we have 2k+3 terms

    memcpy(out,m_mp_states[0],sizeof(G4double)*fnvar);
    // add remaining terms: x += a_1 theta + a2 theta^2 + ... + a_{2k} theta^{2k}
    G4double theta_pow = theta;
    for(G4int i=0; i<=2*m_k_final+1; ++i)
    {
        for (G4int ivar = 0; ivar < fnvar; ++ivar)
        {
            out[ivar] += m_diffs[i][0][ivar]*theta_pow;
        }
        theta_pow *= theta;
    }
}





