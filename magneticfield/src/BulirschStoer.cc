#include "BulirschStoer.hh"

BulirschStoer::BulirschStoer(G4EquationOfMotion* equation,
                             G4int nvar,
                             G4double eps_rel,
                             G4double max_dt):
    fequation(equation),
    fnvar(nvar),
    m_eps_rel(eps_rel),
    m_midpoint(equation,nvar),
    m_last_step_rejected(false), m_first(true),
    m_max_dt(max_dt),
    STEPFAC1( 0.65 ), STEPFAC2( 0.94 ) , STEPFAC3( 0.02 ) , STEPFAC4( 4.0 ) , KFAC1( 0.8 ) , KFAC2( 0.9 )

{
    /* initialize sequence of stage numbers and work */
    for(G4int i = 0; i < m_k_max+1; ++i)
    {
        m_interval_sequence[i] = 2 * (i+1);
        if(i == 0)
            m_cost[i] = m_interval_sequence[i];
        else
            m_cost[i] = m_cost[i-1] + m_interval_sequence[i];
        for(G4int k = 0; k < i; ++k)
        {
            const G4double r = static_cast< G4double >( m_interval_sequence[i] ) / static_cast< G4double >( m_interval_sequence[k] );
            m_coeff[i][k] = 1.0 / (r*r - 1.0); // coefficients for extrapolation
        }

        // crude estimate of optimal order

        m_current_k_opt = 4;
        // no calculation because log10 might not exist for value_type!

        //const G4double logfact =  -log10( std::max(eps_rel , 1.0e-12 ) ) * 0.6 + 0.5;
        //m_current_k_opt = std::max(1., std::min(static_cast<G4double>(m_k_max-1), logfact));

    }
}

step_result BulirschStoer::try_step(const G4double in[] , const G4double dxdt[] ,
                               G4double& t , G4double out[] , G4double& dt)
{
    if(m_max_dt < dt)
    {
        // given step size is bigger then max_dt
        // set limit and return fail
        dt = m_max_dt;
        return fail;
    }

    if(dt != m_dt_last)
    {
        reset(); // step size changed from outside -> reset
    }

    bool reject =  true;

    G4double new_h = dt;

    /* m_current_k_opt is the estimated current optimal stage number */
    for(G4int k = 0; k <= m_current_k_opt+1; ++k)
    {
        /* the stage counts are stored in m_interval_sequence */
        m_midpoint.set_steps(m_interval_sequence[k]);
        if( k == 0 )
        {
            m_midpoint.do_step(in, dxdt, out, dt);
            /* the first step, nothing more to do */
        }
        else
        {
            m_midpoint.do_step(in, dxdt, m_table[k-1], dt);
            extrapolate(k, out);
            // get error estimate
            for (G4int i = 0; i < fnvar; ++i){
                m_err[i] = out[i] - m_table[0][i];
            }
            const G4double error = check_error(out, m_err, dt, m_eps_rel);
            h_opt[k] = calc_h_opt(dt, error, k);
            work[k] = static_cast<G4double>( m_cost[k] ) / h_opt[k];

            if( (k == m_current_k_opt-1) || m_first )
            { // convergence before k_opt ?
                if( error < 1.0 )
                {
                    //convergence
                    reject = false;
                    if( (work[k] < KFAC2*work[k-1]) || (m_current_k_opt <= 2) )
                    {
                        // leave order as is (except we were in first round)
                        m_current_k_opt = std::min(m_k_max-1 , std::max( 2 , k+1 ) );
                        new_h = h_opt[k];
                        new_h *= static_cast<G4double>( m_cost[k+1] ) / static_cast<G4double>( m_cost[k] );
                    }
                    else
                    {
                       m_current_k_opt = std::min(m_k_max-1, std::max(2, k ) );
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
                        m_current_k_opt = std::max( 2 , m_current_k_opt-1 );
                        new_h = h_opt[m_current_k_opt];
                    }
                    else if( (work[k] < KFAC2*work[k-1]) && !m_last_step_rejected )
                    {
                        m_current_k_opt = std::min(m_k_max-1, m_current_k_opt+1);
                        new_h = h_opt[k];
                        new_h *= static_cast<G4double>(m_cost[m_current_k_opt])/static_cast<G4double>(m_cost[k]);
                        //new_h *= (m_cost[m_current_k_opt])/(m_cost[k]);
                        //were there a error!?
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
            if(k == m_current_k_opt+1)
            { // convergence at k_opt+1 ?
                if(error < 1.0)
                {   //convergence
                    reject = false;
                    if(work[k-2] < KFAC2*work[k-1])
                    {
                        m_current_k_opt = std::max(2 , m_current_k_opt-1);
                    }
                    if((work[k] < KFAC2*work[m_current_k_opt]) && !m_last_step_rejected)
                    {
                        m_current_k_opt = std::min(m_k_max-1 , k);
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
        t += dt;
    }

    if(!m_last_step_rejected || new_h < dt)
    {
        // limit step size
        new_h = std::min(m_max_dt, new_h);
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



void BulirschStoer::reset()
{
    m_first = true;
    m_last_step_rejected = false;
}

void BulirschStoer::extrapolate(size_t k , G4double xest[])
/* polynomial extrapolation, see http://www.nr.com/webnotes/nr3web21.pdf
 * uses the obtained intermediate results to extrapolate to dt->0
 */
{
    for(G4int j=k-1 ; j>0 ; --j)
    {
        for (G4int i = 0; i < fnvar; ++i){
            m_table[j-1][i] = m_table[j][i]*(1. + m_coeff[k][j]) - m_table[j-1][i]*m_coeff[k][j];
        }
    }
    for (G4int i = 0; i < fnvar; ++i){
        xest[i] = m_table[0][i]*(1. + m_coeff[k][0]) - xest[i]*m_coeff[k][0];
    }
}

G4double BulirschStoer::calc_h_opt(G4double h , G4double error , size_t k ) const
/* calculates the optimal step size for a given error and stage number */
{
    G4double expo( 1.0/(2*k+1) );
    G4double facmin = std::pow( STEPFAC3 , expo );

    G4double fac;
    if (error == 0.0)
        fac=1.0/facmin;
    else
    {
        fac = STEPFAC2 / std::pow( error / STEPFAC1 , expo );
        fac = std::max(facmin/STEPFAC4, std::min(1.0/facmin, fac));
    }

    return h*fac;
}

//why is not used!!??
G4bool BulirschStoer::set_k_opt( size_t k, G4double &dt )
/* calculates the optimal stage number */
{
    if( k == 1 )
    {
        m_current_k_opt = 2;
        return true;
    }
    if( (work[k-1] < KFAC1*work[k]) || (k == m_k_max) )
    {   // order decrease
        m_current_k_opt = k-1;
        dt = h_opt[ m_current_k_opt ];
        return true;
    }

    else if( (work[k] < KFAC2*work[k-1]) || m_last_step_rejected || (k == m_k_max-1) )
    {   // same order - also do this if last step got rejected
        m_current_k_opt = k;
        dt = h_opt[ m_current_k_opt ];
        return true;
    }

    else
    {   // order increase - only if last step was not rejected
        m_current_k_opt = k+1;
        dt = h_opt[ m_current_k_opt-1 ] * m_cost[ m_current_k_opt ] / m_cost[ m_current_k_opt-1 ] ;
        return true;
    }

}

G4bool BulirschStoer::in_convergence_window(G4int k ) const
{
    if( (k == m_current_k_opt-1) && !m_last_step_rejected )
        return true; // decrease stepsize only if last step was not rejected
    return ( (k == m_current_k_opt) || (k == m_current_k_opt+1) );
}


G4bool BulirschStoer::should_reject(G4double error , G4int k) const
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
        const G4double d = m_interval_sequence[m_current_k_opt] / m_interval_sequence[0];
        return ( error > d*d );
    }
    else
        return error > 1.0;
}

