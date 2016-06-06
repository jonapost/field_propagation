#ifndef ModifiedMidpoint_HH
#define ModifiedMidpoint_HH

#include <cstring>
#include <vector>
#include <array>

#include "G4FieldTrack.hh"

typedef std::array<double,G4FieldTrack::ncompSVEC> state_type;
typedef std::vector<state_type> state_vector_type;

class ModifiedMidpoint
{
public :
    ModifiedMidpoint(unsigned int steps = 2, unsigned int nvar = 6):fsteps(steps),fnvar(nvar)
    {
        fnvar = std::min(fnvar,state_type.size());
    }
    ~ModifiedMidpoint(){
    }

    template<class System>
    void do_step_impl(System system, const state_type& yIn, const state_type& dydxIn, state_type& yOut, double hstep){

        const double h = hstep / fsteps;
        const double h2 = 2 * h;

        // y1 = yIn + h*dydx
        for (int i = 0; i < fnvar; ++i){
            y1[i] = yIn + h*dydxIn[i];
        }

        system(y1, dydx);

        y0 = yIn;

        // general step
        //tmp = y1; y1 = y0 + h2*dydx; y0 = tmp
        for (int i = 0; i < fsteps; ++i){
            tmp = y1;
            for (int j = 0; j < fnvar; ++j){
                y1[j] = y0[j] + h2*dydx[j];
            }
            y0 = tmp;

            th += h;
            system(y1 , dydx);
        }


        // last step
        // yOut = 0.5*(y0 + y1 + h*dydx)
        for (int i = 0; i < fnvar; ++i){
            yOut[i] = 0.5*(y0[i] + y1[i] + h*dydx[i]);
        }

    }


    void set_steps(unsigned int steps)
    {   fsteps = steps; }


    unsigned int steps() const
    {   return fsteps; }

private:

    unsigned int fsteps;
    unsigned int fnvar;

    state_type y0;
    state_type y1;
    state_type dydx;
    state_type tmp;

};




/* Modified midpoint which stores derivatives and state at dt/2
 * in some external storage for later usage in dense output calculation
 * This Stepper is for use in Bulirsch Stoer only.
 */

class ModifiedMidpointDenseOut
{
public :
    ModifiedMidpointDenseOut(unsigned int steps = 2 , unsigned int nvar = 6)
    : fsteps(steps), fnvar(nvar)
    {
        fnvar = std::min(fnvar,state_type.size());
    }
    ~ModifiedMidpoint(){
    }

    /*
     * performs a modified midpoint step with fsteps intermediate points
     * stores approximation for x(t+dt/2) in x_mp and all evaluated function results in derivs
     *
     */

    template< class System>
    void do_step(System system, const state_type& yIn, const state_type& dydxIn,
                 state_type& yOut,double hstep, state_type y_mp , state_vector_type& derivs)
    {

        const double h = hstep / fsteps;
        const double h2 = 2 * h;

        // y1 = yIn + h*dydx
        for (int i = 0; i < fnvar; ++i){
            y1[i] = yIn + h*dydxIn[i];
        }

        // result of first step already gives approximation at the center of the interval
        if(fsteps == 2)
            y_mp = y1;

        system(y1, derivs[0]);

        y0 = yIn;

        // general step
        //tmp = m_x1; m_x1 = m_x0 + h2*m_dxdt; m_x0 = tmp
        for (int i = 1; i <= fsteps; ++i){
            tmp = y1;
            for (int j = 0; j < fnvar; ++j){
                y1[j] = y0[j] + h2*dydx[j];
            }
            y0 = tmp;

            // save approximation at the center of the interval
            if(i == ((fsteps/2) - 1) )
                y_mp = y1;

            system(y1, derivs[i]);
        }

        // last step
        // yOut = 0.5*(y0 + y1 + h*dydx)
        for (int i = 0; i < fnvar; ++i){
            yOut[i] = 0.5*(y0[i] + y1[i] + h*derivs[fsteps-1]);
        }

     }


    void set_steps( unsigned short steps )
    {   fsteps = steps; }


    unsigned short steps() const
    {   return fsteps; }

private:

    unsigned int fsteps;
    unsigned int fnvar;

    state_type y0;
    state_type y1;
    state_type dydx;
    state_type tmp;

};

#endif
