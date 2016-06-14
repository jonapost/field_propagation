//
// classes ModifiedMidpoint and ModifiedMidpointDenseOut
//
// Class description:
//
// This is an independant copy of modified_midpoint.hpp
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

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
    }
    ~ModifiedMidpoint(){
    }

    template<class System>
    void do_step(System system, const state_type& xIn, const state_type& dxdtIn, double t ,state_type& xOut, double dt){

        const double h = dt / fsteps;
        const double h2 = 2 * h;

        double th = t + h;

        // x1 = xIn + h*dxdt
        for (unsigned int i = 0; i < fnvar; ++i){
            x1[i] = xIn[i] + h*dxdtIn[i];
        }

        system(x1, dxdt, th);

        x0 = xIn;

        // general step
        //tmp = x1; x1 = x0 + h2*dxdt; x0 = tmp
        for (unsigned int i = 1; i < fsteps; ++i){
            tmp = x1;
            for (unsigned int j = 0; j < fnvar; ++j){
                x1[j] = x0[j] + h2*dxdt[j];
            }
            th += h;
            x0 = tmp;

            system(x1 , dxdt, th);
        }


        // last step
        // xOut = 0.5*(x0 + x1 + h*dxdt)
        for (unsigned int i = 0; i < fnvar; ++i){
            xOut[i] = 0.5*(x0[i] + x1[i] + h*dxdt[i]);
        }

    }


    void set_steps(unsigned int steps)
    {   fsteps = steps; }


    unsigned int steps() const
    {   return fsteps; }

private:

    unsigned int fsteps;
    unsigned int fnvar;

    state_type x0;
    state_type x1;
    state_type dxdt;
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
    }
    ~ModifiedMidpointDenseOut(){
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
            y1[i] = yIn[i] + h*dydxIn[i];
        }

        // result of first step already gives approximation at the center of the interval
        if(fsteps == 2)
            y_mp = y1;

        system(y1, derivs[0]);

        y0 = yIn;

        // general step
        //tmp = m_x1; m_x1 = m_x0 + h2*m_dxdt; m_x0 = tmp
        for (int i = 1; i < fsteps; ++i){
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
            yOut[i] = 0.5*(y0[i] + y1[i] + h*derivs[fsteps-1][i]);
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
