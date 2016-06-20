#include "ModifiedMidpointDenseOut.hh"
#include <cstring>

ModifiedMidpointDenseOut::ModifiedMidpointDenseOut(G4EquationOfMotion* Equation,
                                                   unsigned int nvar,
                                                   unsigned int steps):

    Equation_Rhs(Equation),
    fnvar(nvar),
    fsteps(steps)

{
}

ModifiedMidpointDenseOut::~ModifiedMidpointDenseOut()
{
}


void ModifiedMidpointDenseOut::do_step(const G4double xIn[], const G4double dxdtIn[],
                                       G4double xOut[], G4double dt,
                                       G4double x_mp[], G4double* derivs[])
{

    const G4double h = dt / fsteps;
    const G4double h2 = 2 * h;

    // x1 = xIn + h*dxdt
    for (G4int i = 0; i < fnvar; ++i){
        x1[i] = xIn[i] + h*dxdtIn[i];
    }

    // result of first step already gives approximation at the center of the interval
    if(fsteps == 2)
        memcpy(x_mp,x1,sizeof(G4double)*fnvar);

    Equation_Rhs->RightHandSide(x1, derivs[0]);

    memcpy(x0,xIn,sizeof(G4double)*fnvar);

    // general step
    //tmp = x1; x1 = x0 + h2*dxdt; x0 = tmp
    for (G4int i = 1; i < fsteps; ++i){
        memcpy(tmp,x1,sizeof(G4double)*fnvar);
        for (G4int j = 0; j < fnvar; ++j){
            x1[j] = x0[j] + h2*dxdt[j];
        }
        memcpy(x0,tmp,sizeof(G4double)*fnvar);

        // save approximation at the center of the interval
        if(i == fsteps/2 - 1 )
            memcpy(x_mp,x1,sizeof(G4double)*fnvar);

        Equation_Rhs->RightHandSide(x1, derivs[i]);
    }

    // last step
    // xOut = 0.5*(x0 + x1 + h*dxdt)
    for (G4int i = 0; i < fnvar; ++i){
        xOut[i] = 0.5*(x0[i] + x1[i] + h*derivs[fsteps-1][i]);
    }
}

