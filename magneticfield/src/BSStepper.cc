#include "BSStepper.hh"
#include "G4LineSection.hh"


BSStepper::BSStepper(G4EquationOfMotion *Equation,
                     G4int  numIntegrationVariables,
                     G4int  numStateVariables):
    G4MagIntegratorStepper(Equation,numIntegrationVariables,numStateVariables),
    denseMidpoint(Equation),
    interval_sequence{2,4},
    fcoeff(1./(sqr(G4double(interval_sequence[1])/G4double(interval_sequence[0]))-1.))
{
}

BSStepper::~BSStepper()
{
}


void  BSStepper::Stepper(const G4double y[],
                         const G4double dydx[],
                         G4double hstep,
                         G4double yout[],
                         G4double yError[])
{
    //G4cout<<"BSStepper::Stepper should never been called! \n";
    //throw("wrong call!");

    memcpy(yIn,y,sizeof(G4double)*GetNumberOfVariables());

    denseMidpoint.set_steps(interval_sequence[0]);
    denseMidpoint.do_step(yIn, dydx, yOut[0], hstep, yMid[0], derivs[0]);

    denseMidpoint.set_steps(interval_sequence[1]);
    denseMidpoint.do_step(yIn, dydx, yOut[1], hstep, yMid[1], derivs[1]);

    //extrapolation
    for (G4int i = 0; i < GetNumberOfVariables(); ++i)
    {
        yOut[0][i] = yOut[1][i] + (yOut[1][i] - yOut[0][i])*fcoeff;
        yMid[0][i] = yMid[1][i] + (yMid[1][i] - yMid[0][i])*fcoeff;
    }

    //calc error
    for (G4int i = 0; i < GetNumberOfVariables(); ++i)
    {
        yError[i] = yOut[0][i] - yOut[1][i];
    }

    memcpy(yout,yOut[0],sizeof(G4double)*GetNumberOfVariables());
}


G4double  BSStepper::DistChord() const
{
    //G4cout<<"BSStepper::DistChord should never been called! \n";
    //throw("wrong call!");
    //return -1;
    return  G4LineSection::Distline(G4ThreeVector(yMid[0][0],yMid[0][1],yMid[0][2]),
                                    G4ThreeVector(yIn[0],yIn[1],yIn[2]),
                                    G4ThreeVector(yOut[0][0],yOut[0][1],yOut[0][2]));
}



G4int BSStepper::IntegratorOrder() const
{
    //G4cout<<"BSStepper::IntegratorOrder should never been called! \n";
    //throw("wrong call!");
    //return -1;
    return 2;
}
