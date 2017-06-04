// The Butcher table of the Higham & Hall 5(4)7 method is:
//
//   0  |
//  2/9 |      2/9
//  1/3 |      1/12	     1/4
//  1/2 |      1/8       0          3/8
//  3/5 |      91/500    -27/100    78/125    8/125
//   1  |      -11/20    27/20      12/5      -36/5    5
//   1  |      1/12      0          27/32     -4/3     125/96    5/48
//-------------------------------------------------------------------------------------------------------------------
//             1/12      0          27/32     -4/3     125/96    5/48    0
//             2/15      0          27/80     -2/15    25/48     1/24    1/10

#include "RK547FEq1.hh"
#include "G4LineSection.hh"
#include "Utils.hh"


RK547FEq1::RK547FEq1(G4EquationOfMotion *EqRhs, G4int integrationVariables)
   : G4MagIntegratorStepper(EqRhs, integrationVariables)
{
}

void RK547FEq1::saveStep(const G4double yInput[],
                         const G4double dydx[],
                         const G4double hstep,
                         const G4double yOutput[])
{

    memcpy (yIn_, yInput, sizeof(G4double) * GetNumberOfStateVariables());
    memcpy (dydx_, dydx, sizeof(G4double) * GetNumberOfStateVariables());
    hstep_ = hstep;
    memcpy (yOut_, yOutput, sizeof(G4double) * GetNumberOfStateVariables());
    memcpy (dydxNext_, ak7_, sizeof(G4double) * GetNumberOfStateVariables());
}

void RK547FEq1::makeStep(const G4double yInput[],
                        const G4double dydx[],
                        const G4double hstep,
                        G4double yOutput[],
                        G4double yError[]) const
{
    //TODO copy time array[7], copy yInput

    const G4double
       b21 = 2./9.,
       b31 = 1./12., b32 = 1./4.,
       b41 = 1./8., b42 = 0., b43 = 3./8.,
       b51 = 91./500., b52 = -27./100., b53 = 78./125., b54 = 8./125.,

       b61 = -11./20., b62 = 27./20., b63 = 12./5.,
           b64 = -36./5., b65 = 5.,

       b71 = 1./12.,    b72 = 0., b73 = 27./32.,
            b74 = -4./3., b75 = 125./96., b76 = 5./48.;

    const G4double
       dc1 = b71 - 2./15.,
       dc2 = b72 - 0.,
       dc3 = b73 - 27./80.,
       dc4 = b74 + 2./15.,
       dc5 = b75 - 25./48.,
       dc6 = b76 - 1./24.,
       dc7 = 0. - 1./10.;

    //RightHandSide(yInput, dydx);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp_[i] = yInput[i] + hstep * b21 * dydx[i];

    //TODO declare RightHandSide - const!
    G4MagIntegratorStepper* super = const_cast<RK547FEq1*>(this);

    super->RightHandSide(yTemp_, ak2_);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp_[i] = yInput[i] + hstep * (b31 * dydx[i] + b32 * ak2_[i]);

    super->RightHandSide(yTemp_, ak3_);
    for(int i = 0;i < GetNumberOfVariables(); ++i)
        yTemp_[i] = yInput[i] + hstep * (b41 * dydx[i] + b42 * ak2_[i] +
                                        b43 * ak3_[i]);

    super->RightHandSide(yTemp_, ak4_);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp_[i] = yInput[i] + hstep * (b51 * dydx[i] + b52 * ak2_[i] +
                                        b53 * ak3_[i] + b54 * ak4_[i]);

    super->RightHandSide(yTemp_, ak5_);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yTemp_[i] = yInput[i] + hstep * (b61 * dydx[i] + b62 * ak2_[i] +
                                        b63 * ak3_[i] + b64 * ak4_[i] +
                                        b65 * ak5_[i]);

    super->RightHandSide(yTemp_, ak6_);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yOutput[i] = yInput[i] + hstep * (b71 * dydx[i] + b72 * ak2_[i] +
                                          b73 * ak3_[i] + b74 * ak4_[i] +
                                          b75 * ak5_[i] + b76 * ak6_[i]);

    super->RightHandSide(yOutput, ak7_);
    for(int i = 0; i < GetNumberOfVariables(); ++i)
        yError[i] = hstep * (dc1 * dydx[i] + dc2 * ak2_[i] + dc3 * ak3_[i] +
                             dc4 * ak4_[i] + dc5 * ak5_[i] + dc6 * ak6_[i] +
                             dc7 * ak7_[i]);

    //TODO use ak7_ as dydx for next step (FSAL property)
}

void RK547FEq1::Stepper(const G4double yInput[],
                        const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[])
{
    makeStep(yInput, dydx, hstep, yOutput, yError);
    saveStep(yInput, dydx, hstep, yOutput);
}

G4double RK547FEq1::DistChord() const
{
    const G4ThreeVector begin = makeVector(yIn_, Type::Position);
    const G4ThreeVector end = makeVector(yOut_, Type::Position);

    makeStep(yIn_, dydx_, hstep_ / 2., yMid_, yErr_);

    const G4ThreeVector mid = makeVector(yMid_, Type::Position);

    return G4LineSection::Distline(mid, begin, end);
}
