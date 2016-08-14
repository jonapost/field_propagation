#include "G4VIntegrationDriver.hh"

G4VIntegrationDriver::G4VIntegrationDriver(G4double hminimum,
                       G4EquationOfMotion* equation,
                       G4int numberOfComponents,
                       G4int  VerboseLevel):
    fpStepper(nullptr),
    fpEquation(equation),
    fMinimumStep(hminimum),
    fNoIntegrationVariables(numberOfComponents),
    fVerboseLevel(VerboseLevel),
    fInterpolationInterval(std::make_pair(0.0,0.0))
{
}

G4VIntegrationDriver::G4VIntegrationDriver(G4double hminimum,
                       G4MagIntegratorStepper *pItsStepper,
                       G4int numberOfComponents,
                       G4int VerboseLevel):
    fpStepper(pItsStepper),
    fpEquation(nullptr),
    fMinimumStep(hminimum),
    fNoIntegrationVariables(numberOfComponents),
    fVerboseLevel(VerboseLevel),
    fInterpolationInterval(std::make_pair(0.,0.))
{
}

G4VIntegrationDriver::~G4VIntegrationDriver()
{
}

void G4VIntegrationDriver::DoStep(G4FieldTrack& /*track*/, G4double /*hstep*/, G4double /*eps*/)
{
    G4Exception("G4VIntegrationDriver::DoStep","FieldGeom005",
                FatalException, "DoStep is not implemented for this Driver!");
}

void G4VIntegrationDriver::DoInterpolation(G4FieldTrack& /*track*/, G4double /*hstep*/, G4double /*eps*/)
{
    G4Exception("G4VIntegrationDriver::DoInterpolation","FieldGeom005",
                FatalException, "DoInterpolation is not implemented for this Driver!");
}
void G4VIntegrationDriver::Reset()
{

    //reset interploation interval
    fInterpolationInterval.first = fInterpolationInterval.second = 0;
}
