#include "G4VIntegrationDriver.hh"

G4VIntegrationDriver::G4VIntegrationDriver(G4double hminimum,
                       G4EquationOfMotion* equation,
                       G4int numberOfComponents,
                       G4int  statisticsVerbosity):
    // fStepper(nullptr),
    fequation(equation),
    fMinimumStep(hminimum),
    fnvar(numberOfComponents),
    fverb(statisticsVerbosity)
{
}

G4VIntegrationDriver::G4VIntegrationDriver(G4double hminimum,
                                           G4MagIntegratorStepper *pStepper,
                                           G4int numberOfComponents,
                                           G4int statisticsVerbosity):
    // fStepper(pStepper),
    fequation(pStepper->GetEquationOfMotion()),
    fMinimumStep(hminimum),
    fnvar(numberOfComponents),
    fverb(statisticsVerbosity)
{

}

G4VIntegrationDriver::~G4VIntegrationDriver(){

}
