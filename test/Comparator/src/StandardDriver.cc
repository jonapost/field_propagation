#include "StandardDriver.hh"
#include "DriverUtils.hh"

#include "G4FieldTrack.hh"

static const int NCOMP = G4FieldTrack::ncompSVEC;

StandardDriver::StandardDriver(
    G4double minimumStep, G4MagIntegratorStepper* stepper):
    fstepper(stepper),
    ftrials(0),
    fmaxtrials(100),
    fsafety(0.9),
    fminimumStep(minimumStep),
    fpShrink(-1.0 / fstepper->IntegratorOrder()),
    fpGrow(-1.0 / (1.0 + fstepper->IntegratorOrder())),
    fmaxStepIncrease(5),
    fErrcon(std::pow(fmaxStepIncrease / fsafety, 1.0 / fpGrow))
{
}

void StandardDriver::OneGoodStep(G4double y[],
                                 const G4double dydx[],
                                 G4double& trackLength,
                                 G4double htry,
                                 G4double errorTolerance,
                                 G4double& hdid,
                                 G4double& hnext)
{
    G4double yerr[NCOMP], ytemp[NCOMP];
    memset(yerr, 0, sizeof(G4double)* NCOMP);

    // Set stepsize to the initial trial value
    G4double hstep = htry;

    G4double error;
    G4int iter = 0;
    for (; iter < fmaxtrials; ++iter) {
        ++ftrials;

        hstep = std::max(hstep, fminimumStep);
        fstepper->Stepper(y, dydx, hstep, ytemp, yerr);

        error = relativeError(y, yerr, hstep) / errorTolerance;

        // Step succeeded
        if (error <= 1.0) {
            break;
        }

        // Step failed. Compute the size of retrial step
        hstep = shrinkStep(error, hstep);
    }

    if (iter == fmaxtrials) {
        G4Exception("GustafssonDriver::OneGoodStep", "007", FatalException,
                    "Driver could not converge!");
    }

    // Compute size of next step
    hnext = growStep(error, hstep);

    trackLength += (hdid = hstep);

    memcpy(y, ytemp, sizeof(G4double) * NCOMP);
}

G4double StandardDriver::shrinkStep(G4double error, G4double hstep)
{
    G4double htemp = fsafety * hstep * std::pow(error, fpShrink);
    return htemp >= 0.1 * hstep ? htemp : hstep * 0.1;
}

G4double StandardDriver::growStep(G4double error, G4double hstep)
{
    return error > fErrcon ? fsafety * hstep * std::pow(error, fpGrow) :
                             fmaxStepIncrease * hstep;
}
