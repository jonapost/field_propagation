#include "StandardDriver.hh"
#include "G4FieldTrack.hh"

StandardDriver::StandardDriver(
    G4double minimumStep, G4MagIntegratorStepper* stepper):
    fstepper(stepper),
    ftrials(0),
    fsafety(0.9),
    fminimumStep(minimumStep),
    fpShrink(-1.0 / fstepper->IntegratorOrder()),
    fpGrow(-1.0 / (1.0 + fstepper->IntegratorOrder())),
    fmaxStepIncrease(5),
    fErrcon(std::pow(fmaxStepIncrease/fsafety,1.0/fpGrow))
{
}

void StandardDriver::OneGoodStep(G4double y[],
                                 const G4double dydx[],
                                 G4double& trackLength,
                                 G4double htry,
                                 G4double epsRelMax,
                                 G4double& hdid,
                                 G4double& hnext)
{
    const int NCOMP = G4FieldTrack::ncompSVEC;

    G4double yerr[NCOMP], ytemp[NCOMP];
    memset(yerr, 0, sizeof(G4double)* NCOMP);

    G4double hstep = htry; // Set stepsize to the initial trial value


    const G4int MAXTRIALS = 100;
    G4double errmax2;
    G4int iter = 0;
    for (; iter < MAXTRIALS; ++iter) {
        ++ftrials;

        hstep = std::min(hstep, fminimumStep);
        fstepper->Stepper(y, dydx, hstep, ytemp, yerr);

        // square of displacement error
        G4double positionError2 =  sqr(yerr[0]) + sqr(yerr[1]) + sqr(yerr[2]);
        positionError2 /= sqr(hstep);

        // square of momentum vector difference
        G4double momentumError2 =  sqr(yerr[3]) + sqr(yerr[4]) + sqr(yerr[5]);
        G4double momentum2 =  sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
        momentumError2 /= momentum2;

        // Square of maximum error
        errmax2 = std::max(positionError2, momentumError2) / sqr(epsRelMax);

        // Step succeeded.
        if (errmax2 <= 1.0) {
            break;
        }

        // Step failed; compute the size of retrial step
        G4double htemp = fsafety * hstep * std::pow(errmax2, 0.5 * fpShrink);

        if (htemp >= 0.1 * hstep) { // Truncation error too large,
            hstep = htemp;
        } else {  // reduce stepsize, but no more than a factor of 10
            hstep *= 0.1;
        }
    }

    if (iter == MAXTRIALS) {
        G4cerr <<"can not converge! \n";
    }

    // Compute size of next step
    if (errmax2 > sqr(fErrcon)) {
        hnext = fsafety * hstep * std::pow(errmax2, 0.5 * fpGrow);
    } else {
        hnext = fmaxStepIncrease * hstep;
    }

    trackLength += (hdid = hstep);

    memcpy(y, ytemp, sizeof(G4double)* NCOMP);
}
