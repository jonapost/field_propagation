#ifndef G4INTEGRATION_DRIVER_HH
#define G4INTEGRATION_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4IntegrationObserver.hh"

class G4IntegrationDriver: public G4VIntegrationDriver {
public:
    G4IntegrationDriver(G4double hminimum,
                            G4MagIntegratorStepper* pItsStepper,
                            G4int numberOfComponents = 6,
                            G4int statisticsVerbosity = 1);

    virtual ~G4IntegrationDriver() override;

    G4IntegrationDriver(const G4IntegrationDriver&) = delete;
    const G4IntegrationDriver& operator = (const G4IntegrationDriver&) = delete;

    // Drivers for integrator (Runge-Kutta) with stepsize control.
    // Integrates ODE starting values y_current
    // from current s (s=s0) to s=s0+h with accuracy eps.
    // On output ystart is replaced by value at end of interval.
    // The concept is similar to the odeint routine from NRC p.721-722.
    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps,            // Requested y_err/hstep
                                   G4double hinitial = 0) override;  // Suggested 1st interval

    // QuickAdvance just tries one Step - it does not ensure accuracy.
    virtual G4bool QuickAdvance(G4FieldTrack& track,      // INOUT
                                const G4double dydx[],
                                G4double hstep,
                                G4double& dchord_step,
                                G4double& dyerr) override;

    virtual void GetDerivatives(const G4FieldTrack &track,
                                G4double dydx[]) const override;

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;
    virtual G4EquationOfMotion* GetEquationOfMotion() override;

    // Taking the last step's normalised error, calculate
    // a step size for the next step.
    // Do not limit the next step's size within a factor of the
    // current one.
    virtual G4double ComputeNewStepSize(G4double errMaxNorm,  // normalised error
                                        G4double hstepCurrent) override; // current step size

    virtual void SetVerboseLevel(G4int level) override;
    virtual G4int GetVerboseLevel() const override;


    G4bool QuickAdvance(G4FieldTrack& y_posvel,        // INOUT
                        G4double dydx[],
                        G4double hstep,           // IN
                        G4double& dchord_step,
                        G4double& dyerr_pos_sq,
                        G4double& dyerr_mom_rel_sq );
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )
       // TODO remove it?



     inline G4double GetHmin() const;
     inline G4double Hmin() const;     // Obsolete
     inline G4double GetSafety() const;
     inline G4double GetPshrnk() const;
     inline G4double GetPgrow() const;
     inline G4double GetErrcon() const;
     inline G4int GetNoTotalSteps() const;	//Only for debug purposes


     inline void RenewStepperAndAdjust(G4MagIntegratorStepper *pItsStepper);
        // Sets a new stepper pItsStepper for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(G4double new_safety = 0.9);
        //  i) sets the exponents (pgrow & pshrnk), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "errcon" according to the above values.

     inline void SetSafety(G4double valS);
     inline void SetPshrnk(G4double valPs);
     inline void SetPgrow (G4double valPg);
     inline void SetErrcon(G4double valEc);
        // When setting safety or pgrow, errcon will be set to a 
        // compatible value.

     inline G4double ComputeAndSetErrcon();

     inline const G4MagIntegratorStepper* GetStepper() const;
     inline G4MagIntegratorStepper* GetStepper();

     void  OneGoodStep(G4double  ystart[], // Like old RKF45step()
                       const G4double dydx[],
                       G4double& x,
                       G4double htry,
                       G4double  eps,      //  memb variables ?
                       G4double& hdid,
                       G4double& hnext);
        // This takes one Step that is as large as possible while 
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|


     G4double ComputeNewStepSize_WithinLimits(
                          G4double  errMaxNorm,    // normalised error
                          G4double  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Limit the next step's size within a range around the current one.

     inline G4int    GetMaxNoSteps() const;
     inline void     SetMaxNoSteps( G4int val); 
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).
    
//---------------------------------------------------------------------
//The following has been introduced by [hackabot] for testing purposes only
    inline G4int GetTotalNoStepperCalls() const;
//---------------------------------------------------------------------

public:  // without description
     inline void SetHmin(G4double newval);

     inline G4double GetSmallestFraction() const; 
     void     SetSmallestFraction( G4double val );

protected:  // without description

     void PrintStatisticsReport();
       //  Report on the number of steps, maximum errors etc.

private:
     G4double GrowStep(G4double hstep, G4double error);
     G4double ShrinkStep(G4double hstep, G4double error);
     // ---------------------------------------------------------------
     //  INVARIANTS 

     G4double  fMinimumStep;
        // Minimum Step allowed in a Step (in absolute units)
     G4double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
        // Smallest fraction of (existing) curve length - in relative units
        //  below this fraction the current step will be the last 

     const G4int  fNoIntegrationVariables;  // Number of Variables in integration
     const G4int  fMinNoVars;               // Minimum number for FieldTrack
     const G4int  fNoVars;                  // Full number of variable

     G4int   fMaxNoSteps;
     static const G4int  fMaxStepBase;  

     G4double safety;
     G4double pshrnk;   //  exponent for shrinking
     G4double pgrow;    //  exponent for growth
     G4double errcon;
        // Parameters used to grow and shrink trial stepsize.

     static const G4double max_stepping_increase;
     static const G4double max_stepping_decrease;
        // Maximum stepsize increase/decrease factors.

     G4int    fStatisticsVerboseLevel;

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     G4MagIntegratorStepper *pIntStepper;

     // ---------------------------------------------------------------
     //  STATE

     G4int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps;
     G4double fDyerr_max, fDyerr_mx2;
     G4double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
     G4double fSumH_sm, fSumH_lg; 
        // Step Statistics 

     G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues
    
    //For Test Purposes :-
    G4int TotalNoStepperCalls;

    G4IntegrationObserver fObserver;
};

#include "G4IntegrationDriver.icc"

#endif /* G4IntegrationDriver_Def */
