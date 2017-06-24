#ifndef G4MAGINT_DRIVER_HH
#define G4MAGINT_DRIVER_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4VIntegrationDriver.hh"
#include "G4MagIntegratorStepper.hh"


class G4MagInt_Driver: public G4VIntegrationDriver {
public:
    G4MagInt_Driver(
        G4double hminimum,
        G4MagIntegratorStepper* pItsStepper,
        G4int numberOfComponents = 6,
        G4int statisticsVerbosity = 1);

    virtual ~G4MagInt_Driver() override;

    G4MagInt_Driver(const G4MagInt_Driver&) = delete;
    const G4MagInt_Driver& operator =(const G4MagInt_Driver&) = delete;

    // Integrates ODE from current s (s=s0) to s=s0+h with accuracy eps.
    // On output track is replaced by value at end of interval.
    // The concept is similar to the odeint routine from NRC p.721-722.
    virtual G4bool AccurateAdvance(
        G4FieldTrack& track,
        G4double hstep,
        G4double eps,            // Requested y_err/hstep
        G4double hinitial = 0) override;  // Suggested 1st interval


    // QuickAdvance just tries one Step - it does not ensure accuracy.
    virtual G4bool QuickAdvance(
        G4FieldTrack& fieldTrack,      // INOUT
        const G4double dydx[],
        G4double hstep,       // IN
        G4double& dchord_step,
        G4double& dyerr) override;

    virtual void GetDerivatives(
         const G4FieldTrack &track,
         G4double dydx[]) const override;

    // Taking the last step's normalised error, calculate
    // a step size for the next step.
    // Do not limit the next step's size within a factor of the
    // current one.
    virtual G4double ComputeNewStepSize(
        G4double errMaxNorm,    // normalised error
        G4double hstepCurrent) override; // current step size

    inline virtual void SetVerboseLevel(G4int newLevel) override;
    inline virtual G4int GetVerboseLevel() const override;

    virtual G4EquationOfMotion* GetEquationOfMotion() override
    {
        return pIntStepper->GetEquationOfMotion();
    }

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override
    {
        pIntStepper->SetEquationOfMotion(equation);
    }


     inline G4double GetHmin() const;
     inline G4double Hmin() const;     // Obsolete
     inline G4double GetSafety() const;
     inline G4double GetPshrnk() const;
     inline G4double GetPgrow() const;
     inline G4double GetErrcon() const;

        // Accessors.

     inline void RenewStepperAndAdjust(G4MagIntegratorStepper *pItsStepper);
        // Sets a new stepper pItsStepper for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(G4double new_safety= 0.9 );
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

    // This takes one Step that is as large as possible while
    // satisfying the accuracy criterion of:
    // yerr < eps * |y_end-y_start|
    void OneGoodStep(
        G4double ystart[],
        const G4double dydx[],
        G4double& curveLength,
        G4double htry,
        G4double eps,      //  memb variables ?
        G4double& hdid,
        G4double& hnext);


     inline G4int GetMaxNoSteps() const;
     inline void SetMaxNoSteps( G4int val);
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).

     inline void SetHmin(G4double newval);

     inline G4double GetSmallestFraction() const; 
     void     SetSmallestFraction( G4double val ); 

private:
     G4double ShrinkStepSize(G4double h, G4double error) const;
     G4double GrowStepSize(G4double h, G4double error) const;

     // Minimum Step allowed in a Step (in absolute units)
     G4double  fMinimumStep;

     // Smallest fraction of (existing) curve length - in relative units
     // below this fraction the current step will be the last
     // Expected range 1e-12 to 5e-15;
     G4double  fSmallestFraction;

     G4int   fMaxNoSteps;
     static const G4int  fMaxStepBase;  

     // Parameters used to grow and shrink trial stepsize.
     G4double safety;
     G4double pshrnk;   //  exponent for shrinking
     G4double pgrow;    //  exponent for growth
     G4double errcon;

     // Maximum stepsize increase/decrease factors.
     static const G4double max_stepping_increase;
     static const G4double max_stepping_decrease;

     G4int    fStatisticsVerboseLevel;

     G4MagIntegratorStepper *pIntStepper;

     // Step Statistics
     unsigned long fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps;
     unsigned long fNoCalls; 



     G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues
};

#include "G4MagIntegratorDriver.icc"

#endif /* G4MagInt_Driver_Def */
