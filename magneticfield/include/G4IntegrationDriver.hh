#ifndef G4MAGINT_DRIVER_HH
#define G4MAGINT_DRIVER_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4VIntegrationDriver.hh"
#include "G4MagIntegratorStepper.hh"

template <class T>
class G4MagInt_Driver: public G4VIntegrationDriver {
public:
    G4MagInt_Driver(
        G4double hminimum,
        T* stepper,
        G4int numberOfComponents = 6,
        G4int statisticsVerbosity = 1);

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

    virtual void SetVerboseLevel(G4int newLevel) override;
    virtual G4int GetVerboseLevel() const override;

    virtual G4EquationOfMotion* GetEquationOfMotion() override;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    // Accessors.
    G4double GetMinimumStep() const;
    G4double GetSafety() const;
    G4double GetPshrnk() const;
    G4double GetPgrow() const;
    G4double GetErrcon() const;



    void RenewStepperAndAdjust(G4MagIntegratorStepper *pItsStepper);
        // Sets a new stepper pItsStepper for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

    void ReSetParameters(G4double new_safety= 0.9 );
        //  i) sets the exponents (pgrow & pshrnk), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "errcon" according to the above values.

    void SetSafety(G4double valS);
    void SetPshrnk(G4double valPs);
    void SetPgrow (G4double valPg);
    void SetErrcon(G4double valEc);
        // When setting safety or pgrow, errcon will be set to a 
        // compatible value.

    G4double ComputeAndSetErrcon();

    const T* GetStepper() const;
    T* GetStepper();

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


     G4int GetMaxNoSteps() const;
     void SetMaxNoSteps( G4int val);
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).

     void SetMinimumStep(G4double newval);

     G4double GetSmallestFraction() const;
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

     // The (default) maximum number of steps is Base
     // divided by the order of Stepper
     static constexpr G4int  fMaxStepBase = 250;

     // Parameters used to grow and shrink trial stepsize.
     G4double safety;
     G4double pshrnk;   //  exponent for shrinking
     G4double pgrow;    //  exponent for growth
     G4double errcon;

     // Maximum stepsize increase/decrease factors.
     static constexpr G4double max_stepping_increase = 5;
     static constexpr G4double max_stepping_decrease = 0.1;

     G4int    fStatisticsVerboseLevel;

     T* pIntStepper;

     // Step Statistics
     unsigned long fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps;
     unsigned long fNoCalls; 

     G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues
};

#include "G4MagIntegratorDriver.icc"

#endif /* G4MagInt_Driver_Def */
