#ifndef TMagHelicalStepper_hh
#define TMagHelicalStepper_hh

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

template
<class T_Stepper, class T_Equation, int N>
class TMagHelicalStepper : public G4MagIntegratorStepper
{
    const G4double fUnitConstant = 0.299792458*(GeV/(tesla*m));

    public:  // with description

    TMagHelicalStepper(T_Equation *EqRhs)
        : // integrate over 6 variables only !!
            G4MagIntegratorStepper(EqRhs, 6),     
            // position & velocity
            fAngCurve(0.), 
            frCurve(0.), 
            frHelix(0.),
            fEquation_Rhs(EqRhs){;}

    virtual ~TMagHelicalStepper(){;}


    void
        AdvanceHelix( const G4double  yIn[],
                G4ThreeVector   Bfld,    
                G4double  h,
                G4double  yHelix[],
                G4double  yHelix2[] )
        {
            // const G4int    nvar = 6;

            // OLD  const G4double approc_limit = 0.05;
            // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
            // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

            const G4double approc_limit = 0.005;
            G4ThreeVector  Bnorm, B_x_P, vperp, vpar;

            G4double B_d_P;
            G4double B_v_P;
            G4double Theta;
            G4double R_1;
            G4double R_Helix;
            G4double CosT2, SinT2, CosT, SinT;
            G4ThreeVector positionMove, endTangent;

            G4double Bmag = Bfld.mag();
            const G4double *pIn = yIn+3;
            G4ThreeVector initVelocity= G4ThreeVector( pIn[0], pIn[1], pIn[2]);
            G4double      velocityVal = initVelocity.mag();
            G4ThreeVector initTangent = (1.0/velocityVal) * initVelocity;

            R_1=GetInverseCurve(velocityVal,Bmag);

            // for too small magnetic fields there is no curvature
            // (include momentum here) FIXME

            if( (std::fabs(R_1) < 1e-10)||(Bmag<1e-12) )
            {
                LinearStep( yIn, h, yHelix );

                // Store and/or calculate parameters for chord distance

                SetAngCurve(1.);     
                SetCurve(h);
                SetRadHelix(0.);
            }
            else
            {
                Bnorm = (1.0/Bmag)*Bfld;

                // calculate the direction of the force

                B_x_P = Bnorm.cross(initTangent);

                // parallel and perp vectors

                B_d_P = Bnorm.dot(initTangent); // this is the fraction of P parallel to B

                vpar = B_d_P * Bnorm;       // the component parallel      to B
                vperp= initTangent - vpar;  // the component perpendicular to B

                B_v_P  = std::sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B

                // calculate  the stepping angle

                Theta   = R_1 * h; // * B_v_P;

                // Trigonometrix

                if( std::fabs(Theta) > approc_limit )
                {
                    SinT     = std::sin(Theta);
                    CosT     = std::cos(Theta);
                }
                else
                {
                    G4double Theta2 = Theta*Theta;
                    G4double Theta3 = Theta2 * Theta;
                    G4double Theta4 = Theta2 * Theta2;
                    SinT     = Theta - 1.0/6.0 * Theta3;
                    CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;
                }

                // the actual "rotation"

                G4double R = 1.0 / R_1;

                positionMove  = R * ( SinT * vperp + (1-CosT) * B_x_P) + h * vpar;
                endTangent    = CosT * vperp + SinT * B_x_P + vpar;

                // Store the resulting position and tangent

                yHelix[0]   = yIn[0] + positionMove.x(); 
                yHelix[1]   = yIn[1] + positionMove.y(); 
                yHelix[2]   = yIn[2] + positionMove.z();
                yHelix[3] = velocityVal * endTangent.x();
                yHelix[4] = velocityVal * endTangent.y();
                yHelix[5] = velocityVal * endTangent.z();

                // Store 2*h step Helix if exist

                if(yHelix2)
                {
                    SinT2     = 2.0 * SinT * CosT;
                    CosT2     = 1.0 - 2.0 * SinT * SinT;
                    endTangent    = (CosT2 * vperp + SinT2 * B_x_P + vpar);
                    positionMove  = R * ( SinT2 * vperp + (1-CosT2) * B_x_P) + h*2 * vpar;

                    yHelix2[0]   = yIn[0] + positionMove.x(); 
                    yHelix2[1]   = yIn[1] + positionMove.y(); 
                    yHelix2[2]   = yIn[2] + positionMove.z(); 
                    yHelix2[3] = velocityVal * endTangent.x();
                    yHelix2[4] = velocityVal * endTangent.y();
                    yHelix2[5] = velocityVal * endTangent.z();
                }

                // Store and/or calculate parameters for chord distance

                G4double ptan=velocityVal*B_v_P;

                G4double particleCharge = fEquation_Rhs->FCof() / (eplus*c_light); 
                R_Helix =std::abs( ptan/(fUnitConstant  * particleCharge*Bmag));

                SetAngCurve(std::abs(Theta));
                SetCurve(std::abs(R));
                SetRadHelix(R_Helix);
            }
        }



    // The stepper for the Runge Kutta integration.
    // The stepsize is fixed, equal to h.
    // Integrates ODE starting values y[0 to 6]
    // Outputs yout[] and its estimated error yerr[].

    void
        Stepper( const G4double yInput[],  // VIRTUAL for ExactHelix - temporary
                const G4double*,
                G4double hstep,
                G4double yOut[],
                G4double yErr[]  )
        {  
            const G4int nvar = 6;

            G4int i;

            // correction for Richardson Extrapolation.
            // G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );

            G4double      yTemp[7], yIn[7] ;
            G4ThreeVector Bfld_initial, Bfld_midpoint;

            //  Saving yInput because yInput and yOut can be aliases for same array

            for(i=0;i<nvar;i++) { yIn[i]=yInput[i]; }

            G4double h = hstep * 0.5; 

            MagFieldEvaluate(yIn, Bfld_initial) ;      

            // Do two half steps

            DumbStepper(yIn,   Bfld_initial,  h, yTemp);
            MagFieldEvaluate(yTemp, Bfld_midpoint) ;     
            DumbStepper(yTemp, Bfld_midpoint, h, yOut); 

            // Do a full Step

            h = hstep ;
            DumbStepper(yIn, Bfld_initial, h, yTemp);

            // Error estimation

            for(i=0;i<nvar;i++)
            {
                yErr[i] = yOut[i] - yTemp[i] ;
            }

            return;
        }


    // Estimate maximum distance of curved solution and chord ... 
    G4double
        DistChord() const 
        {
            // Check whether h/R >  pi  !!
            // Method DistLine is good only for <  pi

            G4double Ang=GetAngCurve();
            if(Ang<=pi)
            {
                return GetRadHelix()*(1-std::cos(0.5*Ang));
            }
            else
            {
                if(Ang<twopi)
                {
                    return GetRadHelix()*(1+std::cos(0.5*(twopi-Ang)));
                }
                else  // return Diameter of projected circle
                {
                    return 2*GetRadHelix();
                }
            }
        }


    protected:  // with description

    // A linear Step in regions without magnetic field.
    inline void
        LinearStep( const G4double  yIn[],
                G4double  h,
                G4double  yLinear[]) const
        {
            G4double  momentum_val = std::sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]) ;
            G4double  inv_momentum = 1.0 / momentum_val ;
            G4double  yDir[3];
            // G4double  h_div_momentum = 1.0 / momentum_val ;

            for( G4int i = 0; i < 3; i++ ) {
                yDir[i]   = inv_momentum * yIn[i+3];
                yLinear[i]   = yIn[i] + h * yDir[i];
                // yLinear[i]   = yIn[i] + h_div_momentum * yIn[i+3];
                yLinear[i+3] = yIn[i+3];
            }
        }

    void AdvanceHelix( const G4double  yIn[],
            G4ThreeVector   Bfld,
            G4double  h,
            G4double  yHelix[],G4double yHelix2[]=0);    // output 
    // A first order Step along a helix inside the field.

    // Evaluate the field at a certain point.
    inline void
        MagFieldEvaluate(const G4double y[],
                G4ThreeVector& Bfield )   
        {
            G4double B[3];
            GetEquationOfMotion()->  GetFieldValue(y, B);
            Bfield= G4ThreeVector( B[0], B[1], B[2] );
        }

    // Evaluate Inverse of Curvature of Track
    inline G4double
        GetInverseCurve(const G4double Momentum,
                const G4double Bmag)   
        {

            G4double  inv_momentum = 1.0 / Momentum ;
            G4double particleCharge = fEquation_Rhs->FCof() / (CLHEP::eplus*CLHEP::c_light); 
            G4double fCoefficient = -fUnitConstant  * particleCharge*inv_momentum;

            return  fCoefficient*Bmag;
        }

    // Store and use the parameters of track : 
    // Radius of curve, Stepping angle, Radius of projected helix
    inline void  SetAngCurve(const G4double Ang)
    {                                                
        fAngCurve=Ang; 

    }

    inline G4double  GetAngCurve() const 
    {
        return fAngCurve;  

    }

    inline void  SetCurve(const G4double Curve)
    {
        frCurve=Curve;  
    }

    inline G4double  GetCurve() const 
    {
        return frCurve;  

    }

    inline void  SetRadHelix(const G4double Rad)
    {
        frHelix=Rad;  
    }

    inline G4double  GetRadHelix() const 
    {
        return frHelix;  

    }

    protected:  // without description

    // void MagFieldEvaluate( const G4double y[], G4double B[] )   
    //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

    private:

    TMagHelicalStepper(const TMagHelicalStepper&);
    TMagHelicalStepper& operator=(const TMagHelicalStepper&);
    // Private copy constructor and assignment operator.

    static const G4double fUnitConstant;   //  As in G4Mag_EqRhs.hh/cc where it is not used.
    private:

    T_Equation* fEquation_Rhs;

    // Data stored in order to find the chord.
    G4double fAngCurve;
    G4double frCurve;
    G4double frHelix;
    // Data stored in order to find the chord.
    G4ThreeVector yInitial, yMidPoint, yFinal;


};

#include  "TMagHelicalStepper.icc"

#endif  /* TMagHelicalStepper_hh */
