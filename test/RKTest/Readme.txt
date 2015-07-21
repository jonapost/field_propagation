//RKTest is a basic test for the RK steppers.
//It is the combination of various previous tests written 
	as a part of this project.
~~~~~~~~~~~~~~
	 MAN
~~~~~~~~~~~~~~

**********************
*	INTERFACE :
*
	BASIC
	--------------------------------------------------------
	1. Create an object of the class RKTest
		//Code------------
		RKTest mytest:
		//End Code
	2. Call a member function pertaining to a certain test 
		//Code------------
		myTest.(function_Name_with_proper_syntax)
		//End Code--------

	------DONE-------

	DETAILS
	---------------------------------------------------------
	1. Test A stepper using a fixed step size in a Uniform 
		Magnetic Field
		
		//Code------------
		#include "RKTest.hh" //RKTest-suite
		#include "G4CashKarpRKF45.hh" //Stepper to test

		RKTest mytest;
		int columns[] = {1,1,0,0,0,0}
		mytest.testSteppersFixedUMF<G4CashKarpRKF45>(columns);
		//End Code--------

		Here we call the function testSteppersFixedUMF( ) with
		its default parameters and pass the columns array.
		This array specifies the coordinates we want to see in
		the output (with the standard 6-cooridnate convention 
		of Geant4. For each coordinate, we have 3 columns : 
		yOut[i], yErr[i], and yOut[i]-yOut[i].
		The function declaration reads :

		template < class STEPPER >
		testSteppersFixedUMF(int columns[6],
                         	 G4double factor = 1.
                             G4double step_len = 25.0*CLHEP::mm,
                             int no_of_steps = 100))

        factor specifies the field strength; pass factor=10 for
        a 10 times stronger field.

    2. Test a G4 Stepper's performance
    	
    	//Code - continued
    	G4double myFieldFac = 10.
       	myTest.testAnyG4Stepper<G4CashKarpRKF45>("umf", myFieldFac);
		//End Code--------

		The string "umf" specifies that we want to use Uniform
		Magnetic Field. One can also use "qmf".

	3. Test an FSAL Stepper (using FSALDriver)
		Here the code is necessarily the same with the only
		difference being you should instantiate using only
		an FSAL stepper.

		//Code 
		myTest.testFSALStepper<FDormandPrince745>("qmf",0.1);
		//End Code

	4. Test interpolation of a capable stepper

		//Code -------
		myTest.testStepperInterpolant<VernerRK78, G4CashKarpRKF45>(columns);
		//End Code ---------

		Here the first class is the one that is under test and the 
		other one is the reference stepper class whose reduced 
		step_len steps would be used as a reference to compare the
		accuracy of interpolation.

		The declaration :

	    template<class STEPPER, class REF_STEPPER>
	    void testStepperInterpolant(int columns[6], 
	                                std::string field_code = "umf", 
	                                G4double factor = 1., 
	                                G4double step_len_pi_divisor = 6.0, 
	                                G4double maxAngle = 2.0*CLHEP::pi)












	


