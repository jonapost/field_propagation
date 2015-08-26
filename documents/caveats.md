# Caveats

In order to use ```error_by_stepper``` you have to compile with the flag 
```#define INTENDED_FOR_ERROR_BY_STEPPER_PROGRAM``` in the header file ```isTracking.hh```, which is a global place to set this flag and/or the ```#define TRACKING``` flag. 

To use the functionality of ```StepTracker``` to record step history, the ```Tracking``` flag must be defined in the file ```isTracking.hh```, and the nystrom-dev versions of ```G4ChordFinder```, ```G4MagIntegratorDriver```, ```G4MagIntegratorStepper```, ```G4Mag_EqRhs``` and ```G4MagEquationOfMotion``` must be used. The first two classes have been modified to call methods of ```StepTracker``` in the appropriate places. The other classes have been modified to contain a member variable, a pointer to the ```StepTracker``` instance. You must make sure that each object in the stack has this pointer initialized to the same ```StepTracker``` instance. In ```testPropagateMagField``` this has to be explicitly done through the respective ```setTracker()``` methods.

The ```IntegratorOrder()``` of Classical RK4, which is derived from ```G4MagErrorStepper```, has been changed to reflect that Richardson extrapolation promotes the stepper up an order. I may or may not be correct about this.

It is only necessary to wrap regular RK steppers within ```MagIntegratorStepper_byArcLength``` if you want to compare performance to a Nystrom stepper. This template class is included only for the purposes of supplying an interface to the ```StepTracker``` instance.

```MuruaRKN6459``` is faulty because of an error in the coefficients. I have inspected them over and over again from the coefficients published in Murua's original paper, but I cannot find a mistake. I conjecture that there is a transcription error in the published coefficients. One obvious way to check would be to plug them all into the equations they are required to satisfy and see if they do in fact satisfy them.

```ErrorComputer``` and ```ErrorComputer2``` are previous attempts to slickly do what ```interpolation_error``` does in a dumb way. If there is a need for efficiency then these approaches might want to be revisited.

```testPropagateMagField``` now takes multiple arguments, and there are supplied default values:

- An integer giving the stepper type
- ```init_data``` followed by 
    - step_length, 
    - number of calls to ComputeStep() to do
    - max step size to take which becomes the argument to ```pMagFieldPropagator -> SetLargestAcceptableStep(...)```. The default value is 2000.0 (which now that I think about it should be multiplied by mm, but in our case mm is 1.0). If this is set to -1.0, then ```SetLargestAcceptableStep(...)``` is not called within ```testPropagateMagField```.
- Mag Field type (Uniform or quadropole)
    - If Uniform, must be followed by coordinates of the Uniform Mag Field Vector. The entries will be multiplied by tesla's.
    - If quadropole, must be followed by a number which gets multiplied by ```10.*tesla/(50.*cm``` in the constructor of the quadropole field.
 - cached on/off, and cache distance in cm if value is "on"
 - geometry on/off (in effect skips the section of ```BuildGeometry()``` which places the volumes). (default is on)
 - initial_pos/mom followed by 7 doubles
    - 3 position coordinates
    - 3 momentum coordinates (used to construct a normalized vector)
    - 1 momentum value
- file_store_info followed by 6 strings
    - output_filename, where the binary step data is stored
    - a meta_file name, a text file which has two entries
        - the number of stored regular steps
        - the number of stored overshoot steps (found when an intersection point is encountered)
    - an accumulated function calls filename, stores as a text file the number of used function calls up until each step (if there are N stored regular steps, there will be N entries in this file)
    - a similar file for storing the accumulated function calls used up until each overshoot step
    - intersection indices file name, a text file listing the indices of the stored regular steps, immediately after which there is a stored overshoot step (currently this is of limited use).
    - an overshoot steps filename, a binary file to store the step data of the overshoot segments.

These argument patterns can be supplied in any order except that the first entry should be the stepper type integer. Not all arguments need to be supplied. At a minimum a stepper type integer should be given. Here is an example:

```./testPropagateMagField 1 init_data 1000.0 15 1500.0 initial_pos/mom 0.0 0.0 100.0 0.2 0.6 0.9 0.5 quadropole 1.0 cached_on 1.0 geometry_on```



The storage format for the binary files is a line for each completed step, with each line containing either 22 doubles (w/o the flag ```INTENDED_FOR_ERROR_BY_STEPPER_PROGRAM``` defined in ```isTracking.hh```) or 28 (if the flag is defined):
- 2 doubles representing accumulated time and arc length (in that order) at the initial point of the step
- 6 doubles representing position and velocity at the initial point of the step (not momentum, since interpolation for Nystrom methods is done with respect to position and velocity coordinates)
- 3 doubles representing a Right hand side evaluation at the initial point of the step. (This is the derivative of the velocity, not the derivative of the momentum)
- 2 doubles representing the accumulated time/arc length at the endpoint of the step
- 6 doubles ... at the end point of the step
- ... and so on (same as for initial point)

If the ```INTENDED_FOR_ERROR_BY_STEPPER_PROGRAM``` flag is defined, then the 3 doubles used to store a Right hand side evaluation are expanded to 6 doubles, which are used to store the right hand side evaluation of a regular Runge-Kutta stepper (in this case G4ClassicalRK4).



