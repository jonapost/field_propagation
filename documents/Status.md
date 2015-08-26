# Status of Nystrom-dev (Jason Suagee)

## Easier things to do
- ```MuruaRKN5459``` needs an interpolant, and it should be used in the ```DistChord()``` method instead of an auxilary stepper.

- ```MuruaRKN5459``` needs to have its defining coefficients made class variables.

- The inverse of ```FMass()``` in ```MagIntegratorStepper_byTime``` should be precomputed so that we do not repeatedly have to perform the division in the ```Stepper()``` and ```ComputeRightHandSide()``` methods. 

-   ```MagEqRhs_byTime_storeB``` is a class used by the Murua stepper to store field evaluations in mid ```Stepper()``` call. It is needed (or something like it) because for 3 of the ```ComputeRightHandSide()``` calls inside of ```MuruaRKN5459::Stepper()``` the value of the magnetic field is reused. This functionality should either be refactored into the ```ComputeRightHandSide()``` method (maybe piggybacking off an FSAL implementation), or borrowed from the ```G4CachedMagneticField``` class, which stores the previous right hand side evaluation already. This will prevent needless recopying of the field evaluations. Although if we rely on ```G4CachedMagneticField``` for this functionality we bind ourselves into always using ```G4CachedMagneticField``` with the Murua stepper.

## NTST demo program

- Remember to initialize the StepTracker pointer member variable in each stack object to the StepTracker object (if using StepTracker). See the caveat.md file for a more full explanation.


## Design decisions

- ```Interpolant``` should be made an abstract class, and particular interpolants for particular steppers derived from it. For instance we could derive a regular 5th order interpolant based on Hermite interpolation. 



## More General TODO
- All of this code can be refactored to work with the FSAL implementation: ```FSALMagIntegratorDriver```.

- Perhaps the full stack should be rewritten (only partially) to accommodate Nystrom steppers. This means including the option to use position and velocity as the coordinates instead of position and momentum. Checks would have to be implemented to ensure that velocity does not exceed the speed of light. Code would also have to be written to simulate integration by arclength instead of time. Also a field's ```EvaluateRhsGivenB()``` method should be rewritten to not store anything in the first 3 coordinates of ```dydx``` since this is unused by a Nystrom stepper.


## StepTracker related
- ```StepTracker::last_velocity()``` is used to update the the accumulated arclength within a ```G4ChordFinder::AdvanceChordLimited()``` call. Currently ```last_velocity()``` uses the velocity of the last stored point to calculate this, and if the velocity is not constant this might present a problem. At the beginning of each call to ```AdvanceChordLimited()``` the accumulated arclength in the ```StepTracker``` instance is updated with the arc length value in the passed ```G4FieldTrack``` object, so the accumulated arc length in ```StepTracker``` probably never gets too out of sync with the rest of the stack.

