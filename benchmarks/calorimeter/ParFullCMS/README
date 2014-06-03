2013. New example migrated for Geant4 Version 10.0
Andrea Dotti
Example still needs cleanup. For this reason for the moment we
keep it separate from sequential version of the code


01.25.2009 Xin Dong: This example came from the original sequential
program FullCMS. The original program is changed here to support parallel
computing with multiple threads. All events are assigned to each worker
thread in a round robin fashion. All threads share most detector data                        
including physics table and physics vector for some physics processes.                       
The master process initializes the data in a regular way. However, worker                    
threads initialize thread private data only.

We use the previous trick to introduce parallelism by implementating a
new G4RunaManager subclass. However, both Geant4 kernel and CLHEP is
changed accordingly. The original README is attached. Compile this
example just as the original FullCMS. The executable file is named as
ParmainApplication. One more argument is needed to give the number of
workers. For example, use:

$G4BIN/$G4SYSTEM/ParmainApplication bench1.g4 8

to run this program. The third argument "8" is the number of worker
threads that will be created. So the total number of threads for this
application is 9.


The original README:

-------------------------------------------------------------------

$Id: README,v 1.1 2007/10/24 12:38:34 gcosmo Exp $
-------------------------------------------------------------------

 Full CMS Benchmark
 ------------------

 In this directory you can find a CPU benchmark test of Geant4
 based on the full CMS detector, imported via a .gdml file.

 To select a Physics List you have to define one of the following
 environmental variables:
   - LHEP :      for the  LHEP       Physics List;
   - QGSP :      for the  QGSP       Physics List;
   - QGSP_EMV :  for the  QGSP_EMV   Physics List;
   - QGSC :      for the  QGSC       Physics List;
   - FTFP :      for the  FTFP       Physics List;
   - QGSP_BIC :  for the  QGSP_BIC   Physics List;
   - QGSP_BERT : for the  QGSP_BERT  Physics List.
 For example, if you want to use QGSP_EMV Physics List you can
 do:  
                    export QGSP_EMV=1
 or, equivalently:
                    make QGSP_EMV=1

 To build the application, first setup your environmental variables
 (the Bash-shell setup file,  setup.sh , shows an example), and then
 do:
                                make XXX=1
 where "XXX" is the name of the Physics List, or, equivalently:
                                export XXX=1 ; make
 and you get the executable:
                                $G4BIN/$G4SYSTEM/mainApplication
 and to run it:
                                $G4BIN/$G4SYSTEM/mainApplication bench1.g4

 You can run this application with the following macro file:
   -  bench1.g4 : 4000 events, each consisting of a beam particle
	          shot into the full CMS detector, with a uniform
	          magnetic field of 4 Tesla along the Z-axis.
	          The beam particle has the following characteristics:
                    o  random particle type 
	               (draw with equal probability between:
                        mu-, mu+, e-, e+, gamma, pi-, pi+, kaon-,
                        kaon+, kaon0L, neutron, proton, anti_neutron,
                        anti_proton, deuteron, triton, alpha, lambda,
                        sigma+, sigma-, xi-, xi0, anti_lambda,
                        anti_sigma+, anti_sigma-, anti_xi-, anti_xi0,
                        omega-, anti_omega- )
                    o  random kinetic energy
	               (draw uniformily in the interval:  1 - 100 GeV )
                    o  starting at the origin (0,0,0)
                    o  with initial random direction
                       (draw uniformily in 4*pi).
                  NB) You can change any of the above choices
	              (for instance shooting always 50 GeV pi-
	               in a given, fixed direction)
	              by modifying the file:
                         src/MyPrimaryGeneratorAction.cc .

 The CPU time for this test can be obtained in two ways
 (which should be, more or less, in agreement):
   - Look at the value "User=..." at the end of the running,
     after the line "Run Summary": this is the total time,
     in seconds, for all (4000) events, excluding the 
     initialization.
   - Use:
            time $G4BIN/$G4SYSTEM/mainApplication bench1.g4
     when launching the program: you would get, at the end
     of the program, the value: "user ..."  which is the 
     total time for all (4000) events, including the 
     initialization.

 Finally, notice that the macro file starts with the same seed number
 (taken from the file  start.rndm ), so if you run twice in the same
 machine you should get the same result, although the time can vary 
 slightly due to the different condition of the machine.
