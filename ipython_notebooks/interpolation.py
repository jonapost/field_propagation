print '''
This is intended to be a module for an ipython notebook.
If you intend to use the 3D plotting capabilities of the module
you should start ipython with the command:

ipython notebook --gui=wx

or whatever may be more appropriate on your system.

The 3D plotting capabilities in this module depend on the following

mayavi and tvtk
see: http://docs.enthought.com/mayavi/mayavi/
and  http://docs.enthought.com/mayavi/tvtk/

for more information.

      '''


import os

import numpy as N
import matplotlib.pyplot as plt
import math
import array
from matplotlib import colors
import subprocess
from subprocess import check_call

# comment out the following 3 import statements if you do not want to enable 3d plots
# or don't want to go through the trouble of installing these modules
import mayavi.mlab as mylab
import tvtk
from tvtk.tools import visual

import matplotlib.patches as mpatches

class Globals:
   pass # a container object for all globals

G = Globals()


# Set these by calling set_DIRs( binaries, data):
G.Binary_DIR = None
G.Data_DIR = './'
G.Plot_Graphics_DIR = './'

G.BUFFER_COLUMN_LEN = 22 # Other possiblity is 28



G.tPMF_program = None 
G.interpolation_error_program = None
G.tPMF_output_prefix = 'tPMF_' # Prefixed to all non-baseline tPMF outputs

# Data regarding the current baseline solution:
G.baseline_set = False
G.baseline_number_full_steps = None
G.baseline_number_overshoot_segments = None

# This is a dictionary:
G.default_baseline_file_store_info = None

G.default_pos = [0.0, 0.0, 0.0]
G.default_mom_dir = [0.2, 0.6, 0.8]
G.default_initial_momentum = 0.5

G.default_ComputeStep_len = 2000.0
G.default_No_ComputeSteps = 15

G.default_field = 'uniform'
G.default_uniform_field_Bvalue = [0.0, 0.0, -1.0]
G.default_quadropole_field_strength = 1.0 # is multiplied by 10.*tesla/(50.*cm

# This has to do with whether we will use the functionality of G4CachedMagneticField:
G.default_is_cached = False
G.default_cache_distance = 0.0

# Whether geometry is on (intesection calculator):
G.default_is_geometry_on = True

G.mass = 1049.020001226783279 # this is proton_mass_c2 (change as necessary)

G.all_stepper_choices = [ (0, 'Murua5459'), (1, 'Fine45'), (4, 'CashKarp45'), (5, 'ClassicalRK4'), (7, 'DormandPrince745')]


def set_DIRs( binaries, data = './', plot_graphics = './'):
   G.Binary_DIR = binaries
   G.Data_DIR = data
   G.Plot_Graphics_DIR = plot_graphics
   G.tPMF_program = G.Binary_DIR + 'testTemplated/testPropagateMagField'
   G.interpolation_error_program = G.Binary_DIR + 'interpolation_error/interpolation_error'

def print_all_default_values():
   print 'Binary_DIR = ' + str(G.Binary_DIR)
   print 'Data_DIR = ' + str(G.Data_DIR)
   print 'Plot_Graphics_DIR = ' + str(G.Plot_Graphics_DIR)
   print 'BUFFER_COLUMN_LEN = ' + str(G.BUFFER_COLUMN_LEN)
   print 'tPMF_program = ' + str(G.tPMF_program)
   print 'interpolation_error_program = ' + str(G.interpolation_error_program)
   print 'tPMF_output_prefix = ' + str(G.tPMF_output_prefix)
   print 'baseline_set = ' + str(G.baseline_set)
   print 'baseline_number_full_steps = ' + str(G.baseline_number_full_steps)
   print 'baseline_number_overshoot_segments = ' + str(G.baseline_number_overshoot_segments)
   print 'default_baseline_file_store_info = ' + str(G.default_baseline_file_store_info)
   print 'default_pos = ' + str(G.default_pos)
   print 'default_mom_dir = ' + str(G.default_mom_dir)
   print 'default_initial_momentum = ' + str(G.default_initial_momentum)
   print 'default_ComputeStep_len = ' + str(G.default_ComputeStep_len)
   print 'default_No_ComputeSteps = ' + str(G.default_No_ComputeSteps)
   print 'default_field = ' + str(G.default_field)
   print 'default_uniform_field_Bvalue = ' + str(G.default_uniform_field_Bvalue)
   print 'default_quadropole_field_strength = ' + str(G.default_quadropole_field_strength)
   print 'default_is_cached = ' + str(G.default_is_cached)
   print 'default_cache_distance = ' + str(G.default_cache_distance)
   print 'default_is_geometry_on = ' + str(G.default_is_geometry_on)
   print 'mass = ' + str(G.mass)
   print 'all_stepper_choices = ' + str(G.all_stepper_choices)



def compute_new_baseline_with_tPMF(file_store_info, propagator_init_data = None):
   
   #
   # Default values:
   if propagator_init_data != None:
      propagator_init_data = {}
      propagator_init_data['ComputeStep_length'] = 0.1
      propagator_init_data['No_ComputeSteps'] = 100000
      propagator_init_data['max_step_size'] = 0.1
   #
   #stepper_type = propagator_init_data['stepper_type']
   stepper_type = 1 #Fine45 with 4th order interpolation
   ComputeStep_size = propagator_init_data['ComputeStep_length']
   No_ComputeSteps = propagator_init_data['No_ComputeSteps']
   #start_pos_mom = propagator_init_data['start_pos_mom']
   max_step_size = propagator_init_data['max_step_size']
   #
   cmd = [G.tPMF_program]
   cmd += [str(stepper_type), 'init_data'] + map(str, [ComputeStep_size, No_ComputeSteps, max_step_size] )
   cmd += ['initial_pos/mom'] + map(str, G.default_pos + G.default_mom_dir + [G.default_initial_momentum])
   cmd += [G.default_field]
   if G.default_field == 'uniform':
      cmd += map(str, G.default_uniform_field_Bvalue)
   else: # was a Quadropole Field
      cmd += map(str, [G.default_quadropole_field_strength])
   if G.default_is_cached == True:
      cmd += [ 'cached_on' ]
      cmd += map(str, [G.default_cache_distance])
   else:
      cmd += [ 'cached_off' ]
   if G.default_is_geometry_on:
      cmd += [ 'geometry_on' ]
   else:
      cmd += [ 'geometry_off' ]
   cmd += ['file_store_info']
   cmd += [G.Data_DIR + file_store_info['output_filename']]
   cmd += [G.Data_DIR + file_store_info['meta_filename']]
   cmd += [G.Data_DIR + file_store_info['no_function_calls_output_filename']]
   cmd += [G.Data_DIR + file_store_info['no_function_calls_overshoot_filename']]
   cmd += [G.Data_DIR + file_store_info['intersection_indices_filename']]
   cmd += [G.Data_DIR + file_store_info['overshoot_segments_filename']]
   print ''.join( map( lambda x: x + ' ', cmd ) )
   check_call( cmd )
   f = open(G.Data_DIR + file_store_info['meta_filename'])
   G.baseline_number_full_steps = int(f.readline())
   G.baseline_number_overshoot_segments = int(f.readline())
   f.close()
   G.baseline_set = True
      
   
def compute_new_baseline(propagator_init_data, file_store_info = None, use_tPMF = True):
   #global use_tPMF
   if use_tPMF == True:
      if file_store_info == None:
         file_store_info = {}
         file_store_info['output_filename'] = 'baseline_out'
         file_store_info['meta_filename'] = 'baseline_meta'
         file_store_info['no_function_calls_output_filename'] = 'baseline_function_calls'
         file_store_info['no_function_calls_overshoot_filename'] = 'baseline_overshoot_function_calls'
         file_store_info['intersection_indices_filename'] = 'baseline_intersection_indices'
         file_store_info['overshoot_segments_filename'] = 'baseline_overshoot_segments'
      #file_store_info = run_info['file_store_info']
      #propagator_init_data = run_info['propagator_init_data']
      compute_new_baseline_with_tPMF(file_store_info, propagator_init_data)
      G.default_baseline_file_store_info = file_store_info
      
   

def set_tPMF_output_prefix( prefix ):
   G.tPMF_output_prefix = prefix
   
def set_initial_position_momentum( pos, mom_dir, initial_momentum ):
   G.default_pos = pos
   G.default_mom_dir = mom_dir
   G.default_initial_momentum = initial_momentum
   
   
def set_default_ComputeStep_No_and_Length(ComputeStep_len, no_ComputeSteps ):
   G.default_ComputeStep_len = ComputeStep_len
   G.default_No_ComputeSteps = no_ComputeSteps
   
def set_default_is_geometry_on( geometry_is_on = True ):
   G.default_is_geometry_on = geometry_is_on

def set_default_field_data(field_type, field_data, cache_distance = 1.0):
   G.default_field = field_type
   if field_type == 'uniform':
      G.default_uniform_field_Bvalue = field_data
   else: 
      if field_type == 'quadropole':
         G.default_quadropole_field_strength = field_data
      else:
         print 'field type (argument 1) must be either "uniform" or "quadropole".'
         print 'If uniform, then 2nd argument is the B field vector (values will be multiplied by teslas).'
         print 'If quadropole, then 2nd argumemnt must be a float which will determine the field strength??.'
   if cache_distance == 0.0:
      G.default_is_cached = False
   else:
      G.default_is_cached = True
   G.default_cache_distance = cache_distance
   
def set_particle_mass( particle_mass ):
   G.mass = particle_mass
   
def get_trajectory( stepper_no, ComputeStep_len = 2000.0, No_ComputeSteps = 10.0, largest_ComputeStep = -1.0 ):
   if not G.baseline_set:
      print 'Please call "compute_new_baseline()" first.'
      return
   cmd = [G.tPMF_program]
   cmd += [str(stepper_no), 'init_data'] + map(str, [ComputeStep_len, No_ComputeSteps, largest_ComputeStep])
   cmd += ['initial_pos/mom']
   cmd += map(str, G.default_pos) + map(str, G.default_mom_dir) + [str(G.default_initial_momentum)]
   if G.default_field == 'uniform':
      cmd += ['uniform'] + map(str, G.default_uniform_field_Bvalue)
   else:
      if G.default_field == 'quadropole':
         cmd += [ 'quadropole', str(G.default_quadropole_field_strength) ]
      else:
         print 'Did you forget to call set_default_field_data(field_type,..)'
         return
   if G.default_is_cached:
      cmd += ['cached_on', str(G.default_cache_distance)]
   else:
      cmd += ['cached_off']
   if G.default_is_geometry_on:
      cmd += ['geometry_on']
   else:
      cmd += ['geometry_off']
   cmd += ['file_store_info']
   files_list = map(lambda x: G.tPMF_output_prefix + x, ['output', 'meta_output', 'no_function_calls', 'no_function_calls_overshoot'])
   files_list += map(lambda x: G.tPMF_output_prefix + x, ['intersection_indices', 'overshoot_segments'])
   cmd += map(lambda x: G.Data_DIR + x, files_list)
   #cmd += files_list
   #print ''.join( map( lambda x: x + ' ', cmd) )
   # Temp debug:
   #cmd = [ G.tPMF_program, '1']
   check_call( cmd )


def get_interpolated_error( stepper_no, ComputeStep_len = 2000.0, \
                  No_ComputeSteps = 10.0, with_trajectory_data = False, baseline_filename_list = None ):
   get_trajectory( stepper_no, ComputeStep_len, No_ComputeSteps )
   #
   f = open( G.Data_DIR + G.tPMF_output_prefix + 'meta_output' )
   number_full_steps = int(f.readline())
   number_overshoot_segments = int(f.readline())
   f.close()
   cmd = [G.interpolation_error_program]
   cmd += map(str, [G.baseline_number_full_steps, G.baseline_number_overshoot_segments])
   cmd += [str(number_full_steps)]
   if baseline_filename_list != None:
      cmd += [G.Data_DIR + baseline_filename_list['output_filename']]
      cmd += [G.Data_DIR + G.tPMF_output_prefix + 'output']
      cmd += [G.Data_DIR + 'error']
      cmd += [G.Data_DIR + baseline_filename_list['overshoot_segments_filename']]
   else:
      cmd += [G.Data_DIR + G.default_baseline_file_store_info['output_filename']]
      cmd += [G.Data_DIR + G.tPMF_output_prefix + 'output']
      cmd += [G.Data_DIR + 'error']
      cmd += [G.Data_DIR + G.default_baseline_file_store_info['overshoot_segments_filename']]
   check_call( cmd )
   #
   f = open(G.Data_DIR + 'interpolation_error_metafile')
   no_interpolated = int(f.readline())
   f.close()
   f = open( G.Data_DIR + 'error', mode = 'rb')
   binvals = array.array('d')
   binvals.read(f, no_interpolated * (8))
   data = N.array(binvals, dtype = float)
   f.close()
   #
   interpolation_errors = [ [None for j in range(8)] for i in range(no_interpolated)]
   counter = 0
   for i in range(no_interpolated):
      for j in range(8):
         interpolation_errors[i][j] = data[counter]
         counter += 1
   #
   function_evals = []
   f = open(G.Data_DIR + G.tPMF_output_prefix + 'no_function_calls')
   for i in range(no_interpolated):
      function_evals += [int(f.readline())]
   f.close()
   #
   if with_trajectory_data:
      f = open(G.Data_DIR + G.tPMF_output_prefix + 'output', mode = 'rb')
      binvals = array.array('d')
      binvals.read(f, number_full_steps * (G.BUFFER_COLUMN_LEN))
      data = N.array(binvals, dtype = float)
      f.close()
      #
      tPMF_trajectory = [ [None for j in range(G.BUFFER_COLUMN_LEN)] for i in range(number_full_steps)]
      counter = 0
      for i in range(number_full_steps):
         for j in range(G.BUFFER_COLUMN_LEN):
            tPMF_trajectory[i][j] = data[counter]
            counter += 1
      return (interpolation_errors, function_evals, tPMF_trajectory)
   else:
      return (interpolation_errors, function_evals)



def plot_errors(stepper_list, with_function_evaluation_count = False, save_file_name = None, tPMF_baseline = True):
   colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
   if len(stepper_list) > len(colors):
      print "Not enough colors to plot all steppers. (There are 8 colors)"
      return
   stepper_choices = []
   for stepper_no, stepper_name in G.all_stepper_choices:
      if (stepper_no in stepper_list) or (stepper_name in stepper_list):
         stepper_choices += [(stepper_no, stepper_name)]
   #
   steppers_big_data_list = []
   for i, name in stepper_choices:
      steppers_big_data_list += [get_interpolated_error(i, G.default_ComputeStep_len, G.default_No_ComputeSteps)]
   cutoff_list = map( lambda m: len(m[0]), steppers_big_data_list )
   #
   f1 = plt.figure()
   counter = 0
   for A, t in steppers_big_data_list:
       P = [ math.sqrt( sum([A[i][2 + k]**2 for k in range(3)]) ) for i in range(cutoff_list[counter])]
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], P, colors[counter]  )
       counter += 1
   plt.title('Position Error vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_position_abs_error.png')
   plt.show()
   #
   f2 = plt.figure()
   counter = 0
   for A, t in steppers_big_data_list:
       V = [ G.mass*math.sqrt( sum([A[i][2 + k]**2 for k in range(3,6)]) ) for i in range(cutoff_list[counter])]
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], V, colors[counter]  )
       counter += 1
   plt.title('Momentum Error vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_momentum_abs_error.png')
   plt.show()
   #
   if with_function_evaluation_count == False:
      return
   f3 = plt.figure()
   counter = 0
   for A, function_evals in steppers_big_data_list:
       colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], function_evals[:cutoff_list[counter]], colors[counter] )
       counter += 1
   plt.title('No. Function Calls vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_function_calls.png')
   plt.show()
   

# Some Helper functions:
def pos_mag(v):
   return math.sqrt( v[2]**2 + v[3]**2 + v[4]**2 )
def mom_mag(v):
   return math.sqrt( v[5]**2 + v[6]**2 + v[7]**2 )

def pos_relative_error_entry(A, A_tPMF, i):
   #A_tPMF is the approximating solution, A is the component wise difference from the "exact" solution
   # thus, p is the components of the "exact" solution:
   p = (A_tPMF[i][2] - A[i][2], A_tPMF[i][3] - A[i][3], A_tPMF[i][4] - A[i][4])
   return pos_mag(A[i]) / math.sqrt( p[0]**2 + p[1]**2 + p[2]**2 )

def mom_relative_error_entry(A, A_tPMF, i):
   #A_tPMF is the approximating solution, A is the component wise difference from the "exact" solution
   #at the intepolating time.
   # thus, p is the components of the "exact" solution:
   p = (A_tPMF[i][5] - A[i][5], A_tPMF[i][6] - A[i][6], A_tPMF[i][7] - A[i][7])
   return mom_mag(A[i]) / math.sqrt( p[0]**2 + p[1]**2 + p[2]**2 ) # mass scaling on both sides of the '/' so cancel out.
# End Helper functions.

def plot_relative_errors(stepper_list, with_function_evaluation_count = False, save_file_name = None, tPMF_baseline = True):
   colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
   if len(stepper_list) > len(colors):
      print "Not enough colors to plot all steppers. (There are 8 colors)"
      return
   stepper_choices = []
   for stepper_no, stepper_name in G.all_stepper_choices:
      if (stepper_no in stepper_list) or (stepper_name in stepper_list):
         stepper_choices += [(stepper_no, stepper_name)]
   #
   steppers_big_data_list = []
   for i, name in stepper_choices:
      steppers_big_data_list += [get_interpolated_error(i, G.default_ComputeStep_len, G.default_No_ComputeSteps, True)]
   cutoff_list = map( lambda m: len(m[0]), steppers_big_data_list )
   #
   plt.figure()
   counter = 0
   for A, function_evals, trajectory in steppers_big_data_list:
       #P = [ math.sqrt( sum([A[i][2 + k]**2 for k in range(3)]) ) for i in range(cutoff_list[counter])]
       P = [pos_relative_error_entry(A, trajectory, i) for i in range(cutoff_list[counter])]
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], P, colors[counter]  )
       counter += 1
   plt.title('Position Error vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   plt.show()
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_position_rel_error.png')
   plt.show()
   #
   plt.figure()
   counter = 0
   for A, function_evals, trajectory in steppers_big_data_list:
       #V = [ G.mass*math.sqrt( sum([A[i][2 + k]**2 for k in range(3,6)]) ) for i in range(cutoff_list[counter])]
       V = [mom_relative_error_entry(A, trajectory, i) for i in range(cutoff_list[counter])]
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], V, colors[counter]  )
       counter += 1
   plt.title('Momentum Error vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_momentum_rel_error.png')
   plt.show()
   #
   if with_function_evaluation_count == False:
      return
   plt.figure()
   counter = 0
   for A, function_evals, trajectory in steppers_big_data_list:
       colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
       plt.plot( [A[i][1] for i in range(cutoff_list[counter])], function_evals[:cutoff_list[counter]], colors[counter] )
       counter += 1
   plt.title('No. Function Calls vs. ArcLength, initial step dist.: ' + str(G.default_ComputeStep_len))
   patches = []
   for i in range(len(stepper_choices)):
       patches += [mpatches.Patch(color = colors[i], label = stepper_choices[i][1])]
   plt.legend( patches, [stepper_choices[i][1] for i in range(len(stepper_choices))], loc = 2)
   if save_file_name != None:
      plt.savefig(G.Plot_Graphics_DIR + save_file_name + '_function_calls.png')
   plt.show()
   
   
def plot_3dtrajectories( stepper_list, ComputeStep_len, No_steps, tube_rad = 6, with_baseline = False, background = True ):
   colors = [(0.,1.,0.), (1.,0.,0.), (0.,1.,1.), (1.,0.,1.), (1.,1.,0.)]
   if len(stepper_list) > len(colors):
      print "Not enough colors to plot all steppers. (There are 5 colors)"
      return
   tPMF_trajectories = []
   for stepper_no in stepper_list:
      get_trajectory(stepper_no, ComputeStep_len, No_steps)
      #
      f = open(G.Data_DIR + G.tPMF_output_prefix + 'meta_output')
      no_steps_tPMF = int( f.readline() )
      f.close()
      f = open(G.Data_DIR + G.tPMF_output_prefix + 'output', mode = 'rb')
      binvals = array.array('d')
      binvals.read(f, no_steps_tPMF * (G.BUFFER_COLUMN_LEN))
      data = N.array(binvals, dtype = float)
      f.close()
      #
      tPMF_trajectory = [ [None for j in range(G.BUFFER_COLUMN_LEN)] for i in range(no_steps_tPMF)]
      counter = 0
      for i in range(no_steps_tPMF):
         for j in range(G.BUFFER_COLUMN_LEN):
            tPMF_trajectory[i][j] = data[counter]
            counter += 1
      tPMF_trajectories += [( tPMF_trajectory, len(tPMF_trajectory) )]
   #
   f = open(G.Data_DIR + G.default_baseline_file_store_info['output_filename'], mode = 'rb')
   binvals = array.array('d')
   binvals.read(f, G.baseline_number_full_steps * (G.BUFFER_COLUMN_LEN))
   data = N.array(binvals, dtype = float)
   f.close()
   #
   baseline_trajectory = [ [None for j in range(G.BUFFER_COLUMN_LEN)] for i in range(G.baseline_number_full_steps)]
   counter = 0
   for i in range(G.baseline_number_full_steps):
      for j in range(G.BUFFER_COLUMN_LEN):
         baseline_trajectory[i][j] = data[counter]
         counter += 1
   f = mylab.figure()
   if background:
      visual.set_viewer(f)
      visual.box(pos = (0*1000,0,-15*1000), size = (10*1000,10*1000,10*1000))
      visual.box(pos = (0*1000,0,15*1000), size = (10*1000,10*1000,10*1000))
      visual.box(pos = (7.5*1000,0,0), size = (2.5*1000,2.5*1000,2.5*1000))
      visual.box(pos = (-7.5*1000,0,0), size = (2.5*1000,2.5*1000,2.5*1000))
      visual.box(pos = (0,7.5*1000,0), size = (2.5*1000,2.5*1000,2.5*1000))
      visual.box(pos = (0,-7.5*1000,0), size = (2.5*1000,2.5*1000,2.5*1000))
      #
      visual.box(pos = (.3*1000,.3*1000,.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (.3*1000,.3*1000,-.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (.3*1000,-.3*1000,.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (.3*1000,-.3*1000,-.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (-.3*1000,.3*1000,.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (-.3*1000,.3*1000,-.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (-.3*1000,-.3*1000,.3*1000), size = (.25*1000,.25*1000,.25*1000) )
      visual.box(pos = (-.3*1000,-.3*1000,-.3*1000), size = (.25*1000,.25*1000,.25*1000) )
   #
   counter = 0
   for tPMF_trajectory, l in tPMF_trajectories:
      X = [tPMF_trajectory[i][2] for i in range(l)]
      Y = [tPMF_trajectory[i][3] for i in range(l)]
      Z = [tPMF_trajectory[i][4] for i in range(l)]
      mylab.plot3d(X,Y,Z, tube_radius = tube_rad, color = colors[counter])
      counter += 1
   #
   X = [baseline_trajectory[i][2] for i in range(G.baseline_number_full_steps)]
   Y = [baseline_trajectory[i][3] for i in range(G.baseline_number_full_steps)]
   Z = [baseline_trajectory[i][4] for i in range(G.baseline_number_full_steps)]
   mylab.plot3d(X,Y,Z, tube_radius = tube_rad, color = (0.,0.,1.))
   #
   mylab.show()
   
def help():
   print '''
Main functions are:

get_interpolated_error(with arguments stepper_no, ComputeStep_len = 2000.0,
   No_ComputeSteps = 10.0, with_trajectory_data = False, baseline_filename_list = None)

Use this function to get raw interpolation error or trajectory data of a stepper.
Returns (interpolation error matrix, cumulative function evaluation data,
and trajectory data (if the argument with_trajectory_data is set to True).

plot_errors(stepper_list, with_function_evaluation_count = False, 
                       save_file_name = None, tPMF_baseline = True).

This plots absolute errors. stepper_list can be a list of integers corresponding to stepper
entries in the default global all_stepper_choices. You can also specify steppers in stepper_list
by specifying the named string (e.g. [1, 4, 'ClassicalRK4'] is a valid input).
Set with_function_evaluation_count = True if you want to
also output a plot of the cumulative function calls. Specify a save_file_name if you
would like to save the figures as .png's. save_file_name will be prefixed to a descriptive
affix.

plot_relative_errors()

Same as plot_errors(), but plots relative errors. Same argument structure.

plot_3dtrajectories( stepper_list, ComputeStep_len, No_steps, 
                     tube_rad = 6, with_baseline = False, background = True )

Plot's 3D trajectories of a list of steppers (up to 5, because only 5 colors are supplied -
see line 440). If with_baseline is set to True, it will draw the baseline solution also.
Note that this could take a lot of time/memory since the baseline solution will most likely 
have a lot more points to plot.
If background is left as True, the boxes (geometric volumes) that are specified in 
testPropagateMagField will be rendered.

Before calling any of the main functions, you first have to initialize the baseline trajectory.
Do that by calling 

compute_new_baseline(propagator_init_data = None)

where propagator_init_data is a dictionary of the following form 
(which acts as the default value):

propagator_init_data = {}
propagator_init_data['ComputeStep_length'] = 0.1
propagator_init_data['No_ComputeSteps'] = 100000
propagator_init_data['max_step_size'] = 0.1

In addition there are the following get/set methods. Call list_get_set() to see them.
Some of them have to be called to intialize global variables before you can call the 
above main functions.
         '''

def list_get_set():
   print '''
set_DIRs( binaries, data = './', plot_graphics = './')

set_tPMF_output_prefix( prefix ) -- Only useful if you are going to be creating a lot of plots
                                    and have to keep them organized.
                                    
set_default_is_geometry_on( geometry_is_on = True ) -- Turn on or off intersection location
   
set_default_field_data(field_type, field_data, cache_distance = 1.0)
      -- field_type must be "uniform" or "quadropole"
      -- if uniform then field_data must be a list of the form [B_x, B_y, B_z]
      -- if quadropole then field_data must be a double which will be multiplied by 
         10.*tesla/(50.*cm) in the constructor of G4QuadropoleMagField
      -- cache_distance can be set to zero.
            
set_particle_mass(mass) -- The particle mass is used to properly scale the absolute error
                           in velocity back to momentum coordinates. The default is currently
                           proton_mass_c2.

set_default_ComputeStep_No_and_Length(ComputeStep_len, no_ComputeSteps )
                        -- this function must be called prior to all other main functions
                           except compute_new_baseline().
                           
set_initial_position_momentum( pos, mom_dir, initial_momentum )
      -- This function must be called prior to all main functions.
         pos is a list of the form [x0, y0, z0] (they will be multiplied by mm units)
         mom_dir is a list of the form [p_x, p_y, p_z] (this represents a unit momentum vector).
         initial_momentum is a double, which will be multiplied by proton_mass_c2.

call see_it_in_action() to see an example of how this is all put together.

         '''
         
def see_it_in_action():
   print '''
A sample run:
   
inter.set_DIRs('../field_propagation-build/test/', './', [Your Graphics Directory/])
inter.set_initial_position_momentum([-100.0,180.0,500.0],[0.2,0.6,0.8], 0.5)
inter.set_default_field_data('quadropole', 1.0, 0.0)
inter.set_default_ComputeStep_No_and_Length(1000.0, [Number of 
                      G4PropagatorInField::ComputeStep()'s you want the steppers to take])

propagator_init_data = {}
propagator_init_data['ComputeStep_length'] = 0.1
propagator_init_data['No_ComputeSteps'] = 100000
propagator_init_data['max_step_size'] = 0.1

inter.compute_new_baseline(propagator_init_data)

inter.plot_errors([0,1,4,5,7],False, 'YourNameHere')
inter.plot_relative_errors([0,1,4,5],False, 'YourNameHere')
inter.plot_3dtrajectories([0,1,4,5], 100.0, 20)

Now copy/paste and check it out.
         '''


   
# Print out a welcome message on module load
print 'Please set the Binary and Data Directories by calling set_DIRs( binaries_dir, data_dir)'
print 'Please append a forward slash to directory names (e.g. "my_dir/")'
print 'These can be relative paths. The default data_dir is the current directory (so be careful).'
print 'There is no default Binary directory. It should be wherever you put field_propagation/test/'
print ''
print 'Call the function help() to see more information.'
print ''

# Print out all default values
print 'Current default values are (called with print_all_default_values():'
print ''
print_all_default_values()

