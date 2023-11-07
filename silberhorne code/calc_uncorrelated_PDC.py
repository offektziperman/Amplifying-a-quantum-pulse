# fc_pdc_calc: Calculates the time-ordered solutions of ultrafast frequency
# conversion / parametric down-conversion and evaluates the process parameters.
#     Copyright (C) 2013  Andreas Christ
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.


#! /usr/bin/env python
""" This script calculates the analytic and rigorous solution of a 
uncorrelated parametric down-conversion process as given in the corresponding
paper. The parameters can easily be adapted to simulate a variety of
different processes"""

import calc_PDC

###################################################
## Parametric down-conversion process parameters ##
###################################################

# Coupling value giving the overall efficiency of the conversion process.
# It includes the pump power and the nonlinearity of the medium.
coupling = 0.75
# Width of the Gaussian pump field amplitude (not intensity) in sigma. 
pump_width = 0.96231155
# Inverse group velocity of the pump beam
A = 3.0
# Inverse group velocity of the A beam
B = 4.5
# Inverse group velocity of the C beam
C = 1.5

# z_start gives the beginning of the crystal and z_stop its end
# (z_start should remain at 0 since the analytic formulas in the code
# assume a waveguide starting at 0.)
z_start = 0 
z_stop = 2



###########################
## Evaluation parameters ##
###########################

# Wavelength range to consider where the process is centered at 0. The values
# have to cover all excited frequencies. 
# (In some analysis scripts it is assumed  that w_start and w_stop are
# symmetric about zero.)
w_start = -5
w_stop = 5
# Sampling steps for the frequency degree of freedom
w_steps = 500

# Sampling steps for the propagation over the length of the crystal
# (Warning: A too small sampling value will lead to grave errors in the
# z-integration)
z_steps = 500


# Directory to save the data in
save_directory = "Data"


################################################################
## Calculate the PDC process using the above stated parameters ##
## and save the end result into a hdf5 file.                  ##
################################################################

calc_PDC.calc_PDC(coupling, w_start, w_stop, w_steps, z_start, z_stop, z_steps, A, B, C, pump_width, save_directory)

# Iterate over a range of coupling values
#coupling_range = np.arange(1.4, 2.1, 0.1)
#print "coupling_range:", coupling_range
#for coupling in coupling_range:
#    calc_PDC.calc_PDC(coupling, w_start, w_stop, w_steps, z_start, z_stop, z_steps, A, B, C, pump_width, save_directory)
