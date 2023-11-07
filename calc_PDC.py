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


# ! /usr/bin/env python
""" This file contains the function calc_PDC with calculates
the analytic and rigorous solution of a given parametric dwon-conversion
process as specified in the corresponding paper. It saves the data
in a hdf5 file for further analysis.

Details on the math can be found in the corresponding paper."""
import matplotlib.pyplot as plt
import numpy as np
import tables


############################
# Definitions and formulas #
############################

## Define PDC functions

# x = omega
# y = omega'

# Phasematching function
# as defined for E_a: x = omega, y = omega', k_p = A, k_s = B, k_i = C
def dk(x, y, A, B, C):
    """ Return the k-vector mismatch between the three waves. """
    return (A * (x + y) - B * (x) - C * (y))


# Pump distribution
# Note that the the pump distribution is not normalized to replicate
# the formula in the paper.
def alpha(x, sigma_x):
    """ Return the amplitude of the pump distribution alpha center at 0
    with given width sigma and frequency x"""
    # Gaussian pump
    dx = x[1] - x[0]
    alpha = np.exp(-x ** 2 / (2 * sigma_x ** 2))
    return alpha/np.sum(alpha*dx)



# Analytical phasematching function
def pm(x, y, z_start, z_stop, A, B, C):
    """ Analytic phasematching function for a waveguide of length
    z_stop - z_start. It is however assumed that the waveguide
    starts at 0 and ends at z_stop - z_start. This enables us to
    write the solution as a sinc which circumvents various numerical
    issues. """

    return (z_stop - z_start) * np.sinc((z_stop - z_start) / 2 * dk(x, y, A, B, C) / np.pi) \
           * np.exp(-1j * dk(x, y, A, B, C) * (z_stop - z_start) / 2)


# Analytical first-order solution for Va (A beam)
# (I cannot name it f like in the paper since this is occupied by the f
# function from the differential equation)
def ra_analytical(x, y, z_start, z_stop, coupling, A, B, C, pump_width):
    """ Analytical solution for the A beam of first-order parametric
    down-conversion."""
    return -1j * np.conjugate(coupling) * np.conjugate(alpha(x + y, pump_width)) \
           * pm(x, y, z_start, z_stop, A, B, C)


# Analytical first-order solution for Vb (B beam)
# (I cannot name it f like in the paper since this is occupied by the f
# function from the differential equation)
def rb_analytical(x, y, z_start, z_stop, coupling, A, B, C, pump_width):
    """ Analytical solution for the B beam of first-order PDC.
        Note: rb_analytical = x-y switdch of ra_analytical"""
    return -1j * np.conjugate(coupling) * np.conjugate(alpha(y + x, pump_width)) \
           * pm(y, x, z_start, z_stop, A, B, C)


def phi_num(x, y, z, A, B, C):
    """ Return the value of the phi function in the differential equation"""
    return np.exp(-1j * dk(x, y, A, B, C) * z)

    # Additional hypergrating to get rid of sinc side peaks:
    # (Currently not used)
    # return grating(z) * np.exp(-1j * dk(x, y, A, B, C)*z)


# f function with prefactor 
def f_a(x, y, z, coupling, A, B, C, pump_width):
    """ Body of the differential equation for A 
    (on the right of the differential equation)"""
    return -1j * coupling * alpha(x + y, pump_width) \
           * phi_num(x, y, z, A, B, C)


def f_b(x, y, z, coupling, A, B, C, pump_width):
    """ Body of the differential equation for C 
    (on the right of the differential equation)"""
    # This function is just the x-y switched version of f_a
    return f_a(y, x, z, coupling, A, B, C, pump_width)


# The big function creation the analytic and rigorous solution for PDC
def calc_PDC(coupling, w_start, w_stop, w_steps, z_start, z_stop, z_steps, \
             A, B, C, pump_width, save_directory):
    """ Solve PDC in three different way:
    1) Analytical first-order solution
    2) Numerical first-order solution (identical to analytic solution, but 
       features a numerical z-integration)
    3) Rigorous solution Vaa iterative Euler approach
    
    All data is saved in a hdf5 file which can get very big."""

    ## Calculate vectors for the iteration process
    wRange = np.linspace(w_start, w_stop, w_steps)
    zRange = np.linspace(z_start, z_stop, z_steps)



    # Calculate the step sizes
    w_step = wRange[1] - wRange[0]
    z_step = zRange[1] - zRange[0]


    ## Generate Initial Matrices to store the results for Ua, Ub, Va and Vb:
    #   In all generated matrices "dtype=complex" has to be explicitly stated,
    #   otherwise the applied matrix operations cannot handle complex numbers
    # print("# Generating U and V initial #")

    Vinitial = np.zeros((wRange.size, wRange.size), dtype=complex)
    # This is the definition of a discretized delta function (note the /w_step)
    Uinitial = np.identity(wRange.size, dtype=complex) / w_step

    # Create and open Hdf5 File
    ## Attention do not change the naming, it's exact nature is used in the 
    ## analysis scripts.
    directory_name = save_directory + "/PDC_calc"
    evaluation_parameters = "---PAR--pow_" + "%.3f" % coupling \
                            + "--w_" + "%.1f" % w_start + "_" + "%.1f" % w_stop + "_" \
                            + "%.4d" % w_steps + \
                            "--z_" + "%.1f" % z_start + "_" + "%.1f" % z_stop + "_" \
                            + "%.4d" % z_steps
    pdc_parameters = "---PDC--A_" + "%.1f" % A + "--B_" + "%.1f" % B \
                     + "--C_" + "%.1f" % C + "--pumpWidth_" + "%.5f" % pump_width

    filename = directory_name + evaluation_parameters + pdc_parameters + ".h5"
    h5file = tables.open_file(filename, mode="w", title="PDC-Bog-Prop")

    # Create subgroup for Va, Vb, Ua and Ub
    groupV = h5file.create_group("/", "V", "z-Position")
    groupU = h5file.create_group("/", "U", "z-Position")

    # Generate initial Va and Ua Arrays

    # For some reason pytables does not want identifiers with a "."
    # consequently the zPos is identified by  0,1,2,3 and the actual position
    # has to be recalculated with help of the given z_start and z_step 

    for i in np.arange(zRange.size):
        h5file.create_array(groupV, "zPos" + str(i), Vinitial, "z-Position")

    for i in np.arange(zRange.size):
        h5file.create_array(groupU, "zPos" + str(i), Uinitial, "z-Position")


    ## Create f functions

    # Note the flipped Y, X in order to cope with numpys  array definition
    # The matrix multiplication used in our solution will only work with this
    # exact notation
    # Please refer to the memo for the in's and out's of this procedure
    Y, X = np.meshgrid(wRange, wRange)

    groupFa = h5file.create_group("/", "fa", "z Position")
    for i in np.arange(zRange.size):
        famatrix = f_a(X, Y, zRange[i], coupling, A, B, C, pump_width)
        h5file.create_array(groupFa, "zPos" + str(i), famatrix, "zPos-Set")


    # Actually fb is just a transposed fa
    groupFb = h5file.create_group("/", "fb", "z Position")
    for i in np.arange(zRange.size):
        fbmatrix = f_b(X, Y, zRange[i], coupling, A, B, C, pump_width)
        h5file.create_array(groupFb, "zPos" + str(i), fbmatrix, "zPos-Set")

    ## Calculate rigorous solution


    # This value should never be reached
    maxIterations = 100
    for iteration in np.arange(maxIterations):
        # print(iteration)
        VdiffNormMean = 0

        for i in np.arange(zRange.size):
            if (i == 0 or i == 1):
                if i == 0:
                    # At the first data point we cannot perform the trapezoid routine
                    # and instead opt for a pseudo midstep method. Just to be clear 
                    # this isn't a midstep method but gives sufficient results for 
                    # the first step 

                    # Load old Vb data
                    VTmpNode = h5file.get_node('/V', 'zPos' + str(i))
                    VTmp = VTmpNode.read()

                    fbtmpNode_0 = h5file.get_node('/fb', 'zPos' + str(i))
                    UTmpNode_0 = h5file.get_node('/U', 'zPos' + str(i))
                    fbtmp_0 = fbtmpNode_0.read()
                    UTmp_0 = UTmpNode_0.read()

                    VNewTmp_0 = z_step * np.dot(fbtmp_0, np.conjugate(UTmp_0)) * w_step
                    VNewTmp = VNewTmp_0

                    # Calculate difference between old and new Vb
                    VdiffNorm = np.linalg.norm(VNewTmp - VTmp) / np.linalg.norm(VNewTmp)
                    VdiffNormMean += VdiffNorm
                    # Save new Vb
                    VTmpNode[:] = VNewTmp

                if i == 1:
                    ## Trapezoid integration routine for the second data point
                    VTmpNode = h5file.get_node('/V', 'zPos' + str(i))
                    VTmp = VTmpNode.read()

                    fbtmpNode_0 = h5file.get_node('/fb', 'zPos' + str(i - 1))
                    fbtmpNode_1 = h5file.get_node('/fb', 'zPos' + str(i))
                    UTmpNode_0 = h5file.get_node('/U', 'zPos' + str(i - 1))
                    UTmpNode_1 = h5file.get_node('/U', 'zPos' + str(i))
                    fbtmp_0 = fbtmpNode_0.read()
                    fbtmp_1 = fbtmpNode_1.read()
                    UTmp_0 = UTmpNode_0.read()
                    UTmp_1 = UTmpNode_1.read()

                    VNewTmp_0 = z_step * np.dot(fbtmp_0, np.conjugate(UTmp_0)) * w_step
                    VNewTmp_1 = z_step * np.dot(fbtmp_1, np.conjugate(UTmp_1)) * w_step
                    VNewTmp = 0.5 * (VNewTmp_0 + VNewTmp_1)

                    # Calculate difference between old and new Vb
                    VdiffNorm = np.linalg.norm(VNewTmp - VTmp) / np.linalg.norm(VNewTmp)
                    VdiffNormMean += VdiffNorm

                    # Save new Vb
                    VTmpNode[:] = VNewTmp

            else:
                ## Trapezoid integration routine
                VNewTmpNode = h5file.get_node('/V', 'zPos' + str(i - 1))
                VTmpNode = h5file.get_node('/V', 'zPos' + str(i))
                VNewTmp = VNewTmpNode.read()
                VTmp = VTmpNode.read()

                fbtmpNode_0 = h5file.get_node('/fb', 'zPos' + str(i - 1))
                fbtmpNode_1 = h5file.get_node('/fb', 'zPos' + str(i))
                UTmpNode_0 = h5file.get_node('/U', 'zPos' + str(i - 1))
                UTmpNode_1 = h5file.get_node('/U', 'zPos' + str(i))
                fbtmp_0 = fbtmpNode_0.read()
                fbtmp_1 = fbtmpNode_1.read()
                UTmp_0 = UTmpNode_0.read()
                UTmp_1 = UTmpNode_1.read()

                VNewTmp_0 = z_step * np.dot(fbtmp_0, np.conjugate(UTmp_0)) * w_step
                VNewTmp_1 = z_step * np.dot(fbtmp_1, np.conjugate(UTmp_1)) * w_step
                VNewTmp += 0.5 * (VNewTmp_0 + VNewTmp_1)

                # Calculate difference between old and new Vb
                VdiffNorm = np.linalg.norm(VNewTmp - VTmp) / np.linalg.norm(VNewTmp)
                VdiffNormMean += VdiffNorm

                # Save new Vb
                VTmpNode[:] = VNewTmp

        VdiffNormMean = VdiffNormMean / zRange.size

        UdiffNormMean = 0
        UNewTmp = Uinitial.copy()
        for i in np.arange(zRange.size):
            if (i == 0):
                if i == 0:
                    # At the first data point we cannot perform the trapezoid routine
                    # and instead opt for a pseudo midstep method. Just to be clear 
                    # this isn't a midstep method but gives sufficient results for 
                    # the first step 

                    UTmpNode = h5file.get_node('/U', 'zPos' + str(i))
                    UTmp = UTmpNode.read()

                    fatmpNode_0 = h5file.get_node('/fa', 'zPos' + str(i))
                    VTmpNode_0 = h5file.get_node('/V', 'zPos' + str(i))
                    fatmp_0 = fatmpNode_0.read()
                    VTmp_0 = VTmpNode_0.read()

                    UNewTmp_0 = Uinitial + z_step * np.dot(fatmp_0, np.conjugate(VTmp_0)) * w_step
                    UNewTmp = UNewTmp_0

                    # Calculate difference between old and new Ua
                    UdiffNorm = np.linalg.norm(UNewTmp - UTmp) / np.linalg.norm(UNewTmp - Uinitial)
                    UdiffNormMean += UdiffNorm
                    # Save new Ua
                    UTmpNode[:] = UNewTmp

                if i == 1:
                    ## Trapezoid integration routine for the second data point
                    UTmpNode = h5file.get_node('/U', 'zPos' + str(i))
                    UTmp = UTmpNode.read()

                    fatmpNode_0 = h5file.get_node('/fa', 'zPos' + str(i - 1))
                    fatmpNode_1 = h5file.get_node('/fa', 'zPos' + str(i))
                    VTmpNode_0 = h5file.get_node('/V', 'zPos' + str(i - 1))
                    VTmpNode_1 = h5file.get_node('/V', 'zPos' + str(i))
                    fatmp_0 = fatmpNode_0.read()
                    fatmp_1 = fatmpNode_1.read()
                    VTmp_0 = VTmpNode_0.read()
                    VTmp_1 = VTmpNode_1.read()

                    UNewTmp_0 = z_step * np.dot(fatmp_0, np.conjugate(VTmp_0)) * w_step
                    UNewTmp_1 = z_step * np.dot(fatmp_1, np.conjugate(VTmp_1)) * w_step
                    UNewTmp = Uinitial + 0.5 * (UNewTmp_0 + UNewTmp_1)

                    # Calculate difference between old and new Ua
                    UdiffNorm = np.linalg.norm(UNewTmp - UTmp) / np.linalg.norm(UNewTmp - Uinitial)
                    UdiffNormMean += UdiffNorm

                    # Save new Ua
                    UTmpNode[:] = UNewTmp

            else:
                ## Trapezoid integration routine 
                UNewTmpNode = h5file.get_node('/U', 'zPos' + str(i - 1))
                UTmpNode = h5file.get_node('/U', 'zPos' + str(i))
                UNewTmp = UNewTmpNode.read()
                UTmp = UTmpNode.read()

                fatmpNode_0 = h5file.get_node('/fa', 'zPos' + str(i - 1))
                fatmpNode_1 = h5file.get_node('/fa', 'zPos' + str(i))
                VTmpNode_0 = h5file.get_node('/V', 'zPos' + str(i - 1))
                VTmpNode_1 = h5file.get_node('/V', 'zPos' + str(i))
                fatmp_0 = fatmpNode_0.read()
                fatmp_1 = fatmpNode_1.read()
                VTmp_0 = VTmpNode_0.read()
                VTmp_1 = VTmpNode_1.read()

                UNewTmp_0 = z_step * np.dot(fatmp_0, np.conjugate(VTmp_0)) * w_step
                UNewTmp_1 = z_step * np.dot(fatmp_1, np.conjugate(VTmp_1)) * w_step
                UNewTmp = UNewTmp + 0.5 * (UNewTmp_0 + UNewTmp_1)

                # Calculate Difference between old and new Ua
                UdiffNorm = np.linalg.norm(UNewTmp - UTmp) / np.linalg.norm(UNewTmp - Uinitial)
                UdiffNormMean += UdiffNorm
                # Save new Ua
                UTmpNode[:] = UNewTmp

        UdiffNormMean = UdiffNormMean / zRange.size

        if UdiffNormMean <= 1e-6 and VdiffNormMean <= 1e-6:
            break
    U = h5file.get_node('/U', 'zPos' + str(zRange.size - 1))
    V = h5file.get_node('/V', 'zPos' + str(zRange.size - 1))
    # f, (ax1, ax2) = plt.subplots(1, 2)
    # ax1.pcolormesh(wRange, wRange, np.abs(U), cmap='inferno')
    # ax2.pcolormesh(wRange, wRange, np.abs(V), cmap='inferno')
    #
    # plt.show()
    U,V = np.array(U), np.array(V)
    h5file.close()
    return U,V
