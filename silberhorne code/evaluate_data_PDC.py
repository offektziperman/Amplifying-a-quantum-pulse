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
""" This script copies all files from a directory containing various
analysed PDC states into a new folder structure where the individual 
files are grouped together. This is very useful to examine how 
changes in the power or grid sizes affect the solution. 

One has to supply the path to the directory holding all the analysis files
and the path to the new directory where the grouped results should be stored.
 """

import sys, os
import shutil


# Read the parameters from command line
parameters = sys.argv[1:]

input_directory= parameters[0]
output_directory= parameters[1]

# Get list of processes which have to be analyzed
data_directories = os.listdir(input_directory)


## Create folders for evaluation
# U and V matrices
os.mkdir(output_directory + "/U_and_V_matrices")
os.mkdir(output_directory + "/U_and_V_matrices/analytical_solution")
os.mkdir(output_directory + "/U_and_V_matrices/first-order_numerical_solution")
os.mkdir(output_directory + "/U_and_V_matrices/rigorous_solution")
os.mkdir(output_directory + "/U_and_V_matrices/reconstructed_rigorous_solution")

# SVD modes and amplitudes
os.mkdir(output_directory + "/SVD_modes_and_amplitudes")
os.mkdir(output_directory + "/SVD_modes_and_amplitudes/svd_amplitudes")
os.mkdir(output_directory + "/SVD_modes_and_amplitudes/svd_modes_ana")
os.mkdir(output_directory + "/SVD_modes_and_amplitudes/svd_modes_num")
os.mkdir(output_directory + "/SVD_modes_and_amplitudes/svd_modes_rig")
os.mkdir(output_directory + "/SVD_modes_and_amplitudes/svd_modes_rec")

# Check SVD
os.mkdir(output_directory + "/check_SVD")
os.mkdir(output_directory + "/check_SVD/schmidt_amplitudes_test")

os.mkdir(output_directory + "/check_SVD/ana")
os.mkdir(output_directory + "/check_SVD/ana/schmidt_modes_test_ana_0")
os.mkdir(output_directory + "/check_SVD/ana/schmidt_modes_test_ana_1")
os.mkdir(output_directory + "/check_SVD/ana/schmidt_modes_test_ana_2")
os.mkdir(output_directory + "/check_SVD/ana/schmidt_modes_test_ana_3")
os.mkdir(output_directory + "/check_SVD/ana/schmidt_modes_test_ana_4")

os.mkdir(output_directory + "/check_SVD/num")
os.mkdir(output_directory + "/check_SVD/num/schmidt_modes_test_num_0")
os.mkdir(output_directory + "/check_SVD/num/schmidt_modes_test_num_1")
os.mkdir(output_directory + "/check_SVD/num/schmidt_modes_test_num_2")
os.mkdir(output_directory + "/check_SVD/num/schmidt_modes_test_num_3")
os.mkdir(output_directory + "/check_SVD/num/schmidt_modes_test_num_4")

os.mkdir(output_directory + "/check_SVD/rig")
os.mkdir(output_directory + "/check_SVD/rig/schmidt_modes_test_rig_0")
os.mkdir(output_directory + "/check_SVD/rig/schmidt_modes_test_rig_1")
os.mkdir(output_directory + "/check_SVD/rig/schmidt_modes_test_rig_2")
os.mkdir(output_directory + "/check_SVD/rig/schmidt_modes_test_rig_3")
os.mkdir(output_directory + "/check_SVD/rig/schmidt_modes_test_rig_4")

os.mkdir(output_directory + "/check_SVD/rec")
os.mkdir(output_directory + "/check_SVD/rec/schmidt_modes_test_rec_0")
os.mkdir(output_directory + "/check_SVD/rec/schmidt_modes_test_rec_1")
os.mkdir(output_directory + "/check_SVD/rec/schmidt_modes_test_rec_2")
os.mkdir(output_directory + "/check_SVD/rec/schmidt_modes_test_rec_3")
os.mkdir(output_directory + "/check_SVD/rec/schmidt_modes_test_rec_4")

os.mkdir(output_directory + "/check_SVD/rig_utu")
os.mkdir(output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_0")
os.mkdir(output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_1")
os.mkdir(output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_2")
os.mkdir(output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_3")
os.mkdir(output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_4")

os.mkdir(output_directory + "/check_SVD/rig_uut")
os.mkdir(output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_0")
os.mkdir(output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_1")
os.mkdir(output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_2")
os.mkdir(output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_3")
os.mkdir(output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_4")

# Paper plot
os.mkdir(output_directory + "/paper_plot/")
os.mkdir(output_directory + "/paper_plot/compare_schmidt_modes_ana_rig")
os.mkdir(output_directory + "/paper_plot/squeezing_amplitudes_ana_rig")
os.mkdir(output_directory + "/paper_plot/va_abs")



# Copy figures from the individual check_data sections into evaluation folders
for data_directory in data_directories:

    ## U and V matrices
    shutil.copy(input_directory + data_directory \
            + "/U_and_V_matrices/pic/analytical_solution.png", \
            output_directory + "/U_and_V_matrices/analytical_solution/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/U_and_V_matrices/pic/first-order_numerical_solution.png", \
            output_directory + "/U_and_V_matrices/first-order_numerical_solution/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/U_and_V_matrices/pic/rigorous_solution.png", \
            output_directory + "/U_and_V_matrices/rigorous_solution/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/U_and_V_matrices/pic/reconstructed_rigorous_solution.png", \
            output_directory + "/U_and_V_matrices/reconstructed_rigorous_solution/" \
            + data_directory + ".png")


    ## SVD modes and amplitudes
    shutil.copy(input_directory + data_directory \
            + "/SVD_modes_and_amplitudes/pic/svd_amplitudes.png", \
            output_directory + "/SVD_modes_and_amplitudes/svd_amplitudes/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/SVD_modes_and_amplitudes/pic/svd_modes_ana.png", \
            output_directory + "/SVD_modes_and_amplitudes/svd_modes_ana/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/SVD_modes_and_amplitudes/pic/svd_modes_num.png", \
            output_directory + "/SVD_modes_and_amplitudes/svd_modes_num/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/SVD_modes_and_amplitudes/pic/svd_modes_rig.png", \
            output_directory + "/SVD_modes_and_amplitudes/svd_modes_rig/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/SVD_modes_and_amplitudes/pic/svd_modes_rec.png", \
            output_directory + "/SVD_modes_and_amplitudes/svd_modes_rec/" \
            + data_directory + ".png")

    # Check SVD 
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/schmidt_amplitudes_test.png", \
            output_directory + "/check_SVD/schmidt_amplitudes_test/" \
            + data_directory + ".png")

    # ana
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/ana/schmidt_modes_test_ana_0.png", \
            output_directory + "/check_SVD/ana/schmidt_modes_test_ana_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/ana/schmidt_modes_test_ana_1.png", \
            output_directory + "/check_SVD/ana/schmidt_modes_test_ana_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/ana/schmidt_modes_test_ana_2.png", \
            output_directory + "/check_SVD/ana/schmidt_modes_test_ana_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/ana/schmidt_modes_test_ana_3.png", \
            output_directory + "/check_SVD/ana/schmidt_modes_test_ana_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/ana/schmidt_modes_test_ana_4.png", \
            output_directory + "/check_SVD/ana/schmidt_modes_test_ana_4/" \
            + data_directory + ".png")


    # num
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/num/schmidt_modes_test_num_0.png", \
            output_directory + "/check_SVD/num/schmidt_modes_test_num_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/num/schmidt_modes_test_num_1.png", \
            output_directory + "/check_SVD/num/schmidt_modes_test_num_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/num/schmidt_modes_test_num_2.png", \
            output_directory + "/check_SVD/num/schmidt_modes_test_num_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/num/schmidt_modes_test_num_3.png", \
            output_directory + "/check_SVD/num/schmidt_modes_test_num_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/num/schmidt_modes_test_num_4.png", \
            output_directory + "/check_SVD/num/schmidt_modes_test_num_4/" \
            + data_directory + ".png")


    # rig
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig/schmidt_modes_test_rig_0.png", \
            output_directory + "/check_SVD/rig/schmidt_modes_test_rig_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig/schmidt_modes_test_rig_1.png", \
            output_directory + "/check_SVD/rig/schmidt_modes_test_rig_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig/schmidt_modes_test_rig_2.png", \
            output_directory + "/check_SVD/rig/schmidt_modes_test_rig_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig/schmidt_modes_test_rig_3.png", \
            output_directory + "/check_SVD/rig/schmidt_modes_test_rig_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig/schmidt_modes_test_rig_4.png", \
            output_directory + "/check_SVD/rig/schmidt_modes_test_rig_4/" \
            + data_directory + ".png")


    # rec
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rec/schmidt_modes_test_rec_0.png", \
            output_directory + "/check_SVD/rec/schmidt_modes_test_rec_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rec/schmidt_modes_test_rec_1.png", \
            output_directory + "/check_SVD/rec/schmidt_modes_test_rec_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rec/schmidt_modes_test_rec_2.png", \
            output_directory + "/check_SVD/rec/schmidt_modes_test_rec_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rec/schmidt_modes_test_rec_3.png", \
            output_directory + "/check_SVD/rec/schmidt_modes_test_rec_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rec/schmidt_modes_test_rec_4.png", \
            output_directory + "/check_SVD/rec/schmidt_modes_test_rec_4/" \
            + data_directory + ".png")


    # rig_utu
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_utu/schmidt_modes_test_rig_utu_0.png", \
            output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_utu/schmidt_modes_test_rig_utu_1.png", \
            output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_utu/schmidt_modes_test_rig_utu_2.png", \
            output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_utu/schmidt_modes_test_rig_utu_3.png", \
            output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_utu/schmidt_modes_test_rig_utu_4.png", \
            output_directory + "/check_SVD/rig_utu/schmidt_modes_test_rig_utu_4/" \
            + data_directory + ".png")


    # rig_uut
    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_uut/schmidt_modes_test_rig_uut_0.png", \
            output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_0/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_uut/schmidt_modes_test_rig_uut_1.png", \
            output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_1/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_uut/schmidt_modes_test_rig_uut_2.png", \
            output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_2/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_uut/schmidt_modes_test_rig_uut_3.png", \
            output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_3/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/check_SVD/pic/rig_uut/schmidt_modes_test_rig_uut_4.png", \
            output_directory + "/check_SVD/rig_uut/schmidt_modes_test_rig_uut_4/" \
            + data_directory + ".png")


    ## Paper plot
    shutil.copy(input_directory + data_directory \
            + "/paper_plot/compare_schmidt_modes_ana_rig.png", \
            output_directory + "/paper_plot/compare_schmidt_modes_ana_rig/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/paper_plot/squeezing_amplitudes_ana_rig.png", \
            output_directory + "/paper_plot/squeezing_amplitudes_ana_rig/" \
            + data_directory + ".png")

    shutil.copy(input_directory + data_directory \
            + "/paper_plot/va_abs.png", \
            output_directory + "/paper_plot/va_abs/" \
            + data_directory + ".png")
