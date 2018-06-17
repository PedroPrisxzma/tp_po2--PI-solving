import numpy as np
import copy
from math import floor
import scipy
import sys
import fractions
import json
import pl
import simplex

#Funtion for Cutting Planes Solution
def cut_planes(progL_methods, file_f):
    c = 0
    while(check_end(progL_methods) == 0):
        for i in range(len(progL_methods.FPI_b)):
            check_floor = floor(progL_methods.FPI_b[i])
            if(check_floor != progL_methods.FPI_b[i]):
                for j in range(len(progL_methods.i_variables)):
                    if(progL_methods.base[i] == progL_methods.i_variables[j]):

                        #montando a nova PL com a nova restrição
                        new_lines = floor_line(progL_methods, i)
                        size = np.shape(progL_methods.FPI_A)
                        new_progL_A = pl.make_frac_matrix(np.zeros((size[0]+1, size[1]+1)))
                        new_progL_op_matrix = pl.make_frac_matrix(np.zeros((size[0]+1, size[0]+1)))
                        new_progL_b = pl.make_frac_matrix(np.zeros((size[0]+1, 1)))
                        new_progL_c = pl.make_frac_matrix(np.zeros((1, size[1]+1)))
                        new_progL_y = pl.make_frac_matrix(np.zeros((1, size[0]+1)))

                        new_progL_A[:-1,:-1] = progL_methods.FPI_A
                        new_progL_op_matrix[:-1,:-1] = progL_methods.FPI_op_matrix
                        new_progL_b[:-1] = progL_methods.FPI_b
                        new_progL_c[:, :-1] = progL_methods.FPI_c
                        new_progL_y[:, :-1] = progL_methods.FPI_y

                        new_progL_A[size[0], :-1] = new_lines[0]
                        new_progL_A[size[0], size[1]] = 1
                        new_progL_op_matrix[size[0], size[0]] = 1
                        new_progL_b[size[0]] = new_lines[1]

                        progL_methods.FPI_A = new_progL_A
                        progL_methods.FPI_op_matrix = new_progL_op_matrix
                        progL_methods.FPI_b = new_progL_b
                        progL_methods.FPI_c = new_progL_c
                        progL_methods.FPI_y = new_progL_y

                        add_correct_lines(progL_methods)
                        progL_methods.base[size[0]] = size[1]

                        result = simplex.simplex_d(progL_methods, file_f)
#checks for Cutting Planes end
def check_end(progL):
    for i in range(len(progL.i_variables)):
        floor_value = floor(progL.FPI_b[i])
        if(floor_value != progL.FPI_b[i]):
        #    print('go')
            return 0

#applies floor Funtion to a line
def floor_line(progL_methods, line):
    new_A = copy.deepcopy(progL_methods.FPI_A[line])
    new_b = copy.deepcopy(progL_methods.FPI_b[line])
    size = np.shape(progL_methods.FPI_A)#gets matrix dimensions in a tuple as such: (rows, columns)
    for i in range(size[1]):
        new_A[0, i] = floor(new_A[0, i])
    new_b[0] = floor(new_b[0])
    return(new_A, new_b)

#adds the correct lines to the PL, for pivo
def add_correct_lines(progL_methods):
    size = np.shape(progL_methods.FPI_A)
    for i in range(len(progL_methods.base)):
        if(progL_methods.FPI_A[size[0]-1, progL_methods.base[i]] != 0):
            sub_0 = progL_methods.FPI_A[size[0]-1, progL_methods.base[i]] / progL_methods.FPI_A[i, progL_methods.base[i]]
            progL_methods.add_lines(size[0]-1, i, sub_0 * -1 ,0)

#For branch and branch bound
def branch_bound(progL, file_f):
    return 0
