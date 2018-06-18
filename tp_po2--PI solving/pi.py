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
def branch_bound(progL, file_f, case):
    vo_floor = 0
    vo_celling = 0
    if(case == 0):#do all floors first
        vo_floor = do_floor(progL, file_f)
        vo_celing = do_celling(progL, file_f)
    elif(case == 1):#do all cellings first
        vo_celing = do_celling(progL, file_f)
        vo_floor = do_floor(progL, file_f)

    #return the highest int solution found and the corresponding solution vector
    return (vo_floor if vo_floor[0] > vo_celing[0] else vo_celing)

#do floor
def do_floor(progL, file_f):
    size = np.shape(progL.FPI_A)
    branch_floor_PL = copy.deepcopy(progL)
    if(check_end_bb(branch_floor_PL)):
        index = choose_var(branch_floor_PL)
        bb_restrict(branch_floor_PL, 0, index, size) #PL para restrição var < int
        add_correct_lines(branch_floor_PL)
        branch_floor_PL.base[size[0]] = size[1]

        if(check_b_negative(branch_floor_PL, size)):
            result_floor = simplex.simplex_d(branch_floor_PL, file_f)
        else:#case cant even apply simplex
            result_floor = 41
            vo_f = -1
            sol_int = np.zeros(size[1])
            for i in range(len(progL.base)):
                sol_int[progL.base[i]] = progL.FPI_b[i]

        if(result_floor == 42):#case simplex viable
            vo_floor = branch_bound(branch_floor_PL, file_f, 0)
            vo_f = vo_floor[0]
            sol_int = vo_floor[1]
        elif(result_floor != 42 and result_floor != 41):#case simplex not viable
            vo_f = -1
            sol_int = np.zeros(size[1])
            for i in range(len(progL.base)):
                sol_int[progL.base[i]] = progL.FPI_b[i]

    else:
        vo_f = progL.vo
        sol_int = np.zeros(size[1])
        for i in range(len(progL.base)):
            sol_int[progL.base[i]] = progL.FPI_b[i]

    return (vo_f, sol_int)

#do celling
def do_celling(progL, file_f):
    size = np.shape(progL.FPI_A)
    branch_celling_PL = copy.deepcopy(progL)
    if(check_end_bb(branch_celling_PL)):
        index = choose_var(branch_celling_PL)
        bb_restrict(branch_celling_PL, 1, index, size) #PL para restrição var > int
        branch_celling_PL.multiply_line(size[0], -1)
        add_correct_lines(branch_celling_PL)

        branch_celling_PL.base[size[0]] = size[1]

        if(check_b_negative(branch_celling_PL, size)):
            result_celing = simplex.simplex_d(branch_celling_PL, file_f)
        else:#case cant even apply simplex
            result_celing = -1
            vo_c = -1
            sol_int = np.zeros(size[1])
            for i in range(len(progL.base)):
                sol_int[progL.base[i]] = progL.FPI_b[i]

        if(result_celing == 42):#case simplex viable
            vo_celing = branch_bound(branch_celling_PL, file_f, 1)
            vo_c = vo_celing[0]
            sol_int = vo_celing[1]
        elif(result_celing != 42 and result_celing != 41):#case simplex not viable
            vo_c = -1
            sol_int = np.zeros(size[1])
            for i in range(len(progL.base)):
                sol_int[progL.base[i]] = progL.FPI_b[i]
    else:
        vo_c = progL.vo
        sol_int = np.zeros(size[1])
        for i in range(len(progL.base)):
            sol_int[progL.base[i]] = progL.FPI_b[i]

    return (vo_c, sol_int)

#checks end of b&b
def check_end_bb(progL):
    for i in range(len(progL.i_variables)):
        floor_value = floor(progL.FPI_b[i])
        if(floor_value != progL.FPI_b[i]):
            return 1

#checks if simplex dual is applicable
def check_b_negative(progL, size):
     for i in range(size[0]+1):
        if(progL.FPI_b[i] < 0):
            return 1
     return 0

#chooses var to make b&b restrictions
def choose_var(progL):
    for i in range(len(progL.i_variables)):
        floor_value = floor(progL.FPI_b[i])
        if(floor_value != progL.FPI_b[i]):
            return i


#Funtion used for the B&B extra restrictions
def bb_restrict(progL, case, line, size):
    floor_value = floor(progL.FPI_b[line])
    new_line_A = pl.make_frac_matrix(np.zeros((1, size[1]+1))) #nova linha da restrição
    if(case == 0): #case Função chão
        new_line_A[0 ,line] = 1
        new_line_A[0, size[1]] = 1
        new_line_b = floor_value
    else: #case função teto, ou chão + 1
        new_line_A[0, line] = 1
        new_line_A[0, size[1]] = -1
        new_line_b = floor_value+1

    #monta nova PL
    new_progL_A = pl.make_frac_matrix(np.zeros((size[0]+1, size[1]+1)))
    new_progL_op_matrix = pl.make_frac_matrix(np.zeros((size[0]+1, size[0]+1)))
    new_progL_b = pl.make_frac_matrix(np.zeros((size[0]+1, 1)))
    new_progL_c = pl.make_frac_matrix(np.zeros((1, size[1]+1)))
    new_progL_y = pl.make_frac_matrix(np.zeros((1, size[0]+1)))

    new_progL_A[:-1,:-1] = progL.FPI_A
    new_progL_op_matrix[:-1,:-1] = progL.FPI_op_matrix
    new_progL_b[:-1] = progL.FPI_b
    new_progL_c[:, :-1] = progL.FPI_c
    new_progL_y[:, :-1] = progL.FPI_y

    new_progL_A[size[0], :] = new_line_A
    new_progL_op_matrix[size[0], size[0]] = 1
    new_progL_b[size[0]] = new_line_b

    progL.FPI_A = new_progL_A
    progL.FPI_op_matrix = new_progL_op_matrix
    progL.FPI_b = new_progL_b
    progL.FPI_c = new_progL_c
    progL.FPI_y = new_progL_y
