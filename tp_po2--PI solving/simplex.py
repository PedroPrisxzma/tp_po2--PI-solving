import numpy as np
import copy
from math import floor
import scipy
import sys
import fractions
import json
import pl
import pi

#Função para escolher qual simplex aplicar
def simplex(progL):
    file_f = open('primeiro.txt', 'w')
    file_g = open('conclusao.txt', 'w')
    case = tipo_sol(progL)
    sol = np.zeros(progL.m + progL.n)
    print_step(progL, file_f)

    #primal
    if(case == 1):
        result = simplex_p(progL, file_f)
        progL_methods = copy.deepcopy(progL)
        #print to file
        if(result == 42):
#----------------------------------------------------------------------------------------------------------
            sol_int_done = 0
            if(progL.method == 0):
                pi.cut_planes(progL_methods, file_f)
                size = np.shape(progL_methods.FPI_A)
                sol_int = np.zeros(size[1])
            elif(progL.method == 1):
                vo = pi.branch_bound(progL_methods, file_f, 0)
                progL_methods.vo = vo[0]
                size = np.shape(progL_methods.FPI_A)
                sol_int = vo[1]
                sol_int_done = 1
                s_int_aux = str(sol_int)
                s_int = ''
                for h in range(len(progL_methods.i_variables * 2)):
                    s_int = s_int + s_int_aux[h]
                s_int = s_int + ']'
            #dictionary to build solution
            for i in range(len(progL.base)):
                sol[progL.base[i]] = progL.FPI_b[i]
            if(sol_int_done != 1):
                for i in range(len(progL_methods.base)):
                    sol_int[progL_methods.base[i]] = progL_methods.FPI_b[i]
                s_int = str(sol_int[:-(size[1] - len(progL_methods.i_variables))])

            #print to file
            vo_int = str(progL_methods.vo)
            s = str(sol[:-progL.m])
            vo = str(progL.vo)
            y = str(progL.FPI_y)
            print_conc(s, s_int, vo, vo_int, y, 1, file_g)
        elif(result[0] == 3):
            #case unlimited
            #dictionary to help with certificate
            sol[result[1]] = 1
            for i in range(len(progL.base)):
                sol[progL.base[i]] = progL.FPI_A[i, result[1]] * -1
            s = str(sol[:-progL.m])
            print_conc(s,-1, 0, -1, 0, 3, file_g)
    #dual
    elif(case == 2):
        result = simplex_d(progL, file_f)
        progL_methods = copy.deepcopy(progL)
        #print to file
        if(result == 42):
#----------------------------------------------------------------------------------------------------------
            sol_int_done = 0
            if(progL.method == 0):
                pi.cut_planes(progL_methods, file_f)
                size = np.shape(progL_methods.FPI_A)
                sol_int = np.zeros(size[1])
            elif(progL.method == 1):
                vo = pi.branch_bound(progL_methods, file_f, 0)
                progL_methods.vo = vo[0]
                size = np.shape(progL_methods.FPI_A)
                sol_int = vo[1]
                sol_int_done = 1
                s_int_aux = str(sol_int)
                s_int = ''
                for h in range(len(progL_methods.i_variables * 2)):
                    s_int = s_int + s_int_aux[h]
                s_int = s_int + ']'
            #dictionary to build solution
            for i in range(len(progL.base)):
                sol[progL.base[i]] = progL.FPI_b[i]
            if(sol_int_done != 1):
                for i in range(len(progL_methods.base)):
                    sol_int[progL_methods.base[i]] = progL_methods.FPI_b[i]
                s_int = str(sol_int[:-(size[1] - len(progL_methods.i_variables))])

            #print to file
            vo_int = str(progL_methods.vo)
            s = str(sol[:-progL.m])
            vo = str(progL.vo)
            y = str(progL.FPI_y)
            print_conc(s, s_int, vo, vo_int, y, 1, file_g)

        elif(result[0] == 2):
            #case unfeasible
            #dictionary to help with certificate
            sol = progL.FPI_op_matrix[result[1]]
            s = str(sol)
            print_conc(s, -1, 0, -1, 0, 2, file_g)
    #with aux
    elif(case[0] == 3):
        result = simplex_aux(progL, file_f, case[1], case[2])
        progL_methods = copy.deepcopy(progL)
        #print to file
        if(result == 42):
#----------------------------------------------------------------------------------------------------------
            sol_int_done = 0
            if(progL.method == 0):
                pi.cut_planes(progL_methods, file_f)
                size = np.shape(progL_methods.FPI_A)
                sol_int = np.zeros(size[1])
            elif(progL.method == 1):
                vo = pi.branch_bound(progL_methods, file_f, 0)
                progL_methods.vo = vo[0]
                size = np.shape(progL_methods.FPI_A)
                sol_int = vo[1]
                sol_int_done = 1
                s_int_aux = str(sol_int)
                s_int = ''
                for h in range(len(progL_methods.i_variables * 2)):
                    s_int = s_int + s_int_aux[h]
                s_int = s_int + ']'
            #dictionary to build solution
            for i in range(len(progL.base)):
                sol[progL.base[i]] = progL.FPI_b[i]
            if(sol_int_done != 1):
                for i in range(len(progL_methods.base)):
                    sol_int[progL_methods.base[i]] = progL_methods.FPI_b[i]
                s_int = str(sol_int[:-(size[1] - len(progL_methods.i_variables))])

            #print to file
            vo_int = str(progL_methods.vo)
            s = str(sol[:-progL.m])
            vo = str(progL.vo)
            y = str(progL.FPI_y)
            print_conc(s, s_int, vo, vo_int, y, 1, file_g)

        elif(result[0] == 2):
            #case unfeasible
            #dictionary to help with certificate
            sol = result[1]
            s = str(sol)
            print_conc(s,-1, 0,-1, 0, 3, file_g)
        elif(result[0] == 3):
            #case unlimited
            #dictionary to help with certificate
            sol[result[1]] = 1
            for i in range(len(progL.base)):
                if(progL.base[i] < progL.n + progL.m):
                    sol[progL.base[i]] = progL.FPI_A[i, result[1]] * -1
            s = str(sol[:-progL.m])
            print_conc(s, -1, 0, -1, 0, 3, file_g)
    file_f.close()
    file_g.close()

#Simplex Primal
def simplex_p(progL, file_f):
    while(1):
        state = pivo_primal(progL)
        if(state == 42):
            return state
            break
        elif(state[0] == 3):
            return state
            break
        elif(state[0] == 1):
            progL.divide_line(state[1], progL.FPI_A[state[1], state[2]])
            sub_0 = progL.FPI_c[0, state[2]] / progL.FPI_A[state[1], state[2]]
            progL.add_lines(state[1], state[1], sub_0 * -1, 1)
            for i in range(progL.m):
                if(i != state[1]):
                    sub_others = progL.FPI_A[i, state[2]] / progL.FPI_A[state[1], state[2]]
                    progL.add_lines(i, state[1], sub_others * -1, 0)
            print_step(progL, file_f)


#simplex Dual
def simplex_d(progL, file_f):
    while(1):
        state = pivo_dual(progL)
        if(state == 42):
            return state
            break
        elif(state[0] == 2):
            return state
            break
        elif(state[0] == 1):
            progL.divide_line(state[1], progL.FPI_A[state[1], state[2]])
            sub_0 = progL.FPI_c[0, state[2]] / progL.FPI_A[state[1], state[2]]
            progL.add_lines(state[1], state[1], sub_0 * -1, 1)
            for i in range(progL.m):
                if(i != state[1]):
                    sub_others = progL.FPI_A[i, state[2]] / progL.FPI_A[state[1], state[2]]
                    progL.add_lines(i, state[1], sub_others * -1, 0)
            print_step(progL, file_f)


#Primal with aux
def simplex_aux(progL, file_f, extra_col, aux_base):
    #creating aux matrix
    new_c = pl.make_frac_matrix(np.zeros(progL.n + progL.m + extra_col))
    for i in range(extra_col):
        new_c[0, progL.n + progL.m:] = 1/1

    new_A = pl.make_frac_matrix(np.zeros((progL.m, progL.n + progL.m + extra_col)))
    for i in range(progL.m):
        new_A[i, aux_base[i]] = 1

    for i in range(len(aux_base)):
        new_A[i, aux_base[i]] = 1

    for i in range(progL.m):
        for j in range(progL.n + progL.m):
            new_A[i, j] = progL.FPI_A[i, j]

    base = aux_base

    progL_aux = pl.PL(progL.A, progL.c, progL.FPI_b, progL.n, progL.m, progL.FPI_op_matrix, base, progL.method)
    progL_aux.make_FPI()

    progL_aux.FPI_c = new_c
    progL_aux.FPI_A = new_A

    for i in range(progL.m):
        progL_aux.FPI_c[progL.n + progL.m:] = 1
    print_step(progL_aux, file_f)
    for i in range(progL.m):
        progL_aux.add_lines(i, i, -1, 1)
    print_step(progL_aux, file_f)

    aux_result = simplex_p(progL_aux, file_f)

    check = progL_aux.vo[0] *-1
    if(progL_aux.vo[0] < 0 and check < 0.0000001):
            progL_aux.vo[0] = 0
    if(progL_aux.vo[0] < 0.0000001 and check < 0):
            progL_aux.vo[0] = 0

    if(progL_aux.vo < 0):
        return (2, progL_aux.FPI_y) #case unfeasible

    #case use aux on original PL
    elif(progL_aux.vo == 0):
        progL.FPI_A[:,:] = progL_aux.FPI_A[:, :-extra_col]
        progL.base = progL_aux.base
        print_step(progL, file_f)

        for i in range(len(progL_aux.base)):
            col = progL_aux.base[i]
            progL.divide_line(i, progL.FPI_A[i, col])
            sub_0 = progL.FPI_c[0, col] / progL.FPI_A[i, col]
            progL.add_lines(i, i, sub_0 * -1, 1)
            for j in range(progL.m):
                if(j != i):
                    sub_others = progL.FPI_A[j, col] / progL.FPI_A[i, col]
                    progL.add_lines(j, i, sub_others * -1, 0)


        aux_result = simplex_p(progL, file_f)

        return aux_result


#Define o tipo de simplex a ser aplicado
def tipo_sol(progL):
    count  = 0
    aux_base = progL.base.copy()
    if(all(i >= 0 for i in progL.FPI_b)):
        if(check_a(progL.FPI_c)):
        #    print("Case B & C positives\n")
            return 42
        else:
        #    print("Case Primal\n")
            return 1
    else:
        if(check_a(progL.FPI_c)):
        #    print("Case Dual")
            return 2
        else:
        #    print("Case Primal with aux, with b also negative")
            for x in range(len(progL.FPI_b)):
                if(progL.FPI_b[x] < 0):
                    progL.multiply_line(x, -1)
            for x in range(len(progL.FPI_b)):
                aux_base[x] = (progL.m + progL.n + count) #marca o novo local do pivo na base auxiliar
                count += 1
            return (3, progL.m, aux_base)

#Funções que escolhem o pivo
def pivo_primal(progL):
    min = 100000000000 #min between b / (possible values of pivo)
    aux_min = 0 #current lowest value
    index_min = 0 #the chosen row index
    positivo = 1 #flag for no positive values in column
    column = choose_pivo_primal(progL)
    size = np.shape(progL.FPI_A)
    if(column >= 0):
        for i in range(size[0]):
            if(progL.FPI_A[i, column] > 0):
                positivo = 1
                break
            elif(progL.FPI_A[i, column] <= 0):
                positivo = -1
        if(positivo == -1):
                return (3, column, progL.FPI_c[0, column]) #caso ilimitada
    if(column >= 0):
        for i in range(progL.m):
            if(progL.FPI_A[i, column] > 0):
                aux_min = progL.FPI_A[i, column]
                aux_min = progL.FPI_b[i] / aux_min #finding lowest b / (possible values of pivo)
                if(aux_min < min):
                    min = aux_min
                    index_min = i
        att_base(progL, index_min, column)
        return (1, index_min, column) #continua simplex
    else:
        return 42 #end of simplex

def pivo_dual(progL):
    min = 100000000000 #min between c / (possible values of pivo)
    aux_min = 0 #current lowest value
    index_min = 0 #the chosen column index
    inviavel = -1 #Viability flag
    line = choose_pivo_dual(progL)
    size = np.shape(progL.FPI_A)
    if(line >= 0):
        for i in range(size[1]):
            if(progL.FPI_A[line, i] < 0):
                aux_min = progL.FPI_A[line, i] * -1
                aux_min = progL.FPI_c[0, i] / aux_min #finding lowest c / (possible values of pivo)
                inviavel = 1
                if(aux_min < min):
                    min = aux_min
                    index_min = i
        if(inviavel > 0):
            att_base(progL, line, index_min)
            return (1, line, index_min) #continue simplex
        else: #case infeasible
            return (2, line)
    else: #line < 0, no negative
        return 42 #end of simplex


#Functions that choose the row(dual) or column(primal) of the pivo
def choose_pivo_primal(progL):
    for i in range(progL.FPI_c.shape[1]):
        if(progL.FPI_c[0, i] < 0 ):
            #column of pivo
            return i #index of coluna
    return -1 #no negatives

def choose_pivo_dual(progL):
    for i in range(progL.FPI_b.shape[0]):
        if(progL.FPI_b[i] < 0):
            #row of pivo
            return i #indice da linha
    return -1 #no negatives

#Base updating function
def att_base(progL, row, col):
    progL.base[row] = col


#array value checking function
def check_a(array):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if(array[i, j] < 0):
                return False

#Printing Function/
def print_step(progL, f):
    #prints each step on a file (primeiro.txt)
    f.write('{0}{1}{2}'.format(progL.FPI_y, progL.FPI_c, progL.vo))
    f.write("\n")
    for i in range(progL.m):
        f.write('{0}{1}{2}'.format(progL.FPI_op_matrix[i], progL.FPI_A[i], progL.FPI_b[i]))
        f.write("\n")
    f.write("\n")

#Conclusion printing
def print_conc(s, s_int, vo, vo_int, y, case, file_g):
    if(case == 1):
        #Optimal
        file_g.write("2\n")
        file_g.write(s_int)
        file_g.write("\n")
        file_g.write(vo_int)
        file_g.write("\n")
        file_g.write(s)
        file_g.write("\n")
        file_g.write(vo)
        file_g.write("\n")
        file_g.write(y)
        file_g.write("\n")
    elif(case == 2):
        #unfeasible
        file_g.write("0\n")
        file_g.write(s)
        file_g.write("\n")
    elif(case == 3):
        #unlimited
        file_g.write("1\n")
        file_g.write(s)
        file_g.write("\n")
