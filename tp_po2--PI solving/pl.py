import numpy as np
import copy
import scipy
import sys
import fractions
import json
import pl
import pi

#The PL class
class PL:
   #Pl, constructor
   def __init__(self, A, c, b, n, m, op, base, method):

       self.A = A #matrix A
       self.c = c #vector c
       self.b = b #vector b

       #metodo de solucao plano de corte = 0, Branch and bound = 1
       self.method = method


       #matrix A dimensions
       self.n = n #columns
       self.m = m #lines

       #FPI elements for the tableu
       self.FPI_A = 0
       self.FPI_b = 0
       self.FPI_c = 0
       self.FPI_y = 0
       self.FPI_op_matrix = op
       self.vo = np.zeros(1)
       self.base = base #dicionario contendo base {'linha': 'coluna'}
       self.i_variables = [x for x in range(m)] #list of important variables (non aux variables)

   #Default prints
   def displayA(self):
       print("Matrix A\n",self.A)

   def displayc(self):
       print("Vector c\n",self.c)

   def displayb(self):
       print("Vector b\n",self.b)

   def displayDimensoes(self):
       print("linhas:",self.m,"colunas:", self.n)

   def display_base(self):
       print("Base Canonica\n", self.base)

   #Format FPI para o tableu
   def make_FPI(self):
       #Faz matriz A da FPI
       self.FPI_A = pl.make_frac_matrix(np.zeros((self.m, self.n + self.m)))
       self.FPI_A[:,:-self.m] = self.A
       self.FPI_A[:,self.n:] = pl.make_frac_matrix(np.identity(self.m))

       #Faz o vetor c
       self.FPI_c = pl.make_frac_matrix(np.zeros(self.n + self.m))
       self.FPI_c[0, :-self.m] = self.c * (-1)

       self.FPI_y = pl.make_frac_matrix(np.zeros(self.m))
       self.FPI_b = self.b



   def display_FPI(self):
       print("FPI A\n", self.FPI_A)
       print("FPI b\n", self.FPI_b)
       print("FPI c\n", self.FPI_c)
       print("FPI y\n", self.FPI_y)
       print("FPI Op\n", self.FPI_op_matrix)
       print("Valor Objetivo\n", self.vo)

   #operations on the PL
   def multiply_line(self, line, value):
       self.FPI_A[line] = self.FPI_A[line] * value
       self.FPI_b[line] = self.FPI_b[line] * value
       self.FPI_op_matrix[line] = self.FPI_op_matrix[line] * value

   def divide_line(self, line, value):
       self.FPI_A[line] = self.FPI_A[line] / value
       self.FPI_b[line] = self.FPI_b[line] / value
       self.FPI_op_matrix[line] = self.FPI_op_matrix[line] / value

   def add_lines(self, line1, line2, value, case):
       if(case == 0):
           self.FPI_A[line1] = self.FPI_A[line1] + self.FPI_A[line2] * value
           self.FPI_b[line1] = self.FPI_b[line1] + self.FPI_b[line2] * value
           self.FPI_op_matrix[line1] = self.FPI_op_matrix[line1] + self.FPI_op_matrix[line2] * value
       else:
           self.FPI_c = self.FPI_c + self.FPI_A[line2] * value
           self.FPI_y = self.FPI_y + self.FPI_op_matrix[line2] * value
           self.vo = self.vo + self.FPI_b[line2] * value

   #Makes original A, b, c matrixes equal to their FPI versions
   def make_equal(self):
       self.FPI_A = self.A
       self.FPI_b = self.b
       self.FPI_c = self.c



#Returns a tuple that has matrix A and vector b
def split_b(matrix):
    A = [line[:-1] for line in matrix]
    b = [line[-1] for line in matrix]
    return A, b

#Test function, prints all
def print_test(progL):
    progL.displayA()
    progL.displayb()
    progL.displayc()
    progL.displayDimensoes()
    progL.display_FPI()
    progL.display_base()

#make fraction matrix
def make_frac_matrix(input_pl):
    input_pl = np.matrix(input_pl)
    input_pl = input_pl.astype('object')
    for i in range(input_pl.shape[0]):
        for j in range(input_pl.shape[1]):
            input_pl[i, j] = fractions.Fraction(input_pl[i, j])
    return input_pl
