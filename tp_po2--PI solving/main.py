import numpy as np
import scipy
import sys
import fractions
import json
import pl
import simplex
import pi

#Found at https://stackoverflow.com/questions/42209365/numpy-convert-decimals-to-fractions
#Changes decimals into fractions on numpy arrays (for screen printing)
np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})

name = sys.argv[1]
file = open(name, "r")
method = int(file.readline()) # Guarda o tipo de solução a ser usada
m = file.readline()
n = file.readline()
hue = file.readline()

input_pl = np.array(json.loads(hue), dtype=float)

#converts to fraction matrix
input_pl = np.matrix(input_pl)
input_pl = input_pl.astype('object')
for i in range(input_pl.shape[0]):
    for j in range(input_pl.shape[1]):
        input_pl[i, j] = fractions.Fraction(input_pl[i, j])

#constructs matrix A, in a numpy array format, separating from c vector
Amatrix =  input_pl[1:, :-1]
Bmatrix = input_pl[1:, int(n):]
Cmatrix = input_pl[0, :-1]
base = {}

#base dictionary contruction
j = 1
for i in range(int(m), 0, -1):
    base[i-1] = (int(n) + int(m) - j)
    j += 1

op_matrix = pl.make_frac_matrix(np.identity(int(m)))

#PL contruction
progL = pl.PL(Amatrix, Cmatrix, Bmatrix, int(n), int(m), op_matrix, base, method)

#puts it in FPI format
progL.make_FPI()

#simplex call
simplex.simplex(progL)
