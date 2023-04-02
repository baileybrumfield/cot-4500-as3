# Name: Bailey Brumfield
# Date: 4/1/2023
# Professor: Juan Parra
# Asignment: Assignment 3


# import required libraries

import numpy as np


# Question 1 Function: Euler's Method

def given_function(t, y):
    return (t - y**2)

def eulers_method(t, y, iterations, x):
    h = (x - t) / iterations

    for unused_variable in range(iterations):
        y = y + (h * given_function(t, y))
        t = t + h
   
    print("%.5f" % y, "\n")


# Question 2 Function: Runge-Kutta

def func(t, y):
    return (t - y**2)

def runge_kutta(t, y, iterations, x):
    h = (x - t) / iterations
   
    for another_unused_variable in range(iterations):
        k_1 = h * func(t, y)
        k_2 = h * func((t + (h / 2)), (y + (k_1 / 2)))
        k_3 = h * func((t + (h / 2)), (y + (k_2 / 2)))
        k_4 = h * func((t + h), (y + k_3))

        y = y + (1 / 6) * (k_1 + (2 * k_2) + (2 * k_3) + k_4)

        t = t + h

    print("%.5f" % y, "\n")


# Question 3 Function: Gaussian Elimination Function

def gaussian_elimination(gaussian_matrix):
    size = gaussian_matrix.shape[0]

    for i in range(size):
        pivot = i
        while gaussian_matrix[pivot, i] == 0:
            pivot += 1
   
        gaussian_matrix[[i, pivot]] = gaussian_matrix[[pivot, i]]

        for j in range(i + 1, size):
            factor = gaussian_matrix[j, i] / gaussian_matrix[i, i]
            gaussian_matrix[j, i:] = gaussian_matrix[j, i:] - factor * gaussian_matrix[i, i:]

    inputs = np.zeros(size)

    for i in range(size - 1, -1, -1):
        inputs[i] = (gaussian_matrix[i, -1] - np.dot(gaussian_matrix[i, i: -1], inputs[i:])) / gaussian_matrix[i, i]
   
    final_answer = np.array([int(inputs[0]), int(inputs[1]), int(inputs[2])], dtype=np.double)
    print(final_answer, "\n")


# Question 4: LU Factorization Function
#           a) print out the matrix determinant
#           b) print out the L matrix
#           c) print out the U matrix

def lu_factorization(lu_matrix):
    size = lu_matrix.shape[0]

    l_factor = np.eye(size)
    u_factor = np.zeros_like(lu_matrix)

    for i in range(size):
        for j in range(i, size):
            u_factor[i, j] = (lu_matrix[i, j] - np.dot(l_factor[i, :i], u_factor[:i, j]))
   
        for j in range(i + 1, size):
            l_factor[j, i] = (lu_matrix[j, i] - np.dot(l_factor[j, :i], u_factor[:i, i])) / u_factor[i, i]
   
    determinant = np.linalg.det(lu_matrix)

    print("%.5f" % determinant, "\n")
    print(l_factor, "\n")
    print(u_factor, "\n")


# Question 5: Function to determine if the following matrix is diagonally dominate (true/false)

def diagonally_dominant(dd_matrix, n):

    for i in range(0, n):
        total = 0
        for j in range(0, n):
            total = total + abs(dd_matrix[i][j])
       
        total = total - abs(dd_matrix[i][i])
   
    if abs(dd_matrix[i][i]) < total:
        print("False\n")
    else:
        print("True\n")


