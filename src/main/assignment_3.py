# Name: Bailey Brumfield
# Date: 4/1/2023
# Professor: Juan Parra
# Asignment: Assignment 3


# import libraries
import numpy as np

## Question 1 - Euler Method with the following details: function - t - y^2; range - 0 < t < 2; iterations - 10;
##              initial point - f(0) = 1

def function_given(t, y):
    return (t - y**2)

def eulers_method(t, y, iterations, x):
    h = (x - t) / iterations

    for unused_variable in range(iterations):
        y = y + (h * function_given(t, y))
        t = t + h
   
    print(round(y, ndigits = 5), "\n")

## Question 2 - Runge-Kutta with the following details: function - t - y^2; range - 0 < t < 2; iterations - 10;
##              initial point - f(0) = 1

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

    print(round(y, ndigits = 5), "\n")

