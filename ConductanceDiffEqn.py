# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:08:55 2017

@author: Joshua Hibbard
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
from scipy.optimize import curve_fit
#import bigfloat as bf

#Leaf Constants
X = 1
m = 2
lambda_g = 0.01
lambda_e = 0.02
w_o = 657959000
T_w = 4982.85
sigma = 0.05
R = 8.314
c_g = 3088
v_w = 1.8054197671114737e-05
const = 2156
const2 =0.5e6
c_e = 0.4


#Data
data = pd.read_csv('Conductance1.csv')
data1 = np.array(data)
time = data1[:,[0]]
time = time - time[0]
conductance = data1[:,[1]]
#plot.plot(time,conductance)
#plot.figure()


plant_data = pd.read_csv('Data1.csv')
plant_data = np.array(plant_data)
plant_time = plant_data[:,0]
plant_time = plant_time - plant_time[0]
plant_w_a = plant_data[:,1]
plant_T_e = plant_data[:,2]


std = 0.5*np.random.random(size=plant_time.size)
#initial_guess = np.array([0.0,0.0,0.0,0.0,2100,3088])
initial_guess = 0., 0., 0., 0., 2100., 3088.

def model_g_s(independent_variables,lambda_g,lambda_e,const2,c_e,const,c_g):
    t, w_a, T_e = independent_variables
    return X*(P_g(t,w_a,T_e) - m*P_e(t,w_a,T_e))
    
#print(curve_fit(model_g_s, (plant_time, plant_w_a, plant_T_e), conductance, initial_guess, std))

print(g_s(plant_time,plant_w_a,plant_T_e))

ww=g_s(plant_time,plant_w_a,plant_T_e)


#Theoretical Conductance
def g_s(t,w_a,T_e):
    g_cell_pressure = P_g(t,w_a,T_e)
    e_cell_pressure = P_e(t,w_a,T_e)
    return X*(g_cell_pressure - m*e_cell_pressure)

#Saturated mole water fraction for a given leaf temperature. Units: mmol air/mol water
def w_i(T_e):
    return w_o*np.exp(-T_w/(T_e+273))

#Mole water fraction in the pore. Units: mmol air/mol water
def w_p(w_a,T_e):
    return sigma*w_a + (1-sigma)*w_i(T_e)
    
#Pressure in the guard cell. Units: MPa
def P_g(t,w_a,T_e):
    return (((R*(T_e+273))/v_w)*np.log((w_p(w_a,T_e))/(w_i(T_e))) + c_g*R*(T_e+273) + (const)*np.exp(-lambda_g*t))*(10e-6)


#Functions for epidermal cell pressure function
def d(T_e):
    return c_e*R*T_e
def b(w_a,T_e):
    return ((R*(T_e+273))/v_w)*np.log((sigma*w_a+(1-sigma)*w_o*np.exp(-T_w/(T_e+273)))/(w_o*np.exp(-T_w/(T_e+273))))
def c(T_e):
    return c_g*R*T_e
def f(w_a,T_e):
    return R*X*(w_a-w_o*np.exp(-T_w/(T_e+273)))
def l(w_a,T_e):
    return (1+R*X*m*(w_a-w_o*np.exp(-T_w/(T_e+273))))

#Pressure in the epidermal cells. Units: MPa
def P_e(t,w_a,T_e):
    #return (np.exp(l(w_a,T_e)*lambda_e*t))*(-np.exp(-t*(lambda_g + l(w_a,T_e)*lambda_e))*lambda_e*((-np.exp(lambda_g*t)*(d(T_e) + (b(w_a,T_e) + c(T_e))*f(w_a,T_e)))/(l(w_a,T_e)*lambda_e) + ((-5e6 + b(w_a,T_e) + c(T_e))*f(w_a,T_e))/(lambda_g + l(w_a,T_e)*lambda_e)) + const2) 
    return ((np.exp(l(w_a,T_e)*lambda_e*t))*(-np.exp(-t*(lambda_g + l(w_a,T_e)*lambda_e))*lambda_e*((-np.exp(lambda_g*t)*(d(T_e) + (b(w_a,T_e) + c(T_e))*f(w_a,T_e)))/(l(w_a,T_e)*lambda_e) + ((-5e6 + b(w_a,T_e) + c(T_e))*f(w_a,T_e))/(lambda_g + l(w_a,T_e)*lambda_e)) + const2))*(10e-6)
