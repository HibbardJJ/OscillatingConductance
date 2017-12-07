# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:08:55 2017

@author: Joshua Hibbard
"""
import numpy as np

#Leaf Constants
X = 1
m = 2
lambda_g = 50
lambda_e = 60
w_o = 657959000
T_w = 4982.85
sigma = 0.05
R = 8.314
c_g = 3088
v_w = 1.8054197671114737e-05
const = 2156
const2 =0.5e6
c_e = 0.4
c_g_values=[]

#Theoretical Conductance
def g_s(t,w_a,T_e):
    g_cell_pressure = P_g(t,w_a,T_e)
    e_cell_pressure = P_e(t,w_a,T_e)
    return X*(g_cell_pressure - m*e_cell_pressure)
#Saturated mole water fraction for a given leaf temperature
"""def w_i(T_e):
    return w_o*np.exp(-T_w/(T_e+271))"""
#Mole water fraction in the pore
"""def w_p(w_a,w_i,*T_e):
    return sigma*w_a + (1-sigma)*w_i(*T_e)"""
#Pressure in the guard cell
def P_g(t,w_a,T_e):
    global c_g
    p_g = ((R*(T_e+273))/v_w)*np.log((sigma*w_a+(1-sigma)*w_o*np.exp(-T_w/(T_e+273)))/(w_o*np.exp(-T_w/(T_e+273)))) + c_g*R*(T_e+273) + (const)*np.exp(-lambda_g*t)
    if p_g < 0:
        c_g+=1
        P_g(t,w_a,T_e)
    else:
        return p_g
        c_g_values.append(c_g)

#Pressure in the epidermal cells
def P_e(t,w_a,T_e):
    return (c_e*R*(T_e+273))/(1-R*X*m*(w_a - w_o*np.exp(-T_w/(T_e+273)))) - (lambda_g*np.exp(1)/lambda_e)*(1/(1-R*X*m*(w_a - w_o*np.exp(-T_w/(T_e+273)))))*((R*(T_e+273)/v_w)*np.log((sigma*w_a+(1-sigma)*w_o*np.exp(-T_w/(T_e+273)))/(w_o*np.exp(-T_w/(T_e+273)))) - c_g*R*(T_e+273)) + (const*lambda_e/(1-R*X*m*(w_a - w_o*np.exp(-T_w/(T_e+273)))*lambda_e*lambda_g-lambda_g**2))*np.exp((1-lambda_g*t)) + const2*np.exp(-(1-R*X*m*(w_a - w_o*np.exp(-T_w/(T_e+273)))*(lambda_e*t)))

def P_g1(t,w_a,T_e):
    return ((R*(T_e+273))/v_w)*np.log((sigma*w_a+(1-sigma)*w_o*np.exp(-T_w/(T_e+273)))/(w_o*np.exp(-T_w/(T_e+273))))


x=P_g1(1,20,23)
y=P_g(1,19,24)
print(y)
