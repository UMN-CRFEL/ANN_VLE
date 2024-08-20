import numpy as np
from SpeciesData import *

def CorrectOrder(x,y,Sp_i):
    x_copy = x
    y_copy = y
    indices = np.argsort(Sp_i)

    for i in range(len(Sp_i)):
        x[i] = x_copy[indices[i]]
        y[i] = y_copy[indices[i]]
    
    return x,y

def DescendingOrder(Ki,Zi,Sp_i,Sp_n):
    #find indices of largest to smallest Ki values
    indices = np.argsort(-Ki)   #argsort returns indices in ascending, using -Ki gives descending as smallest -ve is largest +ve
    K = np.zeros(len(Sp_i))
    Z = np.zeros(len(Sp_i))
    Sp = np.zeros(len(Sp_i))
    Spn = Sp_n

    for i in range(len(Ki)):
        K[i] = Ki[indices[i]]
        Z[i] = Zi[indices[i]]
        Sp[i] = Sp_i[indices[i]]
        Spn[i] = Sp_n[indices[i]]
    
    return K,Z,Sp,Spn

# def u_VT_mix(x,y,vaporfrac,T,P):

# def cv_VT_mix(x,y,vaporfrac,T,P):

def L(u,u_mix):
    return u - u_mix

def V_PR(T,P,x,Species,Species_names):
    R = 8.314

    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']
    Pc = []
    Tc = []
    omega = []

    for i in range(len(Species)):
        Pc.append(SpeciesData[Species_names[i]]['Pc'])
        Tc.append(SpeciesData[Species_names[i]]['Tc'])
        omega.append(SpeciesData[Species_names[i]]['omega'])

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    Ai = ai*P/(R*T)**2
    Bi = bi*P/(R*T)

    Amix = 0
    Bmix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            Amix += x[i]*x[j]*(1-beta[i][j])*np.sqrt(Ai[i]*Ai[j])
        Bmix += x[i]*Bi[i]
    
    coeffs = [1, -(1-Bmix), (Amix - 2*Bmix -3*Bmix**2), -(Amix*Bmix - Bmix**2 - Bmix**3)]
    Z = np.roots(coeffs)
    v = Z*R*T/P
    for i in range(len(v)):
        if(np.imag(v[i]) == 0):
            vol = np.real(v[i])
    #print(vol)
    return vol

def P_PR(T,v,x,Species,Species_names):
    R = 8.314

    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']
    Pc = []
    Tc = []
    omega = []

    for i in range(len(Species)):
        Pc.append(SpeciesData[Species_names[i]]['Pc'])
        Tc.append(SpeciesData[Species_names[i]]['Tc'])
        omega.append(SpeciesData[Species_names[i]]['omega'])

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    #Ai = ai*P/(R*T)**2
    #Bi = bi*P/(R*T)

    amix = 0
    bmix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            amix += x[i]*x[j]*(1-beta[Species[i]][Species[j]])*np.sqrt(ai[i]*ai[j])
        bmix += x[i]*bi[i]
    
    # amix = 2.8066
    # bmix = 0.000102
    # x[0] = 0.07576
    # x[1] = 1-x[0]
    
    P = (R*T)/(v-bmix) - (amix)/(v**2 + 2*bmix*v - bmix**2)
    return P
