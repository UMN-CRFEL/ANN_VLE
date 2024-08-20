import numpy as np 
from utilities import *
from SpeciesData import *

def VTFlash(v,T,Z,Pres,Species,Species_names):
    P = Pres
    Ki = wilsonEquation(T,P,Species,Species_names)
    Sp_i = Species
    Sp_n = Species_names
    #theta = 0.7 #initial guess
    er = 100
    iter = 0

    while (er>1E-10):
        # print("Ki = "+str(Ki))
        # print('Species indices = ' + str(Sp_i))
    
            #step 1 - solving Rachford Rice
        if (np.sum(Ki*Z) >= 1 and np.sum(Z/Ki) >= 1 ):    
            error_RR = 10

            count = 0
            
            theta_min = 0
            theta_max = 1
            theta = 0.5*(theta_min + theta_max) #initial guess

            flag = -1
            #bisection limit corrections
            if( np.sum(Z*(Ki - 1)/(1 + theta*(Ki - 1))) > 0):
                flag = 0
                theta_min = theta
            else: 
                flag = 1
                theta_max = theta

            if(flag == 0 or flag == 1):
                while error_RR>1E-10:
                    count += 1
                    g = np.sum(Z*(Ki - 1)/(1 + theta*(Ki - 1)))
                    gd = -np.sum(Z*((Ki - 1)/(1 + theta*(Ki - 1)))**2)
                    if(g > 0):
                        theta_min = theta
                    elif(g < 0):
                        theta_max = theta
                    theta_new = theta - g/gd
                    error_RR = np.abs((theta_new - theta) / theta)
                    theta = theta_new

                    # if(theta_new > theta_min and theta_new < theta_max):
                    #     theta = theta_new
                    # else:
                    #     theta = 0.5*(theta_min + theta_max)
                   
                    if(count == 100):
                        print('theta = ',theta)
                        print("count exceeded")
                        break

                # elif(flag = 1):
                #     while error_RR>1E-10:
                #         count += 1
                #         g = np.sum(Z*(Ki - 1)/(theta_l + theta_l*(Ki - 1)))
                #         gd = -np.sum(Z*((Ki - 1)/(1 + theta*(Ki - 1)))**2)
                #         theta = theta - g/gd
                #         error_RR = np.abs(g/gd)

                # print('sigma = '+str(sigma))
            #print("count = "+str(count))
                # theta = (c[0] + sigma*c[len(Ki)-1]) / (1 + sigma)
            #print('theta = '+str(theta))
            x = Z/(1+theta*(Ki-1))
            y = x*Ki
            #print('x = '+str(x))
            #print('y = '+str(y)) 

            #step 2
            a_L, b_L = Mixture_AB(x,T,P,Sp_i,Sp_n)
            # print('l a = ' + str(a_L) )
            # print('l b = ' + str(b_L) )
            a_V, b_V = Mixture_AB(y,T,P,Sp_i,Sp_n)
            # print('v a = ' + str(a_V) )
            # print('v b = ' + str(b_V) )
            a_m, b_m = Mixture_AB(Z,T,P,Sp_i,Sp_n)
            #print('mixture a = ' + str(a_m) )
            #print('mixture b = ' + str(b_m) )

            #step 3 
            vl, vv = liq_vol(v,T,x,y,a_L,b_L,a_V,b_V,theta,Species,Species_names)
            #step 4
            R = 8.314
            Z_refl,Z_refv = Z_PR(a_m*P/(R*T)**2,b_m*P/(R*T))
            
            Zl = P*vl/(R*T)
            Zv = P*vv/(R*T)
            #print('Z l = ' + str(Zl) )
            #print('Z v = ' + str(Zv) )
            fl = fugacity(x,T,P,Sp_i,Sp_n,Zl)
            fv = fugacity(y,T,P,Sp_i,Sp_n,Zv)
            #print('fugacity liq = ' +str(fl))
            #print('fugacity vap = ' +str(fv))
            ln_phi_l = ln_fugacity_coeff(x,T,P,Sp_i,Sp_n,Zl)
            ln_phi_v = ln_fugacity_coeff(y,T,P,Sp_i,Sp_n,Zv)
            #print('ln fugacity coeff liq = ' +str(ln_phi_l))
            #print('ln fugacity coeff vap = ' +str(ln_phi_v))

            #objective function
            function = ln_phi_l - ln_phi_v - np.log(Ki)
            er = np.abs(np.sum(fl/fv - 1))
            er = np.sum(np.sqrt(function*function))

            #update step
            Ki = np.exp(ln_phi_l - ln_phi_v)

            #print("\n")

        elif(np.sum(Z*Ki) < 1):
            #print("sub cooled liquid")
            theta = 0
            x = Z
            y = Z*Ki
            y = y/np.sum(y)
            vl = v
            vv = v
            break
    
        elif (np.sum(Z/Ki) < 1):
            #print("supercritical gas")
            theta = 1
            x = Z/Ki
            y = Z 
            x = x/np.sum(x)
            vv = v
            vl = v
            break

        iter += 1

        if(iter == 1000):
            break

        # print("error = "+str(er))
        # print("\n")
        
    
    #print('total itr = '+str(iter))
    return x,y,theta,vl,vv

def wilsonEquation(T,P,Species,Species_names):
    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

    Ki = (Pc/P) * np.exp( 5.373 * (1 + omega) * (1 - Tc/T ))
    return Ki

def Mixture_AB(x,T,P,Species,Species_names):    
    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']
    
    R = 8.314
    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]
    #     print("species : "+str(Species_names[i]))
    #     print("Tr = "+str(T/Tc[i]))
    #     print("Pr = "+str(P/Pc[i]))
    #     print("alpha = "+str((1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2))
    
    # print("ai = "+str(ai))
    # print("bi = "+str(bi))

    # print("Ai = "+str(ai*P/((R*T)**2)))
    # print("Bi = "+str(bi*P/(R*T)))

    a = 0
    b = 0
    for i in range(len(Species)):
        for j in range(len(Species)):
            a += x[i]*x[j]*(1-beta[Species[i]][Species[j]])*np.sqrt(ai[i]*ai[j])
        b += x[i]*bi[i]

    return a,b

def liq_vol(v,T,x,y,al,bl,av,bv,theta,Species,Species_names):
    d1 = 1+np.sqrt(2)
    d2 = 1-np.sqrt(2)
    R = 8.314
    Al = al/(R*T)
    Av = av/(R*T)

    # a1_L = b_l*(d1+d2)-a_l/(R*T)
    # a2_L = b_l*(b_l*d1*d2 + a_l/(R*T))
    # a3_L = b_l*(d1*d2 - 1)
    # a4_L = (b_l**2)*(d1+d2 - d1*d2)
    # a5_L = -(b_l**3)*d1*d2

    # a1_V = b_v*(d1+d2)-a_v/(R*T)
    # a2_V = b_v*(b_v*d1*d2 + a_v/(R*T))
    # a3_V = b_v*(d1*d2 - 1)
    # a4_V = (b_v**2)*(d1+d2 - d1*d2)
    # a5_V = -(b_v**3)*d1*d2


    # c0 = (a2_L*a5_V - a2_V*a5_L)*theta**3 - (a5_L*a1_V - a4_V*a2_L)*v*theta**2 + (a2_L*a3_V - a5_L)*theta*v**2 + a2_L*v**3
    # c1 = (a1_L*a5_V - a5_L*a1_V + a4_V*a2_L - a2_V*a4_L)*theta**3 + (a1_L*a3_V + 3*a2_L - a4_L)*v**2*theta + (v*a1_L - 3*a2_L)*v**2 + 2*(a5_L-a2_L*a3_V)*v*theta + (a1_L*a4_V - a1_V*a4_L + 2*a2_L*a3_V - 2*a5_L)*v*theta**2 + (a5_L*a1_V - a4_V*a2_L)*theta**2
    # c2 = (a1_L*a4_V - a1_V*a4_L + a2_L*a3_V - a2_V*a3_L - a5_L + a5_V)*theta**3 + (v**2 - 3*v*a1_L + 3*a2_L)*v + (2*a1_L*a3_V - a1_V*a3_L + 3*a2_L - 2*a4_L + a4_V)*v*theta**2 + (a1_V*a4_L - a1_L*a4_V - 2*a2_L*a3_V + 2*a5_L)*theta**2 + (3*a1_L - a3_L + a3_V)*v**2*theta + 2*(a4_L - a1_L*a3_V - 3*a2_L)*v*theta + (a2_L*a3_V - a5_L)*theta
    # c3 = (a1_L*a3_V - a1_V*a3_L + a2_L - a2_V - a4_L + a4_V)*theta**3 + (3*a1_L - a1_V -2*a3_L + 2*a3_V)*v*theta**2 + (a1_V*a3_L - 2*a1_L*a3_V - a4_V - 3*a2_L + 2*a4_L)*theta**2 + (-6*a1_L + 2*a3_L - 2*a3_V)*v*theta + (2*v**2 + a1_L*a3_V + 3*a2_L - a4_L)*theta - 3*v**2 + 3*v*a1_L - a2_L
    # c4 = ( (a1_L - a1_V - a3_L + a3_V)*theta**2 + (v - 2*a1_L + a3_L - a3_V)*theta - 3*v + a1_L) * (theta - 1)
    # c5 = - (theta - 1)**2

    c0 = bl * ( Al * (v**3 + bv*v**2*theta - 3*bv**2*v*theta**2 + bv**3*theta**3 ) - bl*( v**3 + (bl + bv)*v**2*theta - (Av*bl - 2*bl*bv + 3*bv**2)*v*theta**2 + bv *(Av*bl - bl*bv + bv**2)*theta**3) )
    #m0 = bl * (-Al * (v**3 + bv*v**2*theta - 3*bv**2*v*theta**2 + bv**3*theta**3 ) + bl*( v**3 + (bl + bv)*v**2*theta - (Av*bl - 2*bl*bv + 3*bv**2)*v*theta**2 + bv *(Av*bl - bl*bv + bv**2)*theta**3) )
    c1 = Al * (-v**3 - bv*v**2*theta + 3*bv**2*v*theta**2 - bv**3*theta**3 + bl*(-1 + theta)*(3*v**2 + 2*bv*v*theta - 3*bv**2*theta**2)) + bl*( bl**2*(-1+theta)*theta*(-2*v + (Av - 2*bv)*theta) + 2*(v**3 + bv*v**2*theta - 3*bv**2*v*theta**2 + bv**3*theta**3) + bl*(3*v**2 - 3*bv*theta**2*(bv - Av*theta) + v*theta*(2*bv - 3*Av*theta + 4*bv*theta)) )
    c2 = v**3 + bv*v**2*theta - bl**3*(-1+theta)**2*theta - 3*bv**2*v*theta**2 + bv**3*theta**3 + Al*(-1 + theta)*(-3*v**2 -2*bv*v*theta + 3*bv**2*theta**2 + bl*(-1 + theta)*(3*v + bv*theta)) + bl**2*(-1 + theta)*(3*v*(1+theta) + theta*(bv - 3*Av*theta + 5*bv*theta)) + bl*(v**2*( -6 + 5*theta) + v*theta*(2*bv*(-2+theta)+Av*theta) -bv*theta**2*(Av*theta + bv*(-6 + 5*theta)))
    #m2 = -v**3 - bv*v**2*theta + bl**3*(-1 + theta)**2*theta + 3*bv**2*v*theta**2 - bv**3*theta**3 - Al*(-1 + theta)*(-3*v**2 - 2*bv*v*theta + 3*bv**2*theta**2 + bl*(-1 + theta)*(3*v + bv*theta)) - bl**2*(-1 + theta)*(3*v*(1 + theta) + theta*(bv - 3*Av*theta + 5*bv*theta)) + bl*(v**2*(6 - 5*theta) - v*theta*(2*bv*(-2 + theta) + Av*theta) + bv*theta**2*(Av*theta + bv*(-6 + 5*theta)))
    c3 = -3*v**2 - 2*bv*v*theta + 2*v**2*theta + 3*bv**2*theta**2 + Av*v*theta**2 - Av*bv*theta**3 - 2*bv**2*theta**3 + bl**2*(-1 + theta)**2*( 1 + 2*theta) + Al*(-1 + theta)**2*(-3*v + bl*(-1 + theta) - bv*theta) + bl*(-1 + theta)*(-6*v - 2*bv*theta + 4*v*theta + Av*theta**2)
    c4 = -( (-1+theta)*(3*v + Al*(-1 + theta)**2 + bv*theta - v*theta - Av*theta**2 + bv*theta**2 - bl*(2 - 3*theta + theta**2)))
    c5 = -(-1+theta)**2

    #using newton iteration to find vl
    vl = 0.99*v #initial guess
    error = 10
    count = 0

    while error > 1E-7:
        
        fx = c5*vl**5 + c4*vl**4 + c3*vl**3 + c2*vl**2 + c1*vl + c0
        fdx = 5*c5*vl**4 + 4*c4*vl**3 + 3*c3*vl**2 + 2*c2*vl + c1
        # print('Al ',Al)
        # print('Av ',Av)
        # print('bl ',bl)
        # print('bv ',bv)
        # print('v ',v)
        # print('vl ',vl)
        # print('theta ',theta)
        
        # print('c5 ',c5)
        # print('c4 ',c4)
        # print('c3 ',c3)
        # print('c2 ',c2)
        # print('m2 ',m2)
        # print('c1 ',c1)
        # print('c0 ',c0)
        # print('m0 ',m0)
        vl = vl - fx/fdx
        error = np.abs(fx/fdx)

        count += 1
        if (count == 1000):
            print("vol not converged")
            break
    
    
    
    #print('vol iterations =' +str(count))
    
    # print(error)

    #update vv from vl  
    vv = (v - (1-theta)*vl)/theta
    return vl, vv

def fugacity(x,T,P,Species,Species_names,Z):
    R = 8.314

    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']

    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    Ai = ai*P/((R*T)**2)
    Bi = bi*P/(R*T)

    Amix = 0
    Bmix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            Amix += x[i]*x[j]*(1-beta[Species[i]][Species[j]])*np.sqrt(Ai[i]*Ai[j])
        Bmix += x[i]*Bi[i]
    
    #Amix = Amix*P/((R*T)**2)
    #Bmix = Bmix*P/(R*T)
    
    #print(Z-Bmix)
    #print((Z+(1+np.sqrt(2))*Bmix)/(Z+(1-np.sqrt(2))*Bmix))
    xjAj = np.zeros(len(Species))

    for i in range (len(Species)):
        for j in range (len(Species)):
            xjAj[i] = xjAj[i] + x[j]*np.sqrt(Ai[i]*Ai[j])*(1-beta[Species[i]][Species[j]])
    
    f = P*x*np.exp( (Bi/Bmix)*(Z-1) - np.log(Z-Bmix) - ( Amix / (2*np.sqrt(2)*Bmix) ) * (2*xjAj/Amix - Bi/Bmix ) * np.log( (Z+(1+np.sqrt(2))*Bmix)/(Z+(1-np.sqrt(2))*Bmix)) )

    return f

def ln_fugacity_coeff(x,T,P,Species,Species_names,Z):
    R = 8.314

    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']

    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    Ai = ai*P/((R*T)**2)
    Bi = bi*P/(R*T)

    Amix = 0
    Bmix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            Amix += x[i]*x[j]*(1-beta[Species[i]][Species[j]])*np.sqrt(Ai[i]*Ai[j])
        Bmix += x[i]*Bi[i]
    
    #Amix = Amix*P/((R*T)**2)
    #Bmix = Bmix*P/(R*T)

    xjAj = np.zeros(len(Species))

    for i in range (len(Species)):
        for j in range (len(Species)):
            xjAj[i] = xjAj[i] + x[j]*np.sqrt(Ai[i]*Ai[j])*(1-beta[Species[i]][Species[j]])

    ln_phi =  (Bi/Bmix)*(Z-1) - np.log(Z-Bmix) - ( Amix / (2*np.sqrt(2)*Bmix) ) * ( 2*xjAj/Amix - Bi/Bmix ) * np.log( (Z+(1+np.sqrt(2))*Bmix)/(Z+(1-np.sqrt(2))*Bmix)) 

    return ln_phi

def Z_PR(a_m,b_m):
    coeff = [ 1, b_m - 1, a_m - 3*b_m**2 - 2*b_m, b_m**3 + b_m**2 - a_m*b_m ]
    Z = np.roots(coeff)
    #print("coeff = " + str(coeff))
    indices = np.argsort(Z)
    #print("Z = " + str(Z))

    return Z[indices[0]],Z[indices[-1]]

def u_cv_mix_ideal(x,T,P,v,Species,Species_names):
    R = 8.314
    Tlow = np.zeros(len(Species))
    Tcommon = np.zeros(len(Species))
    Thigh = np.zeros(len(Species))
    Cp_LowCoeffs = np.zeros((len(Species),7))
    Cp_HighCoeffs = np.zeros((len(Species),7))
    coeffs = np.zeros((len(Species),7))
    cp_ideal = np.zeros(len(Species))
    h_ideal = np.zeros(len(Species))
    cv_ideal = np.zeros(len(Species))
    u_ideal = np.zeros(len(Species))

    for i in range(len(Species)):
        Tlow[i] = SpeciesData[Species_names[i]]['thermo']['Tlow']
        Tlow[i] = SpeciesData[Species_names[i]]['thermo']['Tcommon']
        Thigh[i] = SpeciesData[Species_names[i]]['thermo']['Thigh']
        for j in range(7):
            Cp_LowCoeffs[i][j] = SpeciesData[Species_names[i]]['thermo']['lowCpCoeffs'][j]
            Cp_HighCoeffs[i][j] = SpeciesData[Species_names[i]]['thermo']['highCpCoeffs'][j]
        
        
        # if(T > Tlow[i] and T < Tcommon[i]):
        #     coeffs[i,:] = Cp_LowCoeffs[i,:]
        # elif(T >= Tcommon[i] and T < Thigh[i]):
        #     coeffs[i,:] = Cp_HighCoeffs[i,:]
        
        coeffs[i,:] = Cp_LowCoeffs[i,:]

        # print(Species_names[i])
        # print(coeffs[i])
        # print(Cp_LowCoeffs)
        # print(Cp_HighCoeffs)

        
        cp_ideal[i] = R * ( coeffs[i][0] + coeffs[i][1]*T + coeffs[i][2]*T**2 + coeffs[i][3]*T**3 + coeffs[i][4]*T**4 )
        h_ideal[i] = R * T * (coeffs[i][0] + (coeffs[i][1]/2)*T + (coeffs[i][2]/3)*T**2 + (coeffs[i][3]/4)*T**3 + (coeffs[i][4]/5)*T**4 + coeffs[i][5]/T)

    # print('T = ',T)
    # print('h ideal species = ',h_ideal)
    cv_ideal = cp_ideal - R
    u_ideal = h_ideal #- P*v*x
    
    u_ideal_mix = np.sum(h_ideal*x) #- P*v
    cv_ideal_mix = np.sum(cv_ideal*x)

    return u_ideal_mix, cv_ideal_mix

def u_cv_dev_real(x,T,P,v,Species,Species_names):
    a,b = Mixture_AB(x,T,P,Species,Species_names)
    d1 = 1+np.sqrt(2)
    d2 = 1-np.sqrt(2)
    R = 8.314
    z = P*v/(R*T)

    da_dT, d2a_dT2 = a_derivative_mix(x,T,P,Species,Species_names) 

    u_d = -((a - T*da_dT)/((d1-d2)*b)) * np.log((v + d1*b)/(v + d2*b)) + (z-1)*R*T #- P*v
    cv_d = ((T*d2a_dT2)/((d1-d2)*b)) * np.log((v + d1*b)/(v + d2*b))
    #print("inside u dev = ",(da_dT ) )

    return u_d,cv_d

def u_cv_mix_real(x,y,theta,T,P,vl,vv,Species,Species_names):

    u_ideal_l, cv_ideal_l = u_cv_mix_ideal(x,T,P,vl,Species,Species_names) 
    u_std, cv_std = u_cv_mix_ideal(x,298.15,P,vl,Species,Species_names)
    u_ideal_l = u_ideal_l - u_std

    u_ideal_v, cv_ideal_v = u_cv_mix_ideal(y,T,P,vv,Species,Species_names) 
    u_std, cv_std = u_cv_mix_ideal(y,298.15,P,vv,Species,Species_names)
    u_ideal_v = u_ideal_v - u_std

    u_d_l,cv_d_l = u_cv_dev_real(x,T,P,vl,Species,Species_names)
    u_d_v,cv_d_v = u_cv_dev_real(y,T,P,vv,Species,Species_names)
    #print("u dep = ",u_d_l,"\t",u_d_v)

    u_l = u_ideal_l + u_d_l
    u_v = u_ideal_v + u_d_v
    
    cv_l = cv_ideal_l + cv_d_l
    cv_v = cv_ideal_v + cv_d_v

    u_mix = u_v*theta + u_l*(1-theta) - P*(theta*vv + (1-theta)*vl)
    cv_mix = cv_v*theta + cv_l*(1-theta)

    return u_mix, cv_mix

def a_derivative_mix(x,T,P,Species,Species_names):
    ai = np.zeros(len(Species))
    dai_dT = np.zeros(len(Species))
    d2ai_dT2 = np.zeros(len(Species))
    beta = SpeciesData['beta']
    
    R = 8.314
    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        dai_dT[i] = omega_a *(1 + c_omega*(1 - np.sqrt(T/Tc[i]))) * (R*Tc[i])**2/Pc[i] * (-c_omega/np.sqrt(Tc[i]*T))
        d2ai_dT2[i] = omega_a * ( ( 1 + c_omega * (1 - np.sqrt(T/Tc[i])) ) * (c_omega / ( 2 * T * np.sqrt(Tc[i]*T))) + c_omega * c_omega / (2*Tc[i]*T) ) * (R*Tc[i])**2 / Pc[i]

    dadT_mix = 0
    d2adT2_mix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            dadT_mix += x[i]*x[j]*(1-beta[Species[i]][Species[j]])*0.5*(np.sqrt(ai[i]/ai[j])*dai_dT[j] + np.sqrt(ai[j]/ai[i])*dai_dT[i])
            d2adT2_mix += x[i] * x[j] * 0.25 / (ai[i]*ai[j] * np.sqrt(ai[i]*ai[j])) * ( -(ai[j]*dai_dT[i] - ai[i]*dai_dT[j])**2 + 2*ai[i]*ai[j]*( ai[i]*d2ai_dT2[j] + ai[j]*d2ai_dT2[i] ) ) * (1-beta[Species[i]][Species[j]])
    
    return dadT_mix,d2adT2_mix

def T_SinglePhaseGuess(u,Z,P,v,Species,Species_names):

    T = 350 #random guess
    error = 10
    count = 0

    while(error>1E-7):
        u_s,cv_s = u_cv_single_real(Z,T,P,v,Species,Species_names)
        #print(count,u_s,cv_s)
        Tn = T + 0.1*(u - u_s)/cv_s
        count += 1
        error = np.abs(Tn - T) / T
        T = Tn
        #print(count, 'T iter ',T)
        if(count == 1000):
            print("T guess not converged")
            break

    return T

def u_cv_single_real(x,T,P,v,Species,Species_names):

    u_ideal, cv_ideal = u_cv_mix_ideal(x,T,P,v,Species,Species_names)

    u_d,cv_d = u_cv_dev_real(x,T,P,v,Species,Species_names)

    u = u_ideal + u_d
    cv = cv_ideal + cv_d

    return u, cv

# def c(P,T,vl,vv,x,y,Z,theta,Species,Species_names):
#     rho = density(vl,vv,x,y,Z,Species,Species_names)
#     kappaS_ = kappaS(P,T,vl,vv,x,y,rho,Species,Species_names)
#     return np.sqrt(1.0/(kappaS_*rho))

# def density(vl,vv,x,y,Z,Species,Species_names):
#     Wmix_ = Wmix(Z,Species,Species_names)
#     vmix = theta*vg + (1-theta)*vl
#     rho = Wmix_/vmix
#     return rho

# def kappaS(P,T,vl,vv,x,y,Z,rho,Species,Species_names):
#     kappaT_ = kappaT(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     Cp_ = Cp(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     alphaP_ = alphaP(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     return kappaT_ - T*alphaP_*alphaP_/(rho*Cp_)

# def kappaT(P,T,vl,vv,x,y,Z,rho,Species,Species_names):
#     drhodP_ = drhodP(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     return drhodP_/rho

# def drhodP(P,T,vl,vv,x,y,Z,rho,Species,Species_names):
#     R = 8.314
#     Z_M = theta*(P*vg/(R*T)) + (1.0-theta)*(P*vl/(R*T))
#     dZdp_ = dZdP(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     Wmix_ = Wmix(Z,Species,Species_names)
#     return Wmix_/(Z_M*R*T) - P*W_M*dZdp_/(Z_M*Z_M*R*T)

# def Wmix(Z,Species,Species_names):
#     Wi = np.zeros(len(Species))
#     Wm = 0
#     for i in range(len(Species)):
#         #get species molar weights and convert to kg/mol
#         Wi[i]= SpeciesData[Species_names[i]]['W']*1E-3
#         #find phase molar weights
#         Wm += Wi[i]*Z[i]
#     return Wm

# def dZdP(P,T,vl,vv,x,y,Z,rho,Species,Species_names):
#     R = 8.314
#     dvidP_ = dvidP(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     dXdP_L_,dXdP_G_ = dXdP_P(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     dvfdP_ = dvfdP(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     Zl = P*vl/(R*T)
#     Zg = P*vg/(R*T) 
#     dZdp_L_,dZdp_G_ = dZdP_P(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     dZdxi_L_,dZdxi_G_ = dZdxi(P,T,vl,vv,x,y,Z,rho,Species,Species_names)
#     for i in range(len(Species)):
#         dZdp_L_ += dZdxi_L_[i]*dXdP_L_[i]
#         dZdp_G_ += dZdxi_G_[i]*dXdP_G_[i]
#     return dvfdP_*(Zg-Zl) + theta*dZdp_G_ + (1-theta)*dZdp_L_

# def dvidP(P,T,vl,vv,x,y,Z,rho,Species,Species_names):


# def dXdP_P(P,T,vl,vv,x,y,Z,rho,Species,Species_names):


# def dvfdP(P,T,vl,vv,x,y,Z,rho,Species,Species_names,dvidP):
#     dvfdp_ = np.sum(dvidP_)
#     return dvfdP_ 

# def dZdP_P(P,T,vl,vv,x,y,Z,rho,Species,Species_names):


# def dXdxi(P,T,vl,vv,x,y,Z,rho,Species,Species_names):









