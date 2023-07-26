import math
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from mpmath import nsum, exp, inf
from collections import OrderedDict
from scipy.optimize import fsolve
from sys import exit
import sys
def validate_interval(f, x0, x1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8):
    return f(x0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) * f(x1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) < 0

def error_bound(a, b, err):
    n = np.log((b - a) / err) / np.log(2)
    return int(np.ceil(n))
# solve for root using bisection method
def bisection(f, interval, tol, arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8):
    """
    param f: find root for function
    param interval: within range
    param tol: root accuracy tolerance
    """

    # extract interval start and end points
    x0, x1 = interval[0], interval[1]

    # check interval can be used to solve for root
    if not validate_interval(f, x0, x1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8):
        print("wrong interval")
        return

    # iterations required to find the root within a specified error bound
    n = error_bound(x0, x1, tol)

    counter = 1

    # iterate over error bound
    while True:

        # calculate root approximation
        root_approx = x0 + ((x1 - x0) / 2)

        # evaluate y at current estimate
        y = f(root_approx, arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)

        # check tolerance condition
        if -tol < y < tol:
            # return root approximation
            return root_approx

        # check if next segment is left of bisection
        if validate_interval(f, x0, root_approx, arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8):
            x1 = root_approx
        else:
            x0 = root_approx

        # increment counter
        counter += 1
#Global Variables
nm = 1e-9
k=1.3806488*10**-23
e = 1.60217657*10**-19
h=6.62606957e-34
me=9.11e-31
h_=h/(2*np.pi)
AA=1e-10
Tc=10.1
Delta=2*k*Tc
Rsh=280
d=8*nm
w=250*nm
S=w*d
rho0=Rsh*d
N0DFT=0.9/(e*21*AA**3)
k_F=np.pi/(4*AA)
v_F=h_*k_F/me
n=(2./3)*N0DFT*v_F**2*me
k_F=(3*np.pi**2*n)**(1./3)
v_F=h_*k_F/me
E_F=h_**2*k_F**2/(2*me)
b=3*me/(2*e**2*N0DFT*h_*rho0)
Q=(-b+np.sqrt(b**2+4))/2
l=1/(k_F*Q)
G=h_*v_F/l
N0Q=N0DFT*(1-Q**2)
DQ=v_F*1/(3*k_F*Q)
Zeta_3 = 1.202 # constant from article
# Critical temperatur in [K]
T1= 1
T2 =Tc
T = np.linspace(T1,T2,1000) #
t = T / Tc
Delta_zero = 2*k*Tc
tau = (l / v_F) # 1e+30
N0 =N0Q # Density of state on Fermy surface ???????
#Parametrs for simulations#
N_matsubara_freq = int(sys.argv[1])
Amount_of_points = int(sys.argv[2])
name_of_csv = str(sys.argv[3])

#At first let's build Igl(T) which for current  in range of temperatures closed to Tc. We will need this dependece to plot I(T) for arbitary temperatures #
rho = 1 / (2*tau*(math.pi)*k*Tc / h_ )
Summ = 0
for n in range(1000):
    Summ = Summ + ((2*n+1)**(-2))*((2*n + 1 + rho0)**(-1))
print((8 / (7*Zeta_3))*Summ)
Hi = (8 / (7*Zeta_3))*nsum(lambda n: ((2*n+1)**(-2))*((2*n + 1 + rho)**(-1)), [0, inf])
#print(Hi)
Igl =(16 / (9*(7*Zeta_3)**(0.5))) * (Hi**(0.5))* e * N0*v_F * (math.pi)*k*Tc*((1 - T/Tc)**(3/2)) #(16 / (9*(7*Zeta_3)**(0.5))) * (Hi**(0.5)) * e * N0 * v_F * (math.pi)*k*Tc*((1 -
#plt.plot(T,Igl)
#plt.show()
#u = 0.1#Igl[996]/ (e*N0)
print('Critical currents critical temperature',Igl[998])
#print('Condesate velocity', Delta_zero /k_F )
#print('Vc / Vf',  (Delta_zero /k_F) / v_F  )
#print('D / (2Ef)',  (Delta_zero ) / E_F  )
#print('another Vc', ((Delta_zero ) / E_F)*v_F )
#print('condensate velocity',u)
#print('fermi velocity',v_F)
#print("velocity", Delta_zero / k_F)
#sys.exit(0)
T = 0.9*Tc
def X_finding(x_parametr, order_parametr, cp_velocity, Temperature, mfp_time, v_fermi, h_plank , k_boltzman,freq_number):
     x = x_parametr
     D = order_parametr
     u = cp_velocity
     T = Temperature
     tau = mfp_time
     v_F =  v_fermi
     h_ = h_plank
     k = k_boltzman
     p = freq_number
     #We solve dimensionless equation thus cut them
     u_prime = (me*u*v_F) / (math.pi * k * T)
     tau_prime = (tau / h_)*math.pi * T
     Z =u_prime*(1 / ( (2*p+1) / x +0.5*( 1 / tau_prime  )  )) # h_ / tau
     eq =(((u_prime) / 2 ))*(((1-x**2)*(1 + Z**2))**(0.5))*(1 / (Z)) -  D*(1 -   ((tau_prime*u_prime)**(-1))*math.atan(Z)    )**(-1)  #(((v_F*u*me) / 2 )**2)*((((1-x**2)*(1+ Z**2)   ))**(1))*(1 / (Z**2)) - (D**2)*(( 1 / (1 - ( 1 / (tau*v_F*u*me / h_ ) )*(math.atan(Z))    )  )**2)
     return float(eq)
#print('test',X_finding(1, Delta_zero,u,T,tau,v_F,h_,k,100 #))
#u = ((Delta_zero ) / E_F)*v_F
#values = np.linspace(0.01,1,100)
#function =np.array([X_finding(i,Delta_zero,u,T,tau,v_F,h_,k,1) for i in values ])
#plt.plot(values,function)
#plt.show()
#Vt = fsolve(X_finding, 0.1, args=(Delta_zero,u,T,tau,v_F,h_,k,1))
#B = bisection(X_finding, [0.001, 1], 1e-20,Delta_zero,u,T,tau,v_F,h_,k,1)
#print(B)
#exit(0)
#print("test",float(Vt))
def Delta_finding(order_parametr, cp_velocity, Temperature, mfp_time, v_fermi, h_plank , k_boltzman):
    D = order_parametr
    u = cp_velocity
    T = Temperature
    tau = mfp_time
    v_F =  v_fermi
    h_ = h_plank
    k = k_boltzman
    #We solve dimensionless equation thus cut them
    u_prime = (me*u*v_F) / (math.pi * k * T)
    tau_prime = (tau / h_)*math.pi * T
    Sum = 0
    for p in range(N_matsubara_freq):
        x = float(bisection(X_finding, [0.001, 1], 1e-30,Delta_zero,u,T,tau,v_F,h_,k,p)) #fsolve(X_finding, 0.1, args=(D,u,T,tau,v_F,h_,k,p)) #   fsolve(X_finding, 0, args=(D,u)
        z =me*0.5*v_F*u*(((2*p+1)*k*T) / x +0.5*(h_ / tau ) )**(-1)           #(1 / (((2*p+1)*k*T) / (x) +0.5*(h_ / tau )  ))
        y = ((v_F*u*me) / 2 )*(((1-x**2)*(1 + z**2))**(0.5))*(1 / (z))        #(( 1 / (1 - ( 1 / (me*tau*v_F*u / h_ ) )*(math.atan(z))    )  )) # ((v_F*u) / 2 )*(((1-x**2)*(1+(1 / (z**2)))   )**(0.5))*(1 / z)i
        #y = ((v_F*u*me) / 2 )*(((1-x**2)*(1 + z**2))**(0.5))*(1 / (z))
       # print('diff',(D)*(1 -   ((me*tau*v_F*u / h_ )**(-1))*math.atan(z)    )**(-1) - ((v_F*u*me) / 2 )*(((1-x**2)*(1 + z**2))**(0.5))*(1 / (z)))
        Sum = Sum + (D / ((2*p+1)*k*T)) - (2*y / (me*u*v_F) )*math.atan(z)
    #    print((D / ((2*p+1)*k*T)))
   #     print((2*y / (u*v_F) )*math.atan(z))i
    #print('first term',D * math.log(T / Tc))
    #print('second term',2*math.pi*k*T*Sum)
    eq = D * math.log(T / Tc) + 2*math.pi*k*T*float(Sum) # D*np.log(T / Tc )
    return float(eq)
#values = np.linspace(0,Delta_zero,100)
#function =np.array([Delta_finding(i,u,0.1*Tc,tau,v_F,h_,k) for i in values ])
#plt.plot(values,function)
#plt.show()
#exit(0)
#Delta =float( fsolve(Delta_finding, Delta_zero, args=(u,T,tau,v_F,h_,k)))
#print(Delta)
def Current_finding(order_parametr, cp_velocity, Temperature, mfp_time, v_fermi, h_plank , k_boltzman):
    D = order_parametr
    u = cp_velocity
    T = Temperature
    tau = mfp_time
    v_F =  v_fermi
    h_ = h_plank
    k = k_boltzman
    Sum = 0
    for p in range(N_matsubara_freq): #Amount of Macubara freqiencies
        x =float(bisection(X_finding, [0.001, 1], 1e-30,Delta_zero,u,T,tau,v_F,h_,k,p)) #   fsolve(X_finding, 0, args=(D,u)
        z =me* 0.5*v_F*u*(1 / (((2*p+1)*k*T) / (x) +0.5*(h_ / tau )  ))
        y =  (D)*(( 1 / (1 - ( 1 / (me*tau*v_F*u / h_ ) )*(math.atan(z))    )  )) # ((v_F*u) / 2 )*(((1-x**2)*(1+(1 / (z**2)))   )**(0.5))*(1 / z)i
        print('diff',(D)*(1 -   ((me*tau*v_F*u / h_ )**(-1))*math.atan(z)    )**(-1) - ((v_F*u*me) / 2 )*(((1-x**2)*(1 + z**2))**(0.5))*(1 / (z)))
        Sum = Sum + 2 * ((y / (me*u*v_F) )**2)*(math.atan(z)-(z +  z**(-1))**(-1))
    #    print((D / ((2*p+1)*k*T)))
   #     print((2*y / (u*v_F) )*math.atan(z))
    eq = 4*e*N0*v_F*math.pi*T*k*float(Sum) # D*np.log(T / Tc )
    return float(eq)
#condesate_v = np.linspace(1,100000,1000)
#Delta =float( fsolve(Delta_finding, Delta_zero, args=(u,0.9*Tc,tau,v_F,h_,k)))
#print('Delta',Delta)
#exit(0)
#current =np.array([Current_finding(Delta,i,0.9*Tc,tau,v_F,h_,k) for i in condesate_v  ])
#plt.plot(condesate_v,current)
#plt.show()
#exit(0)
Temperatures = np.linspace(0,Tc,Amount_of_points)
Current_temperature = np.zeros(len(Temperatures)) # this we need
condesate_v = np.linspace(1,80000 ,100)
Current_velocity =np.zeros(100)  # This we will maximaze
Current_data = pd.DataFrame(columns = ['T', 'I'])
Current_data.to_csv(name_of_csv+'.csv')
for i in range(1,len(Temperatures)):
    for j in range(1,len(condesate_v)):
        Delta = float( fsolve(Delta_finding, Delta_zero, args=(condesate_v[j],Temperatures[i],tau,v_F,h_,k)))
       # print('Delta')
        Current_velocity[j] = Current_finding(Delta,condesate_v[j],T,tau,v_F,h_,k)
       # print('Current',Current_velocity[j])
       # print('velocity',condesate_v[j])
    Current_temperature[i] = np.amax(Current_velocity)
    I_critical = pd.DataFrame({'T' : [Temperatures[i]], 'I' : [Current_temperature[i]]})
    I_critical.to_csv(name_of_csv+'.csv', mode='a', index=False, header=False)
    #print('Current',Current_temperature[i])
    #print('Temperature',Temperatures[i] )
plt.plot(Temperatures,Current_temperature)
plt.show()
#print(Current_finding(Delta,u,T,tau,v_F,h_,k))
#condesate_v = np.linspace(0.01,200,1000)
#current =np.array([Current_finding(Delta,i,T,tau,v_F,h_,k) for i in condesate_v  ])
#plt.plot(condesate_v,current)
#plt.show()
#print('Calculated',np.amax(current))
#y = (Delta)*(( 1 / (1 - ( 1 / (me*tau*v_F*u / h_ ) )*(math.atan(z))    )  ))
#z =me* 0.5*v_F*u*(1 / (((2*1+1)*k*T) / (x) +0.5*(h_ / tau )  ))
#print(z)
#print(y)