import matplotlib.pyplot as py
import numpy as np
import sys
from scipy.integrate import *


power = np.loadtxt(sys.argv[1]) 
corr = np.loadtxt(sys.argv[2]) 

k     = power[:,0]
pk    = power[:,1]
s_dat     = corr[:,0]
xi   = corr[:,1]

print "Number of k = ", len(k)
dk = (k[-1]-k[0])/(len(k)-1)
kpk = k*pk

f = lambda r: kpk*np.sin(k*r)/r

s = np.linspace(1,200,1000)

        #calculating the integral                                               
def integrand(k,x,pk):                                                  
    bessel = (np.sin(k*x))/(k*x)                    #bessel function
    return ((k**2)/(2*(np.pi)**2)) * (pk) *bessel 

xi_ft2 = np.zeros(len(s))                                                   

for i in range(len(s)):                                                 
    print "going through i:",i                                      
    I = integrand(k,s[i],pk)                                        
#    xi_ft2[i] = romb(I,dk)    


py.figure()
py.plot(k,f(5))
ksi_ft = np.zeros(len(s))
for i in range(len(s)):
    integrated =  cumtrapz(kpk*np.sin(k*s[i])/s[i], k)
    ksi_ft[i] = integrated[-1]-integrated[0]

ksi_ft /= 1./2/np.pi**2

py.figure()
#py.loglog(s,xi_ft)
py.loglog(s_dat,xi_ft2)
py.show()

