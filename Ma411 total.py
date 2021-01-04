import numpy as np
from matplotlib import pyplot as plt
import math as m

#functions and initial conditions

N=100
a=0
b=1
X=np.linspace(a,b,N)
h=(b-a)/N
alpha=8.33*10**-4 

  


def maxf (f,a,b):
    l=np.linspace(a,b,(b-a)*100)
    max=f(a)
    for i in l:
        if f(i)>max:
            max=f(i)
    return max

def f_w(x):
    return (x**2)
def f_o(x):
    return(((1-x)**2)/4)

#f=alpha*fw/(fw+fo)

def f(x):
    return ((4*alpha*(x**2))/((4*(x**2))+(1-x)**2))

def fder(x):
    return((8*alpha*x-8*(x**2))/(((4*(x**2))+(1-x)**2))**2)

#print(m.ceil(max(fder(X))))

#k=(h/m.ceil(max(fder(X))))
"""
We decide finally not to take this form of k, because it goes out of our time bounds
We take a simplified form, that is following:
"""

k=h/(3*alpha) 

T = 500


K=1*10**-4
phi=0.325
grav=9.81
d_eau=1000
d_huile=800

B=(K/phi)*(d_eau-d_huile)*grav


def g(a,b):
    if (-alpha+B*(f_w(a)))>0:
        x=alpha+(B*f_o(b))
        return ((x*f_w(a))/(f_w(a)+f_o(b)))
    else:
        x=alpha+(B*f_o(a))
        return ((x*f_w(a))/(f_w(a)+f_o(a)))

#schemes

#Without gravity

S1=np.zeros((N,T))
S1[0,:]=1
for n in range (0,T-1):
    for j in range (1,N):
        w1=(f_w(S1[j,n]))
        w2=(f_w(S1[j-1,n]))
        o1=(f_o(S1[j,n]))
        o2=(f_o(S1[j-1,n]))
        S1[j,n+1]=(k/h)*(((-alpha*w1)/(w1+o1))+((alpha*w2)/(w2+o2)))+S1[j,n]

#With gravity

S2=np.zeros((N,T))
S2[0,:]=1
for i in range (1,T-1):
    for j in range (1,N-1):
        S2[j,i+1] =(phi*g(S2[j-1,i],S2[j,i])-phi*g(S2[j,i],S2[j+1,i]))*(k/(h*phi)) + S2[j,i]


#Lax-Friedrich

S3=np.zeros((N,T))
S3[0,:]=1
for i in range (0,T-1):
    for j in range (1,N-1):
        S3[j,i+1]=0.5*((S3[j+1,i]+S3[j-1,i]))-0.5*((k/h)*(f(S3[j+1,i])-f(S3[j-1,i])))



#Godunov

S4=np.zeros((N,T))
S4[0,:]=1
V=np.zeros(N)
V[0]=1
for i in range (0,T-1):
    for j in range (1,N-1):
        V[j]=S4[j,i]-(k/h)*(f(S4[j,i])-f(S4[j-1,i]))
        S4[j,i+1]=(0.5*S4[j,i])+(0.5*(V[j]))-(0.5*(k/h)*(f(V[j])-f(V[j-1])))


#Plots

plt.plot(X,S1[:,-1],label='General without gravity')
plt.plot(X,S2[:,-1],label='General with gravity')
plt.plot(X,S3[:,-1],label='Lax-Friedrich')
plt.plot(X,S4[:,-1],color='darkgreen',label='Godunov')
plt.legend()
temp=k*T
plt.title("Water saturation after %e seconds" %(temp))
plt.show()
            
