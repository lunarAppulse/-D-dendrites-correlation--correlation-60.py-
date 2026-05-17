import numpy as np
import matplotlib.pyplot as plt
import math


# Model parameters
EL=-70    #mV
tau_s=16  #ms
gs=1300   #pA
Cs=370    #pF
tau_ws=100  #100ms
b=200       #pA
VT=-50    #mV

tau_d=7   #ms
gd=1200   #pA
Cd=170    #pF
tau_wd=30   #30ms
a=13        #nS=nA/V=pA/mV

Ed=-38    #mV
Dd=6      #mV

Vr=-70    #mV  reset voltage after a spike

neuronNumber=2000   #number of two-compartment neuron models

# Time resolution and duration settings
dt=0.1    #ms
time=600  #ms
N=round(time/dt)
t=np.arange(0,N+1)*dt #t=[0,0.1,0.2,...,N*dt]

# Noise parameter settings
alpha=1

# Noise-free signal input settings
td0=200    #ms
low_ratio_d=0.4
ts0=200   #ms
low_ratio_s=0.5

#f function
f = lambda x: 1/(1+math.exp((Ed-x)/Dd))

# Convolution definition
def Convolution_K_S(S,t):
    if t<=0.5:
        return 0
    if len(S)==0:
        return 0
    count=0
    for i in range(-1,-1*len(S),-1):
        if S[i]<(t-2.5):
            break
        if S[i]<(t-0.5):
            count=count+1
    return count

# Get the correlation coefficient
def get_correlation(phi, sigma, Id0, Is0):

    Id=np.zeros(N+1)
    for i in range(N+1):
        k=((i-td0/dt*phi)%(td0/dt))/(td0/dt)
        if k<=low_ratio_d:
            Id[i]=0
        else:
            Id[i]=Id0
    

    Is=np.zeros(N+1)
    for i in range(N+1):
        k=i%(ts0/dt)/(ts0/dt)
        if k<=low_ratio_s:
            Is[i]=500
        else:
            Is[i]=Is0
    
    correlation=0
    
    for counter in range(0,neuronNumber):
    
    	# Set noise
        z1=np.random.normal(0,sigma,N)
        z2=np.random.normal(0,sigma,N)
    
        noise1=np.zeros(N+1)
        for i in range(N):
    	    noise1[i+1]=noise1[i]-alpha*noise1[i]*dt+z1[i]*np.sqrt(dt)
        noise2=np.zeros(N+1)
        for i in range(N):
    	    noise2[i+1]=noise2[i]-alpha*noise2[i]*dt+z2[i]*np.sqrt(dt)
    
    
    	#Euler method
        Vs=np.zeros(N+1)
        Vs[0]=EL
        Vd=np.zeros(N+1)
        Vd[0]=EL
    
        ws=np.zeros(N+1)
        wd=np.zeros(N+1)
    
        S=[]   #spike train
        
        for i in range(1,N+1,1):
            Vd[i]=(-(Vd[i-1]-EL)/tau_d+(gd*f(Vd[i-1])+Id[i-1]+noise1[i-1]-wd[i-1]+Cd*Convolution_K_S(S,t[i-1]))/Cd)*dt+Vd[i-1]
            wd[i]=(-wd[i-1]+a*(Vd[i-1]-EL))/tau_wd*dt+wd[i-1]
            if Vs[i-1]>=VT:
                Vs[i]=Vr
                S.append(t[i-1])
                ws[i]=-ws[i-1]/tau_ws*dt+b+ws[i-1]
                continue
            Vs[i]=(-(Vs[i-1]-EL)/tau_s+(gs*f(Vd[i-1])+Is[i-1]+noise2[i-1]-ws[i-1])/Cs)*dt+Vs[i-1]
            ws[i]=-ws[i-1]/tau_ws*dt+ws[i-1]
    
        correlation += np.corrcoef(Vs, Vd)[0,1]
        
    correlation = correlation/neuronNumber
    return correlation

phi=np.linspace(0, 1, 50)    #phase shift 0~1  0.4
sigma=np.linspace(200, 450, 50)   #standard deviation for noises  350
Id0=np.linspace(300, 380, 50)   #pA  340
Is0=np.linspace(550, 700, 50)   #pA  600

variable=[phi, sigma, Id0, Is0]
labels=['Phase Shift', 'Noise Amplitude/pA', 'Dentritic Input/pA', 'Somatic Input/pA']
colors=['brown', 'blue', 'black', 'darkorange']
font = {'family' : 'Times New Roman', 'weight': 'normal','size': 50}
fig, axs = plt.subplots(2, 2, figsize=(30,30))

for i in range(4):
    args=[0.4, 350, 340, 600]
    length=len(variable[i])
    correlation=np.zeros(length)
    for j in range(length):
        args[i]=variable[i][j]
        correlation[j]=get_correlation(args[0], args[1], args[2], args[3])
        
    a=i%2
    b=math.floor((i-a)/2)
    axs[b,a].plot(variable[i], correlation, color=colors[i])
    axs[b,a].set_xlabel(labels[i], font)
    axs[b,a].set_ylabel('Correlation', font)
    axs[b,a].tick_params(labelsize=30)
plt.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig('3.png',format='png',dpi = 400)

