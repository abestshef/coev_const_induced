
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from scipy.integrate import solve_ivp

# h evoles for induced (causes gamma-alpha trade-off)
# c evolves for costitutive (trade-off with a)
# p evolves for parasite(causes beta-alpha trade-off)

def adaptivedynamics(t,V):
    
    h=V[0]
    c=V[1]
    p=V[2]

    gh=1-(h1*h1)/h2*(1-np.exp(h2*(h-hs)/h1))
    a=10-(bet1*bet1)/bet2*(1-np.exp(bet2*(c-bets)/bet1))
    bet=1-(p1*p1)/p2*(1-np.exp(p2*(p-ps)/p1))

    dgh=h1*np.exp(h2*(h-hs)/h1)
    da=bet1*np.exp(bet2*(c-bets)/bet1)
    dbet=p1*np.exp(p2*(p-ps)/p1)

    alpha=ALP0+gh*p
    gamma=GAM0+h

    xx=(alpha+B+gamma)/(bet*(BETA-c)+0.5)
    a1=-Q*FF
    b1=a*FF-Q*xx-Q*FF*xx-(bet*(BETA-c)+0.5)*xx+gamma
    c1=a*xx-Q*xx**2-B*xx
    yy=(-b1-np.sqrt(b1**2-4*a1*c1))/(2*a1)
    #yy=(a-Q*xx-B)*xx/(Q*xx+(BETA-c)*xx-gamma)
    
    dc = (a-Q*(xx+yy)-B-(bet*(BETA-c)+0.5)*yy)*(p*dgh+1)+(bet*(BETA-c)+0.5)*yy
    di = (da+bet*yy)*(alpha+B+gamma)-(gamma+FF*(a-Q*(xx+yy)))*bet*yy+da*FF*(bet*(BETA-c)+0.5)*yy
    dp = mut*(dbet*(BETA-c)*xx-gh)
    
    return [dc,di,dp]

X0=[1,1,1]

hstore=[]
cstore=[]
pstore=[]

Xstore=[]
Ystore=[]
Hstore=[]
prevstore=[]

for i in range (51):

    B=0.3+0.1*i
    Q=0.2
    ALP0=1.0
    GAM0=1.0
    FF=0.001
    BETA=2.0
    mut=1.0

    hs=1
    h1=1.5
    h2=2.5

    bets=1
    bet1=-2
    bet2=-0.5 

    ps=1
    p1=0.5
    p2=-0.4

    XX=solve_ivp(adaptivedynamics,[0,100],X0)
    
   # X0=[XX.y[0,-1],XX.y[1,-1],XX.y[2,-1]]
    h_s=XX.y[0,-1]
    c_s=XX.y[1,-1]
    p_s=XX.y[2,-1]
    ghs=1-(h1*h1)/h2*(1-np.exp(h2*(h_s-hs)/h1))
    a_s=10-(bet1*bet1)/bet2*(1-np.exp(bet2*(c_s-bets)/bet1))
    bets=1-(p1*p1)/p2*(1-np.exp(p2*(p_s-ps)/p1))
    alphas=ALP0+ghs*p_s
    gammas=GAM0+h_s
    
    xeq=(alphas+B+gammas)/(bets*(BETA-c_s)+0.5)
    a1=-Q*FF
    b1=a_s*FF-Q*xeq-Q*FF*xeq-(bets*(BETA-c_s)+0.5)*xeq+gammas
    c1=a_s*xeq-Q*xeq**2-B*xeq
    yeq=(-b1-np.sqrt(b1**2-4*a1*c1))/(2*a1)
        
    hstore.append(h_s)
    cstore.append(c_s)
    pstore.append(p_s)
    
    Xstore.append(xeq)
    Ystore.append(yeq)
    Hstore.append(xeq+yeq)
    prevstore.append(yeq/(xeq+yeq))
    
plt.rc('font', size=16) #controls default text size
q=np.linspace(0.3,5.3,51)

f, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(15,5),sharey=True)
ax1.plot(q,hstore)
ax1.set_title('Induced')
ax1.set_ylim([0,2])
ax1.set_xlim([0,5])
ax1.set_xlabel('Mortality, $b$')
ax1.set_ylabel('Co-CSS Investment')
#ax1.axvspan(0, 0.18, color='k', alpha=0.2, lw=0)

ax2.plot(q,cstore)
ax2.set_title('Constitutive')
ax2.set_xlim([0,5])
ax2.set_xlabel('Mortality, $b$')
#ax2.axvspan(0, 0.18, color='k', alpha=0.2, lw=0)

ax3.plot(q,pstore)
ax3.set_title('Parasite')
ax3.set_xlim([0,5])
ax3.set_xlabel('Mortality, $b$')
#ax3.axvspan(0, 0.18, color='k', alpha=0.2, lw=0)

#plt.savefig('compf1_SI.png')

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
ax1.plot(q,Xstore,label='Susceptible',linestyle='-')
ax1.plot(q,Ystore,label='Infected',linestyle=':')
ax1.plot(q,Hstore,label='Total',linestyle='--')
#ax1.set_title('Population densities')
ax1.set_ylim([0,100])
ax1.set_xlim([0,5])
ax1.set_xlabel('Mortality, $b$')
ax1.set_ylabel('Population density')
ax1.legend()
#ax1.axvspan(0, 0.4, color='k', alpha=0.2, lw=0)

ax2.plot(q,prevstore)
#ax2.set_title('Prevalence')
ax2.set_xlim([0,5])
ax2.set_ylim([0,1])
ax2.set_xlabel('Mortality, $b$')
ax2.set_ylabel('Prevalence')
#ax2.axvspan(0, 0.4, color='k', alpha=0.2, lw=0)

#plt.savefig('compf1_pops.png')



