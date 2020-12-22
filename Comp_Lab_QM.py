import numpy as np
import matplotlib.pyplot as plt
ximax=8
ximin=-ximax
Nsteps=100*2*ximax
nu0=6.0
h = (ximax -ximin)/(Nsteps)




Vo = -12

Xi = np.linspace(ximin,ximax,Nsteps)

''' psi(E) takes in a specific energy value and gives the array of values for that function
     as well as the number of nodes it posses.'''
def psi(E):
         #initial starting values
    node_count=0
    f = Vo*(np.cosh(Xi)**(-2))  -E
    psi_vec=np.zeros(Nsteps)
    psi_vec[0]=1
    psi_vec[1]=np.exp(np.sqrt(-E)*h)
    for i in range(1,len(Xi)-1):    
        psi_vec[i+1]=((2 +(h**2)*f[i])*psi_vec[i] - psi_vec[i-1])
    
        if psi_vec[i-1]*psi_vec[i]<0:
            node_count+=1
    return [psi_vec,node_count]




'''Takes in number of nodes we want energy to have and gives us corresponding energy'''
def bisec(n):
    e_min=Vo
    e_max=-0.1
    if psi(e_min)[1]>n or psi(e_max)[1]<n:
        return print('Cannot find energy with required nodes')
    while abs((e_max-e_min)/(e_min+e_max))> 1e-14:
        e=(e_max+e_min)/2
        if psi(e)[1]<n:
            e_min=e
        elif psi(e)[1]>n:
            e_max=e
        elif psi(e)[0][-2]*psi(e)[0][-1] - (np.e**(np.sqrt(-e)*h))*((psi(e)[0][-1])**2) < 0:
            e_min=e
        else:
            e_max=e
    return (e_min+e_max)/2

ö = 0

m = psi(bisec(ö))[0]   
        
print(bisec(ö))

plt.figure(figsize=(5,5))
plt.title("Plot of the wavefunction at its 3rd excited state")
plt.xlabel("Xi")
plt.ylabel("Phi")
plt.plot(Xi, m)
        
        
        
        
        
        
        
        
        
        
        
        