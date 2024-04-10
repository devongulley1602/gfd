"""
ps2_ATOC513.py

Author: Devon Gulley

Completion date: 2024-02-06

The objectives of this work are to illustrate the Rossby adjustment problem  as follows:
    - Compute the Fourier sums each of (u,v,phi) for geostrophic and Poincare plane waves
    - Graph this time evolution using a Hovmoller plot and the tendency to evolve towards geostropohy
    - Decomposes the total energy to determine the contributions from each mode
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


"""Constants and problem setup"""
n_max = 50 #the number of Fourier terms
epsilon = 0.025 #chosen ratio between deformation and characteristic radius
L = 1 #characteristic length scale
L_d = epsilon*L #deformation radius
k_d = np.ones(n_max)/L_d #constant wavevector associated with deformation radius
aspectRatio = 10**3 #L/H
H = L/aspectRatio #taking H to be much smaller than L
x = np.linspace(start = -0.5*L, stop = 0.5*L, num = 1000)
g = 10#arbitrary gravitational constant
n0 = 1.5/epsilon
f = np.sqrt(g*H)/L_d #shallow water f-plane approximation


"""Wave generation variables and functions"""
#For the nth sample in the Fourier series
n = np.ones(n_max)*range(n_max)

#A_n for the Gaussian pulse Fourier coefficients
A_0 = 2*epsilon*np.pi**0.5
A = A_0*np.exp(-1*(epsilon*np.pi*n)**2)

#Wave vector
k = n*2*np.pi/L

#(Positive) Poincare dispersion relationship
omega = np.sqrt(g*k +np.ones(n_max)*f**2)


"""
function wave(x_,t_,dispersion,phase,type)

    Returns the wave function with a given dispersion method at a specified time.

        dispersion:
            +1 or -1 for corresponding Poincare modes, 0 for geostrophy.

        phase:
            Shift in radians which can allows asymetrical wave representations about zero.

        type:
            'u' and 'v' for zonal and meridional velocities respectively
            'phi' for geopotential heights
"""
def wave(x_,t_=0,dispersion = 0,type = 'phi'):

    #Each kind of wave in (u,v,phi) has different properties set here
    phase = 0
    freqScale = np.ones(n_max)
    match type:
        case 'phi':#for this problem this is the most simplistic base case
            phase = 0
            freqScale = np.ones(n_max)

        case 'u':
            phase = 0
            freqScale = omega/(g*H)

        case 'v':
            phase = np.pi/2
            freqScale = np.ones(n_max)*(-f/(g*H))


    t_ = np.ones(len(x_))*t_
    phase = np.ones(len(x_))*phase

    #Geostrophic modes are associated with a K_d^2/(K_n^2 +K_d^2) and zero dispersion
    if(dispersion==0):#The geostrophic case, checking for zero cases is just for performance 
        if(phase[0] ==0):
            return sum((k_d[0]**2/(k[n]**2 + k_d[0]**2))*A[n]*np.cos(k[n]*x_) for n in range(1,n_max))
        return sum((k_d[0]**2/(k[n]**2 + k_d[0]**2))*A[n]*np.cos(k[n]*x_ -phase) for n in range(1,n_max))


    #Poinare Modes are associated with a K_n^2/(K_n^2+K_d^2) coeficients
    return 0.5*sum((k[n]**2/(k[n]**2 + k_d[0]**2))*A[n]*np.cos(k[n]*x_ - dispersion*omega[n]*t_-phase) for n in range(1,n_max))



#In order to see the time evolution of the waves
speed = 1/1000
displayHeight = 2.0
def animate(t):
    t = t*speed
    ax.cla()
    waveType = 'v'

    y_plus = wave(x,t,1,waveType)
    y_minus = wave(x,t,-1, waveType)
    y_g = wave(x,t,0, waveType)
    y_total = y_plus + y_minus + y_g

    plt.plot(x,y_plus)
    plt.plot(x,y_minus)
    plt.plot(x,y_g)
    plt.plot(x,y_total,c='red')

    view =displayHeight
    plt.legend(["+poincare","-poincare","geostrophy","overall"])
    ax.set_ylim(0-view,0+view)

#Displaying animated figures
fig = plt.figure()
ax = fig.add_subplot()
ani = animation.FuncAnimation(fig,animate,frames= 1200,interval = 1)
plt.show()



""" Creation of the Hovmoller data"""
i_max = 75
data = {'phi':[],'u':[],'v':[]}
plotNum = 0
for form in data.keys():
    superposition = np.zeros([i_max,len(x)])
    poincare = np.zeros([i_max,len(x)])
    for dispersion in [-1,0,1]:
        #Each row t of the hovmollerMatrix is the intensity of the wave as a function of x distance

        description = ''
        if dispersion == 1:
            description = "Poincare (+)"
        elif dispersion == -1:
            description = "Poincare (-)"
        else:
            description = "Geostrophic"

        #the choice of t=1.4*i/i_max is just to run the simulation long enough to see the geostrophic adjustment
        #but not long enough to see interference from the periodic boundaries
        hovmollerMatrix = [wave(x,1.4*i/i_max,dispersion,form) for i in range(i_max,0,-1)]
        #Plotting Hovmoller data
        fig, ax = plt.subplots()
        ax.set_xticks([])
        ax.set_yticks([])
        plt.xlabel(form+'(x)')
        plt.ylabel('time')
        ax.imshow(hovmollerMatrix,cmap='seismic',vmin = -0.5,vmax = 0.5)
        plt.title(form.upper() + ' ' + description)
        plt.show()
        if(dispersion != 0):
            poincare += hovmollerMatrix
        superposition += hovmollerMatrix
    data[form] = superposition
    #Plotting Hovmoller data
    description = 'Superposition'
    fig, ax = plt.subplots()
    ax.set_xticks([])
    ax.set_yticks([])
    plt.xlabel(form+'(x)')
    plt.ylabel('time')
    ax.imshow(superposition,cmap='seismic',vmin = -0.5,vmax = 0.5)
    plt.title(form.upper() + ' ' + description)
    plt.show()

i_max = 40
x = np.linspace(start = -0.5*L, stop = 0.5*L, num = 1000)

#This is highly inefficient but readable
#Calculates the energy contribution from the nth Fourier term at time t_i
def energy(n_i,t_i):
    c = np.sqrt(g*H)
    t_i = np.ones(len(x))*t_i
    kn = k[n_i]
    kd = k_d[0]
    An = A[n_i]
    wn = c*np.sqrt(kn**2+f**2)

    #for phi
    phi_minus =0.5 *((kn**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x - wn*t_i*1.4/i_max))
    phi_plus = 0.5 *((kn**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x + wn*t_i*1.4/i_max))
    phi_g = ((kd**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x))
    phi = phi_g + phi_minus + phi_plus

    #for u
    u_minus =(wn/(2*(c**2)))*((kn**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x - wn*t_i*1.4/i_max))
    u_plus = (wn/(2*(c**2)))*((kn**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x + wn*t_i*1.4/i_max))
    u_g = ((wn/(c**2))*(kd**2)/(kn**2+kd**2))*An*np.average(np.cos(kn*x))
    u = u_g + u_minus + u_plus

    #for v
    v_minus =(f/(2*(c**2)))*((kn**2)/(kn**2+kd**2))*An*np.average(np.sin(kn*x - wn*t_i*1.4/i_max))
    v_plus = (f/(2*(c**2)))*((kn**2)/(kn**2+kd**2))*An*np.average(np.sin(kn*x + wn*t_i*1.4/i_max))
    v_g = (f/(c**2))*((kd**2)/(kn**2+kd**2))*An*np.average(np.sin(kn*x))
    v = v_g + v_minus + v_plus

    #Energy calculations
    KE_g = (c**2)* ( u_g**2 + v_g**2 )
    KE_p = (c**2)* ( (u_minus + u_plus)**2 + (v_minus + v_plus)**2 )

    PE_g = phi_g**2
    PE_p = (phi_plus + phi_minus)**2

    PE = phi**2
    KE = (c**2)*(u**2 + v**2)

    E_g = KE_g + PE_g
    E_p = KE_p + PE_p
    E = KE + PE

    #Energy, geostophic kinetic, geostrophic potential, kinetic Poincare, potential Poincare,at time t
    t = t_i[0]*1.4/i_max
    return E,KE_g,PE_g,KE_p,PE_p,t
wavenumbers = n
plotMe = np.zeros([i_max,len(n)])
E    = np.zeros([i_max,len(n)])
KE_g = np.zeros([i_max,len(n)])
PE_g = np.zeros([i_max,len(n)])
KE_p = np.zeros([i_max,len(n)])
PE_p = np.zeros([i_max,len(n)])
t = np.zeros(i_max) 




for t_i in range(i_max-1,0,-1):
    for n_i in range(len( n)) :
        E[t_i][n_i],KE_g[t_i][n_i],PE_g[t_i][n_i],KE_p[t_i][n_i],PE_p[t_i][n_i],t[t_i] = energy(n_i,t_i)
   # print(KE_p[t_i])
   # plt.plot(k,KE_p[t_i])
   # plt.show()
kinds = zip([E,KE_g,PE_g,KE_p,PE_p],["Total","Kinestic Geostrophic","Potential Geostrophic","Kinetic Poincare","PotentialPoincare"])


#Plotting Hovmoller energy data
for (plotMe,name) in kinds:
    fig, ax = plt.subplots()
    plt.xlabel('Wavenumber')
    plt.ylabel('Time')
    im =plt.imshow(plotMe,cmap='hot')#,vmin = -1,vmax = 1)
    ax.set_xticks(n[0::7])
    ax.set_xticklabels(np.round(k[0::7],2))
    ax.set_yticks(range(i_max-1,0,-1)[0::10])
    ax.set_yticklabels(np.round(t,1)[0::10])
    plt.title("Domain Averaged " + name + " Energy")
    cbar = plt.colorbar(im)
    cbar.set_label("Energy")
    plt.show()


