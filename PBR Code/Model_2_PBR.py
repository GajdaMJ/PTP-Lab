import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sps
import scipy.sparse.linalg as spss

x_end = 1 # meter
tend = 4
d = 1e-3
kr = 1.2 #per second

# mew = volumetric flow / cross sectional area
diam = 0.04 # meter
flow = 6.3e-4 #m^3/s
mew = flow / (np.pi*(diam**2)/4)


def implicit(T, nx, nt, x_end, t_end, fv1,fv2, C, V=137):
    '''Models the behavior of the reaction: Water + Acetic Anhydride -> 2 * Acetic acid in an adiabatic CSTR reactor. \n
    Required Arguments: \n
    T = inlet temperature for the reactor given in unit celsius \n
    nx = number of steps \n
    nt = number of time steps \n
    x_end = the length of the column \n
    t_end = final time \n
    C = values [Water, AAH, Acitic Acid, Temperature] \n
    V = volume of the PBR reactor
    '''
    #Chemical constants
    mm_water = 18.01528 # (g/mol)
    rho_water = 0.999842 # (g/ml)
    cw_pure = rho_water/mm_water # (mol/ml)
    mm_AAH = 102.089 # (g/mol)
    rho_AAH = 1.082 # (g/ml)
    caah_pure = rho_AAH/mm_AAH # (mol/ml)

   # Convert flow rates (ml/min to ml/s)
    fv_w_dm3_s = fv1 / 60  # Water flow rate in ml/s
    fv_a_dm3_s = fv2  / 60  # Anhydride flow rate in ml/s

    flow_array = [fv_w_dm3_s, fv_a_dm3_s]
    
    #Thermodynamic constants
    params = {
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T+273.15,
        "flow": flow_array,
        "V": v_pfr,  # Volume in ml
        "k0": 6e5,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 45622.34,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho": 1,            # Density (g/ml)
        "cp": 4.186             # Heat capacity (J/g/K)
    }

    dx = x_end/nx

    v_pfr= V/x #devide the pfr volume by the number of reactors

    dt = tend/nt
    x = np.linspace(0,x_end,nx +1) # x is the number of steps (reactors)

    #initializing array for results
    c = np.zeros((nx+1,4))
    c[:,0] = C[0] #water fills the whole reactore
    c[0,1] = C[1]
    c[0,2] = C[2] 
    c[:,3] = C[3] #temp is the same across the reactor
    
    #t = np.linspace(0,tend,nt+1)
    t = 0

    A = sps.diags([ - mew*dt/(dx), -1 - (-kr)*dt + mew*dt/dx],[0,1], shape = (nx+1,nx+1))
    A = sps.csr_matrix(A) #efficient to solve using ludecomp
    A[0,0] = 1
    A[0,1] = 0
    A[nx,nx] = 1
    A[nx,nx-1] = -1 #-1 for zero gradient # 0 for...
    for i in range(1, nx):
        if i*dx <0.1 or i*dx>0.9:
            A[i,i] = A[i,i] - kr*dt

    for n in range (nt):
        t += dt
        c_old = np.copy(c)#do not store the solution for every time step
    
        b = c_old
        b[nx]=0
        c = sps.linalg.spsolve(A,b)

        iplot +=1
        if (iplot % 10 == 0):
            plt.title(f"Time = {t:5.3f} s")
            line[0].set_ydata(c)
            figure.canvas.draw()
            figure.canvas.flush_events()
            iplot = 0

    iplot = 0
    plt.ion()
    figure, ax = plt.subplots(figsize = (10,5))
    line = []
    line += ax.plot(x,c)
    plt.title(f"Time = {t:5.3f} s")
    plt.xlabel("Axial position")
    plt.ylabel("concentration")
    plt.xlim(0,x_end)
    # plt.ylim(0,max(c_l,c_r))
    plt.grid()
    plt.pause(40)
    plt.show()

implicit(10, 1000)


