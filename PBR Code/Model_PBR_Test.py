import numpy as np
import matplotlib.pyplot as plt

# Getting all the variables
Nx = 100 ## Grid points
Nt = 10000 # Time steps
diam = 2.5 # Diameter (cm)
fv_w = 23.7/60  # Water Flow Rate (ml/min)
fv_a = 5.2/60 # Anhydride Flow Rate (ml/min)

inital_conditions = [1,0.5,0,298]



params = {
    "k0": np.exp(16.25)*1000,          # Reaction rate constant (ml/mol/s)

    # Thermodynamic constants (taken from Asprey et al., 1996)
    "Ea": 45622.34,             # Activation energy (J/mol)
    "R": 8.314,              # Gas constant (J/mol/K)
    "H": -56.6e3,              # Enthalpy change (J/mol)
    "rho": 1,            # Density (g/ml)
    "cp": 4.186             # Heat capacity (J/g/K)
}

t_end = 6 # 4 sec
x_end = 1 # Length




u = (fv_w+fv_a) / (np.pi*(diam/2)**2) 

dt = t_end / Nt
dx  = x_end / Nx


cR=1
# Courant number
Co = u* dt /dx
print(Co)
#Uncomment to see that Co is less than 1

x = np.linspace(0,x_end, Nx+1)
c = np.zeros((Nx+1,4))

c[0,:] = inital_conditions
c[:,3] = inital_conditions[3]
t = 0


def reaction(c):
    return params["k0"]*np.exp(-params["Ea"]/(params["R"]*c[:,3]))*c[:,0]*c[:,1]

plt.ion()
fig_conc, ax_conc = plt.subplots(figsize=(10, 5))
line_c1, = ax_conc.plot(x, c[:, 0], label="Reactant A")
line_c2, = ax_conc.plot(x, c[:, 1], label="Reactant B")
line_c3, = ax_conc.plot(x, c[:, 2], label="Product")
ax_conc.set_title(f"Time: {t:5.3f} s", fontsize=16)
ax_conc.set_xlabel('Axial position', fontsize=14)
ax_conc.set_ylabel('Concentration', fontsize=14)
ax_conc.set_xlim(0, x_end)
ax_conc.set_ylim(0, max(1, cR))
ax_conc.grid(True)
ax_conc.legend()

# Initialize the plot for temperature in a separate window
fig_temp, ax_temp = plt.subplots(figsize=(10, 5))
line_temp, = ax_temp.plot(x, c[:, 3], label="Temperature", color='r')
ax_temp.set_title(f"Time: {t:5.3f} s", fontsize=16)
ax_temp.set_xlabel('Axial position', fontsize=14)
ax_temp.set_ylabel('Temperature (K)', fontsize=14)
ax_temp.set_xlim(0, x_end)
ax_temp.set_ylim(min(c[:, 3]), max(c[:, 3]) + 100)
ax_temp.grid(True)

iplot = 0
for n in range(Nt):
    t += dt
    c_old = np.copy(c) ## doing c_old = c then modification of c_old will also modify c
    c[0,:] = inital_conditions

    rxn = reaction(c_old)

    for i in range(1, Nx):
        c[i,0] = c_old[i,0] - Co*(c_old[i,0]-c_old[i-1,0]) - rxn[i]*dt
        c[i,1] = c_old[i,1] - Co*(c_old[i,1]-c_old[i-1,1]) - rxn[i]*dt
        c[i,2] = c_old[i,2] - Co*(c_old[i,2]-c_old[i-1,2]) + 2* rxn[i]*dt
        c[i,3] = c_old[i,3] - Co*(c_old[i,3]-c_old[i-1,3]) - (params["H"]/(params["rho"]*params["cp"])) * rxn[i]*dt

    c[Nx] = c[Nx-1]

    iplot += 1
    if iplot % 1000 == 0:
        # Update concentration plot
        ax_conc.set_title(f"Time = {t:5.3f} s")
        line_c1.set_ydata(c[:, 0])
        line_c2.set_ydata(c[:, 1])
        line_c3.set_ydata(c[:, 2])
        fig_conc.canvas.draw()
        fig_conc.canvas.flush_events()

        # Update temperature plot
        ax_temp.set_title(f"Time = {t:5.3f} s")
        line_temp.set_ydata(c[:, 3])
        fig_temp.canvas.draw()
        fig_temp.canvas.flush_events()

        iplot = 0


plt.show()


