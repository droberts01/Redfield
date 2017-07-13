import ODE_integrator
import numpy as np
import matplotlib.pyplot as plt
import pylab

dt = 0.05
gamma = 2 - np.sqrt(2)
list_of_t = np.arange(0,1, dt)
for j in range(len(list_of_t)):
	if j % 2 == 1:
		list_of_t[j] = list_of_t[j - 1] + 2 * gamma * dt

# print (len(list_of_t))
# print (len(list_of_t)/25)
# print(list_of_t)

rho_0 = [1,1]
k = -2

list_of_linblads = [k*np.identity(2, dtype = complex)]*len(list_of_t)

print("hello.")
print(list_of_linblads[0])

list_of_rho_TRBDF2 = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "TR-BDF2")
list_of_rho_Explicit_Euler = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "Explicit Euler")
list_of_rho_TR = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "TR")
list_of_rho_BDF2 = ODE_integrator.run_time_evolution(rho_0, list_of_linblads, list_of_t, "BDF2")


pylab.plot(list_of_t, list(map(np.exp, k*list_of_t)), '-k',label = 'Exact Solution')
pylab.plot(list_of_t[:-2], list_of_rho_TRBDF2[:-2,0], '-m',label = 'TR-BDF2')
pylab.plot(list_of_t[:-2], list_of_rho_TR[:-2,0], '-b',label = 'TR')
pylab.plot(list_of_t[:-2], list_of_rho_BDF2[:-2,0], '-r',label = 'BDF2')


pylab.legend(loc='upper left')


pylab.show()