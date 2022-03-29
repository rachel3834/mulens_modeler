# Event simulator based on example by Etienne Bachelet:
# https://github.com/ebachelet/pyLIMA/blob/master/examples/pyLIMA_example_3.py

import numpy as np
import matplotlib.pyplot as plt
import os, sys

from pyLIMA import microlsimulator
from pyLIMA import microltoolbox
from pyLIMA import microlmodels


simEvent = microlsimulator.simulate_a_microlensing_event(name ='Simulated event',
                                                                ra=270.0, dec=-30.0)

# Parameters:
obsStart = 2460766.500000
obsEnd = 2461356.500000
t0 = 2460953.0
u0 = -0.38
tE = 35.867
mu_rel = 7.0
rho = 0.01
log_s = -0.2
log_q = np.log10(0.02/0.62)
alpha = 6.84*np.pi/180.0 # deg -> radians
pi_en = 0.0
pi_ee = 0.0
source_imag = microltoolbox.magnitude_to_flux(22.0)
g_blend = 0.0

# Name = survey, eventObject, location = 'Earth', start_obs =2457465.500000, end_obs = 2457665.500000,
# sampling(hours) = 4, location='Earth', uniform_sampling=False, filter = 'I', altitude = 1000 m,
# longitude = -109.285399,
# latitude = -27.130814, bad_weather_percentage = 10%, moon_windows_avoidance (degree)=30,
# minimum_alt=30)
tel = microlsimulator.simulate_a_telescope('dataset',simEvent, obsStart, obsEnd, 4, 'Earth','I',
                                                  uniform_sampling=True, altitude=1000, longitude = -109.285399, latitude = -27.130,
                                                  bad_weather_percentage=0.0 / 100, moon_windows_avoidance=0,
                                                  minimum_alt=30)
simEvent.telescopes.append(tel)

eventModel = microlsimulator.simulate_a_microlensing_model(simEvent, model='FSBL', parallax=['None', 0.0],
                                             xallarap=['None', 0.0],
                                             orbital_motion=['None', 0.0], source_spots='None')

modelParameters = microlsimulator.simulate_microlensing_model_parameters(eventModel)
modelParameters = [t0, u0, tE, rho, log_s, log_q, alpha]

fluxParameters = microlsimulator.simulate_fluxes_parameters(simEvent.telescopes)
fluxParameters = [source_imag, g_blend]
modelParameters += fluxParameters

# Transform into pyLIMA standards
pyLIMA_parameters = eventModel.compute_pyLIMA_parameters(modelParameters)

# Simulate the lightcurve of the defined event :
microlsimulator.simulate_lightcurve_flux(eventModel, pyLIMA_parameters,  red_noise_apply='No')

#### Plot it!
show_2d_plot = False
if show_2d_plot:
    for telescope in simEvent.telescopes:
        plt.errorbar(telescope.lightcurve_magnitude[:, 0]-2450000, telescope.lightcurve_magnitude[:, 1],
                     yerr=telescope.lightcurve_magnitude[:, 2], fmt='.',label=telescope.name)


    for key in eventModel.model_dictionnary.keys():
        idx = eventModel.model_dictionnary[key]

        print (key, ' = ', modelParameters[idx] )

    plt.gca().invert_yaxis()
    plt.legend(numpoints=1)
    plt.grid(True)
    plt.show()

time = simEvent.telescopes[0].lightcurve_flux[:, 0]
flux_model = eventModel.compute_the_microlensing_model(simEvent.telescopes[0], pyLIMA_parameters)[0]
magnitude = microltoolbox.flux_to_magnitude(flux_model)

fig = plt.figure()
ax = plt.axes(projection='3d')

icolumns = 5
idx = np.where((time >= (t0-2.0*tE)) & (time <= (t0+2.0*tE)))[0]
x = time[idx]
y = np.linspace(0.0,float(icolumns),icolumns)
X,Y = np.meshgrid(x, y)
Z = np.zeros((icolumns,len(x)))
for i in range(0,icolumns,1):
    Z[i,:] = magnitude[idx]
#print(XX.shape, YY.shape, Z.shape)
#print(XX)

#def f(x, y):
#    return np.sin(np.sqrt(x ** 2 + y ** 2))

#x = np.linspace(-10, 20, 30)
#y = np.linspace(-6, 6, 30)

#X, Y = np.meshgrid(x, y)
#Z = f(X, Y)
colormap = 'copper'
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap=colormap, edgecolor='none')
#ax.contour3D(X, Y, Z, 50, cmap=colormap)
#plt.plot(time, magnitude, 'r-')
#plt.gca().invert_yaxis()
ax.axes.set_xlim3d(left=(t0-2.0*tE), right=(t0+2.0*tE))
ax.axes.set_ylim3d(bottom=0, top=10)
ax.axes.set_zlim3d(bottom=Z.max(), top=Z.min())

ax.view_init(8,-79)
ax.set_axis_off()

#plt.show()
plt.savefig('3D_event_lc.png')

#model_1 = microlmodels.create_model('FSBL', simEvent)
#simEvent.fit(model_1,'DE')

#simEvent.fits[0].produce_outputs()

#print(my_own_creation.fits[0].model.model_type,'Chi2_LM :',my_own_creation.fits[0].outputs.fit_parameters.chichi)

#plt.show()
