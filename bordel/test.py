#!/usr/bin/python
import objects.Star as Star
import objects.Binary as Binary
import objects.Orbit as Orbit
import objects.Observer as Observer
import objects.Function as Fn
import objects.Lightcurve as Lc

import numpy as np

atm = ["black-body", "castelli-kurucz-04"]

primary = Star.Star(mass=1.0, synchronicity_parameter=1.0, potential=5.0, effective_temperature=10500.,
                    gravity_darkening=1.0, metallicity=-0.5, albedo=0.0, verbose=True, atmoshpere=atm[1],
                    phi_steps=8, theta_steps=16)


secondary = Star.Star(mass=0.55, synchronicity_parameter=1.0, potential=7.0, effective_temperature=10500.,
                      gravity_darkening=1.0, metallicity=-0.5, albedo=0.0, verbose=True, atmoshpere=atm[1],
                      phi_steps=8, theta_steps=16)

orbit = Orbit.Orbit(orbital_period=2.8, eccentricity=0.0, inclination=np.pi / 2.0,
                    argument_of_periastron=np.radians(90), verbose=True)

motion = orbit.orbital_motion_beta(from_photometric_phase=0.0, to_photometric_phase=1.0, n=30,
                                   adaptive_orbit_part=0.05, adaptive_multiplicator=5.0,
                                   argument_of_periastron=orbit.get_argument_of_periastron(),
                                   eccentricity=orbit.get_eccentricity())

import objects.Plot as Plt
motion_to_plot = [[item[0] * np.cos(item[1]), item[0] * np.sin(item[1])] for item in motion]
# Plt.plot_2d(points=motion_to_plot, grid=True)

import numpy as np
import pandas
import objects.Plot as Plt
lightcurve = pandas.read_csv("curve2.dat", sep = "\t", header = None)
lc = [[lightcurve.iloc[i][0], lightcurve.iloc[i][1]] for i in range(0, len(lightcurve))]

lc_eclipse = Lc.photometric_phase_of_eclipses(pp_primary=orbit.conjuction[0]["true_phase"],
                                              pp_secondary=orbit.conjuction[1]["true_phase"],
                                              lightcurve=lc)




lc_clear = Lc.fix_minima(eclipses=lc_eclipse, lightcurve=lc)
import matplotlib.pyplot as plt
fig = plt.figure()
# ______________________________________________________________________
# ax = fig.add_subplot(111, aspect="auto")
# xs = zip(*lc_clear)[0]
# ys = zip(*lc_clear)[1] / max(zip(*lc_clear)[1])
# ax.scatter(xs, ys, color="r", marker="o", s=10.0)
# ax.plot(xs, ys, color="r")


# lc_interpoalted = Lc.akima_interpolation(from_photometric_phase=min(xs), to_photometric_phase=max(xs),
#                                          lightcurve=lc_clear, mirror=False)
# ______________________________________________________________________
# ax = fig.add_subplot(111, aspect="auto")
# xs = zip(*lc_interpoalted)[0]
# ys = zip(*lc_interpoalted)[1] / max(zip(*lc_interpoalted)[1])
# # ax.scatter(xs, ys, color="b", marker="o", s=10.0)
# ax.plot(xs, ys, color="y")





# ______________________________________________________________________
# lc_gaussina = Lc.gaussian_smooth(lightcurve=lc_interpoalted)
# ax = fig.add_subplot(111, aspect="auto")
# xs = zip(*lc_gaussina)[0]
# ys = zip(*lc_gaussina)[1] / max(zip(*lc_gaussina)[1])
# ax.scatter(xs, ys, color="r", marker="o", s=2.0)
# ax.plot(xs, ys, color="g")











# ax = fig.add_subplot(111, aspect="auto")
# xs = zip(*lc)[0]
# # ys = zip(*lc)[1]
# ys = zip(*lc)[1] / max(zip(*lc)[1])
#
# ax.scatter(xs, ys, color="r", marker="o", s=1.0)
# ax.plot(xs, ys, color="r")



# ______________________________________________________________________
lightcurve = pandas.read_csv("curve.dat", sep = "\t", header = None)
lc = [[lightcurve.iloc[i][0], lightcurve.iloc[i][1]] for i in range(0, len(lightcurve))]
ax = fig.add_subplot(111, aspect="auto")
xs = zip(*lc)[0]
# ys = zip(*lc)[1]
ys = zip(*lc)[1] / max(zip(*lc)[1])
ax.scatter(xs, ys, color="b", marker="o", s=5.0)
ax.plot(xs, ys, color="b")

# plt.savefig("figure_3" + ".png")





ax.grid(True)
plt.show()




