#!/usr/bin/python
import time

elapseStart = time.time()
print('\033[92m' + 'Info: ' + '\033[0m' + 'Loading libraries...')
import objects.Star as Star
import objects.Binary as Binary
import objects.Orbit as Orbit
import objects.Observer as Observer
import objects.Function as Fn

import numpy as np

print('\033[92m' + 'Info: ' + '\033[0m' + 'Loading libraries done in ' + str(
    round(time.time() - elapseStart, 5)) + " sec.")
elapseStartClear = time.time()

atm = ["black-body", "castelli-kurucz-04"]

primary = Star.Star(mass=4.0, synchronicity_parameter=1.0, potential=2.65, effective_temperature=12000.00,
                    gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=True, atmoshpere=atm[1],
                    phi_steps=25, theta_steps=50,
                    cgal={"min_triangle_angle": np.radians([30.0])[0],
                          "max_triangle_size": 4,
                          "surface_aproximation_error": 0.3,
                          "to_average_spacing": 5}
                    )

secondary = Star.Star(mass=2.0, synchronicity_parameter=1.0, potential=2.65, effective_temperature=12000.00,
                      gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=True, atmoshpere=atm[1],
                      phi_steps=25, theta_steps=50)

orbit = Orbit.Orbit(orbital_period=0.8, eccentricity=0.0, inclination=np.pi / 2.0,
                    argument_of_periastron=np.radians(90), verbose=True)

binary_system = Binary.Binary(primary=primary, secondary=secondary, system="eb", orbit=orbit, verbose=True)

observer = Observer.Observer(passband="Generic/Bessell.U", limb_darkening_model="linear", observe=binary_system,
                             limb_darkening_interp_method="nearest", verbose=True)

# motion = orbit.orbital_motion_beta(from_photometric_phase=0.0, to_photometric_phase=0.5, n=30,
#                                    adaptive_orbit_part=0.05, adaptive_multiplicator=5.0,
#                                    argument_of_periastron=orbit.get_argument_of_periastron(),
#                                    eccentricity=orbit.get_eccentricity())

# motion = orbit.orbital_motion(photometric_phase_from=0.0, photometric_phase_to=0.5,
#                               shift=orbit.get_true_phase_of_periastron()[0],
#                               photometric_phase_step=0.05, eccentricity=orbit.get_eccentricity(),
#                               argument_of_periastron=orbit.get_argument_of_periastron())

# import sys
# import objects.Plot as Plt
#
# motion_to_plot = [[item[0] * np.cos(item[1]), item[0] * np.sin(item[1])] for item in motion]
# Plt.plot_2d(points=motion_to_plot, grid=True)
#
# sys.exit()

lc = observer.compute_lightcurve(
    lightcurve_params={"from_photometric_phase": -0.2, "to_photometric_phase": 1.2, "n": 20,
                       "adaptive_orbit_part": 0.05,
                       "adaptive_multiplicator": 5.0},
    starmodel_params={"critical_angle": np.pi / 4.0, "homo": True},
    postprocess_params={"gaussian_smooth": True, "fix_minima": True, "zbuffer": True, "mirroring": True}
    # limb_darkening_params={"temperature_primary": 5000.0, "temperature_secondary": 5000.0,
    #                       "gravity_primary": 10**4.0, "gravity_secondary": 10**4.0,
    #                       "metallicity_primary": 0.0, "metallicity_secondary": 0.0}
)

elapseEnd = time.time()
print(Fn.color_string("info", "Info: ") + 'Computing elapsed time: ' + str(
    round(elapseEnd - elapseStartClear, 5)) + ' sec')

print(Fn.color_string("info", "Info: ") + Fn.color_string("warning", 'Elapsed time: ' + str(
    round(elapseEnd - elapseStart, 5)) + ' sec'))
