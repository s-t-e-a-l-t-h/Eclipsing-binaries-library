#!/usr/bin/python
import objects.Star as Star
import objects.Binary as Binary
import objects.Orbit as Orbit
import objects.Observer as Observer
import objects.Function as Fn
import time

import numpy as np

atm = ["black-body", "castelli-kurucz-04"]

# mass = 2.0, synchronicity_parameter = 1.0, potential = 5.0, effective_temperature = 10500.,
# gravity_darkening = 1.0, metallicity = -0.5, albedo = 0.0, verbose = True, atmoshpere = atm[0],
# phi_steps = 8, theta_steps = 16

c = 10
mass_p = Fn.rand(0.3, 50.0, c, time.time())
mass_s = Fn.rand(0.3, 50.0, c, time.time())

temp_p = Fn.rand(3500.0, 50000.0, c, time.time())
temp_s = Fn.rand(3500.0, 50000.0, c, time.time())

omega_p = Fn.rand(2.5, 30.0, c, time.time())
omega_s = Fn.rand(2.5, 30.0, c, time.time())

period = Fn.rand(0.1, 10.0, c, time.time())

for mp, ms, tp, ts, p, op, os in zip(mass_p, mass_s, temp_p, temp_s, period, omega_p, omega_s):

    primary = Star.Star(mass=mp, synchronicity_parameter=1.0, potential=op, effective_temperature=tp,
                        gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=True, atmoshpere=atm[1],
                        phi_steps=5, theta_steps=10)

    secondary = Star.Star(mass=ms, synchronicity_parameter=1.0, potential=os, effective_temperature=ts,
                          gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=True, atmoshpere=atm[1],
                          phi_steps=5, theta_steps=10)

    orbit = Orbit.Orbit(orbital_period=p, eccentricity=0.0, inclination=np.pi / 2.0,
                        argument_of_periastron=np.radians(90), verbose=True)

    binary_system = Binary.Binary(primary=primary, secondary=secondary, system="eb", orbit=orbit, verbose=True)

    if not binary_system.init or binary_system.binary_morph == "over-contact":
        continue

    observer = Observer.Observer(passband="Generic/Bessell.U", limb_darkening_model="linear", observe=binary_system,
                                 limb_darkening_interp_method="nearest", verbose=True)

    observer.compute_lightcurve(
        lightcurve_params={"from_photometric_phase": -0.2, "to_photometric_phase": 1.2, "n": 20,
                           "adaptive_orbit_part": 0.05,
                           "adaptive_multiplicator": 5.0},
        starmodel_params={"critical_angle": np.pi / 4.0, "homo": True},
        postprocess_params={"gaussian_smooth": False, "fix_minima": True}
        # limb_darkening_params={"temperature_primary": 5000.0, "temperature_secondary": 5000.0,
        #                       "gravity_primary": 10**4.0, "gravity_secondary": 10**4.0,
        #                       "metallicity_primary": 0.0, "metallicity_secondary": 0.0}
    )








