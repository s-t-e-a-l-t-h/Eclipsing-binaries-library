#!/usr/bin/python
import objects.Star as Star
import objects.Binary as Binary
import objects.Orbit as Orbit
import objects.Observer as Observer
import objects.Function as Fn
import numpy as np
import globe.variables as gv
import MySQLdb
import ast

mysql_conn = MySQLdb.connect(host=gv.HOST,  # your host, usually localhost
                             user=gv.USER,  # your username
                             passwd=gv.PWD,  # your password
                             db="elisa")  # name of the data base
c = mysql_conn.cursor()

mass = [0.5, 2]
temp = [4000.0, 7000.0]
potential = [5.0, 10.0]
atm = ["black-body", "castelli-kurucz-04"]

temp_list = np.arange(temp[0], temp[1] + 100.0, 100.0)
mass_list = np.arange(mass[0], mass[1] + 0.25, 0.25)
init_problem = []
i = 0
f = open("testing_range.log", "w")

orbit = Orbit.Orbit(orbital_period=5.0, eccentricity=0.0, inclination=np.pi / 2.0,
                    argument_of_periastron=np.radians(90), verbose=False)

for tp in temp_list:
    for ts in temp_list:
        for mp in mass_list:
            for ms in mass_list:
                f.write("id: " + str(i) + ", temp P: " + str(
                    tp) + ", temp S: " + str(ts) + ", mass P: " + str(mp) + ", mass S :" + str(ms) + "\n")

                meta = [ 'mp', 'ms', 'tp', 'ts', 'pp', 'ps', 'zp', 'zs', 'gdp', 'gds', 'albp', 'albs', 'period', 'incl', 'bandpass']
                dt = [mp, ms, tp, ts, 10.0, 10.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 5.0, 90.0, "Generic/Bessell.U"]
                d = [str(it1) + "=\"" + str(it2) + "\"" for it1, it2 in list(zip(meta, dt))]
                q = "select id from tb_1 where " + " and ".join(d)

                c.execute(q)
                d = c.fetchall()
                if Fn.empty(d):

                    primary = Star.Star(mass=mp, synchronicity_parameter=1.0, potential=10.0, effective_temperature=tp,
                                        gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=False,
                                        atmoshpere=atm[1],
                                        phi_steps=10, theta_steps=20)

                    secondary = Star.Star(mass=ms, synchronicity_parameter=1.0, potential=10.0, effective_temperature=ts,
                                          gravity_darkening=1.0, metallicity=0.0, albedo=0.0, verbose=False,
                                          atmoshpere=atm[1],
                                          phi_steps=10, theta_steps=20)

                    binary_system = Binary.Binary(primary=primary, secondary=secondary, system="eb", orbit=orbit,
                                                  verbose=False)

                    if not binary_system.init:
                        f.write("binary init problem\n")
                        f.write("---------------------------------------------------------------------------\n")
                        continue
                    else:
                        f.write("binary morph: " + str(binary_system.binary_morph) + "\n")
                        if binary_system.binary_morph != "detached":
                            f.write("continue\n")
                            f.write("---------------------------------------------------------------------------\n")
                            continue
                        #     init_problem.append([mp, ms, tp, ts])
                    #
                    #     print(binary_system.get_exception() + i)
                    #

                    observer = Observer.Observer(passband="Generic/Bessell.U", limb_darkening_model="linear",
                                                 observe=binary_system,
                                                 limb_darkening_interp_method="nearest", verbose=False)

                    lc = observer.compute_lightcurve(
                        lightcurve_params={"from_photometric_phase": -0.2, "to_photometric_phase": 1.2, "n": 20,
                                           "adaptive_orbit_part": 0.05,
                                           "adaptive_multiplicator": 5.0},
                        starmodel_params={"critical_angle": np.pi / 4.0, "homo": True},
                        postprocess_params={"gaussian_smooth": False, "fix_minima": True}
                        # limb_darkening_params={"temperature_primary": 5000.0, "temperature_secondary": 5000.0,
                        #                       "gravity_primary": 10**4.0, "gravity_secondary": 10**4.0,
                        #                       "metallicity_primary": 0.0, "metallicity_secondary": 0.0}
                    )
                    if not Fn.empty(primary.get_exception()) or not Fn.empty(secondary.get_exception()) or not Fn.empty(
                            binary_system.get_exception()) or not Fn.empty(observer.get_exception()):
                        f.write("primary exception: " + str(primary.get_exception()) + "\n")
                        f.write("secondary exception: " + str(secondary.get_exception()) + "\n")
                        f.write("binary exception: " + str(binary_system.get_exception()) + "\n")
                        f.write("observer exception: " + str(observer.get_exception())  + "\n")

                    norm = lc[1]
                    lc = lc[0]

                    if len(lc) != 201:
                        f.write("lightcurve len problem: " + str(len(lc)) + "\n")
                    else:
                        lc_flux =  str(list(zip(*lc))[1] / norm)
                        query = "insert into tb_1 (mp, ms, tp, ts, pp, ps, zp, zs, gdp, gds, albp, albs, period, incl, bandpass, lc) " \
                                "values(\"" + str(
                            mp) + "\", \"" + str(ms) + "\", \"" + str(tp) + "\", \"" + str(
                            ts) + "\", \"10.0\", \"10.0\", \"0.0\", \"0.0\", \"1.0\", \"1.0\", \"0.0\", \"0.0\", \"5.0\", \"90.0\", \"Generic/Bessell.U\", \"" + lc_flux + "\")"

                        # print(list(zip(*lc))[0])
                        # import sys
                        # sys.exit()
                        # c.execute(query)
                        mysql_conn.commit()


                    f.write("---------------------------------------------------------------------------\n")

                    if i % 100 == 0: print(i)
                    i += 1

                    f.close()
                    f = open("testing_range.log", "a")

f.close()
mysql_conn.close()
