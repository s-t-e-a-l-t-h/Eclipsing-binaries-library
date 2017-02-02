#!/usr/bin/python
import time

elapseStart = time.time()
print('\033[92m' + 'Info: ' + '\033[0m' + 'Loading libraries...')
import objects.Star as Star
import objects.Binary as Binary
import globe.variables as gv
import objects.Plot as Plt
import objects.Function as Fn
import objects.Geometry as Geo
import objects.Orbit as Orb
import objects.Observer as Obs
import objects.Iostream as Io
import sys
import numpy as np

print('\033[92m' + 'Info: ' + '\033[0m' + 'Loading libraries done in ' + str(
    round(time.time() - elapseStart, 5)) + " sec.")
elapseStartClear = time.time()

# self testovacie funkcie
# tieto funkcie nie su 100-percentne, kvoli priamym osetrovaniam niektorych veci v kode
# t_Self.init_binary_mass_test( rand_values = 100 , potential = 20 )
# t_Self.lagrangian_points_test(rand_values = 100)

# vytvorenie primarnej zlozky
primary = Star.Star(mass=4., synchronicity_parameter=1.0, potential=2.7, effective_temperature=12000.,
                    gravity_darkening=1., metallicity=0., albedo=0.0, verbose=True, atmoshpere="black-body")
# vytvorenie sekundarnej zlozky
secondary = Star.Star(mass=2., synchronicity_parameter=1.0, potential=2.7, effective_temperature=12000.,
                      gravity_darkening=1., metallicity=0., albedo=0.0, verbose=True, atmoshpere="black-body")
# vytvorenie binarneho systemu zo zloziek

if not primary.init or not secondary.init:
    sys.exit("Star init fail.")

orbit = Orb.Orbit(orbital_period=0.8, eccentricity=0.0, inclination=np.pi / 2.0,
                  argument_of_periastron=np.radians(90.0), verbose=True)

bs = Binary.Binary(primary=primary, secondary=secondary, system="eb", orbit=orbit, verbose=True)

observer = Obs.Observer(passband="Generic/Bessell.U", limb_darkening_model="sqrt", verbose=True, observe=bs)

model = bs.get_3d_model_optimized(t_object="both", actual_distance=1.0, critical_angle=np.pi / 4.0, phi_steps=10,
                                  theta_steps=20, zero_point=False, homo=True)


Plt.plot_3d(normals=None, vertices=[model['system']], faces=None, face_color="w", normals_view=False, points_view=True,
                faces_view=False, point_color="r", normal_color="w", point_size=3., verbose=True, face_alpha=1.,
                azim=30, elev=30)

normal_vectors = \
    Geo.normal_estimation(binary_object=bs, actual_distance=1.0, vertices=np.array(model["system"]),
                          t_object="primary",
                          mode="in_point", verbose=True)

triangulation = Geo.cgal_triangulation(normals=normal_vectors, points=model['system'], verbose=True,
                                       min_triangle_angle=np.radians([30.0])[0], max_triangle_size=10,
                                       surface_aproximation_error=0.3, to_average_spacing=2)

Plt.plot_3d(faces=[triangulation[0]], normals_view=False, points_view=False, faces_view=True,
            point_color="r", point_size=5.0, face_color="c", edge_color="k", verbose=True)




elapseEnd = time.time()
print(Fn.color_string("info", "Info: ") + 'Computing elapsed time: ' + str(
    round(elapseEnd - elapseStartClear, 5)) + ' sec')

print(Fn.color_string("info", "Info: ") + Fn.color_string("warning", 'Elapsed time: ' + str(
    round(elapseEnd - elapseStart, 5)) + ' sec'))