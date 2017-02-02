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
# import objects.SelfTest as t_Self
import sys
import numpy as np

# import scipy

print('\033[92m' + 'Info: ' + '\033[0m' + 'Loading libraries done in ' + str(
    round(time.time() - elapseStart, 5)) + " sec.")
elapseStartClear = time.time()

# self testovacie funkcie
# tieto funkcie nie su 100-percentne, kvoli priamym osetrovaniam niektorych veci v kode
# t_Self.init_binary_mass_test( rand_values = 100 , potential = 20 )
# t_Self.lagrangian_points_test(rand_values = 100)

# vytvorenie primarnej zlozky
primary = Star.Star(mass=4., synchronicity_parameter=1.0, potential=4., effective_temperature=12000.,
                    gravity_darkening=1., metallicity=0., albedo=0.0, verbose=True, atmoshpere="black-body")
# vytvorenie sekundarnej zlozky
secondary = Star.Star(mass=2., synchronicity_parameter=1.0, potential=4., effective_temperature=12000.,
                      gravity_darkening=1., metallicity=0., albedo=0.0, verbose=True, atmoshpere="black-body")
# vytvorenie binarneho systemu zo zloziek

if not primary.init or not secondary.init: sys.exit()

orbit = Orb.Orbit(orbital_period=4.8, eccentricity=0.7, inclination=np.pi / 2.0,
                  argument_of_periastron=np.radians(270.0), verbose=True)

# motion = orbit.orbital_motion_beta(from_photometric_phase=0.0, to_photometric_phase=1.0,
#                                    eccentricity=orbit.get_eccentricity(),
#                                    argument_of_periastron=orbit.get_argument_of_periastron(), n=50,
#                                    adaptive_orbit_part=0.05, adaptive_multiplicator=3.0)
# motion_to_plot = [[item[0] * np.cos(item[1]), item[0] * np.sin(item[1])] for item in motion]
# Plt.plot_2d(points=motion_to_plot, grid=True)


# orbit.get_info()

# motion_to_plot = [[item[0] * np.cos(item[1]), item[0] * np.sin(item[1])] for item in motion]
# Plt.plot_2d(points=motion_to_plot, grid=True)

observer = Obs.Observer(passband="johnson_u", limb_darkening_model="sqrt", verbose=True)
bs = Binary.Binary(primary=primary, secondary=secondary, system="eb", orbit=orbit, verbose=True)

observer.compute_lightcurve()
sys.exit()
# bs.get_info()

# polarne gravitacne zrychlenie
# pocita sa aj v __init__ Binary class, pre d = 1, takze pre periastrum
# tu je to len na ukazku, pripadne nejake testovanie
# bs.primary.polar_gravity = bs.compute_polar_gravity(actual_distance=1.0, angular_velocity=orbit.mean_angular_velocity,
#                                             t_object="primary")
# bs.secondary.polar_gravity = bs.compute_polar_gravity(actual_distance=1.0, angular_velocity=orbit.mean_angular_velocity,
#                                               t_object="secondary")

# Hillova krivka
# xy_plane = bs.compute_equipotential_xy(actual_distance = orbit.periastron_distance, step = 0.01)
# Plt.plot_2d(points = xy_plane, point_marker = "o", savefig = False)

# xy_plane = bs.compute_equipotential_xy(actual_distance = 1.0, step = 0.01)
# Plt.plot_2d(points = xy_plane, point_marker = "o", savefig = False)

# Hillova rovina
# xy_plane = bs.compute_hill_plane_prototype(actual_distance=1.0, fsolve_steps=10, fsolve_radius=5, steps=200)
# Plt.plot_2d(points = xy_plane, point_marker = "o", savefig = False)


# Graf priebehu potencialu
# val = bs.get_potential_value(step = 0.001, interval=[-10, 10], actual_distance=1.0)
# Plt.plot_2d(points = val, point_marker = "o")

# vypocet povrchu
# ak sa jedna o dotykovy system, tak sa nepocita nulovy bod, inak ano
zp = False if bs.binary_morph == "over-contact" else True

model = bs.get_3d_model_optimized(t_object="both", actual_distance=1.0, critical_angle=np.pi / 4.0, phi_steps=10,
                                  theta_steps=20, zero_point=zp, homo=True)

# print modelu
# Plt.plot_3d(normals=None, vertices=[model['system']], faces=None, face_color="w", normals_view=False, points_view=True,
#                 faces_view=False, point_color="r", normal_color="w", point_size=3., verbose=True, face_alpha=1.,
#                 azim=30, elev=30)

# ak sa jedna o oddeleny system, tak sa pouzije konvexna
triangulation_primary, triangulation_secondary = None, None
if bs.binary_morph == "detached" or bs.binary_morph == "semi-contact":
    # convex hull triangulation
    triangulation_primary = Geo.convex_hull_triangulation(vertices=model['primary'], verbose=True)
    triangulation_secondary = Geo.convex_hull_triangulation(vertices=model['secondary'], verbose=True)

elif bs.binary_morph == "over-contact":
    # cgal triangulation

    # obdrzanie normalovych vektorov v bodoch povrchu potrebnych pre cgal triangulaciu
    normal_vectors = \
        Geo.normal_estimation(binary_object=bs, actual_distance=1.0, vertices=np.array(model["system"]),
                              t_object="primary",
                              mode="in_point", verbose=True)

    # navratova hodnota triangulacie: [python list [python list]]
    # struktura:
    # [
    # [[x00, y00, z00],
    #  [x01, y01, z01],
    #  [x02, y02, z02]],
    # [[x10, y10, z10],
    #  [x11, y11, z11],
    #  [x12, y12, z12]],
    # ...
    # [[xn0, yn0, zn0],
    #  [xn1, yn1, zn1],
    #  [xn2, yn2, zn2]]
    # ]

    triangulation = Geo.cgal_triangulation(normals=normal_vectors, points=model['system'], verbose=False,
                                           min_triangle_angle=0.349066, max_triangle_size=2,
                                           surface_aproximation_error=0.375, to_average_spacing=1)

    # plotovacia cast
    # bez bodov
    Plt.plot_3d_dev(normals=None, vertices=None, faces=[triangulation], face_color="w", normals_view=False,
                    points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
                    face_alpha=1., azim=30, elev=30)

    # s bodmi
    # Plt.plot_3d(normals=None, vertices=[model['system']], faces=[triangulation], face_color="w", normals_view=False,
    #                 points_view=True, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
    #                 face_alpha=1., azim=30, elev=30)

    sys.exit()

    # wuma_split, navratova hodnota [python list [python list]]
    # struktura:
    # [[primary faces], [secondary faces]]
    triangulation_splited = bs.wuma_split(faces=triangulation, verbose=True)
    triangulation_primary = triangulation_splited[0]
    triangulation_secondary = triangulation_splited[1]

    # plot primary
    # Plt.plot_3d(normals=None, points=False, faces=[triangulation_primary], face_color="w", normals_view=False,
    #                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
    #                 face_alpha=1., azim=30, elev=30)
    # plot secondary
    # Plt.plot_3d(normals=None, points=False, faces=[triangulation_secondary], face_color="w", normals_view=False,
    #                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
    #                 face_alpha=1., azim=30, elev=30)

primary.set_faces(faces=triangulation_primary)
secondary.set_faces(faces=triangulation_secondary)

# plotovacia cast
# to_plot = []
# to_plot.extend([triangulation_primary, triangulation_secondary])
# Plt.plot_3d(normals=None, vertices=None, faces=to_plot, face_color="w", normals_view=False, points_view=False,
#                 faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True, face_alpha=1.,
#                 azim=30, elev=30)

faces_orientation_primary = Geo.face_orientation_beta(faces=primary.get_faces(), binary_object=bs,
                                                           verbose=True, t_object="primary", actual_distance=1.0)
faces_orientation_secondary = Geo.face_orientation_beta(faces=secondary.get_faces(), binary_object=bs,
                                                             verbose=True, t_object="secondary",
                                                             actual_distance=1.0)

primary.set_faces_orientation(faces_orientation=faces_orientation_primary)
secondary.set_faces_orientation(faces_orientation=faces_orientation_secondary)

# zapinanie a vypinanie plotovania normal
if False:
    # TATO CAST KODU JE POTREBNA LEN PRE PLOTOVANIE NORMAL
    # trba ich zmenist a poposuvat na pozicie faziet
    # normalizacia na 1
    unit_vectors_primary = \
        Geo.vector_array_normalisation(vector_arr=faces_orientation_primary, multi=20.0, verbose=True)
    unit_vectors_secondary = \
        Geo.vector_array_normalisation(vector_arr=faces_orientation_secondary, multi=20.0, verbose=True)

    faces_com_primary = Geo.center_of_mass(faces=triangulation_primary, verbose=True)
    faces_com_secondary = Geo.center_of_mass(faces=triangulation_secondary, verbose=True)

    translation_primary = \
        Geo.vector_array_translation(vector_arr=unit_vectors_primary, translation_arr=faces_com_primary,
                                     verbose=True)
    translation_secondary = \
        Geo.vector_array_translation(vector_arr=unit_vectors_secondary, translation_arr=faces_com_secondary,
                                     verbose=True)

    # farby pre normaly
    c_primary, c_secondary = ["#0000ff"] * len(translation_primary), ["#000055"] * len(translation_secondary)

    Plt.plot_3d_dev(normals=[translation_primary, translation_secondary],
                    vertices=[faces_com_primary, faces_com_secondary],
                    faces=[triangulation_primary, triangulation_secondary], face_color="w", normals_view=True,
                    points_view=False, faces_view=True, point_color=[c_primary, c_secondary],
                    normal_color=[c_primary, c_secondary],
                    point_size=3.0, verbose=True, face_alpha=1., azim=30, elev=30)

    # KONIEC CASTI PRE ZOBRAZOVANIE

gradnorm_primary = Geo.gradient_norm(faces=triangulation_primary, verbose=True, binary_object=bs,
                                     actual_distance=1.0, t_object="primary")

gradnorm_secondary = Geo.gradient_norm(faces=triangulation_secondary, verbose=True, binary_object=bs,
                                       actual_distance=1.0, t_object="secondary")

# gradnorm_primary = Fn.gradient_norm(gradient=faces_orientation_primary[1])
# gradnorm_secondary = Fn.gradient_norm(gradient=faces_orientation_secondary[1])

primary.compute_gravity_distribution(gradnorm=gradnorm_primary)
secondary.compute_gravity_distribution(gradnorm=gradnorm_secondary)

primary.gravity_darkening_factor_distribution()
secondary.gravity_darkening_factor_distribution()

primary.compute_polar_temperature()
secondary.compute_polar_temperature()

primary.compute_temperature_distribution()
secondary.compute_temperature_distribution()

observer.set_limb_darkening_model("sqrt")
ldc = observer.limb_darkening_coefficients(limb_darkeing_model=observer.get_limb_darkening_model(),
                                           passband=observer.get_passband(),
                                           metallicity=primary.get_metallicity(),
                                           temperature=primary.get_temperature_distribution(),
                                           gravity=primary.get_gravity_distribution(),
                                           interpolation_method="nearest",
                                           verbose=True)

# plot gravity distribution primary
# rgb_gp = Fn.arr_to_rainbow(arr=primary.get_gravity_distribution(), minimum=min(primary.get_gravity_distribution()),
#                            maximum=max(primary.get_gravity_distribution()))  # primary.get_polar_gravity())
# #
# hex_gp = Fn.rgb_to_hex(color=rgb_gp, sharp=True)
# Plt.plot_3d(normals=None, vertices=None, faces=[triangulation_primary], face_color=[hex_gp], normals_view=False,
#                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
#                 face_alpha=1., azim=30, elev=30)

# print primary.get_gravity_scalling_factor()
# print secondary.get_gravity_scalling_factor()


# plot temperature
# faces = np.concatenate((triangulation_primary, triangulation_secondary), axis=0)
# temperature = np.concatenate((primary.get_temperature_distribution(), secondary.get_temperature_distribution()), axis=0)
#
# rgb = Fn.arr_to_rainbow(arr=temperature,
#                        minimum=min(temperature),
#                        maximum=max(temperature))
# hex_c = Fn.rgb_to_hex(color=rgb, sharp=True)
# Plt.plot_3d(normals=None, vertices=None, faces=[faces], face_color=[hex_c], normals_view=False,
#                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
#                 face_alpha=1., azim=30, elev=30)


# vrati zrotovany system v tvare
# [[normals_primary, normals_secondary], [faces_primary, faces_secondary]]


# rotation_angle = np.radians(15.0)
rotation_angle = np.radians(170.0)
orbit.set_inclination(np.pi / 2.2)
rot_sys = orbit.rotate_system(faces=[primary.get_faces(), secondary.get_faces()],
                              normals=[primary.get_faces_orientation(), secondary.get_faces_orientation()],
                              rotation_angle=rotation_angle, inclination_rotation=True, faces_rotation=True,
                              inclination=orbit.get_inclination())

# zobrazenie zrotovaneho systemu
# Plt.plot_3d(normals=None, vertices=None, faces=[rot_sys[1][0], rot_sys[1][1]], face_color="w", normals_view=False,
#                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
#                 face_alpha=1., azim=0, elev=0)

# geometria
# darksite filter
# darksite_filter_p = Geo.darkside_filter(faces=primary.get_faces(), normals=primary.get_faces_orientation(),
#                                         verbose=True)
# darksite_filter_s = Geo.darkside_filter(faces=secondary.get_faces(), normals=secondary.get_faces_orientation(),
#                                         verbose=True)  # faces, normals, indices

# darkside filter pre zrotovany system
darksite_filter_p = Geo.darkside_filter(faces=rot_sys[0][0], normals=rot_sys[1][0],
                                        verbose=True)
darksite_filter_s = Geo.darkside_filter(faces=rot_sys[0][1], normals=rot_sys[1][1],
                                        verbose=True)  # faces, normals, indices

# vyplotovanie vyfiltrovanych faziet
# Plt.plot_3d(normals=None, vertices=None, faces=[darksite_filter_p[0], darksite_filter_s[0]], face_color="w", normals_view=False,
#                 points_view=False, faces_view=True, point_color="r", normal_color="w", point_size=3., verbose=True,
#                 face_alpha=1., azim=0, elev=0)

eclipse_filter = Geo.eclipse_filter(verbose=True, actual_distance=1.0,
                                    faces=[darksite_filter_p[0], darksite_filter_s[0]],
                                    normals=[darksite_filter_p[1], darksite_filter_s[1]], orbital_angle=rotation_angle,
                                    indices=[darksite_filter_p[2], darksite_filter_s[2]],
                                    inclination=orbit.get_inclination())

rpwr = primary.radiation_power(passband_model=observer.get_passband_model(),
                               passband_range=observer.get_passband_range(), wavelength_step=10.0)

# print Fn.find_surounded(gv.TEMPERATURE_LIST_ATM, 3500)







elapseEnd = time.time()
print(Fn.color_string("info", "Info: ") + 'Computing elapsed time: ' + str(
    round(elapseEnd - elapseStartClear, 5)) + ' sec')

print(Fn.color_string("info", "Info: ") + Fn.color_string("warning", 'Elapsed time: ' + str(
    round(elapseEnd - elapseStart, 5)) + ' sec'))
