#!/usr/bin/python

""" Geometry trieda
-------------------
  Obsahuje "geometricke" funkcie

Inicializacne parametre:
------------------------
  Staticka trieda, ziadne inicializacne parametre

Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
------------------------------------------------------------------------------------------------------

  find_nearest_dist_3d()
  ======================

    Vstupne parametre:
    ------------------
    data : [python list], defaultna hodnota prazdny python list [], pythonovy list 3D bodov, pre urcenie najkratsej
                          vzdialenosti medzi nimi, struktura listu: [[x0, y0, z0], [x1, y1, z1], ..., [xn, yn, zn]]

    Return:
    -------
    [tuple], vrati najkratsiu vzdialenost medzi 3D bodmi na vstupe

    Popis:
    ------
    vrati najkratsiu vzdialenost medzi 3D bodmi na vstupe

    Bug:
    ----
    ziaden znamy bug

  convex_hull_triangulation()
  ===========================

    Vstupne parametre:
    ------------------
    xyz : [numpy array], defautlna hodnota "[]", numpy array bodov v priestore k triangulacii, body z triedy Binary,
                         funkcie get_3d_model_optimized()
    verbose : [bool], defaultna hodnota "False", premenna pre zapinanie a vypinanie zobrazovania priebehu vykonavania
                      kodu

    Return:
    -------
    [numpy array], v tvare [[[x0, y0, z0], [x1, y1, z1], [x2, y2, z2]], ... atd] v priklade zobrazeny 1 trojuholnik

    Popis:
    ---------
    funkcia pre dalaunayovu triangulaciu convexneho point cloudu

    Bug:
    ----
    ziaden znamy bug

  Funkcia()
  =========
    Vstupne parametre:
    ------------------
    name : [type], description

    Return:
    -------
    void

    Popis:
    ------

    Bug:
    ----
    ziaden znamy bug
"""

import numpy as np
import scipy as sp

import globe.variables as gv
import objects.Function as Fn
import warnings
import objects.Iostream as Io
import objects.Orbit as Orb
import os

import matplotlib.path as mpltpath
import intersection.edge_intersection as ei
import intersection.plane_line_3d_intesection as pli

from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from copy import copy


def find_nearest_dist_3d(data=None):
    from scipy.spatial import KDTree
    points = data[:]
    test_points, distances = points[:], []

    for i in range(0, len(test_points) - 1):
        points.remove(test_points[i])
        tree = KDTree(points)
        distance, ndx = tree.query([test_points[i]], k=1)
        distances.append(distance[0])
    return min(distances)


def average_spacing(data=None, verbose=False, neighbours=6):
    if type(data) != type(np.array([])):
        data = np.array(data)

    # ########################## !!!!!!!!!!!!!!!!! IMPORTANT
    # neighbours je premenna zhodna s premennou v cgal funkcii
    from scipy.spatial import distance
    try:
        dist = sp.spatial.distance.cdist(data, data, 'euclidean')
        tot = 0
        for line in dist:
            tot += np.sort(line)[1:1 + neighbours].sum() / (neighbours + 1)
        return tot / dist.shape[0]
    except:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: average_spacing(), line: " + str(
                Fn.lineno()) + ". Error has been occurred.")
        return False


def convex_hull_triangulation(vertices=None, verbose=False):
    if verbose:
        print(Fn.color_string(color="info",
                              string="Info: ") + "Triangulation is running, used function is convex_hull_triangulation().")

    empty_list_err = False
    if Fn.is_numpy_array(vertices):
        if Fn.is_numpy_array_empty(vertices):
            empty_list_err = True

    elif Fn.is_python_list(vertices):
        if not vertices:
            empty_list_err = True
    else:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: convex_hull_triangulation(), line: " + str(
                Fn.lineno()) + ". Variable `vertices` is invalid.")
        return False
    if empty_list_err:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="EmptyVariableError: ") + "In reference: Geometry, function: convex_hull_triangulation(), line: " + str(
                Fn.lineno()) + ". Nothing to triangulate.")
        return False

    try:
        from scipy.spatial import Delaunay
        if len(vertices) >= 3:
            vertices = np.array(vertices)
            triangulation = Delaunay(vertices)
            convex_hull = triangulation.convex_hull
            vertices = vertices[convex_hull]
            # return simplices, vertices
            return convex_hull, np.array(vertices)
    except:
        if verbose:
            print('Triangulation failed' + ' [' + gv.COLOR_ERROR + ' FAIL ' + gv.COLOR_END + ']')
        return False


def cgal_triangulation(
        normals=None,
        points=None,
        verbose=False,
        min_triangle_angle=0.349066,  # a lower bound on the minimum angle in degrees of the surface mesh facets.
        # an upper bound on the radius of surface Delaunay balls. A surface Delaunay ball is a ball circumscribing
        # a facet, centered on the surface and empty of vertices. Such a ball exists for each facet of the current
        # surface mesh. Indeed the current surface mesh is the Delaunay triangulation of the current sampling restricted
        # to the surface which is just the set of facets in the three dimensional Delaunay triangulation of the sampling
        # that have a Delaunay surface ball.
        max_triangle_size=1,
        # an upper bound on the center-center distances of the surface mesh facets. The center-center distance of
        # a surface mesh facet is the distance between the facet circumcenter and the center of its surface Delaunay ball.
        surface_aproximation_error=0.375,
        to_average_spacing=5
        # zasadne pouizvat hodnotu 5 pre average spacing, pretoze to zodpoveda rozumenj vzdialenosti
        # dvoch blizkych bodoch na povrchu

):
    if verbose:
        print(Fn.color_string(color="info",
                              string="Info: ") + "Triangulation is running, used function is cgal_triangulation().")

    # variables test
    if Fn.empty(var=normals):
        if verbose:
            print("EmptyListError: In reference: Geometry, function: cgal_triangulation(), line: " + str(
                Fn.lineno()) + ". Variable normals is empty [ " + gv.COLOR_ERROR + 'ERROR' + gv.COLOR_END + ' ]')
        return False

    if Fn.empty(var=points):
        if verbose:
            print("EmptyListError: In reference: Geometry, function: cgal_triangulation(), line: " + str(
                Fn.lineno()) + ". Variable points is empty [ " + gv.COLOR_ERROR + 'ERROR' + gv.COLOR_END + ' ]')
        return False

    # ak sa neulozi subor s bodmi a normalami, tak funkcia vrati False
    if not Io.save_cgal_pcd_with_normals(filename="input.xyz", filedir="tmp/", verbose=verbose, points=points,
                                         normals=normals):
        return False



    # Min triangle angle
    # Max triangle size w.r.t. point set average spacing;
    # Surface Approximation error w.r.t. point set average spacing
    # to average spacing

    # tato krkolomna konstrukcia ''np.degrees([min_triangle_angle])[0]'' tu je preto, aby pycharm nepapuloval
    if verbose:
        terminal = 'bin/poisson_reconstruction ' + str(np.degrees([min_triangle_angle])[0]) + ' ' \
                   + str(max_triangle_size) + ' ' + str(surface_aproximation_error) + ' ' + str(to_average_spacing)
    else:
        terminal = 'bin/poisson_reconstruction ' + str(np.degrees([min_triangle_angle])[0]) + ' ' \
                   + str(max_triangle_size) + ' ' + str(surface_aproximation_error) + ' ' + str(to_average_spacing) \
                   + ' > /dev/null 2>&1'

    os.system(terminal)  # returns the exit status
    os.remove('tmp/input.xyz')
    if os.path.isfile('tmp/output.off'):
        tri = Io.load_cgal_3d_mesh(filename='output.off', verbose=verbose)
        os.remove('tmp/output.off')
    else:
        if verbose:
            print("Error: In reference: Geometry, function: cgal_triangulation(), line: " + str(
                Fn.lineno()) + ". File tmp/output.off missing in triangulation process [" + gv.COLOR_ERROR,
                  "ERROR" + gv.COLOR_END + " ]")
        return False

    if not Fn.empty(var=tri):
        return tri[2], tri[0], tri[1]
    else:
        if verbose:
            print("Error: In reference: Geometry, function: cgal_triangulation(), line: " + str(
                Fn.lineno()) + ". Error has been occured during cgal triangulation process [" + gv.COLOR_ERROR,
                  "ERROR" + gv.COLOR_END + " ]")
        return False


def center_of_mass(
        faces=None,
        verbose=False
):
    if not Fn.is_numpy_array(faces):
        faces = np.array(faces)
    if Fn.numpy_array_is_empty_verbose(arr=faces, function_name="center_of_mass()", class_name="Geometry",
                                       var_name="faces", verbose=verbose, line=str(Fn.lineno())):
        return False
    return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])


def normal_estimation(
        binary_object=None,
        actual_distance=None,
        vertices=None,
        t_object="primary",
        mode="in_center",  # in_center/in_point
        verbose=False
):
    if verbose:
        print(Fn.color_string(color="info", string="Info: ") + "Normal estimation is running. Object: " + str(
            t_object) + ".")

    warnings.simplefilter(action="ignore", category=FutureWarning)
    # kontrola premennych
    if binary_object is None:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="ValueError: ") + "In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Binary object is set to None.")
        return False

    if actual_distance is None:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="ValueError: ") + "In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Variable `actual_distance` is set to None.")
        return False

    try:
        # ten try tu je z toho dovodu, ze ked actual_distance je nejakym nedopatrenim uzivatela string, tak np.isnan
        # sa zdrbe a hodi error, tak to treba proste obist cez vynimku a hodit si vlastny error
        if np.isnan(actual_distance):
            if verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In reference: Geometry, function normal_estimation(), variable actual_distance is NaN.")
            return False
    except:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="ValueError: ") + "In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Variable `actual_distance` is not a number.")
        return False

    if Fn.numpy_array_is_empty_verbose(arr=vertices, function_name="normal_estimation", class_name="Geometry",
                                       var_name="vertices", verbose=verbose, line=str(Fn.lineno())):
        return False

    # ak je zle nastaveny objekt pre vypocet
    if t_object != "primary" and t_object != "secondary":
        if verbose:
            print(Fn.color_string(color="error",
                                  string="ValueError: ") + "In reference: Geometry, function: normal_estimation(), line: " + str(
                Fn.lineno()) + ". Variable `object` is invalid, use `primary` or `secondary`.")
        return False

    # koniec kontroly premennych
    if t_object == 'primary':
        synchronicity_parameter = binary_object.primary.synchronicity_parameter
    elif t_object == 'secondary':
        synchronicity_parameter = binary_object.secondary.synchronicity_parameter
    mass_ratio, vector, index = binary_object.mass_ratio, [], 0

    # tie ifka su to preto, ze mam zadefinovany parameter mozno pre buduce pouzitie vo funkcii, a ked nie je pouzity,
    # tak ten priblby PyCharm pyskuje
    if mode == "in_point":
        pass
    if mode == "in_center":
        pass

    for point in vertices:
        try:
            x = point[0]
            r_square = x ** 2 + point[1] ** 2 + point[2] ** 2
            rw_square = (actual_distance - x) ** 2 + point[1] ** 2 + point[2] ** 2

            # vypocet derivacie podla premennej x
            # je potrebne rozlisovat primarnu a sekundarnu zlozku
            if t_object == "primary":
                der_x = - (x / r_square ** (3.0 / 2.0)) + (
                    (mass_ratio * (actual_distance - x)) / rw_square ** (
                        3.0 / 2.0)) + synchronicity_parameter ** 2 * (
                    mass_ratio + 1) * x - (mass_ratio / actual_distance ** 2)
            else:
                der_x = - (x / r_square ** (3.0 / 2.0)) + (
                    (mass_ratio * (actual_distance - x)) / rw_square ** (
                        3.0 / 2.0)) - synchronicity_parameter ** 2 * (
                    mass_ratio + 1) * (1. - x) + (1. / actual_distance ** 2)

            # pre derivacie podla premennej y a z ine je potrebne rozlisovat, ci sa pocitaju pre sekudnarnu alebo
            # pre primarnu zlozku
            der_y = - point[1] * ((1.0 / r_square ** (3.0 / 2.0)) + (mass_ratio / rw_square ** (3.0 / 2.0)) - (
                (synchronicity_parameter ** 2) * (mass_ratio + 1)))
            der_z = - point[2] * ((1.0 / r_square ** (3.0 / 2.0)) + (mass_ratio / rw_square ** (3.0 / 2.0)))

            # vector.append([-der_x, -der_y, -der_z, index])
            vector.append([-der_x, -der_y, -der_z])

        except:
            if verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In reference: Geometry, function: normal_estimation(), line: " + str(
                    Fn.lineno()) + ". Cannot estimate normal vector in step " + str(index) + ".")
            return False
        index += 1

    return np.array(vector)


def face_orientation_a(
        face=None,
        t_object="primary",
        actual_distance=None,
        verbose=False
):
    ################################################################################################################
    #                                                                                                              #
    #                                    FUNKCIA JE LEN POLOFUNKCNY PROTOTYP                                       #
    #                                                                                                              #
    ################################################################################################################

    if Fn.empty(actual_distance) and t_object == "secondary":
        if verbose:
            print("Error: In reference: Geometry, function: face_orientation(), line: " + str(
                Fn.lineno()) + ". Variable actual_distance is empty [ " + gv.COLOR_ERROR + "ERROR" + gv.COLOR_END + " ]")
        return False

    com, vec = center_of_mass(faces=face, verbose=verbose), []

    center = np.array([0., 0., 0.])
    if t_object == "secondary":
        center[0] = actual_distance

    if len(com) != len(face):
        if verbose:
            print("Error: In reference: Geometry, function: face_orientation(), line: " + str(
                Fn.lineno()) + ". com variable has a different length as variable face [ " + gv.COLOR_ERROR + "ERROR" + gv.COLOR_END + " ]")
        return False

    for idx in range(0, len(com)):
        vec_a = (face[idx][0] - face[idx][1])
        vec_b = (face[idx][1] - face[idx][2])
        cross_p = np.cross(vec_a, vec_b)

        # ked sa spocita skalarny sucin, nie je mozne (resp. neviem) urcit ci vektor smeruje do stredu alebo smeruje
        # von, cize ci sa jedna o vonkajsiu alebo vnutornu normalu; ja potrebujem pre pracu bud vsetky vonkajsie alebo
        # vsetky vnutorne, co je problem lebo tu moze vzniknut mix; myslienka ako to eliminovat, je porovnat
        # suradnice taziska trojuholnika so suradnicami zratanej normaly; ak x-ova suradnica taziska je v zapornej casti
        # priestoru, tak je jasne, ze aj x-ova suradnica normaly musi smerovat do zapornej casti (normaly maju pociatok
        # v strede suradnicoveho systemu); rovnako tak to plati pre vsetky dalsie suradnice; problem sa moze vyskytnut
        # v pripade ze x == 0, alebo y == 0 alebo z == 0, ci uz pre normalu alebo tazisko; staci to rozhodnut pre jednu
        # suradnicu;
        # problem je treba rozdelit na primarnu a sekundarnu zlozku, pretoze sekundarna sa neda porovnavat s nulou,
        # ale je nutne ju porovnavat s aktualnou vzdialenostou
        # TOTO ZAKOMENTOVANE TU NECHAVAM PRE UPLNOST, MOZNO MI V BUDUCNOSTI NIECO NAPADNE, CO BY TO ROZCHODILO,
        # ALE TAK AKO TO JE, POROVNAVANIE TAZISKA A NORMALY NEFUNGUJE V NIEKTORCYH EXTREMNYCH PRIPADOCH (EXTREMNYCH
        # V ZMYSLE AK SA TLACI NIEKTORA SURADNICA BLIZKO 0.0, HLAVNE x)
        # PRE SEKUNDARNU ZLOZKU NEBOL KOD TESTOVANY
        #
        # neviem, ci nebude treba pridat osetrenie, ak su vsetky suradnice nulove (taky pripad je invalidny)
        # if ((cross_p[0] < 0 < com[idx][0]) or (cross_p[1] < 0 < com[idx][1]) or (
        #                 cross_p[2] < 0 < com[idx][2])) or \
        #                 ((cross_p[0] > 0 > com[idx][0]) or (cross_p[1] > 0 > com[idx][1]) or (
        #                                 cross_p[2] > 0 > com[idx][2])) and \
        #                         t_object is "primary":
        #     cross_p *= -1
        #
        # if ((cross_p[0] < 0 < com[idx][0] - actual_distance) or (cross_p[1] < 0 < com[idx][1]) or (
        #         cross_p[2] < 0 < com[idx][2])) or (
        #         (cross_p[0] > 0 > com[idx][0] - actual_distance) or (cross_p[1] > 0 > com[idx][1]) or (
        #         cross_p[2] > 0 > com[idx][2])) and t_object is "secondary":
        #     cross_p *= -1


        # myslienka zistovania vonkajsej a vnutronej normaly je v principe jednoducha; spravi sa testovaci vektor stale
        # zo stredu hviezdy po tazisko trojuholnika a spravi sa dot produkt tohto vektora so spocitanou normalou
        # ak je dot produkt zaporny, tak vektory zvieraju uhol vacsi ako np.pi / 2.0 a teda normala je vnutrona
        # typ pado, treba zmenit jej znamienko, aby sme dostali vonkajsiu
        #           ^ n
        #           |
        #           |
        #    _______|com____ face
        #           ^
        #          /
        #         /
        #        /
        #
        #      x O
        # ak je uhol medzi normalou vacsi ako np.pi / 2.0 teda, dot produkt je zaporny, potom je to vonkajsia normala
        # a vice-versa

        ################################################################################################################
        # 1) kod testovany na odddelenych systemoch a tam funguje
        # 2) kod nebol zatial otestovany na kontaktnych systemoch

        test_v = com[idx] - center

        if np.dot(test_v, cross_p) < 0:
            cross_p *= -1.0

        vec.append(cross_p)

    return np.array(vec)


def face_orientation_beta(
        faces=None,
        t_object="primary",
        binary_object=None,
        actual_distance=None,
        verbose=False
):
    ################################################################################################################
    #                                                                                                              #
    #                                    FUNKCIA JE LEN POLOFUNKCNY PROTOTYP                                       #
    #                                                                                                              #
    ################################################################################################################
    if verbose:
        print(Fn.color_string(color="info", string="Info: ") + "Resolving faces orientation. Object: " + str(
            t_object) + ".")

    if Fn.empty(actual_distance):
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: face_orientation_prototype(), line: " + str(
                Fn.lineno()) + ". Variable actual_distance is empty.")
        return False

    if binary_object is None:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="ValueError: ") + "In reference: Geometry, function: face_orientation_prototype(), line: " + str(
                Fn.lineno()) + ". Binary object is set to None.")
        return False

    # pocita sa tu com trojuholnikov povrchu a v tom sa budu ratat porovnavacie normaly, nie je to ale nutne a stacil
    # by vrchol, ale pride mi to elegantnejsie aj ked to spomaluje beh programu XD
    com, vec = center_of_mass(faces=faces, verbose=verbose), []

    if len(com) != len(faces):
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: face_orientation_prototype(), line: " + str(
                Fn.lineno()) + ". Variable `com` has a different length as variable `face`.")
        return False

    # vypocet porovnavacich normalovych vektorov
    # premenne t_object, binary_object, actual_distence prichadzaju ako externe argumenty z tejto funkcie
    fn_normals = normal_estimation(binary_object=binary_object, actual_distance=actual_distance, vertices=com,
                                   t_object=t_object, mode="in_point", verbose=verbose)

    # v tomto for cykle sa normalne spocitaju normaly pomocou vektoroveho sucinu
    # ich smerovanie (ci sa jedna o vonkajsiu alebo vnutornu normalu) sa zisti porovnanim
    for idx in range(0, len(com)):
        cross_p = np.cross((faces[idx][0] - faces[idx][1]), (faces[idx][1] - faces[idx][2]))

        # v tejto podmienke sa otestuje, ci sa jedna o vnutornu aleb vonkajsiu normalu
        # ak je dot produkt s vektorom spocitany cez funkciu (derivaciu potencialu) kladny, tak sa jedna o vonkajsiu
        # normalu, vtedy je vsetko vporiadku, inak treba prehodit znamienko normaly spocitanej cez ako scross produkt
        if np.dot(fn_normals[idx], cross_p) < 0: cross_p *= -1
        vec.append(cross_p)

    # if average_norm:
    #     # normaly pre vypocet noriem (cize pre teplotu, grav. zrychlenie atd.) sa zoberu ako priemerne hodnoty normal
    #     # z vrcholov jednotlivych faziet
    #     try:
    #         vertices = faces.reshape(len(faces) * 3, 3)
    #         normals = normal_estimation(binary_object=binary_object, actual_distance=actual_distance,
    #                                     vertices=vertices,
    #                                     t_object=t_object, mode="in_point", verbose=verbose)
    #
    #         fn_normals = np.array([(normals[idx * 3] + normals[idx * 3 + 1] + normals[idx * 3 + 2]) / 3.0 for idx in
    #                                range(0, len(normals) / 3)])
    #     except:
    #         if verbose:
    #             print Fn.color_string(color="error", string="Error: ") + \
    #                   "In reference: Geometry, function: face_orientation_prototype(), line: " + str(
    #                 Fn.lineno()) + ". Error has been occurred during average norm computation."
    #         return False

    return np.array(vec)


def vector_array_normalisation(
        vector_arr=None,
        multi=1.0,
        verbose=False
):
    # variables test
    if Fn.empty(var=vector_arr):
        if verbose:
            print(Fn.color_string(color="error",
                                  string="EmptyVariableError: ") + "In reference: Geometry, function: vector_array_normalisation(), line: " + str(
                Fn.lineno()) + ". Variable `vector_arr` is empty.")
        return False

    # ak sa jedna o dictionary
    if not Fn.is_numpy_array(arr=vector_arr):
        vector_arr = np.array(vector_arr)
    try:
        normalisated = np.array([vector / (np.linalg.norm(x=vector) * multi) for vector in vector_arr])
    except:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="EmptyVariableError: ") + "In reference: Geometry, function: vector_array_normalisation(), line: " + str(
                Fn.lineno()) + ". Error occurred during normalisation process.")
        return False

    return np.array(normalisated)


def vector_array_translation(
        vector_arr=None,
        translation_arr=None,
        verbose=False
):
    # variables test
    if Fn.empty(var=vector_arr):
        if verbose:
            print(Fn.color_string(color="error",
                                  string="EmptyVariableError: ") + "In reference: Geometry, function: vector_array_translation(), line: " + str(
                Fn.lineno()) + ". Variable `vector_arr` is empty.")
        return False

    if Fn.empty(var=translation_arr):
        if verbose:
            print(Fn.color_string(color="error",
                                  string="EmptyVariableError: ") + "In reference: Geometry, function: vector_array_translation(), line: " + str(
                Fn.lineno()) + ". Variable `translation_arr` is empty.")
        return False

    # porovnanie dlzky poli
    # bude treba dopisat funkciu
    try:
        translated = [np.array(vector_arr[idx]) + np.array(translation_arr[idx]) for idx in range(len(vector_arr))]
    except:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: vector_array_translation(), line: " + str(
                Fn.lineno()) + ". Error occured during translation process.")
        return False

    return np.array(translated)


def darkside_filter(
        to_spectator=None,
        normals=None,
        faces=None,
        verbose=False
):
    if verbose:
        print(Fn.color_string("info", "Info: ") + "Darkside filter is running.")

    if to_spectator is None:
        to_spectator = np.array([1, 0, 0])

    if Fn.empty(normals) or Fn.empty(faces):
        if verbose:
            print(Fn.color_string("error",
                                  "EmptyVariableError: ") + "In reference: Geometry, function: darkside_filter(), line: " + str(
                Fn.lineno()) + ". Variables normals or faces is empty.")
        return False

    if len(normals) != len(faces):
        if verbose:
            print(
            Fn.color_string("error", "Error: ") + "In reference: Geometry, function: darkside_filter(), line: " + str(
                Fn.lineno()) + ". Variables normals and faces has not same length.")
        return False

    normals_filtered, faces_filtered, indices_filtered, critical_angle = [], [], [], np.pi / 2.0

    for t_object in range(0, len(normals)):
        t_normals, t_faces, t_indices = [], [], []

        for normal, face, idx in list(zip(normals[t_object], faces[t_object], np.arange(0, len(faces[t_object]), 1))):
            delta = angle(u=normal, v=to_spectator)
            if delta <= critical_angle:
                t_normals.append(normal)
                t_faces.append(face)
                t_indices.append(idx)

        normals_filtered.append(np.array(t_normals))
        faces_filtered.append(np.array(t_faces))
        indices_filtered.append(np.array(t_indices))

    return [np.array(faces_filtered), np.array(normals_filtered), np.array(indices_filtered)]


def eclipse_filter(
        normals=None,
        indices=None,
        vertices=None,
        simplices=None,
        orbital_angle=None,
        verbose=False
):

    if verbose:
        print(Fn.color_string("info", "Info: ") + "Eclispe filter is running.")

    # prototyp navratoveho pola
    return_array = [[[], []], [[], []], [[], []]]
    # 0 - indices of native created triangulation
    # 1 - surface area of those faces
    # 2 - faces (just for visualization)

    # primarna zlozka prakryva sekundarnu
    if np.cos(orbital_angle) < 0:
        f, b = 0, 1
    else:
        f, b = 1, 0

    front = {"faces": np.array(vertices[f])[simplices[f]],
             "faces2d": None,
             "norms": normals[f],
             "indices": indices[f],
             "vertices": vertices[f],
             "vertices2d": None,
             "simplices": simplices[f]}
    behind = {"faces": np.array(vertices[b])[simplices[b]],
              "faces2d": None,
              "norms": normals[b],
              "indices": indices[b],
              "vertices": vertices[b],
              "vertices2d": None,
              "simplices": simplices[b]}

    # projekcia do roviny zy
    front["faces2d"] = [[np.array([v[1], v[2]]) for v in f] for f in front["faces"]]
    behind["faces2d"] = [[np.array([v[1], v[2]]) for v in f] for f in behind["faces"]]


    front["vertices2d"] = [[f[1], f[2]]
                           for f in np.array(front["vertices"])[np.unique(front["simplices"])]]

    behind["vertices2d"] = [[f[1], f[2]]
                           for f in np.array(behind["vertices"])[np.unique(behind["simplices"])]]

    # boundary of fron object
    bound = ConvexHull(front["vertices2d"])
    hull_points = Fn.array_mask(front["vertices2d"], bound.vertices)
    bb_path = mpltpath.Path(np.array(hull_points))

    # save indices and surface area of faces for front element
    return_array[0][f], return_array[1][f] = front["indices"], triangle_surface_area(triangle=front["faces"],
                                                                                     verbose=verbose)
    return_array[2][f] = front["faces"]

    # odstranit elementy (wuma problem), ktore neboli vyfiltrovane pri dakrside filtri na zlozke v popredi





    for f, i, idx in list(zip(behind["faces2d"], range(0, len(behind["faces2d"])), behind["indices"])):
        # test vertices of each face if is inside of bb_path
        p_inside = [bb_path.contains_points([p])[0] for p in f]
        outside = [None for j in p_inside if j == False]

        if len(outside) == 3:
            return_array[0][b].append(idx)
            return_array[1][b].append(triangle_surface_area(triangle=[behind["faces"][i]], verbose=verbose)[0])
            return_array[2][b].append(behind["faces"][i])

        elif len(outside) == 2:
            # len jeden bod je v zakryte

            # identifikacia indexof faziet podla polohy (v zakryte / mimo zakryt)
            inside_idx = Fn.indices(p_inside, True)[0]
            outside_idx = Fn.indices(p_inside, False)

            # segmenty fazety pre hladanei prieniku s 2d hranicou "front" zlozky
            segments = np.array([[f[inside_idx], f[outside_idx[0]]], [f[inside_idx], f[outside_idx[1]]]])

            # test first and second face segment with boundary
            zy = []
            for seg in range(0, 2):
                # iterate huull_points (boundary)
                for x in range(-1, len(hull_points) - 1):
                    b_segment = [hull_points[x], hull_points[x + 1]]
                    intersection = ei.edge_intersection_2d(b_segment[0], b_segment[1],
                                                           segments[seg][0], segments[seg][1])
                    if not intersection[1] != True:
                        zy.append([intersection[2], intersection[3]])
                        break

            # create entire set of points of face which are outside of boundary and intersection points
            zy = np.concatenate((zy, [f[outside_idx[0]]], [f[outside_idx[1]]]), 0)
            delaunay_simplices = Delaunay(zy).simplices
            # delaunay_faces = zy[delaunay_simplices]

            # faces from 2d to 3d (order with respect to zy list, there was added in the same order)
            delaunay_vertices_3d = [None, None, behind["faces"][i][outside_idx[0]], behind["faces"][i][outside_idx[1]]]


            # najdenie 3d projekcie zodpovedajucej 2d bodu v priesecniku boundary so segmentami fazety
            # planeline_intersection berie na vstup
            # (bod_urcujuci_priamku_a, bod_urcujuci_priamku_b, bod_B_urcujuci_rovinu,
            # normalovy_vektor_urcujuci_rovinu_v_bode_B)

            for inter, order in list(zip(zy[0:2], range(0, 2))):
                intersection = pli.planeline_intersection([0.0, inter[0], inter[1]], [1.0, inter[0], inter[1]],
                                                          behind["faces"][i][inside_idx],
                                                          np.cross(np.array(behind["faces"][i][inside_idx]) -
                                                                   np.array(behind["faces"][i][outside_idx[0]]),
                                                                   np.array(behind["faces"][i][inside_idx]) -
                                                                   np.array(behind["faces"][i][outside_idx[1]])))
                # priradenie 3d bodu na spravne miesto
                delaunay_vertices_3d[order] = intersection
            # convert 3d points set back to faces
            delaunay_faces_3d = np.array(delaunay_vertices_3d)[delaunay_simplices]

            return_array[0][b].append(idx)

            return_array[1][b].append(
                sum([triangle_surface_area(triangle=[val], verbose=verbose)[0] for val in delaunay_faces_3d]))

            for t in delaunay_faces_3d:
                return_array[2][b].append(t)

        elif len(outside) == 1:
            inside_idx = Fn.indices(p_inside, True)
            outside_idx = Fn.indices(p_inside, False)[0]
            segments = [[f[outside_idx], f[inside_idx[0]]], [f[outside_idx], f[inside_idx[1]]]]

            zy = []

            for seg in range(0, 2):
                for x in range(-1, len(hull_points) - 1):
                    b_segment = [hull_points[x], hull_points[x + 1]]
                    intersection = ei.edge_intersection_2d(b_segment[0], b_segment[1],
                                                           segments[seg][0], segments[seg][1])
                    if not intersection[1] != True:
                        zy.append([intersection[2], intersection[3]])
                        break

            zy = np.concatenate((zy, [f[outside_idx]]), 0)

            # len jeden index je mimo v takomto pripade a +2 priesecniky automaticky urcuju trojuholnik
            delaunay_faces_3d = [None, None, behind["faces"][i][outside_idx]]

            # najst priesecnik v 3d priestore pre body na hranici
            for inter, order in list(zip(zy[0:2], range(0, 2))):
                intersection = pli.planeline_intersection([0.0, inter[0], inter[1]], [1.0, inter[0], inter[1]],
                                                          behind["faces"][i][outside_idx],
                                                          np.cross(np.array(behind["faces"][i][outside_idx]) -
                                                                   np.array(behind["faces"][i][inside_idx[0]]),
                                                                   np.array(behind["faces"][i][outside_idx]) -
                                                                   np.array(behind["faces"][i][inside_idx[1]])))
                delaunay_faces_3d[order] = intersection

            return_array[0][b].append(idx)
            return_array[1][b].append(triangle_surface_area(triangle=[delaunay_faces_3d], verbose=verbose)[0])
            return_array[2][b].append(delaunay_faces_3d)


    # prekonvertovanie vsetkeho an numpy aray
    return_array[0][0] = np.array(return_array[0][0])
    return_array[0][1] = np.array(return_array[0][1])
    return_array[1][0] = np.array(return_array[1][0])
    return_array[1][1] = np.array(return_array[1][1])
    return_array[2][0] = np.array(return_array[2][0])
    return_array[2][1] = np.array(return_array[2][1])

    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    # fig, ax = plt.subplots()
    #
    # patches = []
    # for f in behind["faces2d"]:
    #     polygon = Polygon(f, True)
    #     patches.append(polygon)
    #
    # p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.3)
    # p.set_color("r")
    # ax.add_collection(p)
    #
    # patches = []
    # for f in front["faces2d"]:
    #     polygon = Polygon(f, True)
    #     patches.append(polygon)
    #
    # p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.3)
    # p.set_color("b")
    # ax.add_collection(p)
    #
    # # colors = 100 * np.random.rand(len(patches))
    # # p.set_array(np.array(colors))
    #
    # ax.add_collection(p)
    # ax.set_xlim([-.5, .5])
    # ax.set_ylim([-.5, .5])
    # plt.gca().set_aspect('equal')
    #
    # plt.show()

    del(front, behind, bb_path, hull_points)
    return return_array




def angle(
        u=None,
        v=None,
        verbose=False
):
    try:
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return np.arccos(c)
    except:
        if verbose:
            print(Fn.color_string("error", "Error: ") + "In reference: Geometry, function: angle(), line: " + str(
                Fn.lineno()) + ".")
        return False


def triangle_surface_area(
        triangle=None,
        verbose=False
):
    if Fn.empty(triangle):
        if verbose:
            print(Fn.color_string("error",
                                  "Error: ") + "In reference: Geometry, function: triangle_surface_area(), line: " + str(
                Fn.lineno()) + ". Variable triangle is empty.")
        return False
    try:
        surface = [np.linalg.norm(np.cross(np.subtract(triangle[idx][0], triangle[idx][1]),
                                           np.subtract(triangle[idx][0], triangle[idx][2]))) / 2.0 for idx in
                   range(len(triangle))]
    except:
        if verbose:
            print(Fn.color_string("error",
                                  "Error: ") + "In reference: Geometry, function: triangle_surface_area(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during computing process.")
        return False

    return np.array(surface)


def gradient_norm(
        faces=None,
        binary_object=None,
        actual_distance=None,
        t_object="primary",
        verbose=False
):
    if Fn.empty(faces):
        if verbose:
            print(Fn.color_string("error",
                                  "EmptyVariableError: ") + "In reference: Geometry, function: gradient_norm(), line: " + str(
                Fn.lineno()) + ". Variable faces is empty or invalid.")
        return False
    try:
        vertices = faces.reshape(len(faces) * 3, 3)
        gradient = normal_estimation(t_object=t_object, actual_distance=actual_distance, binary_object=binary_object,
                                     verbose=verbose, vertices=vertices)
        gradnorm = Fn.gradient_norm(gradient=gradient, verbose=verbose)
        return np.array([(gradnorm[idx * 3] + gradnorm[idx * 3 + 1] + gradnorm[idx * 3 + 2]) / 3.0 for idx in
                         range(int(len(gradnorm) / 3))])
    except:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Geometry, function: gradient_norm(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during gradient norm computation.")
        return False


def cgal_separation(cgal_simplex=None, x_separation=None):
    x_separation = x_separation[0]

    # face orientation
    face_orientation = []
    for face in cgal_simplex[1]:
        left = []
        for v in face:
            a = True if v[0] <= x_separation else False
            left.append(a)

        on_left = Fn.indices(left, True)
        if len(on_left) == 3 or len(on_left) == 2:
            face_orientation.append(face_orientation_a([face], "primary", 1.0)[0])
        else:
            face_orientation.append(face_orientation_a([face], "secondary", 1.0)[0])

    output = {"faces": {"primary": [], "secondary": []},
              "simplices": {"primary": [], "secondary": []},
              "vertices": {"primary": [], "secondary": []}}

    simplex_map = {"primary": {}, "secondary": {}}
    simplex_counter = {"primary": 0, "secondary": 0}

    for t, s, n in list(zip(cgal_simplex[1], cgal_simplex[0], face_orientation)):

        left = []
        for vertex in t:
            a = True if vertex[0] <= x_separation else False
            left.append(a)

        on_right = Fn.indices(left, False)
        on_left = Fn.indices(left, True)

        if len(on_left) == 3:
            # musi sa zmenit cela struktura, pretoze body sa od seba oddelia, takze uz simplex sediet nebude
            # rovnako tak pribudnu nove body potom

            output["faces"]["primary"].append(t)

            c_simplex = []
            for sim in s:
                try:
                    # skusi ci uz je v mape taky index z povodneho simplexu a ak nie je tak sa vykona
                    # exception;
                    # v opacnom pripade sa zapise o ktory simplex sa jedna, ak uz vertex je v poli

                    sm = simplex_map["primary"][sim]
                    c_simplex.append(sm)
                except KeyError:
                    output["vertices"]["primary"].append(cgal_simplex[2][sim])
                    simplex_map["primary"][sim] = simplex_counter["primary"]
                    c_simplex.append(simplex_counter["primary"])

                    simplex_counter["primary"] += 1
            output["simplices"]["primary"].append(c_simplex)

        elif len(on_left) == 2:
            # tu pride magia
            n_angle = angle([0., -1., 0.], n)
            inverse_rotation = False if n[2] > 0 else True
            t_rotated = [Fn.rotate(n_angle, r, inverse_rotation, "x") for r in t]
            t_projection_2d = [[r[0], r[2]] for r in t_rotated]

            on_left_idx = on_left
            on_right_idx = on_right[0]

            segments = [[t_projection_2d[on_right_idx], t_projection_2d[on_left_idx[0]]],
                        [t_projection_2d[on_right_idx], t_projection_2d[on_left_idx[1]]]]

            # point_to_right je premenna, ktora nesie 3d suradnice bodu, ktory by lezal xovou suradnicou
            # priamo na oddelovaci zloziek systemu a teda by tak vznikli len 2 troj. a nie 3
            #                           | C
            #                           o
            #                           |
            #                A  o       |
            #                           |      o B
            # v tomto priklade by niesol suradnice bodu C

            xz, point_to_right = [], [None, None]

            # najprv sa do xz pripoja intersekcne body
            for seg in range(0, 2):
                # ak xova suradnica laveho bodu je rozna od x_separacie
                # pripade, ze je rovna, tak bod lezi priamo na separacnej hranici tak na lavej strane bude len 1 troj.
                # a nie 2;
                # v pripade, ze oba lave body lezia priamo na separacnej hranici, tak cely trojuholnik patri pravej zlozke
                # a vtedy je xz nulovej dlzky po vystupe z tohto foru
                if segments[seg][1][0] != x_separation:
                    e_interscection = ei.edge_intersection_2d([x_separation, -0.5], [x_separation, 0.5],
                                                              segments[seg][0], segments[seg][1])
                    xz.append([x_separation, e_interscection[3]])
                else:
                    point_to_right = [t_rotated[on_left_idx[seg]], s[on_left_idx[seg]]]



            # len tie body, ktore maju naozaj intersekciu, teda inak myslim pod tym to, ze intersekcia nie je priamo
            # jeden vertex testovaneho trojuhlnika
            intersected = copy(xz)

            if len(intersected) == 0:
                # cely troj. patri pravej zlozke, pretoze testovane "lave" body lezia priamo na separacnej lajne
                #               |
                #               o
                #               |     o
                #               o
                # tento pripad z obrazka

                # ano, tento kod je totalna copypasta z else z nizsia, ale bohuzial sa mi teraz nechce pisat
                # na to funkciu
                output["faces"]["secondary"].append(t)
                c_simplex = []
                for sim in s:
                    try:
                        sm = simplex_map["secondary"][sim]
                        c_simplex.append(sm)
                    except KeyError:
                        output["vertices"]["secondary"].append(cgal_simplex[2][sim])
                        simplex_map["secondary"][sim] = simplex_counter["secondary"]
                        c_simplex.append(simplex_counter["secondary"])
                        simplex_counter["secondary"] += 1
                output["simplices"]["secondary"].append(c_simplex)

                continue

            elif len(intersected) == 1:
                # ak je pretnuty len jedna lajna, teda jeden bod lezi rovno na separacnej hranici,
                # tak sa tento bod zo separacnej a ten co je vlavo pripoji do pola na triangulaciu,
                # cize na vystupe tu budu v "xz" 3 body
                #               |X
                #         C o   x       o B
                #               |
                #               o A
                #               |
                # bod "X" je z intersekcie, to je ten preco len(intersected) == 1 a bod "A" a "C" su tie, ktore
                # sa pripoja do pola
                xz = np.concatenate((xz, [t_projection_2d[on_left_idx[0]]], [t_projection_2d[on_left_idx[1]]]), 0)

            elif len(intersected) == 2:
            # v tomto pripade staci pripojit 2 body z lavej strany, teda z primaru + 2 body z intersekcie, teda budu
            # 2 trojuholniky z tychto 4 bodov v "xz"
            #                 |
            #               o |
            #                 x
            #       o         |
            #                 x     o
            #
                for i in range(0, 2):
                    xz.append(t_projection_2d[on_left_idx[i]])

            delaunay_tri = Delaunay(xz)
            delaunay_sim = delaunay_tri.simplices

            # projekcia intersection bodov do 3d
            # treba zistit spatnu projekciu do 3d

            # delaunay_sim_map je na premapovanie delaunay_sim indexov na indexy z povodneho pola vertexov, ktore
            # prislo na vstupe tejto funkcie
            # ak sa jedna o novy bod, ktory bol vytvoreny tu, ako intersekcny, tak potom je v mape None
            t3d_vertices, delaunay_sim_map = [], []
            for p in xz[0:len(intersected)]:
                intersection_3d = pli.planeline_intersection([p[0], 0.0, p[1]], [p[0], 1.0, p[1]],
                                                             t_rotated[on_right_idx],
                                                             np.cross(
                                                                 np.array(t_rotated[on_right_idx]) -
                                                                 np.array(t_rotated[on_left_idx[0]])
                                                                 ,
                                                                 np.array(t_rotated[on_right_idx]) -
                                                                 np.array(t_rotated[on_left_idx[1]])
                                                             ))

                # prvy bod v tomto zozname bude vzdy patrit intersekcii s x_separation lajnou,
                # pri druhom bode to uz tak byt nemusi, lebo ak je dlzka intersected == 1, tak uz to tak nebude
                t3d_vertices.append(intersection_3d)
                delaunay_sim_map.append(None)


            # v kazdom pripade pokracujuce body budu dva lave body (tie co by mali byt na primarke) z povodneho
            # trojuholnika;
            # toto poradie 3d bodov v t3d_vertices musi byt zhodne s poradim bodov v xz poli,
            # lebo inak nebude sediet simplex z delaunay_sim
            for i in range(0, 2):
                t3d_vertices.append(t_rotated[on_left_idx[i]])
                delaunay_sim_map.append(s[on_left_idx[i]])

            t3d_vertices = [Fn.rotate(n_angle, r, not inverse_rotation, "x") for r in t3d_vertices]
            d3 = np.array(t3d_vertices)[delaunay_sim].tolist()

            for td, sd in list(zip(d3, delaunay_sim)):
                output["faces"]["primary"].append(td)

                c_simplex = []
                for sim in sd:
                    try:
                        sm = simplex_map["primary"][delaunay_sim_map[sim]]
                        c_simplex.append(sm)
                        continue
                    except KeyError:
                        if delaunay_sim_map[sim] is None:
                            output["vertices"]["primary"].append(t3d_vertices[sim])
                        else:
                            output["vertices"]["primary"].append(cgal_simplex[2][delaunay_sim_map[sim]])
                            simplex_map["primary"][delaunay_sim_map[sim]] = simplex_counter["primary"]

                        c_simplex.append(simplex_counter["primary"])
                        simplex_counter["primary"] += 1

                output["simplices"]["primary"].append(c_simplex)

            # este tu treba pridat ten pravy trojuholnik a priradit ho spravne k sekundarnej zlozke
            # START ====================================================================================================
            sd = [0, 1, 2]
            if point_to_right[0] is None and point_to_right[1] is None:
                td = [t3d_vertices[0], t3d_vertices[1], t[on_right_idx]]
                delaunay_sim_map = [None, None, s[on_right_idx]]

            else:
                td = [t3d_vertices[0], point_to_right[0], t[on_right_idx]]
                delaunay_sim_map = [None, point_to_right[1], s[on_right_idx]]

            output["faces"]["secondary"].append(td)
            c_simplex = []
            for sim in sd:
                try:
                    sm = simplex_map["secondary"][delaunay_sim_map[sim]]
                    c_simplex.append(sm)
                    continue
                except KeyError:
                    if delaunay_sim_map[sim] is None:
                        output["vertices"]["secondary"].append(t3d_vertices[sim])
                    else:
                        output["vertices"]["secondary"].append(cgal_simplex[2][delaunay_sim_map[sim]])
                        simplex_map["secondary"][delaunay_sim_map[sim]] = simplex_counter["secondary"]

                    c_simplex.append(simplex_counter["secondary"])
                    simplex_counter["secondary"] += 1

            output["simplices"]["secondary"].append(c_simplex)
            # END ======================================================================================================


        elif len(on_left) == 1:
            # aj tu pride magia
            n_angle = angle([0., -1., 0.], n)
            inverse_rotation = False if n[2] > 0 else True
            t_rotated = [Fn.rotate(n_angle, r, inverse_rotation, "x") for r in t]
            t_projection_2d = [[r[0], r[2]] for r in t_rotated]

            on_left_idx = on_left[0]
            on_right_idx = on_right

            # tu treba len osetrit naviac to, ze ak "x" suradnica pre ten lavy bod == x_separation,
            # tak potom cely trojuholnik patri automaticky sekundarnej zlozke a ziadne intersekcie netreba hladat
            # START ====================================================================================================
            if t[on_left_idx][0] == x_separation:
                output["faces"]["secondary"].append(t)
                c_simplex = []
                for sim in s:
                    try:
                        sm = simplex_map["secondary"][sim]
                        c_simplex.append(sm)
                    except KeyError:
                        output["vertices"]["secondary"].append(cgal_simplex[2][sim])
                        simplex_map["secondary"][sim] = simplex_counter["secondary"]
                        c_simplex.append(simplex_counter["secondary"])
                        simplex_counter["secondary"] += 1
                output["simplices"]["secondary"].append(c_simplex)
                continue
            # END ======================================================================================================

            segments = [[t_projection_2d[on_left_idx], t_projection_2d[on_right_idx[0]]],
                        [t_projection_2d[on_left_idx], t_projection_2d[on_right_idx[1]]]]
            xz = []
            # v tomto pripade sa nemoze stat, ze nejaky bod bude priamo lezat na hranici, lebo to uz bolo
            # poriesene pre elif == 2 a v uvode tohto elif == 1
            for seg in range(0, 2):
                e_interscection = ei.edge_intersection_2d([x_separation, -0.5], [x_separation, 0.5],
                                                          segments[seg][0], segments[seg][1])
                xz.append([x_separation, e_interscection[3]])

            for i in range(0, 2):
                xz.append(t_projection_2d[on_right_idx[i]])

            delaunay_tri = Delaunay(xz)
            delaunay_sim = delaunay_tri.simplices

            t3d_vertices, delaunay_sim_map = [], []
            for p in xz[0:2]:
                intersection_3d = pli.planeline_intersection([p[0], 0.0, p[1]], [p[0], 1.0, p[1]],
                                                             t_rotated[on_left_idx],
                                                             np.cross(
                                                                 np.array(t_rotated[on_left_idx]) -
                                                                 np.array(t_rotated[on_right_idx[0]])
                                                                 ,
                                                                 np.array(t_rotated[on_left_idx]) -
                                                                 np.array(t_rotated[on_right_idx[1]])
                                                             ))

                # prvy bod v tomto zozname bude vzdy patrit intersekcii s x_separation lajnou,
                # pri druhom bode to uz tak byt nemusi, lebo ak je dlzka intersected == 1, tak uz to tak nebude
                t3d_vertices.append(intersection_3d)
                delaunay_sim_map.append(None)

            for i in range(0, 2):
                t3d_vertices.append(t_rotated[on_right_idx[i]])
                delaunay_sim_map.append(s[on_right_idx[i]])

            t3d_vertices = [Fn.rotate(n_angle, r, not inverse_rotation, "x") for r in t3d_vertices]
            d3 = np.array(t3d_vertices)[delaunay_sim].tolist()

            for td, sd in list(zip(d3, delaunay_sim)):
                output["faces"]["secondary"].append(td)

                c_simplex = []
                for sim in sd:
                    try:
                        sm = simplex_map["secondary"][delaunay_sim_map[sim]]
                        c_simplex.append(sm)
                        continue
                    except KeyError:
                        if delaunay_sim_map[sim] is None:
                            output["vertices"]["secondary"].append(t3d_vertices[sim])
                        else:
                            output["vertices"]["secondary"].append(cgal_simplex[2][delaunay_sim_map[sim]])
                            simplex_map["secondary"][delaunay_sim_map[sim]] = simplex_counter["secondary"]

                        c_simplex.append(simplex_counter["secondary"])
                        simplex_counter["secondary"] += 1

                output["simplices"]["secondary"].append(c_simplex)


            # a nakoniec este priradit primarnej zlozke jej cast rozdeleneho trojuholnika
            sd = [0, 1, 2]
            td = [t3d_vertices[0], t3d_vertices[1], t[on_left_idx]]
            delaunay_sim_map = [None, None, s[on_left_idx]]

            output["faces"]["primary"].append(td)
            c_simplex = []
            for sim in sd:
                try:
                    sm = simplex_map["primary"][delaunay_sim_map[sim]]
                    c_simplex.append(sm)
                    continue
                except KeyError:
                    if delaunay_sim_map[sim] is None:
                        output["vertices"]["primary"].append(t3d_vertices[sim])
                    else:
                        output["vertices"]["primary"].append(cgal_simplex[2][delaunay_sim_map[sim]])
                        simplex_map["primary"][delaunay_sim_map[sim]] = simplex_counter["primary"]
                    c_simplex.append(simplex_counter["primary"])
                    simplex_counter["primary"] += 1
            output["simplices"]["primary"].append(c_simplex)

        else:
            # rovnake cachre machre ako pre == 1 ale len pre sekundarnu zlozku
            output["faces"]["secondary"].append(t)

            c_simplex = []
            for sim in s:
                try:
                    # skusi ci uz je v mape taky index z povodneho simplexu a ak nie je tak sa vykona
                    # exception;
                    # v opacnom pripade sa zapise o ktory simplex sa jedna, ak uz vertex je v poli
                    sm = simplex_map["secondary"][sim]
                    c_simplex.append(sm)
                except KeyError:
                    output["vertices"]["secondary"].append(cgal_simplex[2][sim])
                    simplex_map["secondary"][sim] = simplex_counter["secondary"]
                    c_simplex.append(simplex_counter["secondary"])

                    simplex_counter["secondary"] += 1
            output["simplices"]["secondary"].append(c_simplex)

    output["simplices"]["primary"] = np.array(output["simplices"]["primary"])
    output["faces"]["primary"] = np.array(output["faces"]["primary"])
    output["vertices"]["primary"] = np.array(output["vertices"]["primary"])

    output["simplices"]["secondary"] = np.array(output["simplices"]["secondary"])
    output["faces"]["secondary"] = np.array(output["faces"]["secondary"])
    output["vertices"]["secondary"] = np.array(output["vertices"]["secondary"])

    return output
