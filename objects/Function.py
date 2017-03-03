#!/usr/bin/python

""" Function trieda
-------------------
  Obsahuje rozne funkcie potrebne pre pracu

Inicializacne parametre:
------------------------
  Staticka trieda, ziadne inicializacne parametre

Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
------------------------------------------------------------------------------------------------------

  cartesian_to_spheric()
  ======================

    Vstupne parametre:
    ------------------
    vector : [numpy array], [numpy list], defaultna hodnota np.arange(3, dtype=np.float), 3D vektor pre konverziu
                                          z kartezianskych suradnic do sferickych
    degrees : [bool], defaultna hodnota "False", premenna urcena pre testovacie ucely, pri hodnote True zobrazuje
                      azimutalny a polarny uhol v stupnoch namiesto radianov

    Return:
    -------
    vracia vektor sferickych suradnic, [numpy array], navratovy tvar [\r, \phi, \theta] === [vzdialenost, azimutalny
    (rovinny) uhol, polarny (elevacny) uhol] !!! hodnoty uhlov su v radianoch v pripade degrees = False !!!

    Popis:
    ------
    transformuje vektor kartezianksych do vektora sferickych suradnic

    Bug:
    ----
    ziaden znamy bug

  spheric_to_cartesian()
  ======================

    Vstupne parametre:
    ------------------
    vector : [numpy array], [numpy list], defaultna hodnota np.arange(3, dtype=np.float), 3D vektor pre konverziu
                                          zo sferickych do kartezianksych suradnic !!! hodnoty uhlov su v radianoch !!!

    Return:
    -------
    vracia vektor kartezianksych suradnic, [numpy array], navratovy tvar [x, y, z]

    Popis:
    ------
    transformuje vektor sferickych na vektor kartezianksych suradnic

    Bug:
    ----
    1) funkcia je schopna zobrat za polarny uhol zapornu hodnotu, co nasledne vedie nepredpokaldanemu vysledku
    dokaze vziat napr. [1, np.pi / 2, -np.pi / 2], teda v principe snahu o protichodne smery polohoveho vektora

  rand()
  ======

    Vstupne parametre:
    ------------------
    mn : [float], minimalna hodnota nahodneho cisla
    mx : [float], maximalna hodnota nahodneho cisla
    c : [integer], pocet nahodnych cisel, ktore ma funkcia vratint
    seed : [float], seed pre algoritmus nahodneho cisla

    Return:
    -------
    [python list], vrati list nahodnych "c" cisel medzi min a max

    Popis:
    ------
    vrati list nahodnych "c" cisel medzi min a max

    Bug:
    ----
    ziaden znamy bug

  def rotate()
  ============

    Vstupne parametre:
    ------------------
    angle : [float], defaultna hodnota "0", uhol pre zrotovanie, premenna musi byt zadavana v radianoch
    vector : [numpy array], [python list], defaultna hodnota np.arange(3, dtype=np.float), 3D vektor kartezianskych
                                          suradnic pre zrotovanie
    inverse : [bool], defaultna hodnota "False", premenna v pripade nastavenia na True zmeni smer rotacie z matematicky
                      kladneho na matematicky zaporny
    axis : [string], defaultna hodnota "x", urcuje roacnu os, moznosti: "x", "y", "z"

    Return:
    -------
    [numpy array], vrati zrotovany vektor o pozadovany uhol, navratovy tvar [x_rot, y_rot, z_rot]

    Popis:
    ------
    funkcian a rotovanie 3D kartezianskeho vektora v smere alebo preti smeru hodinovych ruciciek o zvoleny uhol

    Bug:
    ----
    ziaden znamy bug

  decimal_non_zero()
  ==================

    Vstupne parametre:
    ------------------
    num : [float], vstupna numericka hodnota, defaultna hodnota "nezadana"

    Return:
    -------
    [integer], vrati pocet cislic pred desatinou ciarkou

    Popis:
    ------
    vrati pocet cislic pred desatinou ciarkou

    Bug:
    ----
    ziaden znamy bug
"""

import numpy as np
import random
import warnings
import math as m
import globe.variables as gv
import inspect
import matplotlib.colors as colors


def cylindrical_to_cartesian(vector=np.arange(3, dtype=np.float)):
    return np.array([vector[0] * np.cos(vector[1]), vector[0] * np.sin(vector[1]), vector[2]])


def cartesian_to_spheric(vector=np.arange(3, dtype=np.float), degrees=False):
    spherical = np.arange(3, dtype=np.float)
    spherical[0] = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)  # vypocet vzdialenosti "r"

    np.seterr(divide='raise', invalid='raise')
    # vypocet a osetrenie vypoctu pre azimutalny uhol
    try:
        spherical[1] = np.arcsin(
            vector[1] / (np.sqrt(vector[0] ** 2 + vector[1] ** 2)))  # vypocet azimutalneho (rovinneho) uhla
    except:
        spherical[1] = 0.0
    # vypocet a osetrenie vypoctu pre polarny uhol
    try:
        spherical[2] = np.arccos(vector[2] / (spherical[0]))  # vypocet polarneho (elevacneho) uhla
    except:
        spherical[2] = 0.0

    np.seterr(divide='print', invalid='print')

    if vector[0] < 0:
        spherical[1] = np.pi - spherical[1]

    if not degrees:
        return spherical
    else:
        return np.array([spherical[0], np.degrees(spherical[1]), np.degrees(spherical[2])])


def spheric_to_cartesian(vector=np.arange(3, dtype=np.float)):
    cartesian = np.arange(3, dtype=np.float)
    cartesian[0] = vector[0] * m.cos(vector[1]) * m.sin(vector[2])
    cartesian[1] = vector[0] * m.sin(vector[1]) * m.sin(vector[2])
    cartesian[2] = vector[0] * m.cos(vector[2])
    return cartesian


def rotate(
        angle=0,
        vector=np.arange(3, dtype=np.float),
        inverse=False,
        axis="x"):
    matrix = np.arange(9, dtype=np.float).reshape((3, 3))
    vector = np.array(vector)
    if axis == "x":
        matrix[0][0], matrix[0][1], matrix[0][2] = 1, 0, 0
        matrix[1][0], matrix[1][1], matrix[1][2] = 0, np.cos(angle), - np.sin(angle)
        matrix[2][0], matrix[2][1], matrix[2][2] = 0, np.sin(angle), np.cos(angle)
        if inverse:
            matrix[1][2], matrix[2][1] = np.sin(angle), - np.sin(angle)
    if axis == "y":
        matrix[0][0], matrix[0][1], matrix[0][2] = np.cos(angle), 0, np.sin(angle)
        matrix[1][0], matrix[1][1], matrix[1][2] = 0, 1, 0
        matrix[2][0], matrix[2][1], matrix[2][2] = - np.sin(angle), 0, np.cos(angle)
        if inverse:
            matrix[2][0], matrix[0][2] = + np.sin(angle), - np.sin(angle)
    if axis == "z":
        matrix[0][0], matrix[0][1], matrix[0][2] = np.cos(angle), - np.sin(angle), 0
        matrix[1][0], matrix[1][1], matrix[1][2] = np.sin(angle), np.cos(angle), 0
        matrix[2][0], matrix[2][1], matrix[2][2] = 0, 0, 1
        if inverse:
            matrix[0][1], matrix[1][0] = + np.sin(angle), - np.sin(angle)
    return np.dot(matrix, vector)


def rand(mn, mx, c, seed):
    random.seed(seed)
    return [mn + random.random() * (mx - mn) for _ in range(c)]


def decimal_non_zero(num):
    return int(abs(np.log10(num))) + 1


def numpy_array_is_empty_verbose(arr, function_name, class_name, var_name, verbose, line):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    try:
        if arr is None:
            if verbose:
                print(color_string(color="error", string="ValueError: ") + "In reference: " + str(
                    class_name) + ", function: " + str(function_name) + ", line: " + line + ". Variable `" + str(
                    var_name) + "` is None.")
            return True

        if arr.size == 0:
            if verbose:
                print(color_string(color="error", string="EmptyVariableError: ") + "In class: " + str(
                    class_name) + ", function: ",
                      str(function_name) + ", line: " + line + ". Variable `" + str(var_name) + "` is empty.")
            return True
        else:
            return False
    except:
        if verbose:
            print(color_string(color="error", string="ValueError: ") + "In class: " + str(
                class_name) + ", function: " + str(function_name) + ", line: " + line + ". Variable `" + str(
                var_name) + "` is invalid.")
        return True


def is_numpy_array(arr):
    if type(arr) == type(np.array([])):
        return True
    return False


def is_numpy_array_empty(arr):
    if arr.size == 0:
        return True
    return False


def is_python_list(lst):
    if type(lst) == type([]):
        return True
    return False


def empty(var, debug=False):
    # if var is numpy arr type
    if type(var) == type(np.array([])):
        if debug:
            print("Variable type is <numpy array>")
        if var.size == 0:
            return True
    # if var is python tuple
    elif type(var) == type(()):
        if debug:
            print("Variable type is <python tuple>")
        if len(var) == 0:
            return True
    # if var is python list type
    elif type(var) == type([]):
        if debug:
            print("Variable type is <python list>")
        if np.array(var).size == 0:
            return True
    # if var is dictionary type
    elif type(var) == type({}):
        if debug:
            print("Variable type is <python dict>")
        if var is {}:
            return True
    elif type(var) == type(True):
        if debug:
            print("Variable type is <bool>")
        if not var:
            return True
    elif var is None:
        if debug:
            print("Variable type is <NoneType>")
        return True
    elif type(var) == type("foo"):
        if debug:
            print("Variable type is <string>")
        if var is "" or var == "0":
            return True
    else:
        try:
            if np.isnan(var):
                if debug:
                    print("Variable type is <numpy.nan>")
                return True
            else:
                if debug:
                    print("Variable type is <number>")
                if var == 0:
                    return True
        except:
            print("Variable type is invalid")
            return True
    return False


def dict_to_list(dict_var=None, verbose=False):
    if empty(var=dict_var, debug=verbose):
        return False
    if type(dict_var) != type({}):
        return False

    return [dict_var[idx] for idx in dict_var]


def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno


def random_color(c_format="hex", size=1):
    if c_format == "hex":
        color = ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(size)]
    else:
        return False
    return color


def gradient_norm(gradient=None, verbose=False):
    if verbose: pass
    if empty(gradient):
        return False
    norm = [np.linalg.norm(np.array(gradient[idx])) for idx in range(0, len(gradient))]
    return np.array(norm)


def color_string(color, string):
    if color == "info":
        return str(gv.COLOR_GREEN + string + gv.COLOR_END)
    if color == "error":
        return str(gv.COLOR_ERROR + string + gv.COLOR_END)
    if color == "warning":
        return str(gv.COLOR_WARNING + string + gv.COLOR_END)
    return 0


def rgb(minimum=0.0, maximum=1.0, value=1.0):
    minimum, maximum = float(minimum), float(maximum)

    ratio = 2 * (value - minimum) / (maximum - minimum)
    b = int(max(0, 255 * (1 - ratio)))
    r = int(max(0, 255 * (ratio - 1)))
    g = 255 - b - r
    return r, g, b


def arr_to_rgb(arr=None, minimum=None, maximum=None):
    if empty(arr):
        print(color_string(color="warning",
                           string="EmptyVariableError: ") + "In reference: Function, function: arr_to_rgb(), line: " + str(
            lineno()) + ". Variable `arr` is empty.")
        return False
    return np.array([rgb(minimum=minimum, maximum=maximum, value=val) for val in arr])


def arr_to_rainbow(arr=None, minimum=None, maximum=None):
    if empty(arr):
        print(color_string(color="warning",
                           string="EmptyVariableError: ") + "In reference: Function, function: arr_to_rainbow(), line: " + str(
            lineno()) + ". Variable `arr` is empty.")
        return False

    return np.array([rainbow(minimum=minimum, maximum=maximum, value=val) for val in arr])


def rainbow(minimum=None, maximum=None, value=None):
    minimum, maximum = float(minimum), float(maximum)
    ratio = (value - minimum) / (maximum - minimum)
    r, g, b = None, None, None
    if 1 >= ratio > .8:
        r = 255
        g = int((1 - ((ratio - .8) / 0.2)) * 255.)
        b = 0
    elif .8 >= ratio > .6:
        r = int(((ratio - .6) / 0.2) * 255.)
        g = 255
        b = 0
    elif .6 >= ratio > .4:
        r = 0
        g = 255
        b = int((1 - ((ratio - .4) / 0.2)) * 255.)
    elif .4 >= ratio > .2:
        r = 0
        g = int(((ratio - .2) / 0.2) * 255.)
        b = 255
    elif .2 >= ratio >= 0:
        r = int((1 - (ratio / 0.2)) * 255.)
        g = 0
        b = 255
    return r, g, b


def rgb_to_hex(color=None, sharp=True):
    if type(color) == type(np.array([])):
        if sharp:
            return np.array([colors.rgb2hex((col[0] / 255.0, col[1] / 255.0, col[2] / 255.0)) for col in color])
        else:
            return np.array(
                [colors.rgb2hex((col[0] / 255.0, col[1] / 255.0, col[2] / 255.0)).replace("#", "") for col in color])
    elif type(color) == type(()):
        return colors.rgb2hex((color[0] / 255.0, color[1] / 255.0, color[2] / 255.0))


def kurucz_model_tablename(temperature, gravity, metallicity):
    metallicity = "m" + str(int(abs(metallicity) * 10)).zfill(2) if metallicity < 0 else  "p" + str(
        int(metallicity * 10)).zfill(2)
    gravity = "g" + str(int(gravity * 10)).zfill(2)
    temperature = str(int(temperature))
    return '"kurucz_' + temperature + '_' + gravity + '_' + metallicity + '"'


def find_nearest_value(array, value):
    array = np.array(array)
    value = array[(np.abs(array - value)).argmin()]
    index = np.where(array == value)[0][0]
    # index = array.tolist().index(value)
    return [value, index]


def find_surounded(array, value):
    # find surounded value in passed array
    arr, ret = np.array(array[:]), []
    f_nst = find_nearest_value(arr, value)
    ret.append(f_nst[0])

    new_arr = []
    if f_nst[0] > value:
        for i in range(0, len(arr)):
            if arr[i] < f_nst[0]:
                new_arr.append(arr[i])
    else:
        for i in range(0, len(arr)):
            if arr[i] > f_nst[0]:
                new_arr.append(arr[i])

    arr = new_arr[:]
    del new_arr

    # arr = np.delete(arr, f_nst[1], 0)
    ret.append(find_nearest_value(arr, value)[0])
    ret = sorted(ret)
    # test
    return ret if ret[0] < value < ret[1] else False


def ramping_necessary():
    return None


def array_mask(array, mask):
    return np.array([array[idx] for idx in mask])


def find_repeated_indices(lst):
    # !!! podhodnotu v poli musia byt listy nie np.array, lebo ta funkcia find_repeated_indices porovnava tvrdo
    # pomocou == a numpy potom pluje, ze pouzi np.all alebo np.any !!!
    lst = lst[:]
    result, repeated = [], []

    if type(lst) == type(np.array([])):
        lst = lst.tolist()

    for value in lst:
        if value == -1:
            continue
        else:
            offset = -1
            repeated = []
            while True:
                try:
                    offset = lst.index(value, offset + 1)
                except ValueError:
                    result.append(np.array(repeated))
                    break

                lst[offset] = -1
                repeated.append(offset)
    return np.array(result)


def find_repeated(array):
    array, repeated = np.array(array)[:], []
    for value in array:
        # ta -1 sa tu dostane tak, ze sa nizsie v else zapise na pozicie, ktore su uz vybavene, aby sa neopakovali
        # hodnoty v poli repeated niekolkrat tie iste
        if value == -1:
            continue
        else:
            indices = np.where(array == value)[0]
            repeated.append(indices)
            array[indices] = -1
    return repeated


def angle_normalisation(angle):
    if angle >= 0:
        return 2.0 * np.pi * m.modf(angle / (2.0 * np.pi))[0]
    elif angle < 0:
        return ((m.modf(angle / (2.0 * np.pi))[0]) * (2.0 * np.pi)) + (2.0 * np.pi)


# matice
def chess_matrix(db_name, band):
    import MySQLdb
    try:
        mysql_db = MySQLdb.connect(host=gv.HOST,  # your host, usually localhost
                                   user=gv.USER,  # your username
                                   passwd=gv.PWD,  # your password
                                   db=db_name)  # name of the data base
        mysql_cur = mysql_db.cursor()

    except:
        print(color_string(color="Error",
                           string="DatabaseError: ") + "In reference: Function, function: chess_matrix(), line: " + str(
            lineno()) + ". Cannot connect to DB.")
        return False

    t_arr = gv.TEMPERATURE_LIST_ATM[:]
    g_arr = gv.GRAVITY_LIST_ATM[:]
    matrix = np.arange(0, (len(g_arr) + 2) * (len(t_arr) + 2), 1).reshape(len(t_arr) + 2,
                                                                          len(g_arr) + 2)

    for i, temp in list(zip(range(0, len(t_arr)), t_arr)):
        # metalicita je irelevantna, pretoze v pripade metalicity je spocitana kazda hodnota;
        # pri niektorych bandoch, mozno niekedy v buducnosti moze chybat niektora hodnota,
        # tak tu tam pre istotu davam, napriek tomu, ze nedavam metalicitu;
        # ak by sa ale nieco zmenilo pre hodnoty v metalicite, tak tu bude treba doplnit parameter
        q = "select distinct(gravity) from ck04_intensity where temperature = " + str(temp) + " and metallicity = " + \
            str(0.0) + " and filter = \"" + str(band) + "\""

        mysql_cur.execute(q)
        logg = np.array(list(zip(*mysql_cur.fetchall()))[0])

        for j, k in list(zip(range(len(g_arr) + 1, -1, -1), range(0, len(g_arr) + 2, 1))):
            # j ide od zadu, teda sa dekrementuje, a pre kazdu hodnotu teploty, sa vyberu existujuce hodnoty
            # gravitacneho zrychlenia v DB;
            # matica ma stale rozmer g_arr (--->) x t_arr (^) a viem, ze je stale dopocitana ta vyssia
            #                                             (|)
            # hodnota pre gravitaciu a teda, ze 5ka npr bude vzdy a teda, ked zavolam logg[j], kde "j" je v dlzke
            # g_arr, ktore ide od zadu, a ak zbehne try. priklad:
            #       i: 29: j: k: 11
            #       logg: [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
            # tak po zavolani logg[1] je mozne citat pole, zbehne try, a k = 11 tak sa 11 - 1 hodnota nastavi na 1
            # a tato 11 - 1 logg = 4.5 hodnota prislucha
            #
            # ak ale 71 2 10 (i, j, k)
            # a mame [ 4.5, 5.0], tak po zavolani 2 je problem, lebo sa neda citat pole a teda hodnota od zadu
            # 10 - 1 musi byt nastaven na 0
            # proste v skratke tie indexy (j, k) musia bezat voci sebe tak, aby toto fungovalo a je to sakra zavisle na
            # tom, ze viem, ze tie horne gravitacne zrychlenia maju vzdy existovat, ak bysa daco pomenilo v modeloch,
            # tak tento kod zdochne jak pes
            try:
                # tu sa skusi zavolat hodnota z DB pola a ked tam je, tak potom tento try zbehne, ak nie,
                # tak to raisne exception a po probleme.
                # nechce sa mi tu totiz davat rozne funkcie na testovanei toho, ci je hodnota isset() co by v php
                # bolo jednoduche, ale tu si to asi isto musim napisat sam a ked sa pomylim, tak to zblbne;
                # tak ako to je teraz to zaberie

                logg[j]
                # toto zapise 1 do matice, cize tam hodnota v DB takejto kombinacie je
                matrix[i + 1][k - 1] = 1
            except:
                # toto zapise 0 do matice, cize tam hodnota v DB takejto kombinacie nie je
                matrix[i + 1][k - 1] = 0
    # vonkajsia obruba matice su 9 a to znamena, ze je to mimo "hernej plochy"
    # tie 9-viatky su tu preto, aby som vedel overit pri hodnotach ako 0,0, ci 0,11 v matici, ci nei som mimo, a teda;
    # samozrejme sa mimo dostanem pomocou sachystickeho algoritmu pre hladanie napadnutych poli figurou, v tomto
    # z pozicie, ktora prislucha najblizsej hodnote pre teplotu a gravitacne zrychlenie na ploske elementarnej
    # pre ktoru chcem interpolovat hodnoty intenzity ziarenia
    for i in range(0, len(g_arr) + 2, 1):
        matrix[0][i] = 9
        matrix[len(t_arr) + 1][i] = 9

    for i in range(0, len(t_arr) + 1, 1):
        matrix[i][0] = 9
        matrix[i][len(g_arr) + 1] = 9

    mysql_db.close()

    return matrix


# kovenrzie suradnic v matici definovanej vyssie
def array2matrix(arr_position):
    # +2 pretoze matica ma obrubu tvorenu 9-viatkami
    gla_length = len(gv.GRAVITY_LIST_ATM) + 2
    return arr_position / gla_length, arr_position % gla_length


def matrix2array(matrix_position):
    gla_length = len(gv.GRAVITY_LIST_ATM) + 2
    return int((matrix_position[0] * gla_length) + (matrix_position[1]))


def pos2matrix(position):
    return position[0] + 1, position[1] + 1


def ck_existence(matrix, t, logg):
    # kontroluje sa len t a logg pretoze pre metallicitu su zratane vsetky hodnoty od -2.5
    # t a logg musia byt polia

    existence = []
    for index, temp, g in list(zip(range(0, len(t)), t, logg)):
        # +1 standardne kvoli obrube matice
        row = find_nearest_value(gv.TEMPERATURE_LIST_ATM, temp)[1] + 1
        column = find_nearest_value(gv.GRAVITY_LIST_ATM, g)[1] + 1
        val = False if matrix[row][column] == 0 else True
        existence.append(val)
    return existence


def replace(what, to, string):
    for w, t in list(zip(what, to)):
        string = string.replace(w, str(t))
    return string


def average_spacing(data=None, neighbours=6):
    if type(data) != type(np.array([])):
        data = np.array(data)

    # ########################## !!!!!!!!!!!!!!!!! IMPORTANT
    # neighbours is variable match same variable in v cgal function
    from scipy.spatial import distance
    dis = distance.cdist(data, data, 'euclidean')
    tot = 0
    for line in dis:
        tot += np.sort(line)[1:1 + neighbours].sum() / (neighbours + 1)
    return tot / dis.shape[0]


def indices(lst, element):
    lst = np.array(lst).tolist()
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


def arbitrary_rotate(theta, omega=None, vector=None, degrees=False):
    # omega - lubovolny vektor okolo ktoreho sa ma rotovat
    # theta - uhol o kolko rotovat
    # vector - vector ktory sa bude rotovat okolo omega

    omega = np.array(omega) / np.linalg.norm(np.array(omega))
    if degrees:
        theta = np.radians(theta)

    matrix = np.arange(9, dtype=np.float).reshape((3, 3))

    matrix[0][0] = (np.cos(theta)) + (omega[0] ** 2 * (1. - np.cos(theta)))
    matrix[0][1] = (omega[0] * omega[1] * (1. - np.cos(theta))) - (omega[2] * np.sin(theta))
    matrix[0][2] = (omega[1] * np.sin(theta)) + (omega[0] * omega[2] * (1. - np.cos(theta)))

    matrix[1][0] = (omega[2] * np.sin(theta)) + (omega[0] * omega[1] * (1. - np.cos(theta)))
    matrix[1][1] = (np.cos(theta)) + (omega[1] ** 2 * (1. - np.cos(theta)))
    matrix[1][2] = (- omega[0] * np.sin(theta)) + (omega[1] * omega[2] * (1. - np.cos(theta)))

    matrix[2][0] = (- omega[1] * np.sin(theta)) + (omega[0] * omega[2] * (1. - np.cos(theta)))
    matrix[2][1] = (omega[0] * np.sin(theta)) + (omega[1] * omega[2] * (1. - np.cos(theta)))
    matrix[2][2] = (np.cos(theta)) + (omega[2] ** 2 * (1. - np.cos(theta)))

    return np.dot(matrix, vector)