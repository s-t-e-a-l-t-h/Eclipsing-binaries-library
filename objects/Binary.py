#!/usr/bin/python

""" Binary trieda
-----------------
 Pomocou triedy Binary sa z dvoch hviezd alebo hviezdy a planety vytvori binarsny system

Inicializacne parametre:
------------------------
  primary : [Object.Star], defaultna hodnota "None", primarna zlozka pre binarny system, objekt triedy Star
  secondary : [Object.Star], defaultna hodnota "None", sekundarna zlozka pre binarny system, objekt triedy Star
  system : [string], defaultna hodnota "None", definuje binarny system ako celok, bude system s planetou alebo s dvomi hviezdami: eb/te
  planete : [string], defaultna hodnota "roche", premenna definujuca tvar planety v pripade systemu s planetou, moznoti: roche/sphere
  verbose : [bool], defaultna hodnota "False", premenna pre zapinanie a vypinanie zobrazovania priebehu vykonavania kodu

Inicilizacia:
-------------
  Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):

  get_info()
  ==========

    Vstupne parametre:
    ------------------
    ziadne

    Return:
    ------
    void

    Popis:
    ------
    momentalne len vytvorena definicia

    Bug:
    ----
    ziaden znamy bug

  primary_potential_fn()
  ======================

    Vstupne parametre:
    ------------------
    radius : [float], defaultna hodnota nie je nastavena, polomer vo sferickych suradniciach od stredu
    *args : [tuple], defaultna hodnota nie je nastavena, argumenty pre vypocet polomeru implicitne zadanej funkcie
                     hodnoty nesuca dana tuple premenna actual_distance, phi, theta = args, pricom
    actual_distance : [float], momentalna separacia zloziek
    phi : [float], azimutalny uhol v radianoch
    theta : [float], polarny uhol v radianoch

    Return:
    -------
    vrati implicitne funkciu zovseobecneneho potencialu pre binarny system

    Popis:
    ------
    funkcia urcena pre implicitny vypocet povrchu primarnej zlozky binarneho systemu

    Bug:
    ----
    ziaden znamy bug

  secondary_potential_fn()
  ========================

    Vstupne parametre:
    pozri funkciu primary_potential_fn()

    Return:
    -------
    vrati implicitne funkciu zovseobecneneho potencialu pre binarny system

    Popis:
    ------
    funkcia urcena pre implicitny vypocet povrchu sekundarnej zlozky binarneho systemu musi byt ina ako pre primarnu,
    pretoze sa musi pouzit symetricka inverzia

    Bug:
    ----
    ziaden znamy bug

  polar_potential_fn()
  ====================

    Vstupne parametre:
    radius : [float], defaultna hodnota nie je nastavena, polomer vo sferickych suradniciach od stredu
    *args : [tuple], defaultna hodnota nie je nastavena, argumenty pre vypocet polomeru implicitne zadanej funkcie,
                     hodnoty nesuce dana tuple premenna su actual_distance, t_object = args, pricom:

                     actual_distance : [float], momentalna separacia zloziek
                     t_object : [string], selektor pre binarnu zlozku "primary/secondary"

    Return:
    -------
    vrati implicitne funkciu zovseobecneneho potencialu pre binarny system v polarnom smere

    Popis:
    ------
    funkcia urcena pre implicitny vypocet polarneho polomeru primarnej alebo sekundarnej zlozky binarneho systemu

    Bug:
    ----
    ziaden znamy bug

  get_polar_radius()
  ==================

    Vstupne parametre:
    ------------------
    t_object : [string], defaultna hodnota "primary", selektor pre zlozku, ktorej polarny polomer chceme spocitat
    actual_distance : [float], defaultna hodnota "None", momentalna separacia zloziek

    Return:
    -------
    vrati polarny polomer vybranej zlozky binarneho systemu

    Popis:
    ------
    funkcia, ktora za vyuzitia polar_potential_fn() spocita polarny polomer vybranej zlozky binarnhe systemu

    Bug:
    ----
    ziaden znamy bug

  get_lagrangian_points()
  =======================

    Vstupne parametre:
    ------------------
    solve_step : [float], defaultna hodnota "0.1", hodnota kroku, po ktorom sa bude alguritmus snazit ziskat riesenie
                          rovnice pre lagrangeove body
    actual_distance : [float], defaultna hodnota "None", momentalna separacia zloziek
    synchronicity_parameter : [float], defaultna hodnota "1.0", parameter synchronicity
                                       (??? otazny je tu tento parameter, kedze pri nezosynchronizovanej orbite tazko
                                       poriesit nejaky spolocny lagrangeove body ???)
    t_object : [string], defaultna hodnota "primary", selektor z pohladu ktorej zlozky v pocaitku sa pocitaju lagrangeove
                       body
                       (??? toto je pouzite presne kvoli parametru synchronicity, kedze pri nesyncrhpnnej orbite su tam
                       dve rozne hodnoty, otazka ja opat ako sa to vlastne pocita, to si bude treba nastudovat ???)

    Return:
    -------
    [numpy array], vrati lagrangeove body

    Popis:
    ------
    spocita a vrati lagrangeove body L1, L2 a L3 pre dany system

    Bug:
    ----
    1) zatial je funkcia sice funkcna, ale neviem ako je to s asynchronnymi orbitami, preto je to tam pripravene dvojako
    a dasa vybrat t_object primary/secondary

  get_3d_model_optmized()
  =======================

    Vstupne parametre:
    ------------------
    phi_steps : [integer], defaultna hodnota "10", pocet krokov na 1/4 hviezdy v rotacnom smere okolo osi x
    theta_steps : [integer], defaultna hodnota "10", pocet krokov na 1/2 hviezdy v rotacnom smere okolo osi y
    t_object : [string], defaultna hodnota "both", definuje zlozku, ktorej povrch sa ma vypocitat, primary/secondary/both
    actual_distance : [float], defaultna hodnota "None", momentalna vzdialenost medzi zlozkami
    critical_angle : [float], defaultna hodnota, "np.pi / 4.", definuje uhol, od ktoreho zacne adaptivny zhusteny
                              vypocet povrchu, rozsah (0, np.pi/2.)
    additional_angle : [float], defautlna hodnota "np.pi / 180.", definuje pridavny uhol od prveho korektne zrataneho
                                bodu na krku typov wuma
    zero_point : [bool], defaultna hodnota "True", definuje, ci sa ma vypocitat dodatocne bod v smere phi = 0.0,
                         theta = np.pi/ 2.0
    homo : [bool], defaultna hodnota "True", definuje, ci sa ma povrch oboch zloziek navzajom udrziavat homogenne
                   pokryty

    Return:
    -------
    [dictionary] v tvare {'primary', 'secondary', 'system'}
    vrati vypocitany povrch systemu ako celku a jednotlivych zloziek

    Popis:
    ------
    spocita povrch binarneho systemu

    Bug:
    ----
    ziaden znamy bug
"""

import numpy as np
import scipy.optimize
import globe.variables as gv
import objects.Function as Fn
import objects.Geometry as Geo
import objects.Orbit as Orb
# import objects.Geometry as Geo
# import sys
import math
import warnings


# import inspect


class Binary:
    def __init__(
            self,
            primary=None,
            secondary=None,
            orbit=None,
            system=None,  # eb/te
            planet="roche",  # roche/sphere
            verbose=False):

        self.verbose = verbose
        self.exception = []

        self.planet = planet  # zadefinovane zatial len preto, aby pycharm nepapuloval
        self.primary = primary
        self.secondary = secondary
        self.system = system
        # premenna init sa nastavuje ako kontrolna premenna, ci zbehla vobec inicializacia vytvarania objektu Binary,
        # ak je init po zavolani na False, tak sa niekde stala chyba a nie je mozne pokracovat dalej v praci, preto
        # je potrebne po vytvoreni objektu skontrolovat vzdy premennu init
        self.init = True
        self.binary_morph = None
        self.relative_semimajor_axis = None

        # nemam cas sa hrat s tym, co z toho pasuje aj pre planetu a co pre hviezdu a ifovat to tam;
        # co sa tyka rychlosti programovania, je momentalne (13.02.2017) jednoduchsie nareplikovat cast kodu
        # pre exoplanety a spravit to cez podmienky eb/te
        if self.system == "eb":
            if self.verbose:
                print(Fn.color_string("info", "Info: ") + "Binary initialisation star.")
            try:
                if not type(primary).__name__ == "Star" or not type(secondary).__name__ == "Star":
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="InitError, TypeError: ") + "In class: Binary, function: __init__(), "
                                                                                 "line: " + str(
                            Fn.lineno()) + ". Variable `primary` or `secondary` is an invalid type.")
                    self.exception.append(
                        "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". `primary` or `secondary` "
                                           "is an invalid type.")
                    self.init = False
                    raise Exception()

                # testing required input parameters
                # variables = vars(primary).keys()
                variables = ["mass", "potential", "albedo", "synchronicity_parameter", "effective_temperature",
                             "gravity_darkening", "metallicity"]
                for var in variables:
                    if primary.__dict__[var] is None or secondary.__dict__[var] is None:
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Variable `" + str(var) + "` is required.")
                        self.exception.append(
                            "In class: Binary, function: __init__(), line: " + str(Fn.lineno()) + ". Variable `" + str(
                                var) + "` is required.")

                        self.init = False
                        raise Exception()

                self.mass_ratio = secondary.mass / primary.mass
                self.invert_mass_ratio = 1.0 / self.mass_ratio

                if orbit is not None:
                    self.orbit = orbit
                else:
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Missing orbit class on initialisation.")
                    self.exception.append("InitError: In class: Binary, function: __init__(), line: " + str(
                        Fn.lineno()) + ". Missing orbit class on initialisation.")
                    self.init = False
                    raise Exception()

                # v kazdom pripade sa nastavi kriticky potencial (ked hviezdy vyplnaju svoj Rocheho lalok)
                # vo vseobecnosti by tieto hodnoty mali byt totozne pri synchronnej orbite
                self.primary.critical_potential = self.critical_potential(t_object="primary",
                                                                          actual_distance=orbit.get_periastron_distance())
                self.secondary.critical_potential = self.critical_potential(t_object="secondary",
                                                                            actual_distance=orbit.get_periastron_distance())

                # osetrit asynchronne orbity

                # filling factor and critical potential
                # ak sa jedna o viazanu rotaciu F1 = F2 = 1.0
                if self.primary.synchronicity_parameter == 1.0 and self.secondary.synchronicity_parameter == 1.0 \
                        and orbit.get_eccentricity() == 0.0:

                    if self.verbose:
                        print(Fn.color_string(color="info", string="Info: ") + "Synchronous rotation")

                    # ceil = 100000.0
                    self.primary.lagrangian_points = self.get_lagrangian_points(actual_distance=1.0,
                                                                                synchronicity_parameter=self.primary.synchronicity_parameter,
                                                                                t_object="primary")
                    # toto nastavenie sekundarnych lagrangeovych bodov je tu len pre uplnost, aby tam neostala hodnota None
                    self.secondary.lagrangian_points = self.primary.lagrangian_points

                    # self.secondary.lagrangian_points = self.get_lagrangian_points(actual_distance = 1.0,
                    # synchronicity_parameter = self.secondary.synchronicity_parameter, t_object = "secondary")


                    argsi, argsii = (1.0, self.primary.lagrangian_points[1], 0., np.pi / 2.), ()
                    # podla pomeru hmotnosti sa za L2 berie ten, ktory je za lahsou zlozkou, cize ak q <= 1, tak sa berie
                    # lagrangian_points[2] a inak sa berie lagrangian_points[0]
                    if 1 >= self.mass_ratio > 0:
                        argsii = (1.0, self.primary.lagrangian_points[2], 0., np.pi / 2.)
                    elif self.mass_ratio > 1:
                        argsii = (1.0, abs(self.primary.lagrangian_points[0]), np.pi, np.pi / 2.)

                    potential_inner = abs(self.potential_value(*argsi))
                    potential_outer = abs(self.potential_value(*argsii))

                    # povodne tam boli zaokruhlene hodnoty a mari sa mi, ze na to bol dovod, ked na neho pridem, tak toto
                    # treba odkomentovat a napisat sem ten dovod potential_inner = math.ceil(abs(self.potential_value(*argsi)) * ceil) / ceil
                    # potential_outer	= math.ceil(abs(self.potential_value(*argsii)) * ceil) / ceil

                    df_potential = potential_inner - potential_outer

                    # premenna nastavena len pre self testing, nema realne vyuzitie ako self.
                    self.potential_inner = potential_inner
                    self.potential_outer = potential_outer
                    self.df_potential = df_potential

                    # nastavenie premennych pre filling faktory primarnej a sekundarnej zlozky
                    self.primary.filling_factor, self.secondary.filling_factor = (
                                                                                     potential_inner - self.primary.potential) / (
                                                                                     df_potential), (
                                                                                     potential_inner - self.secondary.potential) / (
                                                                                     df_potential)

                    # otestovanie over-contact systemu, ci je nastaveny spravny potencial, to znamena, ci su oba potencialy
                    # pri tomto systeme totozne ak nie su, tak sa init nastavi na False
                    # !!! TOTO BUDE MOZNO TREBA ZMENIT TAK, ZE AK SA DETEKUJE OVER-CONTACT, TAK AUTOMATICKY NASTAVIT
                    # HODNOTY POTENCIALU NA HODNOTU PRIMARNEJ ZLOZKY [xx.08.2016]
                    # !!! TO JE LEN DETAIL, UVIDI SA AKO SA TO BUDE SPRAVAT A CO EJ VYHODNEJSIE S ODSTUPOM CASU

                    # musi to byt zaukruhlovane, lebo sa stalo, ze jedna hviezda mala byt presne v laloku a samozrejme
                    # numerika je svian a filling_factor bol 1e-15, a toto sa tu vyhodnotilo jak spravne

                    if ((1 > round(self.secondary.filling_factor, 10) > 0) or (
                                    1 > round(self.primary.filling_factor, 10) > 0)) and (
                                round(self.primary.filling_factor, 10) != round(self.secondary.filling_factor, 10)):
                        self.init = False
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="ValueError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Detected over-contact system, but potentials of components don't match.")
                            # ak by som to chcel niekedy v buducnosti osetrovat, tak to treba spravit tak,
                            # ze sa nastavi pre obe zlozky hodnota tej hviezdy, kde je potencial medzi 0 a 1;
                            # v pripade, ze to je pri oboch, tak asi vybrat jednu, ale osobne teraz [03.11.2016] si myslim,
                            # ze treba raisnut hodit exception a vyriesene napriek tomu, ze vyssie pisem, ze treba nastavit
                            # na primarnu, ale to nemusi byt vzdy overcontact system
                        self.exception.append("ValueError: In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Detected over-contact system, but potentials of components don't match.")
                        raise Exception()

                    if self.primary.filling_factor > 1 or self.secondary.filling_factor > 1:
                        # ak je filling factor pre jednu zo zloziek vacsi ako 1 tak je to nefyzikalny system;
                        # [03.11.2016] - treba si overti teoriu, ale kedze tu mam tuto podmienku, tak som sa na to uz asi
                        # pozeral v minulosti a tak to bolo, minimalne si pamatam, ze som to testoval v kode a empiria
                        # ukazala, ze to tak je
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="ValueError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Non-Physical system.")
                        self.exception.append("ValueError In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Non-Physical system: primary.filling_factor > 1 or secondary.filling_factor > 1")
                        self.init = False
                        raise Exception()

                    if self.primary.filling_factor < 0 and self.secondary.filling_factor < 0:
                        self.binary_morph = "detached"
                        if self.verbose:
                            print(Fn.color_string(color="info", string="Info: ") + "Detached binary system.")
                    elif (round(self.primary.filling_factor, 10) == 0 and self.secondary.filling_factor < 0) or (
                                    self.primary.filling_factor < 0 and round(self.secondary.filling_factor, 10) == 0):
                        self.binary_morph = "semi-detached"
                        if self.verbose:
                            print(Fn.color_string(color="info", string="Info: ") + "Semi-detached binary system.")
                    elif 1 > self.primary.filling_factor > 0:
                        self.binary_morph = "over-contact"
                        if self.verbose:
                            print(Fn.color_string(color="info", string="Info: ") + "Over-contact binary system.")
                    elif self.primary.filling_factor > 1 or self.secondary.filling_factor > 1:
                        if self.verbose:
                            print(Fn.color_string(color="info", string="Info: ") + "Open system.")

                else:
                    self.binary_morph = "WARNING, NOT CIRCULAR ORBIT, ELLIPTICAL IS NOT CURRENTLY SUPPORTED"

                    if verbose:
                        print(Fn.color_string(color="warning",
                                              string="Warning: ") + "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Asynchronous and/or not circular orbit currently not supported.")
                        # /orbitalna trieda
                    self.exception.append("Warning: In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Asynchronous and/or not circular orbit currently not supported.")
                    raise Exception()

                # /filling factor or critical potential

                # ak nie je aspon jeden z polarnych polomerov hviezd spocitany, premenna init sa nastavi na False
                # ak totiz nedoslo k ich zrataniu, pravdepodobne je to sposobene zlymi vstupnymi parametrami, su nefyzikalne
                # preto nema zmysel dalej ratat, za prve preto, lebo sa polomer dalej pouziva, za druhe preto lebo sa nic
                # dalej neporata so zlym systemom

                self.periastron_distance = orbit.get_periastron_distance()
                self.primary.polar_radius = self.get_polar_radius(t_object="primary",
                                                                  actual_distance=self.periastron_distance)
                self.secondary.polar_radius = self.get_polar_radius(t_object="secondary",
                                                                    actual_distance=self.periastron_distance)

                if not self.primary.polar_radius or not self.secondary.polar_radius: self.init = False

                self.primary.backward_radius = self.get_backward_radius(t_object="primary",
                                                                        actual_distance=self.periastron_distance)
                self.secondary.backward_radius = self.get_backward_radius(t_object="secondary",
                                                                          actual_distance=self.periastron_distance)

                # orbitalna trieda
                if not isinstance(self.orbit, Orb.Orbit):
                    if verbose:
                        print(Fn.color_string(color="warning",
                                              string="Warning: ") + "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Orbital class is missing.")
                        # /orbitalna trieda
                    self.exception.append("Warning: In class: Binary, function: __init__(), line: " + str(
                        Fn.lineno()) + ". Orbital class is missing.")
                    self.init = False
                else:
                    # relativna dlzka hlavnej polosi sa nastavi pre instanciu Binary aj pre instanciu Orbit
                    self.relative_semimajor_axis = (((((self.orbit.orbital_period * 86400.0) ** 2) * (
                        gv.G_CONSTANT * gv.SOLAR_MASS * (self.primary.mass + self.secondary.mass))) / (
                                                         4.0 * np.pi ** 2)) ** (
                                                        1.0 / 3.0)) / gv.SOLAR_RADIUS  # in gv.SOLAR_RADIUS unit
                    self.orbit.relative_semimajor_axis = self.relative_semimajor_axis

                    # nastavi sa polarne gravitacne zrychlenie pre primarnu a sekundarnu zlozku v periastre,
                    # cize znova treba zabezpecit, aby pri excentrickej drahe bolo zachovane, ze v periastre d = 1
                    self.primary.polar_gravity = self.compute_polar_gravity(actual_distance=self.periastron_distance,
                                                                            angular_velocity=self.orbit.mean_angular_velocity,
                                                                            t_object="primary")
                    self.secondary.polar_gravity = self.compute_polar_gravity(actual_distance=self.periastron_distance,
                                                                              angular_velocity=self.orbit.mean_angular_velocity,
                                                                              t_object="secondary")

                # nastavenie niektorych polarnych parametrov hviezd (pre periastrum)
                self.primary.polar_gradient_norm = self.compute_polar_gradient_norm(
                    actual_distance=self.periastron_distance,
                    t_object="primary",
                    verbose=self.verbose)
                self.secondary.polar_gradient_norm = self.compute_polar_gradient_norm(
                    actual_distance=self.periastron_distance,
                    t_object="secondary",
                    verbose=self.verbose)

                self.primary.gravity_scalling_factor = self.compute_gravity_scalling_factor(t_object="primary")
                self.secondary.gravity_scalling_factor = self.compute_gravity_scalling_factor(t_object="secondary")

            except:
                self.init = False
        ########################################################################################################
        ###################                    TU ZACIANJU PLANETY
        ########################################################################################################
        elif self.system == "te":
            if self.verbose:
                print(Fn.color_string("info", "Info: ") + "Binary initialisation planetary system.")
            try:
                if not type(primary).__name__ == "Star" and not type(secondary).__name__ == "Planete":
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="InitError, TypeError: ") + "In class: Binary, function: __init__(), "
                                                                                 "line: " + str(
                            Fn.lineno()) + ". Variable `primary` or `secondary` is an invalid type.")
                    self.exception.append(
                        "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". `primary` or `secondary` "
                                           "is an invalid type.")
                    self.init = False
                    raise Exception()

                # testing required input parameters
                variables = ["mass", "potential", "effective_temperature", "gravity_darkening", "metallicity",
                             "synchronicity_parameter"] if self.planet == "roche" \
                    else ["mass", "potential", "effective_temperature", "gravity_darkening", "metallicity",
                          "angular_velocity"]

                for var in variables:
                    if not None != primary.__dict__[var]:
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Variable `" + str(var) + "` is required.")
                        self.exception.append(
                            "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Variable `" + str(
                                var) + "` is required.")

                        self.init = False
                        raise Exception()

                variables = ["mass", "potential", "albedo", "effective_temperature",
                             "synchronicity_parameter"] if self.planet == "roche" else ["mass", "potential", "albedo",
                                                                                        "effective_temperature"]
                for var in variables:
                    if not None != secondary.__dict__[var]:
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Variable `" + str(var) + "` is required.")
                        self.exception.append(
                            "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Variable `" + str(
                                var) + "` is required.")

                        self.init = False
                        raise Exception()

                if orbit is not None:
                    self.orbit = orbit
                else:
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Missing orbit class on initialisation.")
                    self.exception.append("InitError: In class: Binary, function: __init__(), line: " + str(
                        Fn.lineno()) + ". Missing orbit class on initialisation.")
                    self.init = False
                    raise Exception()

                # pre oba pripady (spehere aj roche), tak to budem davat sem
                self.periastron_distance = orbit.get_periastron_distance()

                # orbitalna trieda
                if not isinstance(self.orbit, Orb.Orbit):
                    if verbose:
                        print(Fn.color_string(color="warning",
                                              string="Warning: ") + "In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Orbital class is missing.")
                        # /orbitalna trieda
                    self.exception.append("Warning: In class: Binary, function: __init__(), line: " + str(
                        Fn.lineno()) + ". Orbital class is missing.")
                    self.init = False
                    raise Exception()

                # relativna dlzka hlavnej polosi sa nastavi pre instanciu Binary aj pre instanciu Orbit
                # spocita sa z tretieho Keplerovho zakona
                self.relative_semimajor_axis = (((((self.orbit.orbital_period * 86400.0) ** 2) * (
                    gv.G_CONSTANT * gv.SOLAR_MASS * (self.primary.mass + self.secondary.mass))) / (
                                                     4.0 * np.pi ** 2)) ** (
                                                    1.0 / 3.0)) / gv.SOLAR_RADIUS  # in gv.SOLAR_RADIUS unit
                self.orbit.relative_semimajor_axis = self.relative_semimajor_axis

                # rozdelit na dva pripady: co sa bude pouzivat sfereicke priblizenie alebo klasicke dvojhviezdne
                if self.planet == "roche":
                    self.mass_ratio = secondary.mass / primary.mass
                    self.invert_mass_ratio = 1.0 / self.mass_ratio

                    self.primary.critical_potential = self.critical_potential(t_object="primary",
                                                                              actual_distance=orbit.get_periastron_distance())
                    self.secondary.critical_potential = self.critical_potential(t_object="secondary",
                                                                                actual_distance=orbit.get_periastron_distance())

                    # ak je mensia hodnota potencialu ako je kriticka, tak je to pre exoplanety sprostost, nemoze byt
                    # wuma system
                    if self.primary.potential <= self.primary.critical_potential or \
                                    self.secondary.potential <= self.secondary.critical_potential:
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Invalid potential value.")
                        self.exception.append("InitError: In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Invalid potential value.")

                        self.init = False
                        raise Exception()

                    self.primary.polar_radius = self.get_polar_radius(t_object="primary",
                                                                      actual_distance=self.periastron_distance)
                    self.secondary.polar_radius = self.get_polar_radius(t_object="secondary",
                                                                        actual_distance=self.periastron_distance)

                    if not self.primary.polar_radius or not self.secondary.polar_radius:
                        self.init = False

                    self.primary.backward_radius = self.get_backward_radius(t_object="primary",
                                                                            actual_distance=self.periastron_distance)
                    self.secondary.backward_radius = self.get_backward_radius(t_object="secondary",
                                                                              actual_distance=self.periastron_distance)

                    # nastavi sa polarne gravitacne zrychlenie pre primarnu zlozku v periastre,
                    # periastralna vzialenost v jednotkach dlzky hlavnej poloosi a je v orbit.periastron_distance
                    self.primary.polar_gravity = self.compute_polar_gravity(
                        actual_distance=self.periastron_distance,
                        angular_velocity=self.orbit.mean_angular_velocity,
                        t_object="primary")

                    # nastavenie niektorych polarnych parametrov hviezd (pre periastrum)
                    self.primary.polar_gradient_norm = self.compute_polar_gradient_norm(
                        actual_distance=self.periastron_distance,
                        t_object="primary",
                        verbose=self.verbose)

                    self.primary.gravity_scalling_factor = self.compute_gravity_scalling_factor(t_object="primary")


                elif self.planet == "sphere":
                    self.primary.polar_radius = 1.0 / self.primary.potential
                    self.primary.backward_radius = self.primary.single_start_backward_radius()

                    self.secondary.polar_radius = 1.0 / self.secondary.potential
                    self.secondary.backward_radius = self.secondary.polar_radius

                    if not self.primary.polar_radius or not self.secondary.polar_radius or \
                            not self.primary.backward_radius or not self.secondary.backward_radius or \
                                    round(self.primary.backward_radius - self.primary.polar_radius, 10) < 0:

                        self.init = False

                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                                Fn.lineno()) + ". Invalid radius value.")
                        self.exception.append("InitError: In class: Binary, function: __init__(), line: " + str(
                            Fn.lineno()) + ". Invalid radius value.")
                        raise Exception()

                    primary.polar_gravity = \
                        self.primary.single_star_polar_gravity(
                            polar_radius=self.primary.polar_radius * self.orbit.relative_semimajor_axis * gv.SOLAR_RADIUS)

                    # nastavenie niektorych polarnych parametrov hviezd (pre periastrum)
                    self.primary.polar_gradient_norm = self.primary.single_star_polar_gradient_norm()

                    # tu mozno pouzit aj funkciu z tejto triedy (Binary), pretoze je to len delenie dvoch polarnych
                    # hodnot
                    self.primary.gravity_scalling_factor = self.compute_gravity_scalling_factor(t_object="primary")


            except ValueError:
                self.init = False

        if self.init and self.verbose:
            print(Fn.color_string(color="info", string="Info: ") + "Binary initialisation success.")

        if not self.init and self.verbose:
            print(Fn.color_string(color="error",
                                  string="InitError: ") + "In class: Binary, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation.")
            self.exception.append("InitError: In class: Binary, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation.")

            print(self.get_exception())

    def get_info(self):
        print('\n' + str(self.__class__.__name__) + ' -------------------------------------------------')
        print("primary class:\t\t\t\t\t\t\t" + str(self.primary))
        print("secondary class:\t\t\t\t\t\t" + str(self.secondary))
        print("orbit class:\t\t\t\t\t\t\t" + str(self.orbit))
        print("system:\t\t\t\t\t\t\t\t\t" + str(self.system))
        print("mass ratio:\t\t\t\t\t\t\t\t" + str(self.mass_ratio))
        print("inver mass ratio:\t\t\t\t\t\t" + str(self.invert_mass_ratio))
        print("eccentricity:\t\t\t\t\t\t\t" + str(self.orbit.get_eccentricity()))
        print("binary morphology:\t\t\t\t\t\t" + str(self.binary_morph))
        print("relative semimajor axis:\t\t\t\t" + str(self.relative_semimajor_axis))
        print("periastron distance:\t\t\t\t\t" + str(self.periastron_distance))
        print("init:\t\t\t\t\t\t\t\t\t" + str(self.init))
        print('/', self.__class__.__name__, ' -----------------------------------------------')
        return None

    @classmethod
    def general_cartesian_potential_fn(cls, x, y, z, actual_distance, synchronicity_parameter, mass_ratio):
        block_a = 1.0 / np.sqrt(x ** 2 + y ** 2 + z ** 2)
        block_b = mass_ratio / np.sqrt((actual_distance - x) ** 2 + y ** 2 + z ** 2)
        block_c = -0.5 * mass_ratio ** 2 / (actual_distance ** 4 * (mass_ratio + 1.0) * synchronicity_parameter ** 2)
        block_d = 0.5 * (mass_ratio + 1.0) * ((((mass_ratio / (actual_distance ** 2 * (
            mass_ratio + 1.0) * synchronicity_parameter)) - synchronicity_parameter * x) ** 2) + y ** 2 * synchronicity_parameter ** 2)
        return block_a + block_b + block_c + block_d

    @classmethod
    def general_spheric_potential_fn(cls, r, phi, theta, actual_distance, synchronicity_parameter, mass_ratio):
        lambd = np.cos(phi) * np.sin(theta)
        nu = np.cos(theta)
        block_a = 1.0 / r
        block_b = mass_ratio * ((1.0 / np.sqrt(r ** 2 - (2.0 * r * lambd * actual_distance) + actual_distance ** 2)) - (
            r * lambd / actual_distance ** 2))
        block_c = 0.5 * synchronicity_parameter ** 2 * (1.0 + mass_ratio) * (1.0 - nu ** 2) * r ** 2
        return block_a + block_b + block_c

    # def wuma_cylindric_potential_fn(self, radius, *args):
    #     x, phi = args
    #     block_a = 1.0 / np.sqrt(x ** 2 + radius ** 2)
    #     block_b = self.mass_ratio / np.sqrt((1.0 - x) ** 2 + radius ** 2)
    #     block_c = -0.5 * (self.mass_ratio ** 2) / (self.mass_ratio + 1.0)
    #     block_d = 0.5 * (self.mass_ratio + 1.0) * (
    #         (self.mass_ratio / (self.mass_ratio + 1.0)) ** 2 + radius ** 2 * (np.cos(phi)) ** 2)
    #     return block_a + block_b + block_c + block_d - self.primary.potential

    def wc_potential_fn(self, radius_p, *args):
        actual_distance, phi, x = args
        F = 1.0
        potential_omega = self.primary.potential

        block_a = (1.0 / np.sqrt(radius_p ** 2 + x ** 2))
        block_b = self.mass_ratio * ((1.0 / (np.sqrt((x - actual_distance) ** 2 + radius_p ** 2))) -
                                     (x / actual_distance ** 2))
        block_c = 0.5 * (1.0 + self.mass_ratio) * (F ** 2) * (x ** 2 + (radius_p ** 2 * np.cos(phi) ** 2))
        return block_a + block_b + block_c - potential_omega

    def primary_potential_fn(
            self,
            radius,
            *args):
        # variables
        actual_distance, phi, theta = args
        # /variables

        # block of function
        block_a = (1 / radius)
        block_b = (self.mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2. * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.primary.synchronicity_parameter, 2) * (1 + self.mass_ratio) * np.power(radius, 2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        return (block_a + block_b - block_c + block_d) - self.primary.potential

    def secondary_potential_fn(
            self,
            radius,
            *args):
        # variables
        actual_distance, phi, theta = args
        # /variables

        # block of function
        block_a = (1. / radius)
        block_b = (self.invert_mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2 * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.invert_mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.secondary.synchronicity_parameter, 2) * (1 + self.invert_mass_ratio) * np.power(
                radius,
                2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        inverse_potential = (block_a + block_b - block_c + block_d) / self.invert_mass_ratio + (
            0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))
        return inverse_potential - self.secondary.potential

    def polar_potential_fn(
            self,
            radius,
            *args):
        # variables
        actual_distance, t_object = args
        # /variabreturn - (block_a + block_b - block_c + block_d)les
        if t_object == "primary":
            return (1. / radius) + (self.mass_ratio * (
                (1 / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2)))))) - self.primary.potential
        if t_object == "secondary":
            return (((1. / radius) + (
                self.invert_mass_ratio * (1 / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2)))))) / (
                        self.invert_mass_ratio) + (
                        0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))) - self.secondary.potential

    def get_polar_radius(
            self,
            t_object="primary",
            actual_distance=None):
        # premenna scipy_solve_point sluzi ako startovaci bod pre numercky vypocet implicitnej rovnice pre polarny polomer
        # prvykrat pri vyvoji bola nastavena na nieco okolo 0.01 no pri strasne malych zlozkach, ktorych polomer je pdo touto hodnotou
        # numerika neskonvergovala, tak musel byt nastavena na hodnotu vyrazne mensiu, vsetky dalsie vypocty po polarnom polomere (uz ked je znamy)
        # budu musiet byt nastavene na nejaku 1/10 hodnoty polomeru pre danu zlozku systemu
        ret, scipy_solve_point, solution = False, 1e-20, "NaN"
        np.seterr(divide='raise', invalid='raise')
        try:
            args = (actual_distance, t_object)
            solution, info, ier, msg = scipy.optimize.fsolve(self.polar_potential_fn, scipy_solve_point,
                                                             full_output=True, args=args)
            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                # tu sa kontroluje, ci je spocitany polomer vacsi ako 0 a mesni ako 1
                # urcite nemoze byt mensi ako 0, a nepredpoklada sa, ze presiahen hodnotu separacie 1 pri zlozkach, takyto system ais nemoze nijak vzniknut
                # ak je hodnota ina, pojde urcite o chybu numeriky
                if 1 > solution > 0:
                    ret = True
        except:
            ret = False

        np.seterr(divide='print', invalid='print')
        if ret:
            return solution
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Binary, function: get_polar_radius(), line: " + str(
                    Fn.lineno()) + ", star object: " + gv.COLOR_BLUE + str(
                    t_object) + gv.COLOR_END + ". Error has been occurred while polar radius counting, primary mass (" + self.primary.mass + ") secondary mass (" + self.secondary.mass + ").")
            return False

    def get_backward_radius(
            self,
            t_object="primary",
            actual_distance=None
    ):

        args, ret, ier, solution = (actual_distance, np.pi, np.pi / 2.0), False, 10, 0.0
        try:
            if t_object == "primary":
                scipy_solve_point = self.primary.polar_radius / 10.0
                solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                 scipy_solve_point,
                                                                 full_output=True,
                                                                 args=args)

            elif t_object == "secondary":
                scipy_solve_point = self.secondary.polar_radius / 10.0
                solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                 scipy_solve_point,
                                                                 full_output=True, args=args, )

            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                # tu sa kontroluje, ci je spocitany polomer vacsi ako 0 a mesni ako 1
                # urcite nemoze byt mensi ako 0, a nepredpoklada sa, ze presiahen hodnotu separacie 1 pri zlozkach,
                # takyto system ais nemoze nijak vzniknut; ak je hodnota ina, pojde urcite o chybu numeriky
                if 1 > solution > 0:
                    ret = True
        except:
            ret = False

        if ret:
            return solution
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Binary, function: get_backward_radius(), line: " + str(
                    Fn.lineno()) + ", star object: " + gv.COLOR_BLUE + str(
                    t_object) + gv.COLOR_END + ". Error has been occurred while backward radius counting, primary mass (" + str(
                    self.primary.mass) + ") secondary mass (" + str(self.secondary.mass) + ").")
            return False

    def get_lagrangian_points(
            self,
            solve_step=0.1,
            actual_distance=None,
            synchronicity_parameter=1.0,
            t_object='primary'):

        # ta funkcia funguje tak, ze zacne ratat numericky derivaciu potencialu po x - ovej osi
        # problem je v tom, ze to z nejakej polohy zrata L2, z inej L1 a z inej L3 v zmysle ineho startovacieho bodu
        # z toho dovodu musim prebehnut po osi s nejakym krokom a nechat vyriesit numeriku v kazdom bode
        # dalsi problem je v tom, ze kedze sa jedna o numeriku, tie doratane hodnoty nie su uplne presne
        # teda sa nedaju jednoducho porovnat stylom L_(n) == L_(n-1), pretoze by sa nemuseli rovnat aj keby to boli
        # tie iste lagrangeove body, preto to treba nejak rozumen zaokruhil a tak to porovnavat a "presne" hodnoty
        # ukladat inde
        args, end = (actual_distance, synchronicity_parameter), int(((actual_distance * 6.) / solve_step))
        points, lagrange = [], []
        scipy_solve_point, decimal = - (actual_distance * 3.) + solve_step, int(len(str(solve_step).split('.')[1]))
        derivation_value, ier = None, np.inf

        for i in range(0, end):
            # toto sa vykonava aby sa zistilo, ci nedochadza k deleniu nulou pri danej polohe riesenia rovnice, ak ano,
            # tak sa to proste preskoci
            try:
                np.seterr(divide='raise', invalid='raise')
                if t_object == 'primary': self.primary_potential_derivation_x(round(scipy_solve_point, decimal), *args)
                if t_object == 'secondary': self.secondary_potential_derivation_x(round(scipy_solve_point, decimal),
                                                                                  *args)
                np.seterr(divide='print', invalid='print')
                pass
            except:
                scipy_solve_point += solve_step
                continue

            try:
                if t_object == 'primary':
                    solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_derivation_x,
                                                                     scipy_solve_point, full_output=True, args=args)
                if t_object == 'secondary':
                    solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_derivation_x,
                                                                     scipy_solve_point, full_output=True, args=args)

                if ier == 1:
                    if round(solution[0], 5) not in points:
                        try:
                            if t_object == 'primary': derivation_value = abs(
                                round(self.primary_potential_derivation_x(solution[0], *args), 4))
                            if t_object == 'secondary': derivation_value = abs(
                                round(self.secondary_potential_derivation_x(solution[0], *args), 4))
                            if derivation_value == 0:
                                use = True
                            else:
                                use = False
                        except:
                            use = False
                        if use:
                            points.append(round(solution[0], 5))
                            lagrange.append(solution[0])
            except:
                scipy_solve_point += solve_step
                continue

            scipy_solve_point += solve_step
        return np.array(sorted(lagrange))

    def compute_equipotential_xy(
            self,
            step=np.pi / 50,
            actual_distance=None):
        end, phi, points, ier = int((2. * np.pi) / step), 0., [], np.inf
        for obj in range(0, 2):
            for i in range(0, end):
                args, use = (actual_distance, phi, np.pi / 2.), False
                try:
                    if obj == 0:
                        scipy_solve_point = self.primary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    elif obj == 1:
                        scipy_solve_point = self.secondary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    if ier == 1 and not np.isnan(solution[0]):
                        solution = solution[0]
                        if 30 >= solution >= -30: use = True
                except:
                    use = False
                if use:
                    if obj == 0:
                        points.append([solution * np.cos(phi), solution * np.sin(phi)])
                    elif obj == 1:
                        points.append([- (solution * np.cos(phi) - actual_distance), solution * np.sin(phi)])
                phi += step
        return np.array(points)

    def compute_hill_plane_prototype(
            self,
            steps=100,
            actual_distance=None,
            fsolve_radius=3,
            fsolve_steps=5
    ):
        ###############################################################################################################
        # tato funkcia miestami kresli pekne sprostosti, ale je tu pre uplnost, pretoze niekedy mozno v niecom pomoze #
        ###############################################################################################################

        np.seterr(divide='ignore', invalid='ignore')

        phi, phi_step_length, o, ier = 0.0, 2.0 * np.pi / steps, ['primary', 'secondary'], 10
        points, duplicit_test_points = [], []
        fsolve_step_length = fsolve_radius / float(fsolve_steps)
        scipy_solve_point = min([self.primary.polar_radius, self.secondary.polar_radius])
        first_scipy_solve_point = scipy_solve_point

        for t_object in o:
            phi = 0.0
            while phi < 2.0 * np.pi:
                args, use = (actual_distance, phi, np.pi / 2.), False
                scipy_solve_point = first_scipy_solve_point
                while scipy_solve_point < fsolve_radius:
                    # print scipy_solve_point, phi
                    try:
                        if t_object == 'primary':
                            solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                             scipy_solve_point,
                                                                             full_output=True, args=args)
                        elif t_object == 'secondary':
                            solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                             scipy_solve_point,
                                                                             full_output=True, args=args)
                        if ier == 1 and not np.isnan(solution[0]):
                            solution = solution[0]
                            if 30 >= solution >= -30: use = True
                    except:
                        use = False
                    if use:
                        if t_object == 'primary':
                            result = [solution * np.cos(phi), solution * np.sin(phi)]
                            if not np.around(result, 5).tolist in duplicit_test_points:
                                points.append(result)
                                duplicit_test_points.append(np.around(result, 5).tolist)
                        elif t_object == 'secondary':
                            result = [- (solution * np.cos(phi) - actual_distance), solution * np.sin(phi)]
                            if not np.around(result, 5).tolist in duplicit_test_points:
                                points.append(result)
                                duplicit_test_points.append(np.around(result, 5).tolist)
                    scipy_solve_point += fsolve_step_length

                phi += phi_step_length

        np.seterr(divide='print', invalid='print')
        return points

    def get_potential_value(
            self,
            actual_distance=None,
            step=0.05,
            interval=None):
        if interval is None:
            interval = [-3, 3]
        np.seterr(divide='ignore', invalid='ignore')
        # variable init
        phi, potential = None, None

        points, end, runner, radius = [], int((interval[1] - interval[0]) / step), interval[0] + step, 0.0

        for i in range(0, end):
            try:
                if runner < 0:
                    radius, phi = abs(runner), np.pi
                else:
                    radius, phi = runner, 0.0
                args = (actual_distance, radius, phi, np.pi / 2.)
                potential = round(self.potential_value(*args), 4)
                if not np.isnan(potential):
                    if 10 >= potential >= -10:
                        use = True
                    else:
                        use = False
                else:
                    use = False
            except:
                use = False
            if use: points.append([radius * np.cos(phi), potential])
            runner += step
        np.seterr(divide='print', invalid='print')
        return points

    def primary_potential_derivation_x(
            self,
            x,
            *args):
        actual_distance, synchronicity_parameter = args
        r_sqr, rw_sqr = x ** 2, (actual_distance - x) ** 2
        return - (x / r_sqr ** (3. / 2.)) + (
            (self.mass_ratio * (actual_distance - x)) / rw_sqr ** (3. / 2.)) + synchronicity_parameter ** 2 * (
            self.mass_ratio + 1) * x - self.mass_ratio / actual_distance ** 2

    def secondary_potential_derivation_x(
            self,
            x,
            *args):
        actual_distance, synchronicity_parameter = args
        r_sqr, rw_sqr = x ** 2, (actual_distance - x) ** 2
        return - (x / r_sqr ** (3. / 2.)) + (
            (self.mass_ratio * (actual_distance - x)) / rw_sqr ** (3. / 2.)) - synchronicity_parameter ** 2 * (
            self.mass_ratio + 1) * (1 - x) + (1. / actual_distance ** 2)

    def potential_value(
            self,
            *args):
        # variables
        actual_distance, radius, phi, theta = args
        # /variables

        # block of function
        block_a = (1. / radius)
        block_b = (self.mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2. * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.primary.synchronicity_parameter, 2) * (1 + self.mass_ratio) * np.power(radius, 2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        return - (block_a + block_b - block_c + block_d)

    def inverted_potential_value(
            self,
            *args):
        # variables
        actual_distance, radius, phi, theta = args
        # /variables

        # block of function
        block_a = (1. / radius)
        block_b = (self.invert_mass_ratio / (np.sqrt(np.power(actual_distance, 2) + np.power(radius, 2) - (
            2 * radius * (np.cos(phi) * np.sin(theta)) * actual_distance))))
        block_c = ((self.invert_mass_ratio * radius * (np.cos(phi) * np.sin(theta))) / (np.power(actual_distance, 2)))
        block_d = (
            0.5 * np.power(self.secondary.synchronicity_parameter, 2) * (1 + self.invert_mass_ratio) * np.power(
                radius,
                2) * (
                1 - np.power(np.cos(theta), 2)))
        # /block of function
        inverse_potential = (block_a + block_b - block_c + block_d) / self.invert_mass_ratio + (
            0.5 * ((self.invert_mass_ratio - 1) / self.invert_mass_ratio))
        # /block of function
        return - inverse_potential

    def critical_potential(
            self,
            t_object="primary",
            actual_distance=None):
        # variable init
        solution = None

        if t_object == "primary":
            args = (actual_distance, self.primary.synchronicity_parameter)
            solution = scipy.optimize.newton(self.primary_potential_derivation_x, 0.001, args=args)
        if t_object == "secondary":
            args = (actual_distance, self.secondary.synchronicity_parameter)
            solution = scipy.optimize.newton(self.secondary_potential_derivation_x, 0.001, args=args)
        if not np.isnan(solution):
            if t_object == "primary":
                args = (actual_distance, solution, 0.0, np.pi / 2.)
                return abs(self.potential_value(*args))
            if t_object == "secondary":
                args = (actual_distance, actual_distance - solution, 0.0, np.pi / 2.)
                return abs(self.inverted_potential_value(*args))
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Binary, function: critical_potential(), line: " + str(
                    Fn.lineno()) + ". Wrong value has been encoutered.")
            return False

    def wuma_split(self, faces=None, verbose=False):
        # funkcia rozdeli wuma triangulaciu na lagrangeovom bode L1
        # vsetky fazety s taziskami (x-suradnicami) mensimi ako je L1 budu patrit primarnej zlozke a vsetky vacsie
        # sekundarnej

        if Fn.empty(faces):
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="EmptyVariableError: ") + "In class: Binary, function: wuma_split(), line: " + str(
                    Fn.lineno()) + ". Variable `faces` is empty.")
            return False

        com, lag = Geo.center_of_mass(faces=np.array(faces), verbose=verbose), self.primary.lagrangian_points[1]
        tri_primary, tri_secondary = [], []

        for idx in range(0, len(com)):
            if com[idx][0] < lag:
                tri_primary.append(faces[idx])
            else:
                tri_secondary.append(faces[idx])

        return [np.array(tri_primary), np.array(tri_secondary)]

    def compute_polar_gravity(
            self,
            actual_distance=None,
            angular_velocity=None,
            t_object="primary"
    ):
        mass_ratio = self.mass_ratio if t_object == "primary" else self.invert_mass_ratio
        polar_radius = self.primary.polar_radius if t_object == "primary" else self.secondary.polar_radius
        x_com = (mass_ratio * actual_distance) / (1.0 + mass_ratio)

        primary_mass, secondary_mass = self.primary.mass, self.secondary.mass
        if t_object == "secondary":
            primary_mass, secondary_mass = secondary_mass, primary_mass

        # give length of semi major axis in physical units
        semimajor_axis = gv.SOLAR_RADIUS * self.relative_semimajor_axis
        r_vector = np.array([0.0, 0.0, polar_radius * semimajor_axis])
        centrifugal_distance = np.array([x_com * semimajor_axis, 0.0, 0.0])
        actual_distance = np.array([actual_distance * semimajor_axis, 0., 0.])
        h_vector = r_vector - actual_distance

        block_a = - ((gv.G_CONSTANT * primary_mass * gv.SOLAR_MASS) / np.linalg.norm(r_vector) ** 3) * r_vector
        block_b = - ((gv.G_CONSTANT * secondary_mass * gv.SOLAR_MASS) / np.linalg.norm(h_vector) ** 3) * h_vector
        block_c = - (angular_velocity ** 2) * centrifugal_distance

        g = block_a + block_b + block_c

        return np.linalg.norm(g) * 1e2  # return magnitude of polar gravity acceleration in physical CGS units

    def compute_polar_gradient_norm(
            self,
            t_object='primary',
            actual_distance=None,
            verbose=False
    ):

        polar_vector = np.array([0.0, 0.0, self.primary.polar_radius]) if t_object == "primary" else np.array(
            [actual_distance, 0.0, self.secondary.polar_radius])

        normal = Geo.normal_estimation(binary_object=self, actual_distance=1.0, vertices=np.array([polar_vector]),
                                       t_object=t_object, mode="in_point", verbose=verbose)

        return np.linalg.norm(x=normal[0])

    def compute_gravity_scalling_factor(self, t_object="primary"):
        ret_val = self.primary.polar_gravity / self.primary.polar_gradient_norm if t_object == "primary" else \
            self.secondary.polar_gravity / self.secondary.polar_gradient_norm
        return ret_val

    def get_3d_model_optimized(
            self,
            phi_steps=10,
            theta_steps=10,
            t_object="both",
            actual_distance=None,
            critical_angle=np.pi / 4.,
            additional_angle=np.pi / 180.,
            zero_point=True,
            homo=True):

        if self.verbose:
            print(Fn.color_string(color="info", string="Info: ") + "Computing surface of objects.")

        # inicializacie
        primary_phi_steps, primary_theta_steps = None, None
        secondary_phi_steps, secondary_theta_steps = None, None
        use = False
        z_point, rotation_angle, transform = None, None, None

        # phi_steps, pocet azimutalnych krokov, jedna sa o mierne zavadzajuce pomenovanie, pretoze sa jedna o pocet krokov
        # o kolko sa bude rotovat po kruznici okolo x-ovej osi
        phi_steps = phi_steps
        model, partial, total, point_coordinate = {'primary': 'empty', 'secondary': 'empty',
                                                   'system': 'empty', "separation": "empty"}, \
                                                  [], [], np.arange(3, dtype=np.float)
        ppr, spr = self.primary.polar_radius, self.secondary.polar_radius

        if self.verbose:
            if actual_distance is None:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Binary, function: get_3d_model_testing(), line: " + str(
                    Fn.lineno()) + ". Variable `actual_distance` is set to None.")
                return False

        # vypocet bodov na povrchu hviezdy prebieha z polohy phi = np.pi, theta = np.pi/2., to znamena z jej zadnej casti
        # vypocet po zmene polarneho kroku po kruznici okolo x-ovej osi, ktora sa ziska transformaciou
        # critical_angle sluzi pre zhustenie vypoctu na hrdle wuma systemov, ked polarny uhol narazi na hodnotu np.pi + critical_angle
        # (je to np.pi + critical_angle pretoze ak vypocet zacina na np.pi/2 a pridava sa polarny uhol, tak np.pi dosiahne na hranici
        # tretieho a stvrteho kvadrantu z rovine z-x a dalej teda presahuje np.pi, co nie je problem, pretoze cos je parna funkcia
        # a teada nadobuda rovnaku hodnotu ako pre np.pi + nejaky_uhol ako pre np.pi - nejaky uhol)
        # tu sa skontroluje, ci kriticky uhol nie je vacsi ako pi / 2
        # kriticky uhol sa pocita od 180 stupnov polarneho uhla, vid obrazok
        # oA               |\              Bo
        #  o               | \             o
        #    o             |  \          o
        #      o           | K \       o
        #          o       |    \  o
        #               o  |  o
        # K oznacuje velkost kritickeho uhla
        # ak je K = np.pi/2.0, tak sa pocita klasicky polarny cyklus z polohy A az do polohy B
        if critical_angle > np.pi / 2.: critical_angle = np.pi / 2.

        if t_object == "both":
            j_range_bottom, j_range_top = 0, 2

            if homo:  # homogenny model, snahavytvorit model tak, aby sekundarna aj primarna zlozka mali homogenne plosne elementy
                # ak je premenna homo nastavena na True, tak sa pocet krokov v polarnom a azimutalnom smere nastavi nasledovne
                # 1) ak je polarny polomer primarnej zlozky mensi ako sekundarny, tak na tejto mensje zlozke (primarnej) bude
                #    zadany pocet polarnych a azimutalnych krokov zo vstupu
                #    na sekundarnej zlozke sa sa zvysi pocet azimutalnych a polarnych krokov faktorom spocitanym z pomeru ich vzajomnych
                #    polomerov
                # 2) ak je to naopak, tak sa prepocita pocet krokov na primarnej zlozke
                # myslienka je udrzat homogenne pokrytie a to tak, aby bolo co najviac bodov na hviezdach
                if ppr <= spr:
                    primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
                    # zaokruhluje sa pocet krokov hore, radsej hustejsie pokrytie ako redsie
                    secondary_phi_steps, secondary_theta_steps = math.ceil((spr / ppr) * phi_steps), math.ceil(
                        (spr / ppr) * theta_steps)
                elif ppr > spr:
                    secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
                    primary_phi_steps, primary_theta_steps = math.ceil((ppr / spr) * phi_steps), math.ceil(
                        (ppr / spr) * theta_steps)
            else:
                # ak je premenna homo nastavena na False, tak sa neudrziava homogenita na zlozkach, ale na obe zlozky
                # sa nastavi pocet azimutalnych a polarnych krokov zo vstupu, prakticky teda mensia zlozka bude pokryta
                # hustejsie ako zlozka, ktora ma vasi polomer
                primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
                secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
        # v pripade, ze je snaha spocitat len jednu zo zloziek, tak sa nastavia parametre na danu zlozku
        # rovnako tak rozsah vonkajsieho for cyklu sa nastavi na jedno zbehnutie zodpovedajuce danej zlozke
        elif t_object == "primary":
            j_range_bottom, j_range_top = 0, 1
            primary_phi_steps, primary_theta_steps = phi_steps, theta_steps
        elif t_object == "secondary":
            j_range_bottom, j_range_top = 1, 2
            secondary_phi_steps, secondary_theta_steps = phi_steps, theta_steps
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Binary, function: get_3d_model_testing(), line: " + str(
                    Fn.lineno()) + ". Incorrect `t_object` parameter value.")
            return False

        # for cyklus prebehne cez objekty (primarnu a sekundarnu zlozku)
        for obj in range(j_range_bottom, j_range_top):
            # current_phi_steps a current_theta_steps je pocet azimutalnych a polarnych krokov pre danu hodnotu cyklus
            # teda bude pre primarnu alebo pre sekundarnu zlozku
            if obj == 0:
                current_phi_steps, current_theta_steps = primary_phi_steps, primary_theta_steps
            elif obj == 1:
                current_phi_steps, current_theta_steps = secondary_phi_steps, secondary_theta_steps

            # startovaci uhol theta
            # dlzka poalrneho kroku
            # partial je pole, do ktoreho sa budu ukladat kooordinaty bodov
            theta, theta_step_length, partial, equatorial, meridional = np.pi / 2., np.pi / current_theta_steps, [], [], []

            # posun v polarnom uhle cez rovinu x-z
            #   zo
            #	| o
            #   |  o
            #   |    o
            #   |      o
            #   |         o ->
            #   _________________x
            #
            # postupne tak ako je to zobrazene na obrazku

            for rot in range(0, int(current_theta_steps)):
                # vytvotenie prazdneho vektora pre sfericke a pre kartezianske suradnice
                vector_spheric, vector_xyz = np.arange(3, dtype=np.float), np.arange(3, dtype=np.float)
                # naplnenie spferickeho vektora vychdzou hodnotou pre vypocet na dalsej a dalsej kruznici (dana novou theta)
                # toto osetrenie cez if je tu len pre uplnost, pretoze povodne to bolo takto, ze tam bolo len
                # vector_spheric[1], vector_spheric[2] = 1., np.pi, theta
                # tu je problem, pretoze je divne mat ph = 180 a s napr. theta = 190, su to trochu protichodne uhly
                # funkcia pre transformacie sa s tym vysporiadala nejak tak, ze to fungovalo, ale nepacila sa mi korektnost
                # ked je raz vo sferickych, theta e <0, np.pi>, phi e <0, 2* np.pi>, tak nech je to dodrzane
                if theta <= np.pi:
                    vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1., np.pi, theta
                elif theta > np.pi:
                    vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1.0, 0.0, (2.0 * np.pi) - theta

                # nastavenie hlavneho polomeru pre dalsie vypocty podla momentalneho objektu
                if obj == 0:
                    star_radius = ppr
                elif obj == 1:
                    star_radius = spr

                # cast, ked je polarny uhol mensi ako np.pi + critical_angle
                # polarny uhol bezi stale ako np.pi == 90, 90 + current_theta_steps, 90 + 2 * current_theta_steps ... atd
                # keby nebol kriticky uhol, tak dosiahne az po theta = 270 stupnov, cize prejde celu akoby spodnu vetvu roviny x-z
                # ak je nastaveny kriticky uhol na nejaku hodnotu, tak klasicke ekvidistancne krokovanie pojde len pokial plati
                # podmienka theta <= np.pi + critical_angle, potom sa spravia iste kontroly (ak tato podmienka ine je splnena)
                if theta <= np.pi + critical_angle:
                    # vypocet polomeru v theta
                    # z
                    # |o     |      |
                    # |  o   |      |
                    # |    o | r    | star_radius
                    # |       o     |
                    # |          o  |
                    # |____________ o _________ x
                    # neviem preco tu bolo povodne zaokruhlenie
                    # r = round( star_radius * np.sin( theta - np.pi / 2. ), 4 )

                    r = star_radius * np.sin(theta - np.pi / 2.0)
                    # pocet krokov na kruznici v r (resp v theta), prepocitane pomerovo k polomeru danej mensej kruznice
                    # je to opodmeinkovane tak, ze ak sa jedna o prvy bod, ked sa zacina riesit, tak urcite treba jeden bod
                    # ale kedze je tam polomer r = 0, tak by aj pocet bodov bol nula, preto tam treba nastavit na tvrdo "1"
                    phi_points = 1 if theta == np.pi / 2.0 else math.ceil((r / star_radius) * current_phi_steps)
                    # pocet transformacnych krokov po kruznici
                    # opat opodmienakovane, ak sa jedna o prvy bod, tak je pocet transformacnych korkov 1 (keby bola 0, tak sa cyklus nevykona ani raz)
                    # v ostatnych pripadoch polarneho uhla theta treba + 1 transformacny krok
                    transform_steps = int(phi_points) if theta == np.pi / 2. else int(phi_points) + 1
                    # dlzka kroku, ktory sa bude pouzivat pri transformacii po kruznici
                    # musi sa stale menit, lebo sa meni polomer kruznice a stale sa meni phi_points
                    # nepodstatny je v principe prva dlzka, pretoze je to v zadnej casti hviezdy, phi = np.pi, theta = np.pi / 2.0
                    phi_step_length = (np.pi / 2.) / phi_points
                    # rotacny uhol je to iste ako phi_step_length, je to tu len znova priradene, kvoli nazvu premennej, ak by som to nechcel
                    # pre ilustracne ucely a lepsi prehlad v kode, treba to spravit tak, ze "rotation_angle = ( np.pi / 2. ) / phi_points"
                    # a riadok "rotation_angle = phi_step_length" odstranit
                    rotation_angle = phi_step_length



                ##############################################

                #   TU PRIDE TA KOMPLIKOVANA ADAPTIVNA CAST

                ##############################################

                else:
                    # podmineka theta <= np.pi + critical_angle: nie je splnena ale na druhej strane je theta mensia ako 270 stupnov,
                    # teda pokial sa theta uhol neprejde cely np.pi od np.pi / 2.0
                    #                           O 									   |z
                    # - np.pi/2. o               x                o  (3./2.) * np.pi   |
                    #			   o                           o 					  O|_____x
                    #                 o                     o
                    #                     o            o
                    #                          o   o
                    if theta < (3. / 2.) * np.pi:
                        # testovaci polarny uhol zacina z pozicie np.pi / 2.0 (takze z tej, v ktorej moze nastat problem pri wuma systemoch, na krku
                        # respektive v krku, kde ziaden bod povrchu nie je)
                        testing_theta = np.pi / 2.0
                        # podmienka while cyklu zebezpeci beh len po kriticky uhol, aby to nebezalo do nekonecna
                        while True and testing_theta < (np.pi - critical_angle):
                            # argumenty pre vypocet bodu na povrchu hviezdy
                            args = (actual_distance, 0.0, float(testing_theta))
                            # vypocet bodu na povrchu podobne ako v kode vyssie
                            # pre rozsiahlejsi a lepsi vystup pouzita funkcia fsolve
                            try:
                                if obj == 0:
                                    scipy_solve_point = self.primary.polar_radius / 10.0
                                    solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                                     scipy_solve_point,
                                                                                     full_output=True,
                                                                                     args=args)
                                elif obj == 1:
                                    scipy_solve_point = self.secondary.polar_radius / 10.0
                                    solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                                     scipy_solve_point,
                                                                                     full_output=True, args=args, )
                                if ier == 1 and not np.isnan(solution[0]):
                                    # ak sa naslo korektne riesenie, tak sa ukonci tento cyklus
                                    if 1 > solution[0] > 0: break
                                # posun v testovacom uhle o 1 stupen (mozno by bolo lepsie nastavit mensi krok, ale to by zas spomalilo algoritmus)
                                # ak sa nenaslo korektne riesenie, tak cyklus pokracuje
                                testing_theta += np.pi / 180.0
                            except:
                                # ak sa nenaslo korektne riesenie, tak cyklus pokracuje
                                testing_theta += np.pi / 180.0
                                continue

                        # prvy uhol theta pre adaptivny vypocet sa nastavi ten, kde sa prvykrat v cykle vyssie podarilo najst korektne riesenie
                        # samotny adaptivny uhol zacne tam, kde sa skoncilo ratanie klasickym cyklom na kritickom uhle
                        first_adaptive_theta, theta_adaptive, stop = testing_theta, ((2. * np.pi) - theta), False  #
                        # ak sa zisti, ze uz na theta = np.pi / 2.0 ja koretkne riesenie, tak sa zmeni kriticky uhol na np.pi / 2.0, to znamena,
                        # ze kalsicky cyklus sa zastavi az na 270 stupnoch
                        if first_adaptive_theta == np.pi / 2.0: critical_angle = np.pi / 2.0  # change critical_angle if type is EA or EB

                        # cyklus bezi do nastavenia premennej stop na True
                        # este je to opodmienkovane, aby nezbehol, ak sa vyssie nastavi kriticky uhol na np.pi / 2.0

                        # CYKLUS WHILE POLARNEHO UHLA ADAPTIVNEJ CASTI
                        while not stop and first_adaptive_theta != np.pi / 2.0:
                            # spocitanie polomeru v danom uhle theta (theta_adaptive)
                            r = abs(star_radius * np.sin(theta_adaptive - (np.pi / 2.)))

                            # nastavenie premennych pre zrotovanie po kruznici a vypocet bodov na povrchu hviezdy na tejto kruznici
                            phi_points = math.ceil((r / star_radius) * current_phi_steps)
                            transform_steps = int(phi_points) + 1
                            phi_step_length = (np.pi / 2.) / phi_points
                            rotation_angle = phi_step_length

                            # ked v cykle uz vybehne adaptivny uhol mimo (do oblasti, kde sa neda pocitat (hrdlo WUMa),
                            # tak sa este od kritickeho (prveho adaptivneho) uhla spocita jedna hodnota mierne posunuta o maly uhol, vlastne hranica akoby)
                            # je to na to, aby sa co najviac priblizilo vypoctom k hranici hrdla wuma systemov
                            # premenan stop sa nastavi na True, teda po dokonceni posledneho cyklu po kruznici sa uz viac nevykona
                            if theta_adaptive <= first_adaptive_theta: theta_adaptive, stop = first_adaptive_theta + additional_angle, True

                            vector_spheric, vector_xyz = np.arange(3, dtype=np.float), np.arange(3, dtype=np.float)
                            # tu je na rozdiel od klasickeho cyklu poalrneho uhla spravny polarny uhol, nebezi cez 180 stupnov
                            vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1., 0.0, theta_adaptive

                            # klasicky cyklus po kruznici v polohe "r"
                            # jeho popis je nizzsie pri klasickom cykle, ked sa nerata s kritickym uhlom
                            for transform_adaptive in range(0, transform_steps):
                                args, use = (actual_distance, vector_spheric[1], vector_spheric[2]), False
                                try:
                                    # pre istotu je v tejto casti pouzita funkcia pre riesenie korenov fsolve, pretoze newton moze zratat korektnu hodnotu,
                                    # ktora, ale nelezi na hviezde ale niekde vo vonakjsej casti hillovho priestoru
                                    if obj == 0:
                                        scipy_solve_point = self.primary.polar_radius / 10.
                                        solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn,
                                                                                         scipy_solve_point,
                                                                                         full_output=True, args=args)
                                    elif obj == 1:
                                        scipy_solve_point = self.secondary.polar_radius / 10.
                                        solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn,
                                                                                         scipy_solve_point,
                                                                                         full_output=True, args=args)
                                    if ier == 1 and not np.isnan(solution[0]):
                                        solution = solution[0]
                                        if actual_distance >= solution >= 0: use = True
                                except:
                                    use = False
                                if use:
                                    point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, \
                                                                                                    vector_spheric[1], \
                                                                                                    vector_spheric[2]
                                    xyz = Fn.spheric_to_cartesian(point_coordinate)
                                    if obj == 0:
                                        x = xyz[0]
                                    elif obj == 1:
                                        x = - (xyz[0] - actual_distance)
                                    y, z = xyz[1], xyz[2]
                                    if transform_adaptive == transform_steps - 1:
                                        equatorial.append([x, y, z])
                                    elif transform_adaptive == 0:
                                        meridional.append([x, y, z])
                                    else:
                                        partial.append([x, y, z])
                                vector_xyz = Fn.spheric_to_cartesian(vector_spheric)
                                rotate = Fn.rotate(vector=vector_xyz, angle=rotation_angle)
                                vector_spheric = Fn.cartesian_to_spheric(rotate)
                                if vector_spheric[1] < 0: vector_spheric[1] += (2. * np.pi)
                            # theta adaptive sa defaultne prvykrat nastavil na  2.0 * np.pi - theta, cize sa znormalizovala
                            # jeho pozicia na menej ako 180 stupnov, cize ak mal napr. polarny uhol theta hodnotu 225 stupnov
                            # (samozrjeme vsetko v radianoch), tak sa jeho hodnota zmenila na 135 stupnov, v tomto cykle
                            # polarneho uhla (je nim ten while vyssie) sa teda odpocitava z hodnoty uhla, aby sme sa pohybovali s uhlom
                            # v smere tak ako v povodnom for cykle, teda smerom od zapornej casti z osi ku kladnej casti z pohladu
                            # roviny x-z
                            # to preco sa odpocitava taky zvlastny krok??? preto lebo ked som to skusal, tak toto sa dalo pouzit,
                            # mozno by sa dalo este aj nieco logaritmicke, ale toto funguje
                            theta_adaptive -= np.sin(theta_adaptive) * theta_step_length / 1.5

                        # toto je osetrenie, aby v pripade, ze naozaj bolo pouzite adaptivne pocitanie nedoslo k dokonceniu klasickeho for cyklus
                        # pre polarny uhol, pretoze k dokonceniu vypoctu povrchu hviezdnej zlozky doslo adaptvne
                        # tato cast je uz mimo adaptivnu, preto po break-u tu, dojde k ukonceniu polarneho for cyklu
                        if first_adaptive_theta > np.pi / 2.: break

                # cyklus po kruznici v polomere "r" (t.j. v uhle theta)
                for transform in range(0, transform_steps):
                    # argumenty pre funkciu potencialu
                    # premenna use pre rozhodnutie pouzitia vysledku numerickeho vypoctu (nemusi stale skonvergovat)
                    args, use = (actual_distance, vector_spheric[1], vector_spheric[2]), False
                    # pokus o vyriesenie implicitnej rovnice (pokus, preto je to v try)
                    try:
                        # pouzije sa prislusny potencial podla zlozky
                        if obj == 0:
                            scipy_solve_point = self.primary.polar_radius / 10.
                            solution = scipy.optimize.newton(self.primary_potential_fn, scipy_solve_point, args=args)
                        elif obj == 1:
                            scipy_solve_point = self.secondary.polar_radius / 10.
                            solution = scipy.optimize.newton(self.secondary_potential_fn, scipy_solve_point, args=args)
                        # otestuje sa, ci vratena hodnota nie je NaN
                        if not np.isnan(solution):
                            # tu je taky logicky test, hviezda isto nemoze byt vacsi ako je vzdialenost medzi hviezdami
                            # a risenie urcite nemoze byt zaporne, kedze sa jedna o polomer vo sferickych suradniciach
                            # ak je vsetko yhovujuce, tak sa nastavi premenna pre pouzitie daneho riesenia na True
                            if actual_distance >= solution >= 0: use = True
                    except:
                        use = False

                    # ak sa riesenie pouzije
                    if use:
                        # risenie sa ulozi do vektora a pretransformuje prislusnou funkciou zo sferickcyh do kartezianskych suradnic
                        point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, vector_spheric[1], \
                                                                                        vector_spheric[2]
                        xyz = Fn.spheric_to_cartesian(point_coordinate)
                        # ak sa jedna o primarnu zlozku, tak sa priamo pouzije x-ova suradnica
                        if obj == 0:
                            x = xyz[0]
                        # ak sa jedna o sekundarnu zlozku, je potrebne urobit urcite upravy, pretoze body sa pocitali po preneseni
                        # sekundarnej zlozky do pociatku suradnicoveho systemu, treba preto tuto zlozku posunut o vzdialenost hviezd
                        # a zrdkadlovo ju otocit
                        elif obj == 1:
                            x = - (xyz[0] - actual_distance)
                        y, z = xyz[1], xyz[2]
                        if transform == transform_steps - 1:
                            equatorial.append([x, y, z])
                        elif transform == 0:
                            meridional.append([x, y, z])
                        else:
                            partial.append([x, y, z])
                    # pricahdza na rad posunutie sa po kruznici v polomere "r"
                    # vektor sferickych suradnic, ktory ukazuje na dane miesto na kruznici sa pretrasformuje do kartezianskych suradnic
                    vector_xyz = Fn.spheric_to_cartesian(vector_spheric)
                    # nasledne sa karteziansky vektor ukazujuci na dane miesto na kruznici zrotuje o prislusncy uhol (zodpovedajuci
                    # dlzke kroku v danom polomere "r")
                    rotate = Fn.rotate(vector=vector_xyz, angle=rotation_angle)
                    # z kartezianskeho vektora uz prerotovaneho sa spatne ziska sfericky vektor, aby sme vedeli prislusny azimutalny
                    # a polarny uhol pre dalsie zratanie bodu na povrchu hviezdy
                    vector_spheric = Fn.cartesian_to_spheric(rotate)
                    # funkcia pre transformaciu moze vratit zaporny azimutalny uhol, tak sa to tu upravi, aby to bolo stale kladne
                    if vector_spheric[1] < 0: vector_spheric[1] += (2. * np.pi)

                # posun dalej v polarnom uhle
                theta += theta_step_length

            # funkcia ma na vstupe ovladaciu premennu, tzv. zero_point; tato premenna ovlada moznost vypoctu bodu na povrchu
            # hviezdy v smere (phi = 0, theta = np.pi / 2.); je to takto riesene z dvoch dovodov
            # 1) pri krokovani v polarnom smere sa nemusi trafit funkcia presne na np.pi / 2.0
            # 2) ked sa spravne netrafi u napr. zlozky, ktora vyplna svoj rocheho lalok, tak model hviezdy bude dost nekorektny
            # je to v podstate len prekopirovana cast vypoctu z kody vyssie
            if zero_point:
                # try count point with (r, 0, 180)
                args = (actual_distance, 0.0, np.pi / 2.)
                try:
                    # rozdiel oproti kodu vyssie je v tom, ze sa tu pouziva ina funkcia pre numericke riesenie korenov, pretoze fcia
                    # fsolve() disponuje obsiahlejsim vystupom a informaciami o vypocte, resp. an konci vypoctu, ktore je potreben overit
                    # snaha zratat tento bod pri wuam systemoch s tym, ze by sa tam naozaj nejaky bod nedopatrenim doratal (v hrdle),
                    # by viedla k nekorektnemu modelu
                    if obj == 0:
                        scipy_solve_point = self.primary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.primary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    elif obj == 1:
                        scipy_solve_point = self.secondary.polar_radius / 10.
                        solution, info, ier, msg = scipy.optimize.fsolve(self.secondary_potential_fn, scipy_solve_point,
                                                                         full_output=True, args=args)
                    if ier == 1 and not np.isnan(solution[0]):
                        solution = solution[0]
                        if actual_distance >= solution >= 0:
                            use = True
                except:
                    use = False
                if use:
                    point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, 0.0, np.pi / 2.
                    xyz = Fn.spheric_to_cartesian(point_coordinate)
                    if obj == 0:
                        x = xyz[0]
                    elif obj == 1:
                        x = - (xyz[0] - actual_distance)
                    y, z = xyz[1], xyz[2]
                    # riesenie sa ulozi do premennej z_point a nie priamo do pola partial, pretoze s tym sa potom bude
                    # este dalej pracovat tak, ze sa budu body preklapat a tento bod, by sa len duplikoval, co by pripadne
                    # pri pouziti cgal triangulacie sposobilo pad kodu
                    # partial.append([x, y, z])
                    z_point = [x, y, z]

            # prvy bod (phi = np.pi, theta = np.pi / 2.0) sa oddeli od rovnika, aby nedoslo k jeho zduplikovaniu (resp.
            # ku zrkadlovemu preklopeniu cez x os, co povedie k dvom bodom, ktore budu od seba vzdialene radovo 1e-15)
            f_point, equatorial = equatorial[0], equatorial[1:]

            # dopni sa prvy bod (first point), len raz, aby nedoslo k jeho zduplikovaniu
            # doplni sa zero point, aby nedoslo k jeho duplikacii

            full_partial = []
            if Fn.empty(z_point) and zero_point:
                if self.verbose:
                    print(Fn.color_string(color="error",
                                          string="ValueError: ") + "In class: Binary, function: get_3d_model_optimized(), line: " + str(
                        Fn.lineno()) + ". Variable `zero_point` is set to True, but there is not correct solution. Target object: " + str(
                        obj) + ".")
                return False
            elif not Fn.empty(z_point) and zero_point:
                full_partial.append(z_point)

            if not Fn.empty(partial) and not Fn.empty(equatorial) and not Fn.empty(meridional):
                full_partial.append(f_point)
                for point in partial:
                    full_partial.append(point)
                    full_partial.append([point[0], -point[1], point[2]])
                    full_partial.append([point[0], point[1], -point[2]])
                    full_partial.append([point[0], -point[1], -point[2]])

                for point in equatorial:
                    full_partial.append(point)
                    full_partial.append([point[0], -point[1], point[2]])
                for point in meridional:
                    full_partial.append(point)
                    full_partial.append([point[0], point[1], -point[2]])

                total = total + full_partial

                if obj == 0:
                    model['primary'] = full_partial
                elif obj == 1:
                    model['secondary'] = full_partial
            else:
                if self.verbose:
                    print(Fn.color_string(color="error",
                                          string="ValueError: ") + "In class: Binary, function: get_3d_model_optimized(), line: " + str(
                        Fn.lineno()) + ". One of lists (`partial`, `equatorial`, `meridional`) is empty.")
                return False

        # computing of separation neck
        if self.binary_morph == "over-contact":
            # import objects.Plot as Plt
            # Plt.plot_3d(vertices=[total], normals_view=False,
            #             points_view=True, faces_view=False, point_color="r", point_size=3.0,
            #             verbose=self.verbose)

            separation, neck_radius = self.get_separation(actual_distance=actual_distance)
            r = min([self.primary.polar_radius, self.secondary.polar_radius])

            # cylindric phi angle steps and cylindric phi angle
            cps, cp = int(math.ceil(phi_steps * r / neck_radius)) * 2, 0.0
            step_length = 2.0 * np.pi / cps
            plane = []
            while True:
                if cp >= 2.0 * np.pi:
                    break
                args = actual_distance, cp, separation
                use = False
                solution = scipy.optimize.newton(self.wc_potential_fn,
                                                 min([self.primary.polar_radius / 10.0,
                                                      self.secondary.polar_radius / 10.0]),
                                                 args=args)
                # otestuje sa, ci vratena hodnota nie je NaN
                if not np.isnan(solution):
                    if max([self.primary.polar_radius, self.secondary.polar_radius]) >= solution >= 0:
                        use = True
                if use:
                    plane.append([separation, solution * np.cos(cp), solution * np.sin(cp)])
                cp += step_length

            # reorganise vertices (if any vertex of primary are located behind separation value put it to secondary
            # and vice-versa)
            primary, secondary = [], []
            for t in total:
                if t[0] > separation:
                    secondary.append(t)
                elif t[0] < separation:
                    primary.append(t)
            p, s = np.concatenate((primary, plane), 0), np.concatenate((secondary, plane), 0)
            u, v, w = len(p), len(s), len(plane)
            # average spacing
            avsp = [Fn.average_spacing(data=p, neighbours=4), Fn.average_spacing(data=s, neighbours=4)]

            remove, v, o, l = [[], []], [p, s], [None, None], [u, v]

            for t_o in [0, 1]:
                for plane_vertex in plane:
                    a = ([ndx for ndx, x in list(zip(range(0, l[t_o] - w), v[t_o][0:l[t_o] - w]))
                          if np.linalg.norm(np.array(plane_vertex) - np.array(x)) < avsp[t_o] / 2.0])
                    remove[t_o] = np.concatenate((remove[t_o], a), axis=0)

                remove[t_o] = np.array(list(set(remove[t_o])), dtype="int32").tolist()
                o[t_o] = np.array([x for ndx, x in list(zip(range(0, l[t_o] - w), v[t_o][0:l[t_o] - w]))
                                   if ndx not in remove[t_o]]).tolist()

            model["primary"] = np.concatenate((o[0], plane), 0)
            model["secondary"] = np.concatenate((o[1], plane), 0)
            model["system"] = np.concatenate((o[0], o[1], plane), 0)
            model["separation"] = separation, neck_radius

            # import Plot.Plot as Plt
            # Plt.plot_3d(vertices=[model["primary"]], faces=None, normals_view=False, points_view=True,
            #             faces_view=False, point_color="r", point_size=5.0, face_color="c", edge_color="k")

            del (primary, secondary, p, s, o, v, total)
            return model

        if total:
            model["separation"] = None
            model['system'] = np.array(total)
        return model

    def get_orbit(self):
        return self.orbit

    def get_exception(self):
        return self.exception

    def get_separation(self, actual_distance=None):
        linear_steps, x, phi = 0.005, 0.0, np.pi / 2.0
        separation = []
        while True:
            args = actual_distance, phi, x
            try:
                use = False
                solution = scipy.optimize.newton(self.wc_potential_fn,
                                                 min([self.primary.polar_radius / 10.0,
                                                      self.secondary.polar_radius / 10.0]),
                                                 args=args)
                # otestuje sa, ci vratena hodnota nie je NaN
                if not np.isnan(solution):
                    if max([self.primary.polar_radius, self.secondary.polar_radius]) >= solution >= 0:
                        use = True
                if use:
                    separation.append([x, solution])
            except:
                pass

            x += linear_steps
            if x >= 1:
                break

        z = np.polyfit(list(zip(*separation))[0], list(zip(*separation))[1], 17)
        p = np.poly1d(z)
        # find minimum of function
        crit = p.deriv().r
        r_crit = crit[crit.imag == 0].real
        test = p.deriv(2)(r_crit)
        x_min = [d for d in r_crit[test > 0] if 0.0 < d < 1.0][0]

        # xs = np.arange(0., 0.7, 0.01)
        # ys = p(xs)
        # import matplotlib.pyplot as plt
        # plt.scatter(list(zip(*separation))[0], list(zip(*separation))[1], c="b")
        # plt.plot(xs, ys, c="r")
        # plt.axis('equal')
        # plt.show()

        return x_min, p(x_min)

    def create_spots(self, meta=None):

        def set_false(x):
            for t in x:
                t.append(False)

        spots, sizes, centers, alt_centers, boundaries = [], [], [], [], []
        spots_meta = meta
        x_separation = None

        # if wuma system, get separation and make additional test to location of each point (if primary
        # spot doesn't intersect with secondary, if does, then such spot will be skiped completly)
        if self.binary_morph == "over-contact":
            x_separation = self.get_separation(actual_distance=self.orbit.get_periastron_distance())[0]

        for meta in spots_meta:
            lon, lat, diameter = meta["lon"], meta["lat"], meta["diameter"]  # lon -> phi, lat -> theta
            stps_azimuthal, stps_radial = meta["steps_azimuthal"], meta["steps_radial"]

            # angle to move in radius of spot
            mov_angle = float(diameter * 0.5) / float(stps_radial)

            # radial vector (for spot point computation)
            # radial_v = Fn.spheric_to_cartesian(vector=np.array([1.0, lon, lat]))
            radial_v = np.array([1.0, lon, lat])  # unit radial vector to center of current spot

            boundary, spot, solution = [], [], False
            args, use = (self.orbit.get_periastron_distance(), radial_v[1], radial_v[2]), False

            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)
            # radial distance to the center of current spot from begining of coordinate system
            spot_center_r = solution

            # # for False testing
            # if meta["diameter"] == 0.2:
            #     use = False

            if not use:
                # in case of spots, every point should be usefull, otherwise skip current spot computation
                set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                continue

            center = Fn.spheric_to_cartesian(vector=np.array([solution, radial_v[1], radial_v[2]])).tolist()
            spot.append(center)

            # adaptive azimuthal steps in plane of spot
            n0 = int(stps_azimuthal)

            # we have to obtain distance between center and firts point in first circle of spot
            args, use = (self.orbit.get_periastron_distance(), lon, lat + mov_angle), False
            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)

            # # for False testing
            # if meta["diameter"] == 0.2:
            #     use = False
            if not use:
                # in case of spots, every point should be usefull, otherwise skip current spot computation
                set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                continue

            x0 = np.sqrt(spot_center_r ** 2 + solution ** 2 - (2.0 * spot_center_r * solution * np.cos(mov_angle)))

            # if problem will occured in while loop down bellow, we have to break current spot computation
            # while will break and next_step will set to False
            next_step = True

            for i in range(0, stps_radial):
                spheric_v = np.array([1.0, lon, lat + (mov_angle * (float(i) + 1.0))])

                if spheric_v[1] < 0:
                    spheric_v[1] += 2.0 * np.pi

                # adaptive steps in plane of spot
                ni = n0 * (float(i) + 1.0)
                # ni = n0 * ((float(i) + 1.0) * x0) / x0

                rot_angle, phi = 2.0 * np.pi / ni, 0.0

                while phi < (2.0 * np.pi - (rot_angle * 0.5)):
                    args, use = (self.orbit.get_periastron_distance(), spheric_v[1], spheric_v[2]), False
                    solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)

                    # # for False testing
                    # if meta["diameter"] == 0.2:
                    #     use = False

                    if not use:
                        # in case of spots, every point should be usefull, otherwise skip current spot computation
                        set_false(x=[alt_centers, boundaries, spots, centers, sizes])
                        next_step = False
                        break

                    current_point = np.array([solution, spheric_v[1], spheric_v[2]])
                    xyz = Fn.spheric_to_cartesian(vector=current_point)

                    # spot.append([self.x_flip(t_object=meta["t_object"], x=xyz[0],
                    #             star_separation=self.binary.periastron_distance), xyz[1], xyz[2]])

                    spot.append(xyz.tolist())

                    cartesian_v = Fn.spheric_to_cartesian(vector=spheric_v)
                    # vector radial_v is normalized in fucntion arbitrary_rotation
                    rot_v = Fn.arbitrary_rotate(theta=rot_angle, vector=cartesian_v,
                                                omega=Fn.spheric_to_cartesian(radial_v))
                    spheric_v = Fn.cartesian_to_spheric(vector=rot_v, degrees=False)

                    if spheric_v[1] < 0:
                        spheric_v[1] += 2.0 * np.pi
                    phi += rot_angle

                    if i == int(stps_radial) - 1:
                        boundary.append(xyz.tolist())

            if not next_step:
                continue

            # alternative center of spot (obtained from spot boundary)
            # boundary center of mass
            # return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])

            b_com = sum(np.array(boundary)) / len(boundary)
            spheric_v = Fn.cartesian_to_spheric(vector=b_com, degrees=False)

            # i am not testing use, beacuse if points of entire spot exist, then this center point also has to
            args, use, solution = (self.orbit.get_periastron_distance(), spheric_v[1], spheric_v[2]), False, False
            solution, use = self.solve(args, meta["t_object"], x_separation=x_separation)
            current_point = np.array([solution, spheric_v[1], spheric_v[2]])
            alt_center = Fn.spheric_to_cartesian(vector=current_point).tolist()

            spot.insert(0, alt_center)
            # ----------------------------

            sizes.append(max([np.linalg.norm(np.array(alt_center) - np.array(b)) for b in boundary]))
            if meta["t_object"] == "primary":
                alt_centers.append(alt_center)
                boundaries.append(boundary)
                spots.append(spot)
                centers.append(center)
            elif meta["t_object"] == "secondary":
                alt_centers.append([self.x_flip(t_object="secondary", x=alt_center[0],
                                                star_separation=self.orbit.get_periastron_distance()),
                                    -alt_center[1], alt_center[2]])

                centers.append([self.x_flip(t_object="secondary", x=center[0],
                                            star_separation=self.orbit.get_periastron_distance()),
                                -center[1], center[2]])

                b = [[self.x_flip(t_object="secondary", x=x[0],
                                  star_separation=self.orbit.get_periastron_distance()),
                      -x[1], x[2]] for x in boundary]
                boundaries.append(b)

                s = [[self.x_flip(t_object="secondary", x=x[0],
                                  star_separation=self.orbit.get_periastron_distance()),
                      -x[1], x[2]] for x in spot]
                spots.append(s)

            else:
                return False

        spots_d = {"primary": [], "secondary": []}
        for meta, i in list(zip(spots_meta, range(0, len(spots_meta)))):
            if not Fn.empty(spots[i]):
                norms = Geo.normal_estimation(binary_object=self,
                                                   actual_distance=self.orbit.get_periastron_distance(),
                                                   vertices=np.array(spots[i]), t_object=meta["t_object"])

                spots_d[meta["t_object"]].append(
                    {"vertices": spots[i], "center": centers[i], "alt_center": alt_centers[i],
                     "size": sizes[i], "boundary": boundaries[i], "norms": norms, "meta": meta})

        if not Fn.empty(spots_d["primary"]):
            self.primary.spots = spots_d["primary"]
        if not Fn.empty(spots_d["secondary"]):
            self.secondary.spots = spots_d["secondary"]

        return spots_d


    def solve(self, args, t_object=None, x_separation=None):
        # args = actual_distance, phi, theta
        args, use, solution = args, False, False
        try:
            scipy_solve_point = self.primary.polar_radius / 10.0 if t_object == "primary" else self.secondary.polar_radius / 10.0
            potential_fn = self.primary_potential_fn if t_object == "primary" else self.secondary_potential_fn
            solution, _, ier, _ = scipy.optimize.fsolve(potential_fn, scipy_solve_point, full_output=True, args=args)

            if ier == 1 and not np.isnan(solution[0]):
                solution = solution[0]
                if 1 > solution > 0:
                    use = True
        except:
            use = False

        # test if point is not "behind" separation in case of wuma (if not wuma, then x_separation is set to None, and
        # condition is skipped)
        if x_separation is not None and use == True:
            # x value
            x = Fn.spheric_to_cartesian([solution, args[1], args[2]])[0]
            x = x if t_object == "primary" else self.x_flip(t_object="secondary", x=x,
                                                            star_separation=args[0])

            # also equal, i dont care abot bullshit spots
            if (t_object == "primary" and x >= x_separation) or (t_object == "secondary" and x <= x_separation):
                use = False

        return solution, use

    @classmethod
    def x_flip(cls, t_object=None, x=None, star_separation=None):
        if t_object == "primary":
            return x
        if t_object == "secondary":
            return - (x - star_separation)


    def get_binary_morphology(self):
        return self.binary_morph