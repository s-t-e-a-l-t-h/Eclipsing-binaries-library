#!/usr/bin/python

""" SelfTest trieda
-------------------
  Trieda na vnutorne testovanie funkcii celeho balika

Inicializacne parametre:
------------------------
  Staticka trieda, ziadne inicializacne parametre

Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
------------------------------------------------------------------------------------------------------

  init_binary_mass_test()
  =======================

    Vstupne parametre:
    ------------------
    rand_values : [integer], defaultna hodnota "100", premenna urcujuca pocet testov, pocet dvojic hmotnosti pre system
    potential : [float], defaultna hodnota "100", premenna definujuca potencial pre obe zlozky pouzite vo vsetkych
                         testoch
    min_mass : [float], defaultna hodnota "0.2", minimalna hmotnost, ktora sa nahodne generuje
    max_mass : [float], defaultna hodnota "0.2", maximalna hmotnost, ktora sa nahodne generuje

    Return:
    -------
    void

    Popis:
    ------
    funkcia vygeneruje rand_values dvojic hmotnosti a a snazi sa pri zadanom potential vygenerovat binarny system
    v pripade poruchy vypise problem

    Bug:
    ----
    ziaden znamy bug

  lagrangian_points_test()
  ========================

    Vstupne parametre:
    ------------------
    rand_values : [integer], defaultna hodnota "100", premenna urcujuca pocet testov, pocet dvojic hmotnosti pre system
    potential : [float], defaultna hodnota "100", premenna definujuca potencial pre obe zlozky pouzite vo vsetkych
                         testoch
    min_mass : [float], defaultna hodnota "0.2", minimalna hmotnost, ktora sa nahodne generuje
    max_mass : [float], defaultna hodnota "0.2", maximalna hmotnost, ktora sa nahodne generuje

    Return:
    -------
    void

    Popis:
    ------
    funkcia vygeneruje rand_values dvojic hmotnosti a a snazi sa pri zadanom potential spocitat lagrangeove body

    Bug:
    ----
    ziaden znamy bug
"""

import objects.Star as Star
import objects.Binary as Binary
# import globe.variables as gv
# import objects.Plot as Plt
import objects.Function as Fn
# import objects.Geometry as Geo
# import sys

# import numpy as np
import time as time


def init_binary_mass_test(rand_values=100, potential=100, min_mass=0.2, max_mass=20):
    start_time = time.time()

    random_primary_mass = Fn.rand(min_mass, max_mass, rand_values, time.time())
    random_scondary_mass = Fn.rand(min_mass, max_mass, rand_values, random_primary_mass[0])

    stars_combination = []

    for pm, sm in list(zip(random_primary_mass, random_scondary_mass)):
        stars_combination.append([
            Star.Star(mass=pm, synchronicity_parameter=1., potential=potential, effective_temperature=5000.,
                      gravity_darkening=1., metallicity=0., albedo=0.),
            Star.Star(mass=sm, synchronicity_parameter=1., potential=potential, effective_temperature=5000.,
                      gravity_darkening=1., metallicity=0., albedo=0.),
        ])

    for system in stars_combination:
        bin_sys = Binary.Binary(primary=system[0], secondary=system[1], system="eb", verbose=True)
        if not bin_sys.init:
            print("Problem in binary system with parameters:")
            print("primary mass:\t\t" + str(bin_sys.primary.mass))
            print("   potential:\t\t" + str(bin_sys.primary.potential))
            print("secondary mass:\t\t" + str(bin_sys.secondary.mass))
            print("     potential:\t\t" + str(bin_sys.secondary.potential))
            # tato testovacia funkcia je napisana tak, ze sa na obe zlozky pouziva rovnaky potentcial a preto staci vypisovat len jeden filling_factor
            print("Filling factor:\t\t" + str(bin_sys.primary.filling_factor))
            print("Inner potential:\t" + str(bin_sys.potential_inner))
            print("Outer potential:\t" + str(bin_sys.potential_outer))
            print("Potential dif:\t\t" + str(bin_sys.df_potential))
            print("Critical:\t\t" + str(bin_sys.primary.critical_potential))
            print("_________________________________________")
            print("")

    ### hodnoty pri ktorych raz vybehol error na primarnom polarnom polomere a viackrat sa uz neprejavil
    ### primary.potential = secondary.potentcial = 20.0
    ### primary.mass = 1.00116324935
    ### secondary.mass = 4.24889197062


    end_time = time.time()
    print("init_binary_mass_test() finished, elapsed time:" + str(round(end_time - start_time, 5)) + "sec")


def lagrangian_points_test(rand_values=100, potential=100, min_mass=0.2, max_mass=20):
    # tato funkcia je len pre dualne pretestovanie, v skutocnosti sa totiz uz sama odreze pri nicializacii triedy Binary
    # tam je tottiz vypocet kritickych bodov a tej je v try: except: konstrukcii, takze ak tam nastane nejaka chyba,
    # tak sa nastavi init premenna na False a uz sa potom tu netestuje
    start_time = time.time()

    random_primary_mass = Fn.rand(min_mass, max_mass, rand_values, time.time())
    random_scondary_mass = Fn.rand(min_mass, max_mass, rand_values, random_primary_mass[0])

    stars_combination = []

    for pm, sm in list(zip(random_primary_mass, random_scondary_mass)):
        stars_combination.append([
            Star.Star(mass=pm, synchronicity_parameter=1., potential=potential, effective_temperature=5000.,
                      gravity_darkening=1., metallicity=0., albedo=0.),
            Star.Star(mass=sm, synchronicity_parameter=1., potential=potential, effective_temperature=5000.,
                      gravity_darkening=1., metallicity=0., albedo=0.),
        ])

    for system in stars_combination:
        # tu je vypnute verbose, pretoze ak by bolo zpanute, platilo by to pre vsetko v triede Binary a teda aj pre
        # vypocet polarneho polomeru, kde by to tym padom pritovalo chybu a mna teraz nezaujima, ze je nefyzikalna
        # hodnota velicin pre cely system, ten sa proste podla init vyhodi prec
        bin_sys = Binary.Binary(primary=system[0], secondary=system[1], system="eb", verbose=False)
        if bin_sys.init:
            # zmena verbose, aby to zacalo vypisovat chybove hlasky
            bin_sys.verbose = True

            lp = bin_sys.get_lagrangian_points(actual_distance=1.0, t_object="secondary")

            # po dopocitani, to treba opat vypnut, pretoze sa bude vytvarat opat novy system
            bin_sys.verbose = False

            # ak je lagrangeovych bodov menej ako 3 v zozname, tak sa niekde stala chyba a vypise sa to
            if len(lp) != 3:
                print("Lagrangian points")
                print("Problem in binary system with parameters:")
                print("primary mass:\t\t" + str(bin_sys.primary.mass))
                print("   potential:\t\t" + str(bin_sys.primary.potential))
                print("secondary mass:\t\t" + str(bin_sys.secondary.mass))
                print("     potential:\t\t" + str(bin_sys.secondary.potential))
                print("_________________________________________")
                print("")
    end_time = time.time()
    print("lagrangian_points_test() finished, elapsed time:" + str(round(end_time - start_time, 5)) + "sec")
