#!/usr/bin/python

""" Star trieda
---------------
  Vytvori objekt Star, objekt hviezdy s parametrami prisluchajucimi hviezde.

Inicializacne parametre:
------------------------
  mass : [float], defaultna hodnota "None", hmotnost hviezdy v hmotnostiach slnka
  synchronicity_parameter : [float], defaultna hodnota "None", parameter synchronicity hviezdy v binarnom
  systeme, bezrozmerna
  potential : [float], defaultna hodnota "None", hodnota potencialu urcujuca velkost a tvar hviezdy v Rosheovskom
  zmysle, bezrozmenra,
  effective_temperature : [float], defaultna hodnota "None", efektivna teplota hviezdy v Kelvinoch,
  gravity_darkening : [float], defaultna hodnota "None", beta faktor gravitacneho stemnenia (vid Lucy law a pod.),
  bezrozmerna velicina
  metallicity : [float], metallicita hviezdy v jednotkach [M/H], defaultna hodnota "None"
  albedo : [float], bezrozmenray faktor albeda hviezdy, defaultna hodnota "None"

Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
------------------------------------------------------------------------------------------------------

  get_info()
  ==========

    Vstupne parametre:
    ------------------
    ziadne

    Return:
    -------
    void

    Popis:
    ------
    vypise niektore informacie o danej instancii triedy Star

    Bug:
    ----
    ziaden znamy bug
"""

import pymysql
pymysql.install_as_MySQLdb()
import MySQLdb

import numpy as np
import globe.variables as gv
import objects.Function as Fn
import objects.Geometry as Geo
import objects.Atmosphere as Atm
from scipy import interpolate


class Star:
    def __init__(
            self,
            mass=None,  # unit: solar mass //
            synchronicity_parameter=None,  # unit: [-]
            potential=None,  # unit: [-]
            effective_temperature=None,  # unit: [K]
            gravity_darkening=None,  # unit: [0 - 1]
            metallicity=None,  # unit: [M/H]
            albedo=None,  # unit: [0 - 1]
            verbose=False,
            atmoshpere="kurucz",
            phi_steps=10,
            theta_steps=20,
            cgal=None
    ):
        # <variables>
        if cgal is None:
            cgal = {"min_triangle_angle": np.radians([20.0])[0],
                    "max_triangle_size": 10,
                    "surface_aproximation_error": 0.3,
                    "to_average_spacing": 5}

        self.cgal = cgal
        self.verbose = verbose
        self.init = True
        self.exception = []

        self.phi_steps = phi_steps
        self.theta_steps = theta_steps

        if min(gv.METALLICITY_LIST_ATM) <= metallicity <= max(gv.METALLICITY_LIST_ATM):
            self.metallicity = metallicity
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Star, function: __init___(), line: " + str(
                    Fn.lineno()) + ". Variable `metallicity` is out of range. Use value in range [" + str(
                    min(gv.METALLICITY_LIST_ATM)) + "; " + str(max(gv.METALLICITY_LIST_ATM)) + "].")
            self.exception.append("ValueError: In class: Star, function: __init___(), line: " + str(
                Fn.lineno()) + ". Variable `metallicity` is out of range. Use value in range [" + str(
                min(gv.METALLICITY_LIST_ATM)) + "; " + str(max(gv.METALLICITY_LIST_ATM)) + "].")
            self.init = False

        if min(gv.TEMPERATURE_LIST_ATM) <= effective_temperature <= max(gv.TEMPERATURE_LIST_ATM):
            self.effective_temperature = effective_temperature
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Star, function: __init___(), line: " + str(
                    Fn.lineno()) + ". Variable `effective_temperature` is out of range. Use value in range [" + str(
                    min(gv.TEMPERATURE_LIST_ATM)) + "; " + str(max(gv.TEMPERATURE_LIST_ATM)) + "].")
            self.exception.append("ValueError: In class: Star, function: __init___(), line: " + str(
                Fn.lineno()) + ". Variable `effective_temperature` is out of range. Use value in range [" + str(
                min(gv.TEMPERATURE_LIST_ATM)) + "; " + str(max(gv.TEMPERATURE_LIST_ATM)) + "].")
            self.init = False

        self.mass = float(mass)
        self.synchronicity_parameter = synchronicity_parameter
        self.potential = float(potential)
        self.gravity_darkening = gravity_darkening
        self.albedo = albedo
        self.backward_radius = None
        self.filling_factor = None
        self.critical_potential = None
        self.lagrangian_points = None
        self.atmosphere = atmoshpere

        self.polar_radius = None
        self.polar_gravity = None
        self.polar_gradient_norm = None
        self.gravity_scalling_factor = None
        self.polar_tempterature = None

        # <distributions>
        self.local_gravity = None
        self.gradient_norm = None
        self.local_temperature = None
        self.gravity_darkening_factor = None
        self.local_radiation_power = None
        # </distributions>

        # <surface>
        self.faces = None
        self.vertices = None
        self.faces_orientation = None
        self.simplices = None
        # </surface>

        if not self.init and self.verbose:
            print(
            Fn.color_string(color="error", string="InitError: ") + "In class: Star, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation.")
            self.exception.append("InitError: In class: Star, function: __init__(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during initialisation.")

            # </variables>

    def get_info(self, output=False):

        if output:
            info = ["Info:\n",
                    "<" + str(self.__class__.__name__) + "> -------------------------------------------------\n",
                    "mass:\t\t\t\t\t\t\t" + str(self.mass) + "\n",
                    "potential:\t\t\t\t\t\t" + str(self.potential) + "\n",
                    "effective temperature:\t\t\t" + str(self.effective_temperature) + "\n",
                    "metallicity:\t\t\t\t\t" + str(self.metallicity) + "\n",
                    "gravity darkening:\t\t\t\t" + str(self.gravity_darkening) + "\n",
                    "albedo:\t\t\t\t\t\t\t" + str(self.albedo) + "\n",
                    "synchronicity parameter:\t\t" + str(self.synchronicity_parameter) + "\n",
                    "polar radius:\t\t\t\t\t" + str(self.polar_radius) + "\n",
                    "backward radius:\t\t\t\t" + str(self.backward_radius) + "\n",
                    "gravity scalling factor:\t\t" + str(self.gravity_scalling_factor) + "\n",
                    "polar gradient norm:\t\t\t" + str(self.polar_gradient_norm) + "\n",
                    "polar gravity:\t\t\t\t\t" + str(self.polar_gravity) + "\n",
                    "polar temperature:\t\t\t\t" + str(self.polar_tempterature) + "\n",
                    "phi steps:\t\t\t\t\t\t" + str(self.phi_steps) + "\n",
                    "theta steps:\t\t\t\t\t" + str(self.theta_steps) + "\n"]

            if not Fn.empty(self.local_gravity):
                info.append("local gravity:\t\t\t\t\t" + str(
                    str(self.local_gravity[0]) + ", ... ," + str(self.local_gravity[-1])) + "\n")
                # info.append("local gravity:\t\t\t\t\t" + str(", ".join([str(x) for x in self.local_gravity])) + "\n")
            else:
                info.append("local gravity:\t\t\t\t\tNone\n")

            if not Fn.empty(self.gradient_norm):
                # info.append("gradient norm:\t\t\t\t\t" + str(", ".join([str(x) for x in self.gradient_norm])) + "\n")
                info.append("gradient norm:\t\t\t\t\t" + str(
                    str(self.gradient_norm[0]) + ", ... ," + str(self.gradient_norm[-1])) + "\n")

            else:
                info.append("gradient norm:\t\t\t\t\tNone\n")

            if not Fn.empty(self.gravity_darkening_factor):
                # info.append("gravity darkening factor:\t\t" + str(", ".join([str(x) for x in self.gravity_darkening_factor])) + "\n")
                info.append("gravity darkening factor:\t\t" + str(
                    str(self.gravity_darkening_factor[0]) + ", ... ," + str(self.gravity_darkening_factor[-1])) + "\n")
            else:
                info.append("gravity darkening factor:\t\tNone\n")

            if not Fn.empty(self.local_temperature):
                # info.append("local temperature:\t\t\t\t" + str(", ".join([str(x) for x in self.local_temperature])) + "\n")
                info.append("local temperature:\t\t\t\t" + str(
                    str(self.local_temperature[0]) + ", ... ," + str(self.local_temperature[-1])) + "\n")
            else:
                info.append("local temperature:\t\t\t\tNone\n")

            if not Fn.empty(self.local_radiation_power):
                # info.append("local radiation:\t\t\t\t" + str(", ".join([str(x) for x in self.local_radiation_power])) + "\n")
                info.append("local radiation:\t\t\t\t" + str(
                    str(self.local_radiation_power[0]) + ", ... ," + str(self.local_radiation_power[-1])) + "\n")
            else:
                info.append("local radiation:\t\t\t\tNone\n")
            info.append('</' + str(self.__class__.__name__) + '> -----------------------------------------------\n')
            return info
        else:
            print()
            print(Fn.color_string("info", "Info:"))
            print("<" + str(self.__class__.__name__) + "> -------------------------------------------------")
            print("mass:\t\t\t\t\t\t\t" + str(self.mass))
            print("potential:\t\t\t\t\t\t" + str(self.potential))
            print("effective temperature:\t\t\t" + str(self.effective_temperature))
            print("metallicity:\t\t\t\t\t" + str(self.metallicity))
            print("gravity darkening:\t\t\t\t" + str(self.gravity_darkening))
            print("albedo:\t\t\t\t\t\t\t" + str(self.albedo))
            print("synchronicity parameter:\t\t" + str(self.synchronicity_parameter))
            print("polar radius:\t\t\t\t\t" + str(self.polar_radius))
            print("backward radius:\t\t\t\t" + str(self.backward_radius))
            print("gravity scalling factor:\t\t" + str(self.gravity_scalling_factor))
            print("polar gradient norm:\t\t\t" + str(self.polar_gradient_norm))
            print("polar gravity:\t\t\t\t\t" + str(self.polar_gravity))
            print("polar temperature:\t\t\t\t" + str(self.polar_tempterature))
            print("phi steps:\t\t\t\t\t\t" + str(self.phi_steps))
            print("theta steps:\t\t\t\t\t" + str(self.theta_steps))

            print("\nDistributions:")

            if not Fn.empty(self.local_gravity):
                print("local gravity:\t\t\t\t\t" + str(self.local_gravity[0]) + ", ... ," + str(self.local_gravity[-1]))
            else:
                print("local gravity:\t\t\t\t\t" + str(self.local_gravity))

            if not Fn.empty(self.gradient_norm):
                print("gradient norm:\t\t\t\t\t" + str(self.gradient_norm[0]) + ", ... ," + str(self.gradient_norm[-1]))
            else:
                print("gradient norm:\t\t\t\t\t" + str(self.gradient_norm))

            if not Fn.empty(self.gravity_darkening_factor):
                print("gravity darkening factor:\t\t" + str(self.gravity_darkening_factor[0]) + ", ... ," + str(
                    self.gravity_darkening_factor[-1]))
            else:
                print("gravity darkening factor:\t\t" + str(self.gravity_darkening_factor))

            if not Fn.empty(self.local_temperature):
                print("local temperature:\t\t\t\t" + str(self.local_temperature[0]) + ", ... ," + str(
                    self.local_temperature[-1]))
            else:
                print("local temperature:\t\t\t\t" + str(self.local_temperature))

            if not Fn.empty(self.local_radiation_power):
                print("local radiation:\t\t\t\t" + str(self.local_radiation_power[0]) + ", ... ," + str(
                    self.local_radiation_power[-1]))
            else:
                print("local radiation:\t\t\t\t" + str(self.local_radiation_power))

            print('</' + str(self.__class__.__name__) + '> -----------------------------------------------')
            print()

    def compute_polar_temperature(
            self
    ):
        if self.verbose:
            print(Fn.color_string("info", "Info: ") + "Computing polar temperature for " + str(self))

        try:
            surface = Geo.triangle_surface_area(triangle=self.faces, verbose=self.verbose)
            gdf = self.gravity_darkening_factor
            self.polar_tempterature = self.effective_temperature * (surface.sum() / (gdf * surface).sum()) ** 0.25
        except:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In class: Star, function: compute_polar_temperature(), line: " + str(
                    Fn.lineno()) + ". Error has been occurred during computatuion process.")
            self.exception.append("Error: In class: Star, function: compute_polar_temperature(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during computatuion process.")
            return False
        return self.polar_tempterature

    def compute_temperature_distribution(
            self,
            indices=None
    ):
        if self.verbose:
            print(Fn.color_string("info", "Info: ") + "Computing temperature distribution for " + str(self))

        if not Fn.is_numpy_array(arr=self.gravity_darkening_factor) or Fn.empty(var=self.gravity_darkening_factor):
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In class: Star, function: compute_temperature_distribution(), line: " + str(
                    Fn.lineno()) + ". Variable `self.gravity_darkening_factor` is invalid")
            self.exception.append("Error: In class: Star, function: compute_temperature_distribution(), line: " + str(
                Fn.lineno()) + ". Variable `self.gravity_darkening_factor` is invalid")
            return False

        try:
            if Fn.empty(indices): indices = np.arange(0, len(self.get_gravity_darkeninig_factor_dsitribution()), 1)
            gravity_darkening_factor = Fn.array_mask(array=self.get_gravity_darkeninig_factor_dsitribution(),
                                                     mask=indices)
            self.local_temperature = gravity_darkening_factor ** 0.25 * self.polar_tempterature
            return self.local_temperature
        except:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In class: Star, function: compute_temperature_distribution(), line: " + str(
                    Fn.lineno()) + ". Error has been occurred during temperature distribution computing.")
            self.exception.append("Error: In class: Star, function: compute_temperature_distribution(), line: " + str(
                Fn.lineno()) + ". Error has been occurred during temperature distribution computing.")

    def compute_gravity_distribution(
            self,
            gradnorm=None,
            indices=None,
            log=False
    ):
        if self.verbose:
            print(Fn.color_string(color="info", string="Info: ") + "Computing local gravity for " + str(self) + ".")

        if Fn.empty(var=gradnorm):
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="EmptyVariableError: ") + "In class: Star, function: local_gravity(), line: " + str(
                    Fn.lineno()) + ". Variable `gradnorm` is empty.")
            self.exception.append("EmptyVariableError: In class: Star, function: local_gravity(), line: " + str(
                Fn.lineno()) + ". Variable `gradnorm` is empty.")
            return False

        local_gravity = []
        if Fn.empty(indices): indices = np.arange(0, len(gradnorm), 1)

        for idx in indices:
            current_g = gradnorm[idx] * self.gravity_scalling_factor if not log else np.log10(
                gradnorm[idx] * self.gravity_scalling_factor)

            # if current_g > self.get_polar_gravity(): current_g = self.get_polar_gravity()

            local_gravity.append(current_g)

        self.local_gravity = np.array(local_gravity)
        self.gradient_norm = gradnorm

        return self.local_gravity

    def gravity_darkening_factor_distribution(self):
        if self.verbose:
            print(Fn.color_string("info", "Info: ") + "Computing gravity darkening factor for " + str(self))
        if Fn.empty(self.gradient_norm) or Fn.empty(self.polar_gradient_norm) or Fn.empty(self.gravity_darkening):
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="EmptyVariableError: ") + "In class: Star, function: gravity_darkening_factor_distribution(), line: " + str(
                    Fn.lineno()) + ". Variable one of variables (`self.gradient_norm`, `self.polar_gradient_norm`, `self.gravity_darkening`) " "is empty.")
            self.exception.append(
                "EmptyVariableError: In class: Star, function: gravity_darkening_factor_distribution(), line: " + str(
                    Fn.lineno()) + ". Variable one of variables (`self.gradient_norm`, `self.polar_gradient_norm`, `self.gravity_darkening`) is empty.")
            return False
        try:
            gdf = [(self.gradient_norm[idx] / self.polar_gradient_norm) ** self.gravity_darkening for idx in
                   range(0, len(self.gradient_norm))]
        except:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In class: Star, function: gravity_darkening_factor_distribution(), line: " + str(
                    Fn.lineno()) + ". Error has been occurred during computing of gravity darkening factor. Probably wrong input data.")
            self.exception.append(
                "Error: In class: Star, function: gravity_darkening_factor_distribution(), line: " + str(
                    Fn.lineno()) + ". Error has been occurred during computing of gravity darkening factor. Probably wrong input data.")
            return False

        self.gravity_darkening_factor = np.array(gdf)
        return self.gravity_darkening_factor

    def radiation_power(
            self,
            passband_model=None,
            passband_range=None,
            passband=None,  # len pre castelli-kurucz-04
            wavelength_step=10,  # len pre atmosferu black-body
            indices=None
    ):
        if self.verbose:
            print(Fn.color_string("info", "Info: ") + "Computing radiation power.")

        from scipy import integrate

        if self.atmosphere == "black-body":
            temperature_distribution = self.get_temperature_distribution()

            # indexy su dolzeite pre eccentricke drahy, pretoze tam pri kazdej polohe, je zbytocne pocitat tok ziarenia
            # na kazdom trojuholniku; ak nie je fazetu vidno, tak jej tok ziarenia nebude pouzitelny pre dalisu polohu
            # na orbite, pretoze sa zmeni tozlozenie teploty;
            # u velicin ako norma gradientu, orientacia faziet, distribucia faktorov gravitacneho stemnenia, je potrebne
            # ratat ich hodnotu aj na fazetach, ktore nie su viditelne, pretoze vsetupuju do vypoctu teploty aj na fazetach,
            # ktore viditelne su
            if Fn.empty(indices): indices = np.arange(0, len(temperature_distribution), 1)
            temperature_distribution = Fn.array_mask(array=temperature_distribution, mask=indices)
            temperature_distribution = np.around(temperature_distribution)
            flux = np.array([None for _ in range(len(temperature_distribution))])
            # repeated_temperature = Fn.find_repeated(temperature_distribution)
            # tato funkcia je vseobecnejsia ako ta vyssie v zakomentovanom riadku (i ked trebalo by testnut, ktora je
            # rychlejsia, kedze tu v planckovom pripade staci aj fcia find_repeated, ono problem je ak sa jedna o
            # atmosfery, kde vsetupuje aj gravitacne zrychlenie a teda plosny element sa vyznacuje viac ako jednym
            # parametrom a teda je potrebne kontrolovat zhodu v celom poli napr, kolkokrat sa nachadza v poli
            # hodnota, ktora je opat pole [10000.0, 2.0])
            repeated_temperature = Fn.find_repeated_indices(temperature_distribution)

            print(Fn.color_string("info", "Info: ") + "Using black body atmosphere model.")

            # flux = []
            # for temperature in self.get_temperature_distribution():
            #     wavelength_step, wavelength = wavelength_step, passband_range[0]  # nm, nm
            #     x_values, y_values = [], []
            #     while wavelength <= passband_range[1]:
            #         y_values.append(
            #             Atm.planck(wavelength=wavelength, temperature=temperature, passband_fn=passband_model))
            #         x_values.append(wavelength * 1e-9)
            #         wavelength += wavelength_step
            #     flux.append(np.pi * integrate.simps(y_values, x_values))

            # optimalizacia voci zakomentoanemu kodu je postavena na tom, ze sa zaokruhli teplota na cele islo
            # (v principe aj tak nema zmysel hrat sa na desatiny kelvina), zisti sa teplota, ktora sa opakuje v poli,
            # nasledne sa spocita flux len pre jednu z opakujucich sa hodnot a na ostatne sa nakopiruje ta ista hodnota
            # najviac casu zabera prave integracia hodnot na kazdu fazetu, preto to trebalo obist, aby sa v pripade,
            # ze sa ineco opakuje nemusel robit duplicinty vypocet;
            # takto doslo k vyrznemu usetreniu casu

            for repeated in repeated_temperature:
                wavelength_step, wavelength = wavelength_step, passband_range[0]  # nm, nm
                x_values, y_values = [], []
                while wavelength <= passband_range[1]:
                    y_values.append(
                        Atm.planck(wavelength=wavelength, temperature=temperature_distribution[repeated[0]],
                                   passband_fn=passband_model))
                    x_values.append(wavelength * 1e-9)
                    wavelength += wavelength_step
                flx = np.pi * integrate.simps(y_values, x_values)

                for idx in repeated: flux[idx] = flx
            return flux
        elif self.atmosphere == "castelli-kurucz-04":
            if self.verbose:
                print(Fn.color_string("info", "Info: ") + "Using Castelli & Kurucz 2004 - ATLAS9 atmosphere model.")

            # redukcia teploty, gravitacneho zrychlenia a metalicity
            temperature_distribution = np.around(self.get_temperature_distribution())
            gravity_distribution = np.around(np.log10(self.get_gravity_distribution()), 2)

            if Fn.empty(indices):
                if len(gravity_distribution) != len(temperature_distribution):
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="ValueError: ") + "In class: Star, function: radiation_power(), line: " + str(
                            Fn.lineno()) + ". Variables `gravity_distribution` and `temperature_distributin` have a different length.")
                    self.exception.append("ValueError: In class: Star, function: radiation_power(), line: " + str(
                        Fn.lineno()) + ". Variables `gravity_distribution` and `temperature_distributin` have a different length.")
                    return False
                else:
                    indices = np.arange(0, len(temperature_distribution), 1)

            # aplikovanie masky ak nahodou existuje (preto ta poznamka nahodou, lebo ak nie, tak vyssie sa vytvori
            # len na oko, aby nemusela ist dalsia podmienka)
            temperature_distribution = Fn.array_mask(array=temperature_distribution, mask=indices)
            gravity_distribution = Fn.array_mask(array=gravity_distribution, mask=indices)
            metallicity = Fn.find_nearest_value(array=gv.METALLICITY_LIST_ATM, value=self.metallicity)[0]

            # preprocesing dat na interpolaciu
            interpolation_method = np.array(["linear", "nearest"])

            # vytiahnutie vestkych dat z DB pre interpolaciu pre dane M/H
            mysql_conn = MySQLdb.connect(host=gv.HOST,  # your host, usually localhost
                                         user=gv.USER,  # your username
                                         passwd=gv.PWD,  # your password
                                         db="elisa_assets")  # name of the data base
            mysql_cur = mysql_conn.cursor()

            q = "SELECT " \
                "intensity, temperature, gravity " \
                "FROM ck04_intensity " \
                "WHERE metallicity = " + \
                str(metallicity) + " AND filter = \"" + str(passband) + "\" " \
                                                                        "ORDER BY intensity"

            mysql_cur.execute(q)
            db_data = mysql_cur.fetchall()
            mysql_conn.close()

            # vstupne data pre interpolaciu
            points, values = [], []
            for row in db_data:
                points.append([row[1], row[2]])  # [temperature, gravity]
                values.append(row[0])  # [intensity]

            # vytvorenie zoznamu
            distribution = [[t, g] for t, g in zip(temperature_distribution, gravity_distribution)]
            t_interval = [min(temperature_distribution), max(temperature_distribution)]
            g_interval = [min(gravity_distribution), max(gravity_distribution)]

            # ak je min/max hodnota v distribucii mimo maximalne dovolene extrapolacne/interpolacne pole tak program
            # skonci, nema zmysel pocitat s hodnotami, ktore su s velkou pravdepodobnostou mimo reality
            if 0.0 > g_interval[0] or g_interval[1] > 5.0 or 3500.0 > t_interval[0] or t_interval[1] > 50000.0:
                if self.verbose:
                    print(Fn.color_string(color="error",
                                          string="ValueError: ") + "In class: Star, function: radiation_power(), line: " + str(
                        Fn.lineno()) + ". Physical quantities distributien is out of inter/extra - polation range.")
                self.exception.append("ValueError: In class: Star, function: radiation_power(), line: " + str(
                    Fn.lineno()) + ". Physical quantities distributien is out of inter/extra - polation range.")
                return False

            # existencna matica
            np.set_printoptions(threshold=np.nan)
            ext_matrix = Fn.chess_matrix(db_name="elisa_assets", band=passband)

            # uvw - pole s hodnotami pre intra/extra polaciu
            # distribution - nove pole pre distibuciu
            distribution, uvw = [], []

            for t, g in zip(temperature_distribution, gravity_distribution):
                t_nearest, g_nearest = Fn.find_nearest_value(array=gv.TEMPERATURE_LIST_ATM, value=t), \
                                       Fn.find_nearest_value(array=gv.GRAVITY_LIST_ATM, value=g)

                matrix_value = ext_matrix[t_nearest[1] + 1][g_nearest[1] + 1]

                # interpolacia, resp. v priapde, ze teplota je v rozsahu 3250.0 <= t < 3500, < 50000 t < 51000
                # tak extrapolacia; tento rozsah je osetreny v podmienke vyssie
                if matrix_value == 1:
                    distribution.append([t, g])

                # extrapolacia
                # ak najblizsia hodnota v smere teploty a gravitacneho zrychlenia spada mimo ck04;
                # v takom pripade je pripuste doextrapolovat len vtedy, ak je v okoli niekde v existencnej matici
                # hodnota jedna [1] aspon v 4 pripadoch
                elif matrix_value == 0:
                    ext_array_pos = Fn.matrix2array([t_nearest[1] + 1, g_nearest[1] + 1])
                    a2m_pos = [Fn.array2matrix(ext_array_pos + direction) for direction in gv.CHESS_DIRECTION[0]]
                    matrix_around_value = np.array([ext_matrix[val[0]][val[1]] for val in a2m_pos])

                    # extrapolacia je dovolena len ak je jednotka v poli aspon v 4 pripadoch; ak sa tam nachadza
                    # v menej pripadoch, tka bud stoji momentalna hodnota v nulovom poli alebo strasne mimo rozsah
                    the_one = len(np.where(matrix_around_value == 1)[0])

                    # for val in a2m_pos:
                    #     ext_matrix[val[0]][val[1]] = 5
                    # print(ext_matrix)

                    if the_one < 4:
                        if self.verbose:
                            print(Fn.color_string(color="error",
                                                  string="ValueError: ") + "In class: Star, function: radiation_power(), line: " + str(
                                Fn.lineno()) + ". Physical quantities distributien is out of inter/extra - polation range.")
                        self.exception.append("ValueError: In class: Star, function: radiation_power(), line: " + str(
                            Fn.lineno()) + ". Physical quantities distributien is out of inter/extra - polation range.")
                        return False

                    distribution.append([t, g])

            distribution = np.array(distribution)

            # !!! podhodnotu v poli musia byt listy nie np.array, lebo ta funkcia find_repeated_indices porovnava tvrdo
            # pomocou == a numpy potom pluje, ze pouzi np.all alebo np.any !!!
            distribution_repeated = Fn.find_repeated_indices(distribution)
            distribution_mask = [index[0] for index in distribution_repeated]
            uvw = Fn.array_mask(array=distribution, mask=distribution_mask)

            intensity = np.arange(0, len(distribution), 1) * 0.0
            try:
                interpolated_intensity = interpolate.griddata(np.array(points), np.array(values), uvw,
                                                              method=interpolation_method[0])

                if not Fn.empty(np.where(interpolated_intensity == np.nan)[0]):
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="ValueError: ") + "In class: Star, function: radiation_power(), line: " + str(
                            Fn.lineno()) + ". Variable `interpolated_intensity` contain a nan value.")
                    self.exception.append("ValueError: In class: Star, function: radiation_power(), line: " + str(
                        Fn.lineno()) + ". Variable `interpolated_intensity` contain a nan value.")
                    return False

                # hodnoty pre interpolaciu boli v takom poradi ake je poradie blokov opakujucich sa elementov v
                # distribution_repeated, teda ak bol prvy blok s opakujucimi sa hodnotami an indexoch [1, 2, 5],
                # tak sa rovnakym indexom v poli intensity priradi hodnota nainterpolovanej intenzity
                for i, intens in zip(range(0, len(interpolated_intensity)), interpolated_intensity):
                    for idx in distribution_repeated[i]:
                        intensity[idx] = intens

            except:
                if self.verbose:
                    print(Fn.color_string(color="error",
                                          string="ValueError: ") + "In class: Star, function: radiation_power(), line: " + str(
                        Fn.lineno()) + ". Error has been occured during interpolation.")
                self.exception.append("ValueError: In class: Star, function: radiation_power(), line: " + str(
                    Fn.lineno()) + ". Error has been occured during interpolation.")
                return False

            del ext_matrix, interpolated_intensity, uvw, distribution, \
                distribution_repeated, distribution_mask, gravity_distribution, temperature_distribution
            return intensity

    # <SETTERs>
    def set_faces(self, faces=None):
        self.faces = faces

    def set_simplices(self, simplices=None):
        self.simplices = simplices

    def set_vertices(self, vertices=None):
        self.vertices = vertices

    def set_faces_orientation(self, faces_orientation=None):
        self.faces_orientation = faces_orientation

    def set_gradient_norm(self, gradient_norm=None):
        self.gradient_norm = gradient_norm

    def set_radiation_power(self, radiation_power=None):
        self.local_radiation_power = radiation_power

    # < /SETTERs>

    # <GETTERs>
    def get_simplices(self):
        return self.simplices

    def get_gravity_distribution(self):
        return self.local_gravity

    def get_temperature_distribution(self):
        return self.local_temperature

    def get_gradient_distribution(self):
        return self.gradient_norm

    def get_gravity_darkeninig_factor_dsitribution(self):
        return self.gravity_darkening_factor

    def get_polar_gravity(self):
        return self.polar_gravity

    def get_gravity_scalling_factor(self):
        return self.gravity_scalling_factor

    def get_polar_gradient_norm(self):
        return self.polar_gradient_norm

    def get_polar_temperature(self):
        return self.polar_tempterature

    def get_faces_orientation(self):
        return self.faces_orientation

    def get_faces(self):
        return self.faces

    def get_metallicity(self):
        return self.metallicity

    def get_vertices(self):
        return self.vertices

    def get_radiation_power(self):
        return self.local_radiation_power

    def get_exception(self):
        return self.exception
        # < /GETTERs>
