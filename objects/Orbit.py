#!/usr/bin/python

""" Orbit trieda
-----------------
 Zadefinovanie orbitalnych vlastnosti

Inicializacne parametre:
------------------------


Inicilizacia:
-------------
  Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
"""

import numpy as np
import scipy as sp
import math
import objects.Function as Fn


# import scipy.optimize
# import globe.variables as gv
# import objects.Geometry as Geo


class Orbit:
    def __init__(
            self,
            orbital_period=None,
            eccentricity=None,
            inclination=None,
            argument_of_periastron=None,
            verbose=False
    ):
        self.verbose = verbose
        self.init = True
        self.inclination = None
        self.argument_of_periastron = argument_of_periastron

        if self.argument_of_periastron >= np.pi / 2.0:
            self.argument_of_periastron_azimut = argument_of_periastron - np.pi / 2.0
        else:
            self.argument_of_periastron_azimut = 2.0 * np.pi - self.argument_of_periastron

        if inclination < 0 or inclination > np.pi:
            if self.verbose:
                print(Fn.color_string("error",
                                      "ValueError: ") + "In reference: Orbit, function: __init__, line: " + str(
                    Fn.lineno()) + " Variable `inclination` is invalid. Use `inclination` in range [0, pi].")
            self.init = False
        else:
            self.inclination = inclination

        self.relative_semimajor_axis = None
        self.orbital_period = orbital_period
        self.eccentricity = eccentricity
        self.mean_angular_velocity = self.angular_velocity()  # $\rad.sec^{-1}$

        try:

            self.conjuction = self.conjuction(eccentricity=self.eccentricity,
                                              argument_of_periastron=self.argument_of_periastron)

            # hodnoty pre primarny [0] a sekundarny [1] zakryt

            # prava anomalia konjunkcie (teda merana od periastra, od apsidalnej priamky)
            # \nu_{con}
            self.true_anomaly_of_conjuction = [self.conjuction[0]["true_anomaly"], self.conjuction[1]["true_anomaly"]]

            # excentricka anomalia konjunkcie, merana od apsidalnej priamky
            self.eccentric_anomaly_of_conjunction = [self.conjuction[0]["eccentric_anomaly"],
                                                     self.conjuction[1]["eccentric_anomaly"]]

            # stredna anomalia konjunkcie, merana od apsidalnej priamky
            self.mean_anomaly_of_conjunction = [self.conjuction[0]["mean_anomaly"], self.conjuction[1]["mean_anomaly"]]

            # skutocna faza periastra a konjunkcie je totozna, neviem co som tym chcel povedat, ked som to programoval
            self.true_phase_of_conjunction = [self.conjuction[0]["true_phase"], self.conjuction[1]["true_phase"]]

            # \Phi_{per}
            # faza periastra merana od nuly suradnioveho systemu po nulu sytemu orbitalnej drahy, v kladnom smere teda po
            # periastrum toto je vlastne shift na prevod medzi fazou orbitalnou a fotometrickou, kedze orbitalna je merana
            # od apsidalnej priamky

            self.true_phase_of_periastron = [self.conjuction[0]["true_phase"], self.conjuction[1]["true_phase"]]

            # standardne by tam ako shift slo self.true_phase_of_periastron[0] a vod 0.0 do 0.0 fotometrickej fazy
            # by zratalo polohu prave pre 0.0 fotometrickej, lenze ja chcem periastrum, tak to oblenem,
            # ze sa shiftnem o nulu od nuly fotometrickej, teda vypocet bude vychadzat z nuly a M_c = 0 len ak
            # ak sme na priamke apsid a teda v periastre;
            # keby som sa posunul o skutocnu fazu periastra od nuly (mysli sa casovy shift), tak sa posuniem v case
            # prave na M_c konjunkcie a to nechcem
            self.periastron_distance = self.orbital_motion(photometric_phase_from=0.0, photometric_phase_to=0.0,
                                                           shift=0.0,# self.true_phase_of_periastron[0],
                                                           photometric_phase_step=1.0, eccentricity=self.eccentricity,
                                                           argument_of_periastron=self.argument_of_periastron)[0][0]
        except:
            if not self.init and self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Orbit, function: __init__(), line: " + str(
                    Fn.lineno()) + ". Error has been occurred during initialisation.")

    def get_info(self, output=False):

        if output:
            info = ["Info:",
                    '\n<', str(self.__class__.__name__) + '> -------------------------------------------------\n',
                    "orbital period:\t\t\t\t\t\t" + str(self.orbital_period) + "\n",
                    "mean angular velocity:\t\t\t\t" + str(self.mean_angular_velocity) + "\n",
                    "eccentricity:\t\t\t\t\t\t" + str(self.eccentricity) + "\n",
                    "inclination:\t\t\t\t\t\t" + str(self.inclination) + "\n",
                    "argument of periastron:\t\t\t\t" + str(self.argument_of_periastron) +
                    ", [" + str(np.degrees(self.argument_of_periastron)) + "]\n"

                                                                           "relative semimajor axis:\t\t\t" + str(
                        self.relative_semimajor_axis) + "\n",
                    "periastron distance:\t\t\t\t" + str(self.periastron_distance) + "\n",
                    "true phase of periastron:\t\t\t" + str(self.true_phase_of_periastron) + "\n",
                    "true anomaly of conjuction:\t\t\t" + str(self.true_anomaly_of_conjuction) + "\n",
                    "eccentric anomaly of conjuction:\t" + str(self.eccentric_anomaly_of_conjunction) + "\n",
                    "mean anomaly of conjuction:\t\t\t" + str(self.mean_anomaly_of_conjunction) + "\n",
                    "true phase of conjuction:\t\t\t" + str(self.true_phase_of_conjunction) + "\n",
                    '/' + str(self.__class__.__name__) + ' -----------------------------------------------\n']
            return info
        else:
            print('\n' + str(self.__class__.__name__) + ' -------------------------------------------------')
            print("orbital period:\t\t\t\t\t\t" + str(self.orbital_period))
            print("mean angular velocity:\t\t\t\t" + str(self.mean_angular_velocity))
            print("eccentricity:\t\t\t\t\t\t" + str(self.eccentricity))
            print("inclination:\t\t\t\t\t\t" + str(self.inclination))
            print("argument of periastron:\t\t\t\t" + str(self.argument_of_periastron) + str(
                [np.degrees(self.argument_of_periastron)]))
            print("relative semimajor axis:\t\t\t" + str(self.relative_semimajor_axis))
            print("periastron distance:\t\t\t\t" + str(self.periastron_distance))
            print("true phase of periastron:\t\t\t" + str(self.true_phase_of_periastron))
            print("true anomaly of conjuction:\t\t\t" + str(self.true_anomaly_of_conjuction))
            print("eccentric anomaly of conjuction:\t" + str(self.eccentric_anomaly_of_conjunction))
            print("mean anomaly of conjuction:\t\t\t" + str(self.mean_anomaly_of_conjunction))
            print("true phase of conjuction:\t\t\t" + str(self.true_phase_of_conjunction))
            print('/' + str(self.__class__.__name__) + ' -----------------------------------------------')

    def angular_velocity(self):
        actual_distance = 1.0
        return ((2.0 * np.pi) / (self.orbital_period * 86400.0 * (actual_distance ** 2))) * np.sqrt(
            (1.0 - self.eccentricity) * (1.0 + self.eccentricity))

    @classmethod
    def mean_anomaly(cls, phase=None):
        return 2.0 * np.pi * phase

    @classmethod
    def mean_anomaly_fn(cls, eccentric_anomaly, *args):
        mean_anomaly, eccentricity = args
        return (eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)) - mean_anomaly

    @classmethod
    def eccentric_anomaly(cls, eccentricity=None, mean_anomaly=None):
        import scipy.optimize
        args = (mean_anomaly, eccentricity)
        try:
            solution = scipy.optimize.newton(cls.mean_anomaly_fn, 1.0, args=args, tol=1e-10)

            if not np.isnan(solution):
                if solution < 0: solution += 2.0 * np.pi
                return solution
            else:
                return False
        except:
            return False

    @classmethod
    def true_phase(cls, photometric_phase=None, shift=None):
        return photometric_phase + shift

    @classmethod
    def true_anomaly(cls, eccentricity=None, eccentric_anomaly=None):
        anomaly = 2.0 * np.arctan(
            np.sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) * np.tan(eccentric_anomaly / 2.0))
        if anomaly < 0: anomaly += 2.0 * np.pi
        return anomaly

    @classmethod
    def radius(cls, true_anomaly=None, eccentricity=None):
        return (1.0 - eccentricity ** 2) / (1.0 + eccentricity * np.cos(true_anomaly))

    @classmethod
    def true_anomaly_to_azimuth(cls, true_anomaly=None, argument_of_periastron=None):
        azimut = np.float64(true_anomaly + argument_of_periastron - np.pi * 0.5)
        if azimut > (2.0 * np.pi): azimut -= 2.0 * np.pi
        if azimut < 0: azimut += 2.0 * np.pi
        return np.float64(azimut)

    @classmethod
    def conjuction(cls, argument_of_periastron=None, eccentricity=None):
        ret = {}
        for alpha, idx in list(zip([np.pi / 2.0, 3.0 * np.pi / 2.0], [0, 1])):
            # prava anomalia konjunkcie (teda merana od periastra, od apsidalnej priamky)

            true_anomaly_of_conjuction = alpha - argument_of_periastron  # \nu_{con}
            if true_anomaly_of_conjuction < 0: true_anomaly_of_conjuction += 2.0 * np.pi

            # excentricka anomalia konjunkcie, merana od apsidalnej priamky
            eccentric_anomaly_of_conjunction = 2.0 * np.arctan(
                np.sqrt((1.0 - eccentricity) / (1.0 + eccentricity)) * np.tan(true_anomaly_of_conjuction / 2.0))

            if eccentric_anomaly_of_conjunction < 0:
                eccentric_anomaly_of_conjunction += 2.0 * np.pi

            # stredna anomalia konjunkcie, merana od apsidalnej priamky
            mean_anomaly_of_conjunction = eccentric_anomaly_of_conjunction - eccentricity * np.sin(
                eccentric_anomaly_of_conjunction)

            if mean_anomaly_of_conjunction < 0:
                mean_anomaly_of_conjunction += 2.0 * np.pi

            true_phase_of_conjunction = mean_anomaly_of_conjunction / (2.0 * np.pi)

            if true_phase_of_conjunction < 0:
                true_phase_of_conjunction += 1.0

            ret[idx] = {}
            ret[idx]["true_anomaly"] = true_anomaly_of_conjuction
            ret[idx]["eccentric_anomaly"] = eccentric_anomaly_of_conjunction
            ret[idx]["mean_anomaly"] = mean_anomaly_of_conjunction
            ret[idx]["true_phase"] = true_phase_of_conjunction
        return ret

    @classmethod
    def orbital_motion(cls, photometric_phase_from=None, photometric_phase_to=None, shift=None,
                       photometric_phase_step=None, eccentricity=None, argument_of_periastron=None):

        phase, position = photometric_phase_from, []
        while True:
            if phase > photometric_phase_to: break

            true_phase = cls.true_phase(photometric_phase=phase, shift=shift)

            # toto prerobi cislo vacsie ako 1.0 na desatine, napr 1.1 na 0.1
            if abs(true_phase) > 1: true_phase = math.modf(true_phase)[0]
            while true_phase < 0:
                true_phase += 1.0

            mean_anomaly = cls.mean_anomaly(phase=true_phase)
            eccentric_anomaly = cls.eccentric_anomaly(eccentricity=eccentricity, mean_anomaly=mean_anomaly)
            true_anomaly = cls.true_anomaly(eccentric_anomaly=eccentric_anomaly, eccentricity=eccentricity)

            if true_anomaly < 0: true_anomaly = (2.0 * np.pi) - abs(true_anomaly)

            actual_distance = cls.radius(true_anomaly=true_anomaly, eccentricity=eccentricity)
            azimuthal_angle = cls.true_anomaly_to_azimuth(true_anomaly=true_anomaly,
                                                         argument_of_periastron=argument_of_periastron)

            if azimuthal_angle < 0: azimuthal_angle = (2.0 * np.pi) - abs(azimuthal_angle)

            position.append([actual_distance, azimuthal_angle, true_anomaly, phase])

            phase += photometric_phase_step
        return position

    @classmethod
    def rotate_inclination(cls, inclination=None, verbose=False, arr=None):
        if inclination < 0.0 or inclination > np.pi:
            if verbose:
                print(Fn.color_string("error",
                                      "ValueError: ") + "In reference: Orbit, function: rotate_inclination, line: " + str(
                    Fn.lineno()) + " Variable `inclination` is invalid. Use `inclination` in range [0, pi].")
            return False

        if Fn.empty(arr):
            if verbose:
                print(Fn.color_string("error",
                                      "EmptyVariableError: ") + "In reference: Orbit, function: rotate_inclination, "
                                                                "line: " + str(Fn.lineno()) + " Variable `arr` is empty.")
            return False

        return np.array(
            [Fn.rotate(angle=(np.pi / 2.0) - abs(inclination), vector=val, inverse=True, axis="y") for val in arr])

    @classmethod
    def rotate_system(cls, normals=None, faces=None, vertices=None,
                      inclination=None, inclination_rotation=False, faces_rotation=False,
                      rotation_angle=None, verbose=False):

        if verbose:
            print(Fn.color_string("info", "Info: ") + "Rotating in angle " + str(rotation_angle))
        normals_to_return, faces_to_return, vertices_to_return = [], [], []
        for t_object in range(0, len(normals)):
            # local vertices to return
            vtr = [Fn.rotate(angle=rotation_angle, vector=np.array(vertex), inverse=False, axis="z")
                   for vertex in vertices[t_object]]

            if inclination_rotation:
                vtr = cls.rotate_inclination(inclination=inclination, verbose=verbose, arr=vtr)

            # local normals to return
            ntr = [Fn.rotate(angle=rotation_angle, vector=np.array(normal), inverse=False, axis="z")
                   for normal in normals[t_object]]

            if inclination_rotation:
                ntr = cls.rotate_inclination(inclination=inclination, verbose=verbose, arr=ntr)

            # rotacia faziet je volitelny parameter, je to potrebne len pre zobrazovanie, pre vypocet nie je potrebne
            # trojuholnik otocit, staci vediet, o ktory sa jedna
            ftr = []
            if faces_rotation:
                ftr = [np.array([Fn.rotate(angle=rotation_angle, vector=np.array(face[i]), inverse=False, axis="z")
                                 for i in range(0, 3)]) for face in faces[t_object]]

                if inclination_rotation:
                    ftr = cls.rotate_inclination(inclination=inclination, verbose=verbose,
                                                 arr=np.array(ftr).reshape(len(ftr) * 3, 3))
                    ftr = [np.array([ftr[idx * 3], ftr[(idx * 3) + 1], ftr[(idx * 3) + 2]]) for idx in
                           range(0, int(len(ftr) / 3))]
            else:
                ftr = faces[t_object]

            vertices_to_return.append(np.array(vtr))
            normals_to_return.append(np.array(ntr))
            faces_to_return.append(np.array(ftr))

        #     # normals_to_return = np.array([])
        return [np.array(faces_to_return), np.array(normals_to_return), np.array(vertices_to_return)]

    # <GETTERs>
    def get_mean_angular_velocity(self):
        return self.mean_angular_velocity

    def get_true_anomaly_of_conjuction(self):
        return self.true_anomaly_of_conjuction

    def get_true_phase_of_periastron(self):
        return self.true_phase_of_periastron

    def get_eccentric_anomaly_of_conjuction(self):
        return self.eccentric_anomaly_of_conjunction

    def get_mean_anomaly_of_conjuction(self):
        return self.mean_anomaly_of_conjunction

    def get_relative_semimajor_axis(self):
        return self.relative_semimajor_axis

    def get_orbital_period(self):
        return self.orbital_period

    def get_eccentricity(self):
        return self.eccentricity

    def get_periastron_distance(self):
        return self.periastron_distance

    def get_inclination(self):
        return self.inclination

    def get_argument_of_periastron(self):
        return self.argument_of_periastron

    # < /GETTERs>

    # <SETTERs>
    def set_argument_of_periastron(self, argument_of_periastron):
        self.argument_of_periastron = argument_of_periastron

    def set_relative_semimajor_axis(self, relative_semimajor_axis):
        self.relative_semimajor_axis = relative_semimajor_axis

    def set_inclination(self, inclination):
        self.inclination = inclination

    def set_mean_angular_velocity(self, mean_angular_velocity):
        self.mean_angular_velocity = mean_angular_velocity

    def set_orbital_period(self, orbital_period):
        self.orbital_period = orbital_period

    def set_eccentricity(self, eccentricity):
        self.eccentricity = eccentricity

    # < /SETTERs>

    # <testovacky>
    @classmethod
    def ellipse_inc(cls, t, *args):
        s, m = args
        return s - sp.special.ellipeinc(t, m)

    @classmethod
    def ellipse(cls, a, e, from_t, to_t, n):
        b = np.sqrt(a ** 2 * (1.0 - e ** 2))
        from_s = a * sp.special.ellipeinc(from_t, e ** 2)
        to_s = a * sp.special.ellipeinc(to_t, e ** 2)
        # step_length, s = (to_s - from_s) / n, from_s
        ellipse_length = a * sp.special.ellipeinc(2.0 * np.pi, e ** 2)

        step_length, s, x = ellipse_length / n, from_s, []
        # for k in np.arange(0, n + 1, 1):
        while s <= to_s:
            t = sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, e ** 2))
            s += step_length
            r = (a * b) / (np.sqrt((b * np.cos(t)) ** 2 + (a * np.sin(t)) ** 2))
            x.append([r * np.cos(t), r * np.sin(t)])
        return x

    @classmethod
    def ellipse_eqd(cls, a, e, from_t, to_t, n):
        b = np.sqrt(a ** 2 * (1.0 - e ** 2))
        from_s = a * sp.special.ellipeinc(from_t, e ** 2)
        to_s = a * sp.special.ellipeinc(to_t, e ** 2)

        # ked sa pouziva priama transformacia na x, y, to je ten zapis do premennej x, tak to musi byt posunute, aby
        # nula zacinala standardne vo vektore r = [1.0, 0.0]
        # from_s = a * sp.special.ellipeinc(from_t - np.pi / 2.0, e ** 2)
        # to_s = a * sp.special.ellipeinc(to_t - np.pi / 2, e ** 2)

        step_length, s = (to_s - from_s) / n, from_s
        # x = []
        y = []
        for _ in np.arange(0, n + 1, 1):
            t = sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, e ** 2))
            s += step_length

            # x.append([-a * np.sin(t), b * np.cos(t)])
            r = (a * b) / (np.sqrt((b * np.cos(t)) ** 2 + (a * np.sin(t)) ** 2))
            y.append([r * np.cos(t) - (a * e), r * np.sin(t)])

        import objects.Plot as Plt
        Plt.plot_2d(points=y, grid=True)

    # < /testovacky>

    # <BETAs>

    @classmethod
    def true_anomaly_to_center_angle(cls, ni, a, e):
        if e == 0: return [a, ni]
        l = (a * (1.0 - e ** 2)) / (1.0 + e * np.cos(ni))
        r = np.sqrt((a ** 2 * e ** 2) + l ** 2 - (2.0 * l * a * e * np.cos(np.pi - ni)))
        x = ((l ** 2 - r ** 2 - (a ** 2 * e ** 2)) / (2.0 * r * a * e)) * (-1)

        # toto je opodmienkovane, pretoze niekedy sa pri vypocte stane, ze to da hodnotu typu (-1.000000000000000000001)
        # a np.arcos() drbne error, proste blba numerika
        if -((l ** 2 - r ** 2 - (a ** 2 * e ** 2)) / (2.0 * r * a * e)) < -1: x = -1.0
        if -((l ** 2 - r ** 2 - (a ** 2 * e ** 2)) / (2.0 * r * a * e)) > 1: x = 1.0
        phi = np.arccos(x)

        # ta podmienka je tu preto, ze ak sa sem dostane uhol vacsi ako np.pi, tak arccos bude vracat stale len do
        # rozsahu np.pi, lebo tak je definovany, takze to treba obist
        if ni > np.pi: phi = 2.0 * np.pi - phi
        if ni < 0: phi = -phi

        return [r, phi]

    @classmethod
    def center_angle_to_true_anomaly(cls, phi, a, e, r):
        phi = Fn.angle_normalisation(phi)
        if e == 0: return [a, phi]

        l = np.sqrt((r ** 2) + (a ** 2 * e ** 2) - (2.0 * r * a * e * np.cos(phi)))
        x = np.float64(((a * (1.0 - e ** 2)) / (e * l) - (1.0 / e)))
        if x < -1:
            x = -1.0
        elif x > 1:
            x = 1.0
        ni = np.arccos(x)

        if phi > np.pi: ni = np.float64((2.0 * np.pi) - ni)
        return [l, ni]

    @classmethod
    def true_anomaly_to_phase(cls, true_anomaly=None, e=None, shift=0.0):
        eccentric_anomaly = Fn.angle_normalisation(
            2.0 * np.arctan(np.sqrt((1.0 - e) / (1.0 + e)) * np.tan(true_anomaly / 2.0)))
        mean_anomaly = Fn.angle_normalisation(eccentric_anomaly - e * np.sin(eccentric_anomaly))
        phase = (mean_anomaly / (2.0 * np.pi)) + shift
        return phase

    @classmethod
    def orbital_motion_beta(
            cls,
            from_photometric_phase=None,
            to_photometric_phase=None,
            n=None,  # !!! pocet krokov na jednu orbitu
            eccentricity=None,
            argument_of_periastron=None,
            adaptive_orbit_part=0.0,
            adaptive_multiplicator=1.0
    ):
        # funkcia vracia orbitu s ekvidistancnymi krokmi a k nim prisluchajucu fazu

        con = cls.conjuction(argument_of_periastron=argument_of_periastron, eccentricity=eccentricity)
        b, ellipse_length = np.sqrt(1.0 ** 2 * (1.0 - eccentricity ** 2)), 1.0 * sp.special.ellipeinc(2.0 * np.pi,
                                                                                                      eccentricity ** 2)

        step_length, phase = ellipse_length / float(n), from_photometric_phase
        short_step_length = step_length / adaptive_multiplicator
        shift = con[0]["true_phase"]

        position = []
        # <na debugovanie>
        # x, y, z = [], [], []
        # i = 0
        # < /na debugovanie>

        true_phase_start = cls.true_phase(photometric_phase=math.modf(from_photometric_phase)[0], shift=shift)
        mean_anomaly_start = cls.mean_anomaly(phase=true_phase_start)
        eccentric_anomaly_start = cls.eccentric_anomaly(eccentricity=eccentricity, mean_anomaly=mean_anomaly_start)
        true_anomaly_start = cls.true_anomaly(eccentric_anomaly=eccentric_anomaly_start, eccentricity=eccentricity)
        ellipse_start_angle = cls.true_anomaly_to_center_angle(ni=true_anomaly_start, a=1.0, e=eccentricity)[1]

        s = np.float64(1.0 * sp.special.ellipeinc(ellipse_start_angle, eccentricity ** 2))
        t_minus_on_phase = None

        # <adaptivny fazovy krok v oblasti zakrytov>
        # MYSLIENKA:
        # vstupne argumenty pre adaptivnost kroku pri zakrytoch adaptive_orbit_part a adaptive_multiplicator
        # adaptive_orbit_part je drahova cas, ktora sa prirata +/- k orbite od bodu zakrytu;
        # k nej sa dopocitaju prislusne centralne uhly na elipse, ktore sa potom pouzivaju ako adaptive, to znamena,
        # ze v danom rozsahu bude skrateny krok pohybu po orbitalnej drahe;
        # skratenie kroku sa robi pomocou premennej adaptive_multiplikator, ktorym sa vydeli dlzka standardneho
        # orbitalneho kroku;
        # myslienka jednoducha, vypocet v kode trochu krkolomny, ale funguje

        adaptive_length = ellipse_length * adaptive_orbit_part
        eclipse_length = \
            [1.0 * sp.special.ellipeinc(Fn.angle_normalisation(
                cls.true_anomaly_to_center_angle(ni=con[0]["true_anomaly"], a=1.0, e=eccentricity)[1]),
                eccentricity ** 2), 1.0 * sp.special.ellipeinc(Fn.angle_normalisation(
                cls.true_anomaly_to_center_angle(ni=con[1]["true_anomaly"], a=1.0, e=eccentricity)[1]),
                eccentricity ** 2)]
        eclipse_length = [[eclipse_length[0] - adaptive_length, eclipse_length[0] + adaptive_length],
                          [eclipse_length[1] - adaptive_length, eclipse_length[1] + adaptive_length]]

        adaptive_angle = [
            [Fn.angle_normalisation(
                sp.optimize.newton(cls.ellipse_inc, 0.01, args=(eclipse_length[0][0], eccentricity ** 2))),
                Fn.angle_normalisation(
                    sp.optimize.newton(cls.ellipse_inc, 0.01, args=(eclipse_length[0][1], eccentricity ** 2)))],

            [Fn.angle_normalisation(
                sp.optimize.newton(cls.ellipse_inc, 0.01, args=(eclipse_length[1][0], eccentricity ** 2))),
                Fn.angle_normalisation(
                    sp.optimize.newton(cls.ellipse_inc, 0.01, args=(eclipse_length[1][1], eccentricity ** 2)))
            ]
        ]

        # < /adaptivny fazovy krok v oblasti zakrytov>
        s_length, adaptive_part = step_length, True
        while True:
            if phase >= to_photometric_phase: break
            # tak ako je to naprogramovane momentalne, bude uhol t stale kladny
            # uhol t zodpoveda centralnemu uhlu prisluchajucemu momentalnej dlzke `s`, meranemu od hlavnej poloosi
            # v kladnom smere
            t = np.float64(
                Fn.angle_normalisation(sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, eccentricity ** 2))))

            # <adaptivny fazovy krok v oblasti zakrytov>
            if adaptive_angle[0][0] < adaptive_angle[0][1] and adaptive_angle[1][0] < adaptive_angle[1][1]:
                # print("statemnt 1")
                if (adaptive_angle[0][0] <= t <= adaptive_angle[0][1]) or (
                                adaptive_angle[1][0] <= t <= adaptive_angle[1][1]):
                    if adaptive_part is False:
                        a = None
                        # treba porozdelovat vsetky tieto moznosti, pretoze inak sa neda rozlisit, ktory uhol prislucha
                        # ktoremu pociatku adaptivneho useku
                        if adaptive_angle[0][0] <= t <= adaptive_angle[0][1]:
                            a = adaptive_angle[0][0]
                        elif adaptive_angle[1][0] <= t <= adaptive_angle[1][1]:
                            a = adaptive_angle[1][0]

                        s = 1.0 * sp.special.ellipeinc(a, eccentricity ** 2)
                        t = np.float64(Fn.angle_normalisation(
                            sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, eccentricity ** 2))))
                    s_length, adaptive_part = short_step_length, True
                else:
                    if adaptive_part: adaptive_part = False
                    s_length = step_length

            elif adaptive_angle[0][0] < adaptive_angle[0][1] and adaptive_angle[1][0] > adaptive_angle[1][1]:
                # print("statement 2")
                # print(str(np.degrees(t)) + str(np.degrees(adaptive_angle)))
                if (adaptive_angle[0][0] <= t <= adaptive_angle[0][1]) or (
                                    2.0 * np.pi >= t >= adaptive_angle[1][0]) or (0 <= t <= adaptive_angle[1][1]):
                    if adaptive_part is False:
                        a = None
                        # treba porozdelovat vsetky tieto moznosti, pretoze inak sa neda rozlisit, ktory uhol prislucha
                        # ktoremu pociatku
                        if adaptive_angle[0][0] <= t <= adaptive_angle[0][1]:
                            a = adaptive_angle[0][0]
                        elif 2.0 * np.pi >= t >= adaptive_angle[1][0] or 0 <= t <= adaptive_angle[1][1]:
                            a = adaptive_angle[1][0]

                        s = 1.0 * sp.special.ellipeinc(a, eccentricity ** 2)
                        t = Fn.angle_normalisation(
                            sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, eccentricity ** 2)))
                    s_length, adaptive_part = short_step_length, True
                else:
                    if adaptive_part: adaptive_part = False
                    s_length = step_length

            elif adaptive_angle[0][0] > adaptive_angle[0][1] and adaptive_angle[1][0] > adaptive_angle[1][1]:
                # print("statement 3")
                if (2.0 * np.pi >= t >= adaptive_angle[0][0]) or (0 <= t <= adaptive_angle[0][1]) or \
                        (2.0 * np.pi >= t >= adaptive_angle[1][0]) or (0 <= t <= adaptive_angle[1][1]):
                    if adaptive_part is False:
                        a = None
                        # treba porozdelovat vsetky tieto moznosti, pretoze inak sa neda rozlisit, ktory uhol prislucha
                        # ktoremu pociatku
                        if 2.0 * np.pi >= t >= adaptive_angle[0][0] or 0.0 <= t <= adaptive_angle[0][1]:
                            a = adaptive_angle[0][0]
                        elif 2.0 * np.pi >= t >= adaptive_angle[1][0] or 0.0 <= t <= adaptive_angle[1][1]:
                            a = adaptive_angle[0][1]

                        s = 1.0 * sp.special.ellipeinc(a, eccentricity ** 2)
                        t = Fn.angle_normalisation(
                            sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, eccentricity ** 2)))
                    s_length, adaptive_part = short_step_length, True
                else:
                    if adaptive_part: adaptive_part = False
                    s_length = step_length

            elif adaptive_angle[0][0] > adaptive_angle[0][1] and adaptive_angle[1][0] < adaptive_angle[1][1]:
                # print("statement 4")
                if (2.0 * np.pi >= t >= adaptive_angle[0][0]) or (0 <= t <= adaptive_angle[0][1]) or (
                                adaptive_angle[1][0] <= t <= adaptive_angle[1][1]):
                    # tato cas je tu na znizenie uhla, na zaciatok adaptivneho useku, pretoze pri velkom
                    # celkovom kroku, dojde k velkemu posunu voci pociatku adaptivneho useku
                    if True:
                        if adaptive_part is False:
                            a = None
                            # treba porozdelovat vsetky tieto moznosti, pretoze inak sa neda rozlisit, ktory uhol prislucha
                            # ktoremu pociatku
                            if adaptive_angle[1][0] <= t <= adaptive_angle[1][1]:
                                a = adaptive_angle[1][0]
                            elif 2.0 * np.pi >= t >= adaptive_angle[0][0] or 0 <= t <= adaptive_angle[0][1]:
                                a = adaptive_angle[0][0]

                            s = 1.0 * sp.special.ellipeinc(a, eccentricity ** 2)
                            t = Fn.angle_normalisation(
                                sp.optimize.newton(cls.ellipse_inc, 0.01, args=(s, eccentricity ** 2)))

                    s_length, adaptive_part = short_step_length, True
                else:
                    if adaptive_part: adaptive_part = False
                    s_length = step_length

            elif adaptive_angle[0][0] == adaptive_angle[0][1] and adaptive_angle[1][0] == adaptive_angle[1][1]:
                s_length = step_length
            # < /adaptivny fazovy krok v oblasti zakrytov>
            # print(str(np.degrees(t)) + str(adaptive_part) + str(np.degrees(adaptive_angle)))
            # klasicka vzdialenost od centra k danemu bodu na elipse prisluchajucemu `s`
            r = np.float64((1.0 * b) / (np.sqrt((b * np.cos(t)) ** 2 + (1.0 * np.sin(t)) ** 2)))
            # prevod centralnych do fokusalnych suradnic
            foo = cls.center_angle_to_true_anomaly(phi=t, a=1.0, e=eccentricity, r=r)

            # faza prisluchajuca pravej anomalii v mieste `s` na elipse
            t_on_phase = np.float64(cls.true_anomaly_to_phase(true_anomaly=foo[1], e=eccentricity, shift=0.0))

            # !!! toto sa nesnazit testovat funkciou Fn.empty() pretoze ta aj pre 0.0 vrati True,
            # ale t_minus_on_phase = 0.0 je validna hodnota pre pocitanie !!!
            if t_minus_on_phase is None:
                phase_diff, t_minus_on_phase = 0, t_on_phase
                # print "1. phase_diff = ", str(t_minus_on_phase), "- 0 =", str(phase_diff)
            else:
                if t_on_phase < t_minus_on_phase:
                    phase_diff, t_minus_on_phase = abs((1.0 - t_minus_on_phase) + t_on_phase), t_on_phase
                    # print "2. phase_diff = (1.0 -", str(t_minus_on_phase), ") +", str(t_on_phase),"=", str(phase_diff)
                else:
                    phase_diff, t_minus_on_phase = abs(t_on_phase - t_minus_on_phase), t_on_phase
                    # print "3. phase_diff = ", str(t_on_phase), "-", str(t_minus_on_phase),"=", str(phase_diff)

            phase += phase_diff
            s += s_length

            # print phase, cls.true_anomaly_to_phase(true_anomaly=foo[1], e=eccentricity, shift=0.0), shift, np.degrees(t)

            azimuth = cls.true_anomaly_to_azimuth(true_anomaly=foo[1], argument_of_periastron=argument_of_periastron)
            position.append([foo[0], azimuth, foo[1], phase])  # radius, azimuth, true anomaly, phase

            # x.append([r * np.cos(t), r * np.sin(t)])
            # y.append([foo[0] * np.cos(foo[1]), foo[0] * np.sin(foo[1])])
            # z.append([foo[0] * np.cos(azimuth), foo[0] * np.sin(azimuth)])

            # i += 1
            # ak sa fotometricke fazy rovnaju, teda chcem len jeden bod na krivke, tak to hned po prvom kroku breakne
            if from_photometric_phase == to_photometric_phase: break

        # import objects.Plot as Plt
        # Plt.plot_2d(points=x, grid=True)
        return position

        # < /BETAs>
