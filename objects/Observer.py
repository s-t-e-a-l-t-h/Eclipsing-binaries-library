#!/usr/bin/python
# import numpy as np
import pymysql

pymysql.install_as_MySQLdb()
import MySQLdb

import objects.Function as Fn
import numpy as np
# import math as m
import sys
import globe.variables as gv
import os
import objects.Lightcurve as Lc
import objects.Iostream as Io
from objects.Timer import Timer


class Observer:
    def __init__(
            self,
            passband=None,
            limb_darkening_model=None,
            observe=None,  # object to observe, e.g. binary star system creted as Binary class
            limb_darkening_interp_method="nearest",
            verbose=False
    ):
        # <default>
        self.exception = []
        self.passband_model = None
        self.passband_range = None
        self.init = True
        self.passband_list = ['bolometric', 'Generic/Bessell.U', 'Generic/Bessell.B', 'Generic/Bessell.V',
                              'Generic/Bessell.R', 'Generic/Bessell.I', 'SLOAN/SDSS.u',
                              'SLOAN/SDSS.g', 'SLOAN/SDSS.r', 'SLOAN/SDSS.i', 'SLOAN/SDSS.z', 'Generic/Stromgren.u',
                              'Generic/Stromgren.v', 'Generic/Stromgren.b',
                              'Generic/Stromgren.y']
        self.limb_darkening_interp_method = limb_darkening_interp_method
        # < /default>

        self.verbose = verbose
        self.limb_dakening_model = limb_darkening_model

        if passband in self.passband_list:
            self.passband = passband

            # IMPORTANT: toto je dolezite, tu sa do premennej self.passband_model za`se cela funkcia (cely pointer),
            # takze sa potom da volat normalne self.passband_model(var)
            self.set_passband_model()
            self.set_passband_range()
        else:
            if self.verbose:
                print(Fn.color_string("error", "ValueError: ") +
                      "In class Observer, function __init__(), line " + str(
                    Fn.lineno()) + ". Variable `passband` is invalid.")
            self.exception.append("ValueError: In class Observer, function __init__(), line " + str(
                Fn.lineno()) + ". Variable `passband` is invalid.")
            self.init = False

        if observe is not None:
            self.observe = observe
        else:
            if self.verbose:
                print(Fn.color_string("warning", "Warning: ") +
                      "In class Observer, function __init__(), line " + str(
                    Fn.lineno()) + ". Nothing to observer.")
            self.exception.append(
                "Warning: In class Observer, function __init__(), line " + str(Fn.lineno()) + ". Nothing to observer.")
            self.init = False

    @classmethod
    def limb_darkening_linear(cls, gamma, xlin):
        return 1.0 - (xlin * (1. - abs(np.cos(gamma))))

    @classmethod
    def limb_darkening_logarithmic(cls, gamma, xlog, ylog):
        return 1.0 - (xlog * (1.0 - abs(np.cos(gamma)))) - (ylog * abs(np.cos(gamma)) * np.log10(abs(np.cos(gamma))))

    @classmethod
    def limb_darkening_sqrt(cls, gamma, xsqrt, ysqrt):
        return 1.0 - (xsqrt * (1.0 - abs(np.cos(gamma)))) - (ysqrt * (1.0 - np.sqrt(abs(np.cos(gamma)))))

    @classmethod
    def limb_darkening_coefficients(
            cls,
            limb_darkeing_model=None,
            passband=None,
            metallicity=None,
            temperature=None,
            gravity=None,
            interpolation_method="nearest",
            verbose=False
    ):
        if verbose:
            print(Fn.color_string("info", "Info: ") + "Computing limb darkening coefficients.")
        # vrati kompletnu tabulku limb darkening koeficientov potrebnu pre interpolaciu
        ldc_table = np.array(
            cls.limb_darkening_table(verbose=verbose, passband=passband, limb_darkening_model=limb_darkeing_model))

        # body pre interpolaciu
        points, values = [], []
        if limb_darkeing_model == 'linear':
            for item in ldc_table:
                points.append([item[0], item[1], item[2]])  # [temperature, gravity, metallicity]
                values.append(item[3])  # [xlin]

        else:
            values = [[], []]
            for item in ldc_table:
                points.append(np.array([item[0], item[1], item[2]]))  # [temperature, gravity, metallicity]
                values[0].append(item[3])
                values[1].append(item[4])

        from scipy import interpolate

        # hodnoty v ktorych chceme interpolovat
        uvw = np.array(
            [np.array([temp, np.log10(grav), metallicity]) for temp, grav in list(zip(temperature, gravity))])

        if limb_darkeing_model == "linear":
            coefficients = interpolate.griddata(np.array(points), np.array(values), uvw, method=interpolation_method)

        else:
            x = interpolate.griddata(np.array(points), np.array(values[0]), uvw, method=interpolation_method)
            y = interpolate.griddata(np.array(points), np.array(values[1]), uvw, method=interpolation_method)
            coefficients = np.array(list(zip(*[x, y])))

        if len(coefficients) != len(uvw):
            if verbose:
                print(Fn.color_string("error",
                                      "Error: ") + "Error has been occured durring iterpolation. Length of input and output array is not same.")
            return False
        return coefficients

    @classmethod
    def limb_darkening_table(
            cls,
            limb_darkening_model=None,
            passband=None,
            verbose=False
    ):

        # passband translator ------------------------------------------------------------------------------------------
        # kedze neviem s istotou, ktory filter prepisat v tabulke limbdarkening ako Generic/Bessell.U,
        # ci Johnson_U alebo Bessell_UX alebo ktory, tak si tu pre to radsej spravim translator
        inp = np.array(
            ["Generic/Bessell.U", "Generic/Bessell.B", "Generic/Bessell.V", "Generic/Bessell.R", "Generic/Bessell.I"])
        outp = np.array(
            ["Johnson_U", "Johnson_B", "Johnson_V", "Johnson_R", "Johnson_I"])
        index = np.where(inp == passband)[0]
        if not Fn.empty(index):
            passband = outp[index[0]]
        # --------------------------------------------------------------------------------------------------------------

        # dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        # import sqlite3
        # conn = sqlite3.connect(dir_path + '/database/database.db')
        # c = conn.cursor()
        mysql_conn = MySQLdb.connect(host=gv.HOST,  # your host, usually localhost
                                     user=gv.USER,  # your username
                                     passwd=gv.PWD,  # your password
                                     db="elisa_assets")  # name of the data base
        c = mysql_conn.cursor()

        if limb_darkening_model is 'linear':
            # qry = 'SELECT temperature, gravity, metallicity, xlin FROM limbdarkening WHERE filter = "' + str \
            #     (passband) + '" COLLATE NOCASE'
            qry = 'SELECT temperature, gravity, metallicity, xlin FROM limbdarkening WHERE filter = "' + str \
                (passband) + '"'
        elif limb_darkening_model is 'logarithmic':
            # qry = 'SELECT temperature, gravity, metallicity, xlog, ylog FROM limbdarkening WHERE filter = "' + str \
            #     (passband) + '" COLLATE NOCASE'
            qry = 'SELECT temperature, gravity, metallicity, xlog, ylog FROM limbdarkening WHERE filter = "' + str \
                (passband) + '"'
        elif limb_darkening_model is 'sqrt':
            # qry = 'SELECT temperature, gravity, metallicity, xsqrt, ysqrt FROM limbdarkening WHERE filter = "' + str \
            #     (passband) + '" COLLATE NOCASE'
            qry = 'SELECT temperature, gravity, metallicity, xsqrt, ysqrt FROM limbdarkening WHERE filter = "' + str \
                (passband) + '"'
        else:
            if verbose:
                print(Fn.color_string("error", "ValueError: ") +
                      "In class: Observer, function: limb_darkening_table, line: " + str(Fn.lineno()) +
                      ". There is no such limb darkening model. Use `linear`, `logarithmic` or `sqrt`.")
            return False

        c.execute(qry)
        ret_val = np.array(c.fetchall()).tolist()

        if Fn.empty(ret_val):
            if verbose:
                print(Fn.color_string("error", "EmptyVariableError: ") +
                      "In class: Observer, function: limb_darkening_table, line: " + str(Fn.lineno()) +
                      ". Empty list, probably sqlite query problem.")
            mysql_conn.close()
            return False

        mysql_conn.close()
        return ret_val

    @classmethod
    def limb_darkening_factor(
            cls,
            faces_orientation=None,
            limb_darkening_model=None,
            limb_darkening_coefficients=None,
            verbose=False
    ):

        import objects.Geometry as Geo
        gamma = np.array(
            [Geo.angle(u=np.array([1.0, 0.0, 0.0]), v=normal, verbose=verbose) for normal in faces_orientation])
            
        # gamma = np.array([0.0] * len(faces_orientation))

        # tu je to osetrene tak, aby ked sa pouziva len jeden limb darkenink koeficient, tak aby sa nemusel upravovat
        # kod pre vypocet nizsie, tu sa proste zvacsi pole z jednotkoveho na potrebnu dlzku, nakopiruje sa
        if len(limb_darkening_coefficients) == 1:
            if type(limb_darkening_coefficients) == type(np.array([])):
                limb_darkening_coefficients = np.array(len(faces_orientation) * limb_darkening_coefficients.tolist())
            else:
                limb_darkening_coefficients = np.array(len(faces_orientation) * limb_darkening_coefficients)

        if len(limb_darkening_coefficients) != len(faces_orientation):
            if verbose:
                print(Fn.color_string("error",
                                      "Error: ") + "Length of variables `limb_darkeing_coefficients` and `faces_orientation` is not same.")
            return False

        ld_factor = None
        if limb_darkening_model == "linear":
            ld_factor = np.array(
                [cls.limb_darkening_linear(gamma=g, xlin=ldc) for g, ldc in
                 list(zip(gamma, limb_darkening_coefficients))])
        elif limb_darkening_model == "logarithmic":
            ld_factor = np.array(
                [cls.limb_darkening_logarithmic(gamma=g, xlog=ldc[0], ylog=ldc[1]) for g, ldc in
                 list(zip(gamma, limb_darkening_coefficients))])
        elif limb_darkening_model == "sqrt":
            ld_factor = np.array(
                [cls.limb_darkening_sqrt(gamma=g, xsqrt=ldc[0], ysqrt=ldc[1]) for g, ldc in
                 list(zip(gamma, limb_darkening_coefficients))])

        return [ld_factor, gamma]

    @classmethod
    def get_passband_table(
            cls,
            passband=None
    ):
        # dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        # import sqlite3
        # conn = sqlite3.connect(dir_path + '/database/database.db')

        mysql_conn = MySQLdb.connect(host=gv.HOST,  # your host, usually localhost
                                     user=gv.USER,  # your username
                                     passwd=gv.PWD,  # your password
                                     db="elisa_assets")  # name of the data base
        c = mysql_conn.cursor()
        order_by = ' ORDER BY wavelength ASC'
        c.execute('SELECT wavelength, pass FROM passband WHERE filter = "' + passband + '"' + order_by)
        ret_val = np.array(c.fetchall()).tolist()
        mysql_conn.close()
        return ret_val

    @staticmethod
    def passband_interpolation(
            passband=None
    ):
        passband = list(passband)
        from scipy import interpolate
        x, y = np.array(list(zip(*passband))[0]), np.array(list(zip(*passband))[1])
        return interpolate.Akima1DInterpolator(x, y)

    # medota bolometric stale vracia 1.0
    # vsetky ostatne passbandy vracaju interpolacnu funkciu s parametrom vlnovej dlzky, ktora sa potom pouziva pri
    # vypocte ziarenia; ak sa pocita bolometricky tok ziarenia, tak nie je kvazi nastaveny ziaden passband, tak aby to
    # nemuselo byt vo funkciach pre intenzitu ziarenia podmienkovane, tak do rovnic stale vstupuje passband funkcia, ale
    # v pripade bolometric je stale prenasobi len rovnicu 1.0, takze nic sa nedeje
    # NOTE: mozno to bude mierne spomalovat kod, pre bolometric, ale ten sa aj tak casto nepouzia, ale pri nastavenom
    # passband inom ako bolometric sa musi aj tak vyhodnocovat vypoctom, takze to spomaluje este viac, cize je to jedno
    @staticmethod
    def bolometric(w):
        if False: return w  # tento zapis je tu preto, aby nepyskoval PyCharm
        return 1.0

    def set_passband_range(self):
        if self.verbose: print(Fn.color_string("info", "Info: ") + "Setting passband range.")
        self.passband_range = self.passband_model_range(passband=self.passband)

    def set_passband_model(self):
        if self.verbose: print(Fn.color_string("info", "Info: ") + "Setting passband model.")

        table = self.get_passband_table(passband=self.passband)
        self.passband_model = self.passband_interpolation(passband=table)

        try:
            if self.passband == "bolometric":
                self.passband_model = self.bolometric
            else:
                table = self.get_passband_table(passband=self.passband)
                self.passband_model = self.passband_interpolation(passband=table)
            return True
        except:
            if self.verbose:
                print(Fn.color_string("error", "Error: ") +
                      "In class Observer, function set_passband_model(), line " + str(
                    Fn.lineno()) + ". Error has been occurred. Problem occurred probably during passband interpolation.")
            return False

    def get_passband_model(self):
        return self.passband_model

    @classmethod
    def passband_model_range(cls, passband=None):
        pb = cls.get_passband_table(passband=passband)
        pb = list(zip(*cls.get_passband_table(passband=passband)))[0] if passband != "bolometric" else [0.0, 10000.0]
        return [min(pb), max(pb)]

    def get_passband(self):
        return self.passband

    def get_passband_range(self):
        return self.passband_range

    def get_limb_darkening_model(self):
        return self.limb_dakening_model

    def set_limb_darkening_model(self, limb_darkening_model=None):
        self.limb_dakening_model = limb_darkening_model

    def compute_lightcurve(self, lightcurve_params=None, starmodel_params=None, limb_darkening_params=None,
                           postprocess_params=None, verbose=False):
        # bludny riadok
        if verbose: pass

        if Fn.empty(lightcurve_params):
            if self.verbose:
                print(Fn.color_string("error", "Error: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Variable `lightcurve_params` is empty.")
            self.exception.append("Error: In class Observer, function compute_lightcurve(), line " + str(
                Fn.lineno()) + ". Variable `lightcurve_params` is empty.")
            return False

        if Fn.empty(starmodel_params):
            if self.verbose:
                print(Fn.color_string("error", "Error: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Variable `starmodel_params` is empty.")
            self.exception.append("Error: In class Observer, function compute_lightcurve(), line " + str(
                Fn.lineno()) + ". Variable `starmodel_params` is empty.")
            return False

        if Fn.empty(postprocess_params):
            if self.verbose:
                print(Fn.color_string("error", "Error: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Variable `postprocess_params` is empty.")
            self.exception.append("Error: In class Observer, function compute_lightcurve(), line " + str(
                Fn.lineno()) + ". Variable `postprocess_params` is empty.")
            return False

        if not self.observe.init:
            return False

        # <import and definitions>
        import objects.Geometry as Geo
        import objects.Plot as Plt

        binary_system = self.observe
        primary = binary_system.primary
        secondary = binary_system.secondary
        orbit = binary_system.orbit
        lightcurve, norm_value = [], False

        # tu natvrdo zadefinujem spots na False, lebo ked to budem v buducnosti riesit, aby na to bol kod uz
        # opodmienkovany
        spots = True if None != primary.spots_meta or None != secondary.spots_meta else False
        # toto sa tu prekope tak, ze sa do kazdej prislusnej prida dali meta, a to t_object
        # lebo tak bol postaveny kod pre triangulaciu a nechcem vynucovat zadavanie
        # t_object pri inicializacii objektu Star, bolo by to redundantne pre uzivatela
        if not Fn.empty(primary.spots_meta):
            for i in range(0, len(primary.spots_meta)):
                primary.spots_meta[i]["t_object"] = "primary"

        if not Fn.empty(secondary.spots_meta):
            for i in range(0, len(secondary.spots_meta)):
                secondary.spots_meta[i]["t_object"] = "secondary"


        # premenna, ktora definuje, ci su skrvny na zlozkach aj po prepocte (teda napr, ak tam bola jedna
        # skvrna a nedalo sa ju zratat, tak uz tam nebude ziadna a nebude s nou treba uvazovat);
        # zapise sa to do tohto listu po tom, co sa skvrny vytvoria a ztrianguluju, defaultne su obe hodnoty na
        # False
        spots_validation = {"primary": False, "secondary": False}

        # v pripade, ze na hviezde nie su skvrny a kod zatial neuvazuje ziadne dynamicke efekty ako stacanie priamky
        # apsid a pod. tak pre zrychlenie kodu mozno stale ratat krivku len v rozsahu 0.0 az 0.5 fazy a potom ju
        # nainterpolovat
        # interpolated_lightcurve = True

        mirroring = False
        try:
            mirroring = postprocess_params["mirroring"]
        except:
            mirroring = True

        # < /import and definitions>

        # radius, azimuth, true anomaly, phase
        if orbit.get_eccentricity() > 0 or not mirroring:
            ffp, tfp = lightcurve_params["from_photometric_phase"], lightcurve_params["to_photometric_phase"]
        elif spots:
            ffp, tfp = 0.0, 1.0
        else:
            ffp, tfp = 0.0, 0.5

        # zbuffer
        zbuffer = postprocess_params["zbuffer"]

        orbital_motion = orbit.orbital_motion_beta(
            from_photometric_phase=ffp, to_photometric_phase=tfp, n=lightcurve_params["n"],
            adaptive_multiplicator=lightcurve_params["adaptive_multiplicator"],
            # multiplikacia pri zakrytoch v orbitalnom pohybe
            adaptive_orbit_part=lightcurve_params["adaptive_orbit_part"], eccentricity=orbit.eccentricity,
            argument_of_periastron=orbit.argument_of_periastron)

        # motion_to_plot = [[item[0] * np.cos(item[1]), item[0] * np.sin(item[1])] for item in orbital_motion]
        # Plt.plot_2d(points=motion_to_plot, grid=True)

        if binary_system.system == "eb":
            if orbit.eccentricity == 0.0 and primary.synchronicity_parameter == 1.0 and secondary.synchronicity_parameter == 1.0:
                #             zp = False if binary_system.binary_morph == "over-contact" else True

                if binary_system.binary_morph == "detached" or binary_system.binary_morph == "semi-detached":
                    if starmodel_params["homo"]:
                        model = \
                            binary_system.get_3d_model_optimized(t_object="both", actual_distance=1.0,
                                                                 critical_angle=np.pi / 2.0,
                                                                 phi_steps=primary.phi_steps,
                                                                 theta_steps=primary.theta_steps,
                                                                 zero_point=True, homo=True)

                        primary_model, secondary_model = model["primary"], model["secondary"]
                        del model

                    else:
                        primary_model = \
                            binary_system.get_3d_model_optimized(t_object="primary", actual_distance=1.0,
                                                                 critical_angle=np.pi / 2.0,
                                                                 phi_steps=primary.phi_steps,
                                                                 theta_steps=primary.theta_steps,
                                                                 zero_point=True, homo=False)["primary"]
                        secondary_model = \
                            binary_system.get_3d_model_optimized(t_object="secondary", actual_distance=1.0,
                                                                 critical_angle=np.pi / 2.0,
                                                                 phi_steps=secondary.phi_steps,
                                                                 theta_steps=secondary.theta_steps,
                                                                 zero_point=True, homo=False)["secondary"]

                    if spots:
                        primary_spots, secondary_spots = [], []
                        spots_metadata = {"primary": [], "secondary": []}
                        if None != primary.spots_meta:
                            primary_spots = binary_system.create_spots(meta=primary.spots_meta)["primary"]

                            for spot in primary_spots:
                                spots_metadata["primary"].append(spot["meta"])
                        if None != secondary.spots_meta:
                            secondary_spots = binary_system.create_spots(meta=secondary.spots_meta)["secondary"]

                            for spot in secondary_spots:
                                spots_metadata["secondary"].append(spot["meta"])


                        spots = {"primary": primary_spots, "secondary": secondary_spots}

                        norms = \
                            {"primary": Geo.normal_estimation(binary_object=binary_system,
                                                              actual_distance=orbit.get_periastron_distance(),
                                                              vertices=np.array(primary_model), t_object="primary",
                                                              mode="in_point", verbose=True),
                             "secondary": Geo.normal_estimation(binary_object=binary_system,
                                                                actual_distance=orbit.get_periastron_distance(),
                                                                vertices=np.array(secondary_model),
                                                                t_object="secondary", mode="in_point", verbose=True)
                             }
                        triangulation = \
                            Geo.trispot(vertices={"primary": primary_model, "secondary": secondary_model},
                                        norms=norms, spots=spots, binary_morph=binary_system.get_binary_morphology(),
                                        metadata=spots_metadata)

                        primary.set_vertices(vertices=np.array(triangulation[1]["primary"]))
                        secondary.set_vertices(vertices=np.array(triangulation[1]["secondary"]))

                        primary.set_simplices(simplices=triangulation[2]["primary"])
                        secondary.set_simplices(simplices=triangulation[2]["secondary"])

                        primary.set_faces(faces=primary.get_vertices()[primary.get_simplices()])
                        secondary.set_faces(faces=secondary.get_vertices()[secondary.get_simplices()])

                        spots_validation["primary"] = True if not Fn.empty(triangulation[4]["primary"]) else False
                        spots_validation["secondary"] = True if not Fn.empty(triangulation[4]["secondary"]) else False

                        simplex_map, spots_meta = triangulation[3], triangulation[4]


                        # # <kontrolne plotovanie primary>
                        # color, colors = ["r", "g", "b", "c", "y"], []
                        #
                        # for y in simplex_map["primary"]:
                        #     if simplex_map["primary"][y][1] == -1:
                        #         colors.append("w")
                        #     else:
                        #         colors.append(color[simplex_map["primary"][y][1]])
                        #
                        # Plt.plot_3d(vertices=None, faces=[primary.get_faces()], normals_view=False,
                        #             points_view=False, faces_view=True, face_color=[colors], point_size=10.0,
                        #             edge_color="k", elev=45, azim=45, save=False)
                        # # < /kontrolne plotovanie>
                        #
                        # # <kontrolne plotovanie secondary>
                        # color, colors = ["r", "g", "b", "c", "y"], []
                        #
                        # for y in simplex_map["secondary"]:
                        #     if simplex_map["secondary"][y][1] == -1:
                        #         colors.append("w")
                        #     else:
                        #         colors.append(color[simplex_map["secondary"][y][1]])
                        #
                        # Plt.plot_3d(vertices=None, faces=[secondary.get_faces()], normals_view=False,
                        #             points_view=False, faces_view=True, face_color=[colors], point_size=10.0,
                        #             edge_color="k", elev=45, azim=45, save=False)
                        # # < /kontrolne plotovanie>

                    else:
                        primary.set_vertices(vertices=primary_model)
                        secondary.set_vertices(vertices=secondary_model)

                        del primary_model, secondary_model

                        # # <kontrolne plotovanie>
                        # Plt.plot_3d(vertices=[primary.get_vertices(), secondary.get_vertices()], normals_view=False,
                        #                 points_view=True, faces_view=False, point_color="r", point_size=3.0,
                        #                 verbose=self.verbose)
                        # # < /kontrolne plotovanie>

                        # convex hull triangulation
                        triangulation_primary = Geo.convex_hull_triangulation(vertices=primary.get_vertices(),
                                                                              verbose=self.verbose)
                        triangulation_secondary = Geo.convex_hull_triangulation(vertices=secondary.get_vertices(),
                                                                                verbose=self.verbose)

                        primary.set_faces(faces=triangulation_primary[1])
                        secondary.set_faces(faces=triangulation_secondary[1])

                        primary.set_simplices(simplices=triangulation_primary[0])
                        secondary.set_simplices(simplices=triangulation_secondary[0])

                        del triangulation_primary, triangulation_secondary

                        # # <kontrolne plotovanie>
                        # Plt.plot_3d(faces=[primary.get_faces(), secondary.get_faces()], face_color="w",
                        #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                        #                 face_alpha=1.0, azim=30, elev=30, save=False,
                        #             filename="./out/" + str(primary.effective_temperature) + "_" + str(primary.potential) + "_" + str(primary.mass) + "___" +
                        #      str(secondary.effective_temperature) + "_" + str(secondary.potential) + "_" + str(secondary.mass) + "---" + str(orbit.orbital_period))
                        # exit()
                        # # < /kontrolne plotovanie>
                elif binary_system.binary_morph == "over-contact":
                    # doplnit triangulaciu pre over-contact
                    model = \
                        binary_system.get_3d_model_optimized(t_object="both", actual_distance=1.0,
                                                             critical_angle=np.pi / 4.0,
                                                             phi_steps=primary.phi_steps,
                                                             theta_steps=primary.theta_steps,
                                                             zero_point=False, homo=True)
                    system, x_separation = model["system"], model["separation"]
                    primary_model, secondary_model = model["primary"], model["secondary"]
                    del model

                    # # <kontrolne plotovanie>
                    # Plt.plot_3d(vertices=[system], normals_view=False,
                    #                 points_view=True, faces_view=False, point_color="r", point_size=3.0,
                    #                 verbose=self.verbose)
                    # # < /kontrolne plotovanie>

                    if spots:
                        primary_spots, secondary_spots = [], []
                        spots_metadata = {"primary": [], "secondary": []}
                        if None != primary.spots_meta:
                            primary_spots = binary_system.create_spots(meta=primary.spots_meta)["primary"]

                            for spot in primary_spots:
                                spots_metadata["primary"].append(spot["meta"])
                        if None != secondary.spots_meta:
                            secondary_spots = binary_system.create_spots(meta=secondary.spots_meta)["secondary"]

                            for spot in secondary_spots:
                                spots_metadata["secondary"].append(spot["meta"])

                        spots = {"primary": primary_spots, "secondary": secondary_spots}

                        norms = \
                            {"primary": Geo.normal_estimation(binary_object=binary_system,
                                                              actual_distance=orbit.get_periastron_distance(),
                                                              vertices=np.array(primary_model), t_object="primary",
                                                              mode="in_point", verbose=True),
                             "secondary": Geo.normal_estimation(binary_object=binary_system,
                                                                actual_distance=orbit.get_periastron_distance(),
                                                                vertices=np.array(secondary_model),
                                                                t_object="secondary", mode="in_point", verbose=True)
                             }

                        triangulation = \
                            Geo.trispot(vertices={"primary": primary_model, "secondary": secondary_model},
                                        norms=norms, spots=spots, binary_morph=binary_system.get_binary_morphology(),
                                        metadata=spots_metadata)


                        primary.set_vertices(vertices=np.array(triangulation[1]["primary"]))
                        secondary.set_vertices(vertices=np.array(triangulation[1]["secondary"]))

                        primary.set_simplices(simplices=triangulation[2]["primary"])
                        secondary.set_simplices(simplices=triangulation[2]["secondary"])

                        primary.set_faces(faces=primary.get_vertices()[primary.get_simplices()])
                        secondary.set_faces(faces=secondary.get_vertices()[secondary.get_simplices()])

                        spots_validation["primary"] = True if not Fn.empty(triangulation[4]["primary"]) else False
                        spots_validation["secondary"] = True if not Fn.empty(triangulation[4]["secondary"]) else False

                        simplex_map, spots_meta = triangulation[3], triangulation[4]

                        # # <kontrolne plotovanie primary>
                        # color, colors = ["r", "g", "b", "c", "y"], []
                        #
                        # for y in simplex_map["primary"]:
                        #     if simplex_map["primary"][y][1] == -1:
                        #         colors.append("w")
                        #     else:
                        #         colors.append(color[simplex_map["primary"][y][1]])
                        #
                        # Plt.plot_3d(vertices=None, faces=[primary.get_faces()], normals_view=False,
                        #             points_view=False, faces_view=True, face_color=[colors], point_size=10.0,
                        #             edge_color="k", elev=45, azim=45, save=False)
                        # # < /kontrolne plotovanie>

                        # # <kontrolne plotovanie secondary>
                        # color, colors = ["r", "g", "b", "c", "y"], []
                        #
                        # for y in simplex_map["secondary"]:
                        #     if simplex_map["secondary"][y][1] == -1:
                        #         colors.append("w")
                        #     else:
                        #         colors.append(color[simplex_map["secondary"][y][1]])
                        #
                        # Plt.plot_3d(vertices=None, faces=[secondary.get_faces()], normals_view=False,
                        #             points_view=False, faces_view=True, face_color=[colors], point_size=10.0,
                        #             edge_color="k", elev=45, azim=45, save=False)
                        # # < /kontrolne plotovanie>

                    else:
                        # normals necessary for triangulation
                        normal_vectors = \
                            Geo.normal_estimation(binary_object=binary_system, actual_distance=1.0,
                                                  vertices=np.array(system),
                                                  t_object="primary",
                                                  mode="in_point", verbose=False)

                        # Plt.plot_3d(normals=None, vertices=[system], faces=None, face_color="w", normals_view=False, points_view=True,
                        #                 faces_view=False, point_color="r", normal_color="w", point_size=3., verbose=True, face_alpha=1.,
                        #                 azim=30, elev=30)

                        # cgal triangulation
                        # print(os.path.dirname(os.path.realpath(__file__)))
                        # sys.exit()

                        triangulation = Geo.cgal_triangulation(normals=normal_vectors, points=system, verbose=False,
                                                               min_triangle_angle=primary.cgal["min_triangle_angle"],
                                                               max_triangle_size=primary.cgal["max_triangle_size"],
                                                               surface_aproximation_error=primary.cgal[
                                                                   "surface_aproximation_error"],
                                                               to_average_spacing=primary.cgal["to_average_spacing"])

                        # # <kontrolne plotovanie>
                        # Plt.plot_3d(faces=[triangulation[1]], face_color="w",
                        #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                        #                 face_alpha=1.0, azim=30, elev=30)
                        # # < /kontrolne plotovanie>

                        triangulation = Geo.cgal_separation(cgal_simplex=triangulation, x_separation=x_separation)

                        del (system, normal_vectors)

                        primary.set_vertices(vertices=triangulation["vertices"]["primary"])
                        secondary.set_vertices(vertices=triangulation["vertices"]["secondary"])

                        primary.set_faces(faces=triangulation["faces"]["primary"])
                        secondary.set_faces(faces=triangulation["faces"]["secondary"])

                        primary.set_simplices(simplices=triangulation["simplices"]["primary"])
                        secondary.set_simplices(simplices=triangulation["simplices"]["secondary"])

                faces_orientation_primary = Geo.face_orientation_beta(faces=primary.get_faces(),
                                                                      binary_object=binary_system, t_object="primary",
                                                                      actual_distance=1.0, verbose=self.verbose)

                faces_orientation_secondary = Geo.face_orientation_beta(faces=secondary.get_faces(),
                                                                        binary_object=binary_system,
                                                                        t_object="secondary",
                                                                        actual_distance=1.0, verbose=self.verbose)

                primary.set_faces_orientation(faces_orientation=faces_orientation_primary)
                secondary.set_faces_orientation(faces_orientation=faces_orientation_secondary)

                del faces_orientation_primary, faces_orientation_secondary

                # zapinanie a vypinanie plotovania normal
                if False:
                    # TATO CAST KODU JE POTREBNA LEN PRE PLOTOVANIE NORMAL
                    # trba ich zmenist a poposuvat na pozicie faziet
                    # normalizacia na 1
                    unit_vectors_primary = \
                        Geo.vector_array_normalisation(vector_arr=primary.get_faces_orientation(), multi=15.0,
                                                       verbose=self.verbose)
                    unit_vectors_secondary = \
                        Geo.vector_array_normalisation(vector_arr=secondary.get_faces_orientation(), multi=15.0,
                                                       verbose=self.verbose)

                    faces_com_primary = Geo.center_of_mass(faces=primary.get_faces(), verbose=self.verbose)
                    faces_com_secondary = Geo.center_of_mass(faces=secondary.get_faces(), verbose=self.verbose)

                    translation_primary = \
                        Geo.vector_array_translation(vector_arr=unit_vectors_primary, translation_arr=faces_com_primary,
                                                     verbose=self.verbose)
                    translation_secondary = \
                        Geo.vector_array_translation(vector_arr=unit_vectors_secondary,
                                                     translation_arr=faces_com_secondary,
                                                     verbose=self.verbose)

                    # farby pre normaly
                    c_primary, c_secondary = ["#0000ff"] * len(translation_primary), ["#000055"] * len(
                        translation_secondary)
                    Plt.plot_3d(normals=[translation_primary, translation_secondary],
                                vertices=[faces_com_primary, faces_com_secondary],
                                faces=[primary.get_faces(), secondary.get_faces()], face_color="w",
                                normals_view=True, points_view=False, faces_view=True,
                                point_color=[c_primary, c_secondary],
                                normal_color=[c_primary, c_secondary], point_size=3.0, verbose=True, face_alpha=1.0,
                                azim=30, elev=30)

                    # KONIEC CASTI PRE ZOBRAZOVANIE NORMAL

                gradnorm_primary = Geo.gradient_norm(faces=primary.get_faces(), verbose=self.verbose,
                                                     binary_object=binary_system, actual_distance=1.0,
                                                     t_object="primary")
                gradnorm_secondary = Geo.gradient_norm(faces=secondary.get_faces(), verbose=self.verbose,
                                                       binary_object=binary_system, actual_distance=1.0,
                                                       t_object="secondary")

                primary.set_gradient_norm(gradient_norm=gradnorm_primary)
                secondary.set_gradient_norm(gradient_norm=gradnorm_secondary)

                del gradnorm_primary, gradnorm_secondary

                primary.compute_gravity_distribution(gradnorm=primary.get_gradient_distribution())
                secondary.compute_gravity_distribution(gradnorm=secondary.get_gradient_distribution())

                # <kontrolne plotovanie gravity distribution primary>
                # rgb_p = Fn.arr_to_rainbow(arr=primary.get_gravity_distribution(),
                #                           minimum=min(primary.get_gravity_distribution()),
                #                           maximum=max(primary.get_gravity_distribution()))
                #
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                #                 points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)
                # < /kontrolne plotovanie gravity distribution primary>

                primary.gravity_darkening_factor_distribution()
                secondary.gravity_darkening_factor_distribution()

                primary.compute_polar_temperature()
                secondary.compute_polar_temperature()

                if spots_validation["primary"]:
                    primary.compute_temperature_distribution(simplex_map=simplex_map["primary"],
                                                             spots_meta=spots_meta["primary"])
                else:
                    primary.compute_temperature_distribution()

                if spots_validation["secondary"]:
                    secondary.compute_temperature_distribution(simplex_map=simplex_map["secondary"],
                                                             spots_meta=spots_meta["secondary"])
                else:
                    secondary.compute_temperature_distribution()

                # # <kontrolne plotovanie temperature distribution primary>
                # rgb_p = Fn.arr_to_rainbow(arr=primary.get_temperature_distribution(),
                #                        minimum=min(primary.get_temperature_distribution()),
                #                        maximum=max(primary.get_temperature_distribution()))
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                #             points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)
                # # < /kontrolne plotovanie temperature distribution primary>

                # # < kontrolne plotovanie temperature distribution for wuma>
                # faces = np.concatenate((primary.get_faces(), secondary.get_faces()), 0)
                # temperature = np.concatenate((primary.get_temperature_distribution(),
                #                               secondary.get_temperature_distribution()), 0)
                # rgb_p = Fn.arr_to_rainbow(arr=temperature,
                #                           minimum=min(temperature),
                #                           maximum=max(temperature))
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                #
                # Plt.plot_3d(faces=[faces], face_color=[hex_p], normals_view=False,
                #         points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)
                # # < /kontrolne plotovanie temperature distribution for wuma>

                r_pwr = [primary.radiation_power(passband_model=self.get_passband_model(), passband=self.passband,
                                                 passband_range=self.get_passband_range(), wavelength_step=10.0),
                         secondary.radiation_power(passband_model=self.get_passband_model(), passband=self.passband,
                                                   passband_range=self.get_passband_range(), wavelength_step=10.0)
                         ]

                if Fn.empty(r_pwr[0]) or Fn.empty(r_pwr[1]):
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="ValueError: ") + "In class: Observer, function: compute_lightcurve(), line: " + str(
                            Fn.lineno()) + ". Invalid value (`boolean`) encountered during flux computing.")
                    self.exception.append(
                        "ValueError: In class: Observer, function: compute_lightcurve(), line: " + str(
                            Fn.lineno()) + ". Invalid value (`boolean`) encountered during flux computing.")
                    return False

                primary.set_radiation_power(radiation_power=r_pwr[0])
                secondary.set_radiation_power(radiation_power=r_pwr[1])
                del r_pwr

                # import time
                # import sys
                # folder = str(time.time())
                # info = [primary.get_info(output=True), "\n", secondary.get_info(output=True), "\n", orbit.get_info(output=True)]
                # Io.save_csv(filename="info", datafolder="./out/" + folder, data=info, delim="")
                #
                # # <kontrolne plotovanie>
                # Plt.plot_3d(faces=[primary.get_faces(), secondary.get_faces()], face_color="w",
                #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                #                 face_alpha=1.0, azim=30, elev=30, save=True,
                #             filename="./out/" + folder + "/" + "img", dpi=100)
                # # < /kontrolne plotovanie>


                # <kontrolne plotovanie radiation power distribution>
                # rgb_p = Fn.arr_to_rainbow(arr=primary.get_radiation_power(),
                #                        minimum=min(primary.get_radiation_power()),
                #                        maximum=max(primary.get_radiation_power()))
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                #                 points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)

                # rgb_s = Fn.arr_to_rainbow(arr=secondary.get_radiation_power(),
                #                           minimum=min(secondary.get_radiation_power()),
                #                           maximum=max(secondary.get_radiation_power()))
                # hex_s = Fn.rgb_to_hex(color=rgb_s, sharp=True)
                # Plt.plot_3d(faces=[secondary.get_faces()], face_color=[hex_s], normals_view=False,
                #                 points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30,
                #             x_range=[0, 2], y_range=[-1, 1], z_range=[-1, 1])
                # < /kontrolne plotovanie radiation power distribution>


                if not Fn.empty(limb_darkening_params, debug=False):
                    ld_temperature = [[limb_darkening_params["temperature_primary"]],
                                      [limb_darkening_params["temperature_secondary"]]]

                    ld_gravity = [[limb_darkening_params["gravity_primary"]],
                                  [limb_darkening_params["gravity_secondary"]]]

                    ld_metallicity = [limb_darkening_params["metallicity_primary"],
                                      limb_darkening_params["metallicity_secondary"]]
                else:
                    ld_temperature = [primary.get_temperature_distribution(), secondary.get_temperature_distribution()]
                    ld_gravity = [primary.get_gravity_distribution(), secondary.get_gravity_distribution()]
                    ld_metallicity = [primary.get_metallicity(), secondary.get_metallicity()]

                ldc_primary = \
                    self.limb_darkening_coefficients(limb_darkeing_model=self.get_limb_darkening_model(),
                                                     passband=self.get_passband(), metallicity=ld_metallicity[0],
                                                     temperature=ld_temperature[0], gravity=ld_gravity[0],
                                                     interpolation_method=self.limb_darkening_interp_method,
                                                     verbose=self.verbose)
                # <kontrolne plotovanie koeficientov okrajoveho stemnenia>
                # rgb_p = Fn.arr_to_rainbow(arr=ldc_primary,
                #                        minimum=min(ldc_primary),
                #                        maximum=max(ldc_primary))
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                #                 points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)
                # < /kontrolne plotovanie koeficientov okrajoveho stemnenia>

                ldc_secondary = \
                    self.limb_darkening_coefficients(limb_darkeing_model=self.get_limb_darkening_model(),
                                                     passband=self.get_passband(), metallicity=ld_metallicity[1],
                                                     temperature=ld_temperature[1], gravity=ld_gravity[1],
                                                     interpolation_method=self.limb_darkening_interp_method,
                                                     verbose=self.verbose)

                # <kontrolne plotovanie>
                # Plt.plot_3d(faces=[def_position[0][0], def_position[0][1]], face_color="w",
                #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                #                 face_alpha=1.0, azim=0, elev=0)
                # < /kontrolne plotovanie>

                percentage = 0.0
                percentage_step = 100.0 / len(orbital_motion)
                if self.verbose:
                    sys.stdout.write(Fn.color_string("info", "Info: ") + "Making lightcurve... \n")
                    sys.stdout.write("\r\t\t%d%% done\t" % percentage)
                    sys.stdout.flush()

                for orbital_position in orbital_motion:

                    if self.verbose:
                        sys.stdout.write("\r\t\t%d%% done\t" % percentage)
                        sys.stdout.flush()
                        percentage += percentage_step

                    # orbital_position[1] = np.pi - (np.pi / 10)

                    act_position = \
                        orbit.rotate_system(faces=[primary.get_faces(), secondary.get_faces()],
                                            normals=[primary.get_faces_orientation(),
                                                     secondary.get_faces_orientation()],
                                            vertices=[primary.get_vertices(), secondary.get_vertices()],
                                            rotation_angle=orbital_position[1], inclination_rotation=True,
                                            faces_rotation=True, inclination=orbit.get_inclination(), verbose=False)

                    # # <kontrolne plotovanie>
                    # Plt.plot_3d(faces=[act_position[0][0], act_position[0][1]], face_color="w",
                    #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                    #                 face_alpha=1.0, azim=0, elev=0)
                    # exit()
                    # # < /kontrolne plotovanie>


                    # darkside_filter() vrati hodnoty v tvare
                    #  [faces[primary, secondary], normals[primary, secondary], indices[primary, secondary]]
                    # co sa tyka indexov, jedna sa o cislovanie z povodneho pola triangulacie
                    darksite_filter = \
                        Geo.darkside_filter(faces=[act_position[0][0], act_position[0][1]],
                                            normals=[act_position[1][0], act_position[1][1]],
                                            verbose=False)
                    # # <kontrolne plotovanie>
                    # Plt.plot_3d(faces=[darksite_filter[0][0], darksite_filter[0][1]], face_color="w",
                    #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                    #                 face_alpha=1.0, azim=0, elev=0)
                    # # < /kontrolne plotovanie>

                    # eclipse_filter() vrati hodnoty v tvare
                    # [idx of visible faces [primary, secondary], surface are of those faces [primary, secondary]]
                    #
                    eclipse_filter = \
                        Geo.eclipse_filter(indices=[darksite_filter[2][0], darksite_filter[2][1]],
                                           vertices=[act_position[2][0], act_position[2][1]],
                                           simplices=[Fn.array_mask(primary.get_simplices(), darksite_filter[2][0]),
                                                      Fn.array_mask(secondary.get_simplices(), darksite_filter[2][1])],
                                           orbital_angle=orbital_position[1], verbose=False, zbuffer=zbuffer,
                                           resolution=1000)

                    # ak nie je viditelna ziadna zlozka, tak plocha sa nastavi na False, potom sa to dalej vyuzije
                    surface_primary = False if len(eclipse_filter[0][0]) == 0 else eclipse_filter[1][0]
                    surface_secondary = False if len(eclipse_filter[0][1]) == 0 else eclipse_filter[1][1]
                    surface = [surface_primary, surface_secondary]

                    # <kontrolne plotovanie>
                    # Plt.plot_3d(faces=[eclipse_filter[2][0], eclipse_filter[2][1]], face_color="w",
                    #                 normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                    #                 face_alpha=1.0, azim=0, elev=0)
                    # < /kontrolne plotovanie>

                    # faktor okrajoveho stemnenia pre vyfiltrovane elementy
                    # struktura, ktoru vrati limb_darkening_factor() je nasledujuca
                    # [faktor okrajoveho stemnenia pre dany element pri danej pozicii, uhol medzi normalou elementu
                    #  a vektorom smerom k pozorovatelovi teda vektorom (1,0,0)]
                    ldf_primary, ldf_secondary = 0.0, 0.0
                    if type(surface[0]) is not type(True):
                        ldf_primary = \
                            self.limb_darkening_factor(faces_orientation=Fn.array_mask(act_position[1][0],
                                                                                       eclipse_filter[0][0]),
                                                       limb_darkening_model=self.get_limb_darkening_model(),
                                                       limb_darkening_coefficients=Fn.array_mask(array=ldc_primary,
                                                                                                 mask=eclipse_filter[0][
                                                                                                     0]))

                    if type(surface[1]) is not type(True):
                        ldf_secondary = \
                            self.limb_darkening_factor(faces_orientation=Fn.array_mask(act_position[1][1],
                                                                                       eclipse_filter[0][1]),
                                                       limb_darkening_model=self.get_limb_darkening_model(),
                                                       limb_darkening_coefficients=Fn.array_mask(array=ldc_secondary,
                                                                                                 mask=eclipse_filter[0][
                                                                                                     1]))

                    flux = \
                        [0.0 if type(surface[0]) is type(True)
                         else self.compute_flux(e_limb_darkening_factor=ldf_primary, e_surface=surface[0],
                                                e_flux=Fn.array_mask(array=primary.get_radiation_power(),
                                                                     mask=eclipse_filter[0][0])),
                         0.0 if type(surface[1]) is type(True)
                         else self.compute_flux(e_limb_darkening_factor=ldf_secondary, e_surface=surface[1],
                                                e_flux=Fn.array_mask(array=secondary.get_radiation_power(),
                                                                     mask=eclipse_filter[0][1]))]

                    lightcurve.append([orbital_position[3], (flux[0] + flux[1]) * (
                        (orbit.get_relative_semimajor_axis() * gv.SOLAR_RADIUS) ** 2)])

                if self.verbose:
                    sys.stdout.write("\r\t\t100% done")
                    sys.stdout.write("\n")

        elif binary_system.system == "te":
            if binary_system.planet == "roche":
                pass

            elif binary_system.planet == "sphere":
                primary.set_vertices(vertices=primary.get_3d_model(phi_steps=primary.phi_steps,
                                                                   theta_steps=primary.theta_steps))

                secondary.set_vertices(vertices=secondary.get_3d_model(phi_steps=secondary.phi_steps,
                                                                       theta_steps=secondary.theta_steps,
                                                                       radius=secondary.get_polar_radius(),
                                                                       actual_distance=orbit.get_periastron_distance()))

                # <kontrolne plotovanie>
                Plt.plot_3d(vertices=[primary.get_vertices(), secondary.get_vertices()], normals_view=False,
                                points_view=True, faces_view=False, point_color="r", point_size=3.0,
                                verbose=self.verbose)
                # < /kontrolne plotovanie>

                # convex hull triangulation
                triangulation_primary = Geo.convex_hull_triangulation(vertices=primary.get_vertices(),
                                                                      verbose=self.verbose)
                triangulation_secondary = Geo.convex_hull_triangulation(vertices=secondary.get_vertices(),
                                                                        verbose=self.verbose)
                primary.set_faces(faces=triangulation_primary[1])
                secondary.set_faces(faces=triangulation_secondary[1])
                primary.set_simplices(simplices=triangulation_primary[0])
                secondary.set_simplices(simplices=triangulation_secondary[0])
                del triangulation_primary, triangulation_secondary
                faces_orientation_primary = Geo.face_orientation_beta(faces=primary.get_faces(),
                                                                      binary_object=binary_system, t_object="primary",
                                                                      actual_distance=orbit.periastron_distance,
                                                                      verbose=self.verbose)
                faces_orientation_secondary = Geo.face_orientation_beta(faces=secondary.get_faces(),
                                                                        binary_object=binary_system,
                                                                        t_object="secondary",
                                                                        actual_distance=orbit.periastron_distance,
                                                                        verbose=self.verbose)
                primary.set_faces_orientation(faces_orientation=faces_orientation_primary)
                secondary.set_faces_orientation(faces_orientation=faces_orientation_secondary)
                del faces_orientation_primary, faces_orientation_secondary
                # vzidalenost tu nie je podstatna
                gradnorm_primary = Geo.gradient_norm(faces=primary.get_faces(), verbose=self.verbose,
                                                     binary_object=binary_system, actual_distance=None,
                                                     t_object="primary")
                primary.set_gradient_norm(gradient_norm=gradnorm_primary)
                del gradnorm_primary
                primary.compute_gravity_distribution(gradnorm=primary.get_gradient_distribution())

                primary.compute_gravity_distribution(gradnorm=primary.get_gradient_distribution())
                primary.gravity_darkening_factor_distribution()
                primary.compute_polar_temperature()
                primary.compute_temperature_distribution()

                primary.set_temperature_distribution([primary.effective_temperature] * len(primary.get_simplices()))
                primary.local_gravity = [primary.polar_gravity] * len(primary.get_simplices())


                # # <kontrolne plotovanie temperature distribution primary>
                # rgb_p = Fn.arr_to_rainbow(arr=primary.get_temperature_distribution(),
                #                        minimum=min(primary.get_temperature_distribution()),
                #                        maximum=max(primary.get_temperature_distribution()))
                # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                #             points_view=False, faces_view=True, verbose=True, face_alpha=1., azim=30, elev=30)
                # # < /kontrolne plotovanie temperature distribution primary>

                r_pwr = [primary.radiation_power(passband_model=self.get_passband_model(), passband=self.passband,
                                                 passband_range=self.get_passband_range(), wavelength_step=10.0),
                         None]

                if Fn.empty(r_pwr[0]):
                    if self.verbose:
                        print(Fn.color_string(color="error",
                                              string="ValueError: ") + "In class: Observer, function: compute_lightcurve(), line: " + str(
                            Fn.lineno()) + ". Invalid value (`boolean`) encountered during flux computing.")
                    self.exception.append(
                        "ValueError: In class: Observer, function: compute_lightcurve(), line: " + str(
                            Fn.lineno()) + ". Invalid value (`boolean`) encountered during flux computing.")
                    return False

                primary.set_radiation_power(radiation_power=r_pwr[0])
                del r_pwr
                if not Fn.empty(limb_darkening_params, debug=False):
                    ld_temperature = [[limb_darkening_params["temperature_primary"]]]
                    ld_gravity = [[limb_darkening_params["gravity_primary"]]]
                    ld_metallicity = [limb_darkening_params["metallicity_primary"]]
                else:
                    ld_temperature = [primary.get_temperature_distribution()]
                    ld_gravity = [primary.get_gravity_distribution()]
                    ld_metallicity = [primary.get_metallicity()]

                ldc_primary = \
                    self.limb_darkening_coefficients(limb_darkeing_model=self.get_limb_darkening_model(),
                                                     passband=self.get_passband(), metallicity=ld_metallicity[0],
                                                     temperature=ld_temperature[0], gravity=ld_gravity[0],
                                                     interpolation_method=self.limb_darkening_interp_method,
                                                     verbose=self.verbose)

                percentage = 0.0
                percentage_step = 100.0 / len(orbital_motion)
                if self.verbose:
                    sys.stdout.write(Fn.color_string("info", "Info: ") + "Making lightcurve... \n")
                    sys.stdout.write("\r\t\t%d%% done\t" % percentage)
                    sys.stdout.flush()
                for orbital_position in orbital_motion:
                    delta_distance = orbital_position[0] - orbit.get_periastron_distance()
                    delta = np.array([delta_distance, 0.0, 0.0])
                    secondary_vertices = secondary.get_vertices() + delta
                    secondary_faces = secondary_vertices[secondary.get_simplices()]

                    if self.verbose:
                        sys.stdout.write("\r\t\t%d%% done\t" % percentage)
                        sys.stdout.flush()
                        percentage += percentage_step

                    # orbital_position[1] = np.pi - (np.pi / 10)

                    act_position = \
                        orbit.rotate_system(faces=[primary.get_faces(), secondary_faces],
                                            normals=[primary.get_faces_orientation(),
                                                     secondary.get_faces_orientation()],
                                            vertices=[primary.get_vertices(), secondary_vertices],
                                            rotation_angle=orbital_position[1], inclination_rotation=True,
                                            faces_rotation=True, inclination=orbit.get_inclination(), verbose=False)



                    # # <kontrolne plotovanie temperature distribution primary>
                    # rgb_p = Fn.arr_to_rainbow(arr=primary.get_temperature_distribution(),
                    #                           minimum=min(primary.get_temperature_distribution()),
                    #                           maximum=max(primary.get_temperature_distribution()))
                    # hex_p = Fn.rgb_to_hex(color=rgb_p, sharp=True)
                    # Plt.plot_3d(faces=[primary.get_faces()], face_color=[hex_p], normals_view=False,
                    #             points_view=False, faces_view=True, verbose=True, face_alpha=1.,
                    #             azim=-np.degrees(orbital_position[1]), elev=0)
                    # # < /kontrolne plotovanie temperature distribution primary>


                    # # <kontrolne plotovanie>
                    # Plt.plot_3d(faces=[act_position[0][0], act_position[0][1]], face_color="w",
                    #             normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                    #             face_alpha=1.0, azim=-90, elev=90)
                    # # < /kontrolne plotovanie>

                    # from time import time
                    # ddp = [[v[0], v[1]]for v in act_position[2][0]]
                    # dds = [[v[0], v[1]] for v in act_position[2][1]]
                    # v = np.concatenate((ddp, dds), 0)
                    #
                    # from matplotlib import pyplot as plt
                    # plt.scatter(list(zip(*v))[0], list(zip(*v))[1], s=0.5)
                    # plt.axis("equal")
                    # plt.xlim([-2.5, 2.5])
                    # plt.ylim([-2.5, 2.5])
                    # plt.scatter(orbital_position[0] * np.cos(orbital_position[1]),
                    #             orbital_position[0] * np.sin(orbital_position[1]),
                    #                                      s=20, c="r")
                    #
                    # plt.grid(True)
                    # plt.savefig(str(time()) + ".png")
                    # continue

                    # darkside_filter() vrati hodnoty v tvare
                    #  [faces[primary, secondary], normals[primary, secondary], indices[primary, secondary]]
                    # co sa tyka indexov, jedna sa o cislovanie z povodneho pola triangulacie
                    darksite_filter = \
                        Geo.darkside_filter(faces=[act_position[0][0], act_position[0][1]],
                                            normals=[act_position[1][0], act_position[1][1]],
                                            verbose=False)
                    # # <kontrolne plotovanie>
                    # Plt.plot_3d(faces=[darksite_filter[0][0], darksite_filter[0][1]], face_color="w",
                    #             normals_view=False, points_view=False, faces_view=True, verbose=self.verbose,
                    #             face_alpha=1.0, azim=0, elev=0)
                    # # < /kontrolne plotovanie>

                    # eclipse_filter() vrati hodnoty v tvare
                    # [idx of visible faces [primary, secondary], surface are of those faces [primary, secondary]]

                    eclipse_filter = \
                        Geo.eclipse_filter(indices=[darksite_filter[2][0], darksite_filter[2][1]],
                                           vertices=[act_position[2][0], act_position[2][1]],
                                           simplices=[Fn.array_mask(primary.get_simplices(), darksite_filter[2][0]),
                                                      Fn.array_mask(secondary.get_simplices(), darksite_filter[2][1])],
                                           orbital_angle=orbital_position[1], verbose=False, zbuffer=zbuffer,
                                           resolution=500)


                    # eclipse_filter = [[[], []], [[], []], [[], []]]
                    # eclipse_filter[0][0] = darksite_filter[2][0]
                    # eclipse_filter[1][0] = Geo.triangle_surface_area(np.array(act_position[0][0])[darksite_filter[2][0]])
                    # eclipse_filter[2][0] = np.array(act_position[0][0])[darksite_filter[2][0]]



                    # ak nie je viditelna ziadna zlozka, tak plocha sa nastavi na False, potom sa to dalej vyuzije
                    surface_primary = False if len(eclipse_filter[0][0]) == 0 else eclipse_filter[1][0]
                    surface = [surface_primary, None]

                    ldf_primary, ldf_secondary = 0.0, 0.0
                    if type(surface[0]) is not type(True):
                        ldf_primary = \
                            self.limb_darkening_factor(faces_orientation=Fn.array_mask(act_position[1][0],
                                                                                       eclipse_filter[0][0]),
                                                       limb_darkening_model=self.get_limb_darkening_model(),
                                                       limb_darkening_coefficients=Fn.array_mask(array=ldc_primary,
                                                                                                 mask=eclipse_filter[0][
                                                                                                     0]))
                    flux = \
                        np.array([0.0 if type(surface[0]) is type(True)
                         else self.compute_flux(e_limb_darkening_factor=ldf_primary, e_surface=surface[0],
                                                e_flux=Fn.array_mask(array=primary.get_radiation_power(),
                                                                     mask=eclipse_filter[0][0])), None], dtype="float64")

                    lightcurve.append([orbital_position[3], (flux[0]) * (
                        (orbit.get_relative_semimajor_axis() * gv.SOLAR_RADIUS) ** 2)])

                    print(lightcurve[-1][1])

                if self.verbose:
                    sys.stdout.write("\r\t\t100% done")
                    sys.stdout.write("\n")

        try:
            # <interpolacia krivky>
            if mirroring and not spots:
                norm_value = Lc.akima_interpolation(from_photometric_phase=0.0,
                                                    to_photometric_phase=0.5,
                                                    lightcurve=lightcurve, mirror=True, spots=False)

                norm_value = Lc.akima_interpolation(from_photometric_phase=0.25, to_photometric_phase=0.25,
                                                    lightcurve=norm_value, mirror=False, spots=False)

                interp_start, interp_stop = lightcurve_params["from_photometric_phase"], lightcurve_params[
                    "to_photometric_phase"]

                lightcurve = Lc.akima_interpolation(from_photometric_phase=interp_start,
                                                    to_photometric_phase=interp_stop,
                                                    lightcurve=lightcurve, mirror=True)

            elif mirroring and spots:
                norm_value = Lc.akima_interpolation(from_photometric_phase=0.0,
                                                    to_photometric_phase=1.0,
                                                    lightcurve=lightcurve, mirror=True, spots=True)

                norm_value = Lc.akima_interpolation(from_photometric_phase=0.25, to_photometric_phase=0.25,
                                                    lightcurve=norm_value, mirror=False, spots=True)

                interp_start, interp_stop = lightcurve_params["from_photometric_phase"], lightcurve_params[
                    "to_photometric_phase"]

                lightcurve = Lc.akima_interpolation(from_photometric_phase=interp_start,
                                                    to_photometric_phase=interp_stop,
                                                    lightcurve=lightcurve, mirror=True, spots=True)

            # < /interpolacia krivky>
        except:
            if self.verbose:
                print(Fn.color_string("error", "ValueError: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Error has been occured during lightcruve postprocessing")
                print("Primary info:")
                primary.get_info()
                print("Secondary info:")
                secondary.get_info()
                print("Orbit info:")
                orbit.get_info()
            return False

        # Plt.plot_2d(points=lightcurve, aspect="auto", save=True, line=True,
        #             filename="./out/" + folder + "/" + "lcurve_img_raw")
        # Io.save_csv(filename="lc_raw", datafolder="./out/" + folder, data=lightcurve, delim="\t", rewrite=False)

        try:
            # <vyhladenie minim>
            fix_minima = postprocess_params["fix_minima"]
            if fix_minima:
                eclipses = Lc.photometric_phase_of_eclipses(pp_primary=orbit.conjuction[0]["true_phase"],
                                                            pp_secondary=orbit.conjuction[1]["true_phase"],
                                                            lightcurve=lightcurve)
                lightcurve = Lc.fix_minima(lightcurve=lightcurve, eclipses=eclipses)
                # </vyhladenie minim>
        except:
            if self.verbose:
                print(Fn.color_string("error", "ValueError: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Error has been occured during fixing minima process")
                print("Primary info:")
                primary.get_info()
                print("Secondary info:")
                secondary.get_info()
                print("Orbit info:")
                orbit.get_info()
            return False

        # Plt.plot_2d(points=lightcurve, aspect="auto", save=True, line=True,
        #             filename="./out/" + folder + "/" + "lcurve_img_fixed")
        # Io.save_csv(filename="lc_fixed", datafolder="./out/" + folder, data=lightcurve, delim="\t", rewrite=False)

        try:
            # <smoothing krivky>
            smoothed_lightcurve = postprocess_params["gaussian_smooth"]
            if smoothed_lightcurve:
                lightcurve = Lc.gaussian_smooth(lightcurve=lightcurve)
                # < /smoothing krivky>
        except:
            if self.verbose:
                print(Fn.color_string("error", "ValueError: ") +
                      "In class Observer, function compute_lightcurve(), line " + str(
                    Fn.lineno()) + ". Error has been occured during smoothing process")
                print("Primary info:")
                primary.get_info()
                print("Secondary info:")
                secondary.get_info()
                print("Orbit info:")
                orbit.get_info()
            return False

        # v pripade, ze doslo k fixovaniu minim, je potrebne nainterpolovat na novo krivku, aby mala potrebny pocet
        # bodov pre machine learning
        if fix_minima:
            interp_start, interp_stop = lightcurve_params["from_photometric_phase"], lightcurve_params[
                "to_photometric_phase"]
            lightcurve = Lc.akima_interpolation(from_photometric_phase=interp_start,
                                                to_photometric_phase=interp_stop,
                                                lightcurve=lightcurve, mirror=False)


        Plt.plot_2d(points=lightcurve, aspect="auto", save=False, line=True)
        if not Fn.empty(norm_value):
            return [lightcurve, norm_value[0][1]]
        else:
            return [lightcurve, 1.0]

    @classmethod
    def compute_flux(cls, e_limb_darkening_factor=None, e_surface=None, e_flux=None):
        # flux = 0.0
        # for ldf, gamma, s, f in list(zip(e_limb_darkening_factor[0], e_limb_darkening_factor[1], e_surface, e_flux)):
        #     flux += ldf * np.cos(gamma) * f * s
        flux = np.sum(np.array(e_limb_darkening_factor[0]) * np.array(e_surface) * np.array(e_flux) * np.cos(
            np.array(e_limb_darkening_factor[1])))
        return flux

    def get_exception(self):
        return self.exception
