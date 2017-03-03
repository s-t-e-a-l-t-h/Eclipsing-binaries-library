#!/usr/bin/python

import pymysql
pymysql.install_as_MySQLdb()
import MySQLdb

import numpy as np
import globe.variables as gv
import objects.Function as Fn
import objects.Geometry as Geo
import objects.Atmosphere as Atm
from scipy import interpolate

class Planet:
    def __init__(
            self,
            mass=None,  # unit: solar mass //
            synchronicity_parameter=None,  # unit: [-]
            potential=None,  # unit: [-]
            effective_temperature=None,  # unit: [K]
            albedo=None,  # unit: [0 - 1]
            verbose=False,
            phi_steps=10,
            theta_steps=20
    ):
        self.verbose = verbose
        self.init = True
        self.exception = []

        self.phi_steps = phi_steps
        self.theta_steps = theta_steps

        self.effective_temperature = effective_temperature

        self.mass = float(mass)
        self.synchronicity_parameter = synchronicity_parameter
        self.potential = float(potential)  # bude sa pouzivat potencial \Omega = 1/r, takze ak bude nastavena sfera,
                                           # tak bude mat polomer 1/\Omega

        self.albedo = albedo
        self.backward_radius = None
        self.filling_factor = None
        self.critical_potential = None
        self.lagrangian_points = None

        self.polar_radius = None

        # <surface>
        self.faces = None
        self.vertices = None
        self.faces_orientation = None
        self.simplices = None
        # </surface>

        # <distribution>
        self.local_temperature = None
        self.local_radiation_power = None  # tok ziarenia z planet pri teplotach 300K je zanedbatelny, ale fajn, je to
        # tu pre uplnost (intenzita je na urovni 10e-40)
        # </distribution>


    def get_info(self, output=False):
        if output:
            info = ["Info:\n",
                    "<" + str(self.__class__.__name__) + "> -------------------------------------------------\n",
                    "mass:\t\t\t\t\t\t\t" + str(self.mass) + "\n",
                    "potential:\t\t\t\t\t\t" + str(self.potential) + "\n",
                    "effective temperature:\t\t\t" + str(self.effective_temperature) + "\n",
                    "albedo:\t\t\t\t\t\t\t" + str(self.albedo) + "\n",
                    "synchronicity parameter:\t\t" + str(self.synchronicity_parameter) + "\n",
                    "polar radius:\t\t\t\t\t" + str(self.polar_radius) + "\n",
                    "backward radius:\t\t\t\t" + str(self.backward_radius) + "\n",
                    "phi steps:\t\t\t\t\t\t" + str(self.phi_steps) + "\n",
                    "theta steps:\t\t\t\t\t" + str(self.theta_steps) + "\n"]

            return info
        else:
            print()
            print(Fn.color_string("info", "Info:"))
            print("<" + str(self.__class__.__name__) + "> -------------------------------------------------")
            print("mass:\t\t\t\t\t\t\t" + str(self.mass))
            print("potential:\t\t\t\t\t\t" + str(self.potential))
            print("effective temperature:\t\t\t" + str(self.effective_temperature))
            print("albedo:\t\t\t\t\t\t\t" + str(self.albedo))
            print("synchronicity parameter:\t\t" + str(self.synchronicity_parameter))
            print("polar radius:\t\t\t\t\t" + str(self.polar_radius))
            print("backward radius:\t\t\t\t" + str(self.backward_radius))
            print("phi steps:\t\t\t\t\t\t" + str(self.phi_steps))
            print("theta steps:\t\t\t\t\t" + str(self.theta_steps))
            print('</' + str(self.__class__.__name__) + '> -----------------------------------------------')
            print()

    def get_3d_model(self, phi_steps=None, theta_steps=None, radius=None, actual_distance=None):
        use = False
        z_point, rotation_angle, transform = None, None, None
        current_phi_steps, current_theta_steps = phi_steps, theta_steps
        theta, theta_step_length, partial, equatorial, meridional = np.pi / 2., np.pi / current_theta_steps, [], [], []
        point_coordinate = np.arange(3, dtype=np.float)

        for rot in range(0, int(current_theta_steps)):
            vector_spheric, vector_xyz = np.arange(3, dtype=np.float), np.arange(3, dtype=np.float)

            if theta <= np.pi:
                vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1., np.pi, theta
            elif theta > np.pi:
                vector_spheric[0], vector_spheric[1], vector_spheric[2] = 1.0, 0.0, (2.0 * np.pi) - theta

            # if theta <= np.pi:
            r = self.polar_radius * np.sin(theta - np.pi / 2.0)
            phi_points = 1 if theta == np.pi / 2.0 else np.ceil((r / self.polar_radius) * current_phi_steps)
            transform_steps = int(phi_points) if theta == np.pi / 2. else int(phi_points) + 1
            phi_step_length = (np.pi / 2.) / phi_points
            rotation_angle = phi_step_length

            for transform in range(0, transform_steps):
                solution = radius
                # risenie sa ulozi do vektora a pretransformuje prislusnou funkciou zo sferickcyh do kartezianskych suradnic
                point_coordinate[0], point_coordinate[1], point_coordinate[2] = solution, vector_spheric[1], \
                                                                                vector_spheric[2]
                xyz = Fn.spheric_to_cartesian(point_coordinate)
                xyz = [xyz[0] + actual_distance, xyz[1], xyz[2]]
                if transform == transform_steps - 1:
                    equatorial.append(xyz)
                elif transform == 0:
                    meridional.append(xyz)
                else:
                    partial.append(xyz)

                vector_xyz = Fn.spheric_to_cartesian(vector_spheric)
                rotate = Fn.rotate(vector=vector_xyz, angle=rotation_angle)
                vector_spheric = Fn.cartesian_to_spheric(rotate)
                if vector_spheric[1] < 0: vector_spheric[1] += (2. * np.pi)

            theta += theta_step_length

        f_point, equatorial = equatorial[0], equatorial[1:]
        # zero point (flip x coo of f_point)
        z_point = [abs(actual_distance - f_point[0]) + actual_distance, f_point[1], f_point[2]]

        full = []
        if not Fn.empty(partial) and not Fn.empty(equatorial) and not Fn.empty(meridional):
            full.append(f_point)
            full.append(z_point)
            for point in partial:
                full.append(point)
                full.append([point[0], -point[1], point[2]])
                full.append([point[0], point[1], -point[2]])
                full.append([point[0], -point[1], -point[2]])

            for point in equatorial:
                full.append(point)
                full.append([point[0], -point[1], point[2]])
            for point in meridional:
                full.append(point)
                full.append([point[0], point[1], -point[2]])
        else:
            if self.verbose:
                print(Fn.color_string(color="error",
                                      string="ValueError: ") + "In class: Star, function: get_3d_model(), line: " + str(
                    Fn.lineno()) + ". One of lists (`partial`, `equatorial`, `meridional`) is empty.")
            return False

        # # <kontrolne plotovanie>
        # import objects.Plot as Plt
        # Plt.plot_3d(faces=None, face_color="w", vertices=[full],
        #             normals_view=False, points_view=True, faces_view=False, verbose=self.verbose,
        #             face_alpha=1.0, azim=30, elev=30, save=False,
        #             x_range=[0.5, 1.5], y_range=[-0.5, 0.5], z_range=[-0.5,0.5])
        # # < /kontrolne plotovanie>

        return full

    def radiation_power(self,
            passband_model=None,
            passband_range=None,
            passband=None,  # len pre castelli-kurucz-04
            wavelength_step=10,  # len pre atmosferu black-body
            indices=None
    ):
        if self.verbose:
            print(Fn.color_string("info", "Info: ") + "Computing radiation power.")

        from scipy import integrate
        wavelength_step, wavelength = wavelength_step, passband_range[0]  # nm, nm
        x_values, y_values = [], []
        while wavelength <= passband_range[1]:
            y_values.append(
                Atm.planck(wavelength=wavelength, temperature=self.effective_temperature,
                           passband_fn=passband_model))
            x_values.append(wavelength * 1e-9)
            wavelength += wavelength_step
        flx = np.pi * integrate.simps(y_values, x_values)

        return np.array([flx] * len(self.local_temperature))


    # <SETTERs>
    def set_faces(self, faces=None):
        self.faces = faces

    def set_simplices(self, simplices=None):
        self.simplices = simplices

    def set_vertices(self, vertices=None):
        self.vertices = vertices

    def set_faces_orientation(self, faces_orientation=None):
        self.faces_orientation = faces_orientation

    def set_temperature_distribution(self, local_temperature=None):
        self.local_temperature = local_temperature

    def set_radiation_power(self, radiation_power=None):
        self.local_radiation_power = radiation_power
    # </SETTERs>

    # <GETTERs>
    def get_simplices(self):
        return self.simplices

    def get_faces_orientation(self):
        return self.faces_orientation

    def get_faces(self):
        return self.faces

    def get_vertices(self):
        return self.vertices

    def get_polar_radius(self):
        return self.polar_radius

    def get_radiation_power(self):
        return self.local_radiation_power

    def get_temperature_distribution(self):
        return self.local_temperature
    # </GETTERs>
