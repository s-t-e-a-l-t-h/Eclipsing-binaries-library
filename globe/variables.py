#!/usr/bin/python
from scipy import constants

# colors

# global COLOR_WARNING
COLOR_WARNING = '\033[93m'

# global COLOR_ERROR
COLOR_ERROR = '\033[91m'

# global COLOR_GREEN
COLOR_GREEN = '\033[92m'

# global COLOR_BLUE
COLOR_BLUE = '\033[94m'

# global COLOR_END
COLOR_END = '\033[0m'

'''
Cyan = '\033[96m'
White = '\033[97m'
Yellow = '\033[93m'
Magenta = '\033[95m'
Grey = '\033[90m'
Black = '\033[90m'
Default = '\033[99m'
'''

# /colors


# physical constants

# global SOLAR_RADIUS
SOLAR_RADIUS = 6.955E8  # m

# global SOLAR_MASS
SOLAR_MASS = 1.9891E30  # kg

# global G_CONSTANT
G_CONSTANT = constants.G
GRAVITATIONAL_CONSTANT = constants.G

# global BOLTZMAN_CONSTANT
BOLTZMAN_CONSTANT = 1.3806488E-23

# global SPEED_OF_LIGHT
SPEED_OF_LIGHT = 299792458

# global PLANCK_CONSTANT
PLANCK_CONSTANT = 6.62606957E-34

# global SOLAR_ANGULAR_VELOCITY
SOLAR_ANGULAR_VELOCITY = 2.97E-6

# global SOLAR_EFFECTIVE_TEMPERATURE
SOLAR_EFFECTIVE_TEMPERATURE = 5780

# global AU
AU = 149597871000  # m

# /physical constants

# global TEMPERATURE_LIST_LD
TEMPERATURE_LIST_LD = [3500.0, 3750.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0,
                       6000.0, 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0, 8250.0,
                       8500.0, 8750.0, 9000.0, 9250.0, 9500.0, 9750.0, 10000.0, 10250.0, 10500.0, 10750.0,
                       11000.0, 11250.0, 11500.0, 11750.0, 12000.0, 12250.0, 12500.0, 12750.0, 13000.0,
                       14000.0, 15000.0, 16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 21000.0, 22000.0,
                       23000.0, 24000.0, 25000.0, 26000.0, 27000.0, 28000.0, 29000.0, 30000.0, 31000.0,
                       32000.0, 33000.0, 34000.0, 35000.0, 36000.0, 37000.0, 37500.0, 38000.0, 39000.0,
                       40000.0, 41000.0, 42000.0, 42500.0, 43000.0, 44000.0, 45000.0, 46000.0, 47000.0,
                       47500.0, 48000.0, 49000.0, 50000.0]

# global GRAVITY_LIST_LD
GRAVITY_LIST_LD = [0.0, 0.5, 1.0, 1.5, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

# global METALLICITY_LIST_LD
METALLICITY_LIST_LD = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
                       0.1, 0.2, 0.3, 0.5, 1.0]

# global TEMPERATURE_LIST_ATM
TEMPERATURE_LIST_ATM = [3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250,
                        7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750, 11000,
                        11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000, 14000, 15000, 16000, 17000, 18000,
                        19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                        32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000, 40000, 41000, 42000, 43000, 44000,
                        45000, 46000, 47000, 48000, 49000, 50000]  # 76 values

# global GRAVITY_LIST_ATM
GRAVITY_LIST_ATM = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]  # 11 values

# global METALLICITY_LIST_ATM
METALLICITY_LIST_ATM = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]  # 8 values

METALLICITY_LIST_ATM = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2,
                        0.3, 0.5, 1.0]

# global GRAVTEMP_RANGE
GRAVTEMP_RANGE = {0.0: [3500.0, 6000.0],
                  0.5: [3500.0, 7500.0],
                  1.0: [3500.0, 8250.0],
                  1.5: [3500.0, 9000.0],
                  2.0: [3500.0, 14000.0],
                  # pozor, tu su hodnoty len do 11750 a potom je medzera a len jedna hodnota pre 14000
                  2.5: [3500.0, 19000.0],
                  3.0: [3500.0, 26000.0],
                  3.5: [3500.0, 31000.0],
                  4.0: [3500.0, 39000.0],
                  4.5: [3500.0, 49000.0],
                  5.0: [3500.0, 50000.0]}

# chess direction of interpolation on logg and T in matrix
CHESS_DIRECTION = [[-1, +1, -12, +12, -13, +13, -14, +14],
                   [-1, +1, -12, +12, -13, +13, -14, +14, -2, +2, -11, +11, -15, +15, -26, 26, -27, 27]]

# mysql login
HOST = "localhost"
USER = "root"
PWD = "toor"


# atmosphere
ATMOSPHERE = ["black-body", "castelli-kurucz-04"]