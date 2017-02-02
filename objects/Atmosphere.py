#!/usr/bin/python

import globe.variables as gv
import objects.Function as Fn
import numpy as np


# def planck(
#         wavelength=1E-9,
#         verbose=False,
#         *args  # (temperature, passband_fn)
# ):
#     if args == ():
#         # temperature, passband = None, 1.0
#         if verbose:
#             print Fn.color_string(color="error", string="ValueError: ") + \
#                   "In class Atmoshpere, function plack(), line " + str(
#                 Fn.lineno()) + ". Variable `args` is invalid."
#         return 0
#     else:
#         temperature, passband_fn = args[0], args[1]
#
#     default_settings = np.seterr()
#     np.seterr(over='ignore', divide='ignore', invalid='ignore', under='ignore')
#
#     try:
#         # wavelength pre passband funkciu sa nasobi 1e9, pretoze horna a dolna hranica integrovania je nastavena v
#         # metroch a aj tu pride v metroch. passband je ale nainterpolovana na nanometre, takze ak tu dostanem cislo
#         # 200e-9, tak ho pre passband potrebujem previest na nm ako 200e-9 * 1e9.
#
#         # premenna passpand sa tu dostane ako funkcia v argumente, cize premenna passband nie je premenna, ale funkcia
#         result = (((2.0 * gv.PLANCK_CONSTANT * (np.power(gv.SPEED_OF_LIGHT, 2))) / (np.power(wavelength, 5))) * (1.0 / (
#         np.power(np.e, (
#         (gv.PLANCK_CONSTANT * gv.SPEED_OF_LIGHT) / (wavelength * gv.BOLTZMAN_CONSTANT * temperature))) - 1.0))) * abs(
#             passband_fn(wavelength * 1e9))
#     except:
#         return 0
#
#     np.seterr(**default_settings)
#     return result


def planck(
        wavelength=1E-9,
        temperature=None,
        passband_fn=None,
        verbose=False,
        *args  # (temperatur, passband)
):
    if Fn.empty(temperature) or not passband_fn:
        if verbose:
            print(Fn.color_string(color="error", string="ValueError: ") + \
                  "In class Atmoshpere, function plack(), line " + str(
                Fn.lineno()) + ". Variable `args` is invalid.")

    default_settings = np.seterr()
    np.seterr(over='ignore', divide='ignore', invalid='ignore', under='ignore')

    try:
        # wavelength pre planckovu funkciu sa nasobi 1e-9, pretoze argument pride v nanometroch
        # je to sposobene tym, ze passband funkcia je nainterpolovana v nanometroch
        # proste v DB je lepsie drzat hodnotu 300nm ako 3e-7m

        # premenna passpand sa tu dostane ako funkcia v argumente, cize premenna passband nie je premenna, ale funkcia
        result = (
                     ((2.0 * gv.PLANCK_CONSTANT * (np.power(gv.SPEED_OF_LIGHT, 2))) / (np.power(wavelength * 1e-9, 5))) * (
                         1.0 / (np.power(np.e, ((gv.PLANCK_CONSTANT * gv.SPEED_OF_LIGHT) / (
                             wavelength * 1e-9 * gv.BOLTZMAN_CONSTANT * temperature))) - 1.0))) * float(abs(
            passband_fn(wavelength)))
    except:
        return 0

    np.seterr(**default_settings)
    return result
