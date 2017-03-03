import numpy as np
import math as m
import objects.Function as Fn


def photometric_phase_of_eclipses(pp_primary, pp_secondary, lightcurve=None):
    # premenne
    # pp_primary: fotometricka faza primarneho zakrytu (Orbit.conjunction[0]["true_phase"])
    # pp_secondary: fotometricka faza sekundarneho zakrytu (Orbit.conjunction[1]["true_phase"])
    # lightcurve: svetelna krivka v tvare [[phase, flux], ... ]

    # studiom spravania sa faz a toho ako je naprogramovana konjunkcia v tejto triede (teda ze je normalizovana)
    # sa doslo k zaveru, ze krok od primarneho k sekundarnemu zakrytu sa spocita takto a rovnako tak potom
    # krok od sekundarneho k primarnemu ako 1.0 - ...
    phasearc_to_primary = abs(pp_primary - pp_secondary)
    phasearc_to_secondary = 1.0 - phasearc_to_primary

    # najde sa min a max faza na svetelnej krivke
    lc_ziped = list(zip(*lightcurve))
    min_lc_phase, max_lc_phase = min(lc_ziped[0]), max(lc_ziped[0])

    # znamienko minimalnej fazy (moze byt aj zaporna)
    # trosku nestandardne ak je +1 tak nastavujem nulu; to preto aby o riadok dalej som nemusel podmienkovat a len
    # tvrdo priratal
    sign = -1 if min_lc_phase < 0 else 0
    zero_eclipse = np.floor(min_lc_phase) + sign

    # premenne
    current_phase, counter, lc_eclipse = zero_eclipse, 0, []

    # cyklus, ktory prejde cez krivku
    while True:
        add = phasearc_to_secondary if counter % 2 == 0 else phasearc_to_primary
        if min_lc_phase <= current_phase <= max_lc_phase:
            lc_eclipse.append(current_phase)
        # toto elseif je tu preto, ze sa startuje z podstrelenej fazy, teda napr. ak je min faza 0.1 tak sa
        # startuje z 0.0 a teda to nespada do if vyssie a doslo by k brake-u, co je neziaduce, lebo by takto
        # nezbehol cyklus ani raz
        elif zero_eclipse <= current_phase < min_lc_phase:
            pass
        else:
            break
        current_phase += add
        counter += 1

    return lc_eclipse


def fix_minima(lightcurve=None, eclipses=None):
    lc_ziped = list(zip(*lightcurve))
    phase, flux, inverted_mask = lc_ziped[0], lc_ziped[1], []
    # definicia, aby nizsie PyCharm nepyskoval
    index_i, index_j = 0, 0

    # premenna fixed_processing je na to, aby sa vedelo ku ktoremu minimu do stredu eclipsu pridat bod,
    # resp. aby sa vedelo, v ktorom minime sa robila uprava
    fixed_processing = [False for _ in range(0, len(eclipses))]

    for i, eclipse in list(zip(range(0, len(eclipses)), eclipses)):
        # najde najblizsiu dostupnu hodnotu fazy na svetelnej krivke
        nearest = Fn.find_nearest_value(array=phase, value=eclipse)

        # kazde minimum je potrebne prejst dvakrat a to od stredu v smere do prava a rovnako tak do lava
        # smer vpravo je oznaceny ako +1 a smer vlavo ako -1
        for side in range(0, 2):
            index = nearest[1]

            # najblizsia hodnota fazy v krivke moze byt vlavo od stredu alebo vpravo od stredu a podla toho sa
            # zacne smerovost -1 alebo +1;
            # v pripade, ze sa jedna uz o druhy krok cyklu for, tak sa len flipne smer ako *= -1
            if side == 0 and nearest[0] >= eclipse:
                direction = 1
            elif side == 0 and nearest[0] < eclipse:
                direction = -1
            elif side == 1:
                direction *= -1

            # pokial je splnena podmienka na vyhadzovanie bodov, teda ze po sebe naasledujuce fluxy
            # maju tendenciu naberat nespravnu hodnotu (vlavo a vpravo od stredu klesat), tak bezi while
            while True:
                if direction == 1:
                    index_i, index_j = index + 1, index
                elif direction == -1:
                    index_i, index_j = index - 1, index

                # try je tu preto, ze posledna hodnota napr. moze byt prave minimum a uz ziadna dalsia nenasleduje
                # a algoritmus sa bude snazit na nu pristupovat, tak aby sa to nesprasilo za behu, kedze k tomu
                # nie je ziaden dovod
                try:
                    sign = -1 if flux[index_i] - flux[index_j] < 0 else 1
                    if sign == -1:
                        if not fixed_processing[i]:
                            fixed_processing[i] = True
                        if not index in inverted_mask:
                            inverted_mask.append(index)
                    else:
                        break
                except:
                    break
                index += 1 if direction == 1 else -1

    # vytvorenie krivky bez indexov v inverted_mask
    fixed_lightcurve = []
    for i, item in list(zip(range(0, len(lightcurve)), lightcurve)):
        if i in inverted_mask:
            continue
        else:
            fixed_lightcurve.append([item[0], item[1]])

    # pridanie tokov presne v zakryte
    fixed_lightcurve_ziped = list(zip(*fixed_lightcurve))
    for doit, eclipse in list(zip(fixed_processing, eclipses)):
        if doit:
            nearest_index = Fn.find_nearest_value(array=fixed_lightcurve_ziped[0], value=eclipse)[1]
            nearest_flux = fixed_lightcurve_ziped[1][nearest_index]

            fixed_lightcurve.append([eclipse, nearest_flux])
        else:
            continue

    fixed_lightcurve.sort(key=lambda x: x[0])
    return fixed_lightcurve


def gaussian_smooth(lightcurve=None):
    from scipy.signal import gaussian
    from scipy.ndimage import filters

    u, v = np.array(list(zip(*lightcurve))[0]), np.array(list(zip(*lightcurve))[1])
    g = gaussian(3, 1, True)
    w = filters.convolve1d(v, g / g.sum())
    lightcurve = np.array([[u[i], w[i]] for i in range(0, len(u))])
    return lightcurve


def akima_interpolation(from_photometric_phase, to_photometric_phase, lightcurve=None, mirror=False, spots=False):
    u, v = np.array(list(zip(*lightcurve))[0]), np.array(list(zip(*lightcurve))[1])

    from scipy import interpolate
    w = interpolate.Akima1DInterpolator(u, v)

    interp_start, interp_stop = from_photometric_phase, to_photometric_phase
    # preco taky krok? lebo pri nom bolo najlepsie vysledky vyhladzovania minim
    step, current, lightcurve = 0.007, interp_start, []

    if mirror and not spots:
        while True:
            current_decimal = abs(m.modf(current)[0])

            if current_decimal >= 0.5: current_decimal = 1.0 - current_decimal
            current_value = w(current_decimal)
            if not np.isnan(current_value):
                lightcurve.append([current, current_value])
            current += step
            if current > interp_stop: break
    elif mirror and spots:
        while True:
            current_decimal = abs(m.modf(current)[0])

            if current < 0: current_decimal = 1.0 - current_decimal

            current_value = w(current_decimal)
            if not np.isnan(current_value):
                lightcurve.append([current, current_value])
            current += step
            if current > interp_stop: break
    else:
        while True:
            current_value = w(current)
            if not np.isnan(current_value):
                lightcurve.append([current, current_value])
            current += step
            if current > interp_stop: break

    lightcurve = np.array(lightcurve)

    return lightcurve
