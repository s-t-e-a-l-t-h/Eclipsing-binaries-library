import numpy as np
import objects.Function as Fn

def fill_bottom_flat_triangle(face=None, boundary=None, d=None):
    v = [None, None, None]
    # sort by y
    face = sorted(face, key=lambda i: i[1], reverse=True)
    v[0] = face[0]
    v[1], v[2] = tuple(sorted(face[1:], key=lambda i: i[0], reverse=False))

    dx, dy = d
    # prepocet realnej suradnice na diskretnu suradnicu
    dy1, dy2 = int((boundary[3] - v[0][1]) // dy), int((boundary[3] - v[1][1]) // dy)
    # maximalny rozsah pre riadkove naplnanie
    max_x, min_x = int((v[2][0] - boundary[0]) // dx) + 1, int((v[1][0] - boundary[0]) // dx) - 1

    if abs(dy2 - dy1) <= 1 and dy > v[0][1] - v[1][1]:
        return np.array([])

    invslope1 = ((v[1][0] - v[0][0]) / (v[1][1] - v[0][1])) * dy
    invslope2 = ((v[2][0] - v[0][0]) / (v[2][1] - v[0][1])) * dy

    curx1, curx2 = v[0][0], v[0][0]
    curdx1, curdx2 = int((curx1 - boundary[0]) // dx), int((curx2 - boundary[0]) // dx)

    screen_f = []
    # vsetky okrem najnizsieho (flat) riadku
    for y in range(dy1, dy2, 1):
        for x in range(curdx1, curdx2 + 1, 1):
            screen_f.append([y, x])
        curx1 -= invslope1
        curx2 -= invslope2
        curdx1, curdx2 = int((curx1 - boundary[0]) // dx), int((curx2 - boundary[0]) // dx)
    # posledny riadok
    for x in range(curdx1, curdx2 + 1, 1):
        if x < min_x:
            continue
        if x > max_x: break
        screen_f.append([dy2, x])
    return np.array(screen_f)


def fill_top_flat_triangle(face=None, boundary=None, d=None):
    v = [None, None, None]
    # sort by y
    face = sorted(face, key=lambda i: i[1], reverse=True)
    v[2] = face[2]
    v[0], v[1] = tuple(sorted(face[:-1], key=lambda i: i[0], reverse=False))

    dx, dy = d
    # prepocet realnej suradnice na diskretnu suradnicu
    dy1, dy2 = int((boundary[3] - v[2][1]) // dy), int((boundary[3] - v[0][1]) // dy)
    # maximalny rozsah pre riadkove naplnanie
    max_x, min_x = int((v[1][0] - boundary[0]) // dx) + 1, int((v[0][0] - boundary[0]) // dx) - 1

    # otestovanie, ci sa nema vykonat krko v "y" smere, napriek tomu, ze velkost pixela je daleko vacsia
    # ako vzdialenost dvoch bodov v "y" smere
    if abs(dy2 - dy1) <= 1 and dy > v[2][1] - v[0][1]:
        return np.array([])

    invslope1 = ((v[2][0] - v[0][0]) / (v[2][1] - v[0][1])) * dy
    invslope2 = ((v[2][0] - v[1][0]) / (v[2][1] - v[1][1])) * dy
    curx1, curx2 = v[2][0], v[2][0]

    curdx1, curdx2 = int((curx1 - boundary[0]) // dx), int((curx2 - boundary[0]) // dx)

    screen_f = []
    # tu sa vykonaju setky riadky okrem posledneho (teda toho flat)
    for y in range(dy1, dy2, -1):
        for x in range(curdx1, curdx2 + 1, 1):
            screen_f.append([y, x])

        curx1 += invslope1
        curx2 += invslope2
        curdx1, curdx2 = int((curx1 - boundary[0]) // dx), int((curx2 - boundary[0]) // dx)
    # tu sa vykona ten flat riadok (mam to rozbite z dovodu vykonu, aby sa v kazdom loope nemuselo robit if)
    for x in range(curdx1, curdx2 + 1, 1):
        if x < min_x:
            continue
        elif x > max_x: break
        screen_f.append([dy2, x])
    return np.array(screen_f)


def fill_triangle(face=None, boundary=None, d=None):
    # sort by y and then by x
    face.sort(key=lambda i: i[1], reverse=True)

    if face[0][1] == face[1][1]:
        fill_top_flat_triangle(face=face, boundary=boundary, d=d)
    elif face[1][1] == face[2][1]:
        fill_bottom_flat_triangle(face=face, boundary=boundary, d=d)
    else:

        x = face[0][0] + ((face[1][1] - face[0][1]) / (face[2][1] - face[0][1])) * (face[2][0] - face[0][0])
        face_top_flat, face_bottom_flat = [face[1], face[2], [x, face[1][1]]], [face[0], face[1], [x, face[1][1]]]

        a, b = np.array([]), np.array([])
        a = fill_bottom_flat_triangle(face=face_bottom_flat, boundary=boundary, d=d)
        b = fill_top_flat_triangle(face=face_top_flat, boundary=boundary, d=d)

        if Fn.empty(a) and not Fn.empty(b):
            return b
        elif Fn.empty(b) and not Fn.empty(a):
            return a
        return np.concatenate((a, b), 0)
