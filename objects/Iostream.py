#!/usr/bin/python

""" Iostream trieda
-------------------
  Obsahuje funkcie pre vstupy a vystupy, standardne input/output stream

Inicializacne parametre:
------------------------
  Staticka trieda, ziadne inicializacne parametre

Metody a strucny popis (blizsi popis funkcnosti v komentaroch k jednotlivym riadkom kodu v funkciach):
------------------------------------------------------------------------------------------------------

  save_cgal_pcd_with_normals()
  ============================

    Vstupne parametre:
    ------------------
    filedir : [string], defaultna hodnota "tmp/", premenna, urcujuca zlozku, do ktorej sa ma *.pcd ulozit
    filename : [string], defaultna hodnota "3DModel.xyz", urcuje nazov pcd suboru
    points : [python list], defaultna hodnota, prazdny list [], python list 3D bodov pre ulozenie nasledujucej
                            struktury: [[x0, y0, z0], [x1, y1, z1], ..., [xn, yn, zn]]
    normals : [python list], defaultna hodnota, prazdny list [], python list normal indexovo prisluchajuci bodom
                             z premennej points v tvare: [[x0, y0, z0], [x1, y1, z1], ..., [xn, yn, zn]]
    verbose : [bool], defaultna hodnota "False", premenna pre zapinanie a vypinanie zobrazovania priebehu vykonavania
                      kodu

    Return:
    -------
    void

    Popis:
    ------
    ulozi point cloud document (pcd) - body + normaly pre dalsie spracovanie metodami kniznice CGAL

    Bug:
    ----
    ziaden znamy bug

  Funkcia()
  =========
    Vstupne parametre:
    ------------------
    name : [type], description

    Return:
    -------
    void

    Popis:
    ------

    Bug:
    ----
    ziaden znamy bug
"""

# import numpy as np
import globe.variables as gv
import objects.Function as Fn
import numpy as np
import os


def save_cgal_pcd_with_normals(  # point cloud with normals
        filedir='tmp/',
        filename='3DModel.xyz',
        points=None,
        normals=None,
        verbose=False):
    if verbose:
        print(
        Fn.color_string(color="info", string="Info: ") + "Saving point cloud with normals for CGAL triangulation.")

    filename, filedir = str(filename), str(filedir)
    f = open(filedir + filename, 'w')
    if os.path.isfile(filedir + filename):
        if verbose:
            print(Fn.color_string(color="info",
                                  string="Info: ") + "In reference: Iostream, function: save_cgal_pcd_with_normals(), line: " + str(
                Fn.lineno()) + ". File `" + filedir + filename + "` created successfully.")

        for idx in range(0, len(normals)):
            f.write(str(points[idx][0]) + " " + str(points[idx][1]) + " " + str(points[idx][2]) + " " + str(
                normals[idx][0]) + " " + str(normals[idx][1]) + " " + str(normals[idx][2]) + "\n")
    else:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Iostream, function: save_cgal_pcd_with_normals(), line: " + str(
                Fn.lineno()) + " File `" + gv.COLOR_WARNING + filedir + filename + gv.COLOR_END + "` cannot be created.")
        try:
            f.close()
        except:
            pass
        return False
    f.close()
    os.chmod(filedir + filename, 0o755)

    return True


def load_cgal_3d_mesh(
        filename="out.off",
        filedir="tmp/",
        verbose=False
):
    path = str(filedir) + str(filename)

    if os.path.isfile(path):
        vertices, indices, faces = [], [], []  # real points, index of triangle points

        try:
            f = open(path, 'r')
            for _ in range(3): next(f)  # preskoci 3 riadky hlavicky suboru
            for line in f:
                line = str(line).strip()
                if line is not "":
                    # vystupny cgal subor obsahuje samotne body a k nim prisluchajuce indexy trojuholnikov
                    # ak sa riadok zacina cislicou 3, tak sa jedna o znacenie indexov trojuholnikov
                    if line[:2] == "3 ":
                        line = line[3:]
                        num = line.split(" ")
                        indices.append([int(num[0]), int(num[1]), int(num[2])])
                    # inak sa jedna o samostatne body
                    else:
                        num = line.split(" ")
                        vertices.append([float(num[0]), float(num[1]), float(num[2])])
            f.close()

            # ulozi trojuholniky do rovnakeho tvaru ako Delaunayova trianguacia pri convex_hull
            for index in indices:
                faces.append([
                    vertices[index[0]],
                    vertices[index[1]],
                    vertices[index[2]]
                ])

            return np.array(faces), np.array(vertices), np.array(indices)
        except:
            if verbose:
                print(Fn.color_string(color="error",
                                      string="Error: ") + "In reference: Iostream, function: load_cgal_3d_mesh(), line: " + str(
                    Fn.lineno()) + ". Error has been occured during reading file `" + gv.COLOR_WARNING + path + gv.COLOR_END + "`.")
            return False
    else:
        if verbose:
            print(Fn.color_string(color="error",
                                  string="Error: ") + "In reference: Iostream, function: load_cgal_3d_mesh(), line: " + str(
                Fn.lineno()) + ". File `" + gv.COLOR_WARNING + path + gv.COLOR_END + "` does not exist.")
        return False


def save_gnuplot_rgb_surface(
        filename="surface",
        datafolder="surface",
        faces=None,
        rgb=None,
        minimum=None,
        maximum=None
):
    try:
        if os.path.exists(datafolder):
            import shutil
            shutil.rmtree(datafolder)
        os.makedirs(datafolder)
        if os.path.isdir(datafolder):
            filename += ".sh"
            f = open(filename, "w")

            if os.path.isfile(filename): os.chmod(filename, 0o755)

            f.write('#!/bin/bash\necho "set autoscale\nset view equal xyz\nset xlabel \'x/a\'\nset ylabel \'y/a\'\nset '
                    'zlabel \'z/a\'\nunset key\nset view 0, 0, 1.5\nset pm3d ftriangles\nset style fill transparent '
                    'solid 1.0\nset pm3d ftriangles\nset pm3d depthorder hidden3d 7\nset hidden3d\nset cbrange [' +
                    str(minimum) + ':' + str(maximum) + ']\nset palette defined (0 \\"purple\\", 1 \\"blue\\", 2 '
                                                        '\\"cyan\\", 3 \\"green\\", 4 \\"yellow\\", 5 \\"red\\")\nset '
                                                        'terminal wxt enhanced font \'Verdana,14\' persist\n')

            for idx_o in range(0, len(faces)):
                fl = open(datafolder + "/" + str(idx_o), 'w')
                for idx_i in range(0, 3):
                    fl.write(str(faces[idx_o][idx_i][0]) + " " + str(faces[idx_o][idx_i][1]) + " " + str(
                        faces[idx_o][idx_i][2]) + " 0x" + str(rgb[idx_o]) + "\r\n")
                    if idx_i == 0: fl.write("\n")
                fl.close()

                if idx_o == 0:
                    f.write("splot '" + datafolder + "/" + str(idx_o) + "' using 1:2:3:4 lc rgb variable w pm3d, ")
                else:
                    f.write("'" + datafolder + "/" + str(idx_o) + "' using 1:2:3:4 lc rgb variable w pm3d, ")

            f.write("\npause mouse keypress\" | gnuplot -persist")
            f.close()

            print(
            Fn.color_string("info", "Info: ") + "save_gnuplot_rgb_surface(): File `" + filename + "` has been written.")
    except:
        return False

    return True


def save_csv(filename, datafolder, delim, data, rewrite=True):
    try:
        dirext = os.path.exists(datafolder)
        if rewrite and dirext:
            import shutil
            shutil.rmtree(datafolder)

        if not dirext:
            os.makedirs(datafolder)

        if os.path.isdir(datafolder):
            filename += ".csv"
            f = open(datafolder + "/" + filename, "w")
            if os.path.isfile(filename): os.chmod(filename, 0o755)
        else:
            raise Exception()

        try:
            for t_row in data:
                w = []
                for t_col in t_row:
                    w.append(str(t_col) + str(delim))
                f.write(str("".join(w)[:-1] + "\n"))
        except:
            f.close()
    except:
        pass
