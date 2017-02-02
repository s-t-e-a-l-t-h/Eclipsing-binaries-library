import matplotlib.pyplot as plt
import numpy as np


def edge_intersection_2d(pt1_xy, pt2_xy, pt3_xy, pt4_xy):
    # return: tuple
    #       0: intersection_status:
    #               False: parallel
    #               True:  intersection
    #       1: segment intersection:
    #               False:     no intersection
    #               True:      intersection between defined points
    #               numpy.nan: uknown
    #       2: intersection point x value if exists, if not numpy.nan
    #       3: intersection point y value if exists, if not numpy.nan
    #       4: distance if parallel
    # defs:
    #       x1, y1 = pt1_xy + u * (pt2_xy - pt1_xy) = pt1_xy + u * dp1
    #       x2, y2 = pt3_xy + v * (pt4_xy - pt3_xy) = pt3_xy + v * dp2
    #       dp1 = pt2 - pt1 = (pt2_x - pt1_x, pt2_y - pt1_y)
    #       dp2 = pt4 - pt3 = (pt4_x - pt3_x, pt4_y - pt3_y)

    # intersection:
    #       x1, y1 = x2, y2
    #       pt1_xy + u * dp1 = pt3_xy + v * dp2
    #
    #       in coo:
    #       pt1_x + u * dp1_x = pt3_x + v * dp2_x
    #       pt1_y + u * dp1_y = pt3_y + v * dp2_y
    #
    #       variables: u, v
    #       solution:
    #       d = (dp1_x * dp2_y) - (dp1_y * dp2_x)
    #       u = (((pt1_y - pt3_y) * dp2_x) - (dp2_y * (pt1_x - pt3_x))) / d
    #       v = (((pt1_y - pt3_y) * dp1_x) - (dp1_y * (pt1_x - pt3_x))) / d

    # first line
    pt1_x, pt1_y, pt2_x, pt2_y = pt1_xy[0], pt1_xy[1], pt2_xy[0], pt2_xy[1]
    dp1_x, dp1_y = pt2_x - pt1_x, pt2_y - pt1_y
    # parameter for pt1 and pt2

    # second line
    pt3_x, pt3_y, pt4_x, pt4_y = pt3_xy[0], pt3_xy[1], pt4_xy[0], pt4_xy[1]
    dp2_x, dp2_y = pt4_x - pt3_x, pt4_y - pt3_y
    # parameter for pt3 and pt4

    d = (dp1_x * dp2_y) - (dp1_y * dp2_x)

    # test if d < 1e-10
    # testing on zero, but precission should cause a problem
    if abs(d) < 1e-10:
        # test distance between lines
        # if general form is known (ax + by + c1 = 0 and ax + by + c2 = 0),
        # d = abs(c1 - c2) / sqrt(a**2 + b**2)
        # parametric equation in general:
        #   x, y = [pt1_x, pt1_y] + u * [T_x, T_y], where T is tangential vector defined as pt2 - pt1
        # N = (a, b) represent normal vector of line; a, b from general equation of line
        # N = [-Ty, Tx], can be obtained
        # general equation:
        #   -Ty * x + Tx * y + c = 0, then
        # c = Ty * pt1_x - Tx * pt1_y
        # finaly, general equation:
        #   -Ty * x + Tx * y + (Ty * pt1_x - Tx * pt1_y) = 0
        a1, b1, c1 = -dp1_y, dp1_x, (dp1_y * pt1_x) - (dp1_x * pt1_y)
        a2, b2, c2 = -dp2_y, dp2_x, (dp2_y * pt3_x) - (dp2_x * pt3_y)

        # sign of y coefficient (b1, b2)
        # we want same sign at y
        sign1 = +1 if b1 >= 0 else -1
        sign2 = +1 if b2 >= 0 else -1
        c2 *= +1 if sign1 == sign2 else -1

        d = abs(c2 - c1) / (np.sqrt(a1 ** 2 + b1 ** 2))

        int_segment, msg = (True, "OVERLAPPING") if d == 0 else (False, "PARALLEL")
        return int_segment, np.nan, np.nan, np.nan, d, msg

    # +0 because of negative zero (-0.0 is incorrect) formatting on output
    u = ((((pt1_y - pt3_y) * dp2_x) - (dp2_y * (pt1_x - pt3_x))) / d) + 0
    v = ((((pt1_y - pt3_y) * dp1_x) - (dp1_y * (pt1_x - pt3_x))) / d) + 0

    int_x, int_y = pt1_x + (u * dp1_x), pt1_y + (u * dp1_y)
    int_segment = True if 0.0 <= u <= 1.0 and 0.0 <= v <= 1.0 else False
    return True, int_segment, int_x, int_y, np.nan, "INTERSECTING"


# interception example in segment
# pt1, pt2, pt3, pt4 = [1.5, -.5], [1.5, 2.0], [1.5, 1.0], [1.0, 1.0]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# interception example out of segment
# pt1, pt2, pt3, pt4 = [-1.5, -2.5], [1.0, 1.0], [1.0, 0.0], [2.0, 0.0]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# identicall example
# pt1, pt2, pt3, pt4 = [-1.5, 0.0], [1.0, 0.0], [0.5, 0.0], [2.0, 0.0]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# parallel example
# pt1, pt2, pt3, pt4 = [0.0, 0.0], [1.0, 0.0], [0.0, 0.5], [1.0, 0.5]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# parallel example
# pt1, pt2, pt3, pt4 = [0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [2.0, 0.0]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# parallel example
# pt1, pt2, pt3, pt4 = [0.0, 0.5], [1.0, 0.5], [1.0, 0.5], [0.0, 0.5]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# not intersection example
# pt1, pt2, pt3, pt4 = [1.0, 0.5], [0.0, 1.0], [0.0, 0.5], [0.0, 0.0]
# print(edge_intersection_2d(pt1, pt2, pt3, pt4))

# fig = plt.figure()
# ax = fig.add_subplot(111, aspect="auto")
# ax.scatter([pt1[0], pt2[0]], [pt1[1], pt2[1]], color="r", s=200)
# ax.scatter([pt3[0], pt4[0]], [pt3[1], pt4[1]], color="b", s=100)
# ax.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], color="r")
# ax.plot([pt3[0], pt4[0]], [pt3[1], pt4[1]], color="b")
# ax.grid(True)
# plt.show()
