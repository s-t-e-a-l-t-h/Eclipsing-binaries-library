# intersection of two lines, defined by two points each
def intersection(p1, p2, p3, p4):
    dnt = ((p1[0] - p2[0]) * (p3[1] - p4[1])) - ((p1[1] - p2[1]) * (p3[0] - p4[0]))
    if abs(dnt) > 1e-15:
        p_x = ((((p1[0] * p2[1]) - (p1[1] * p2[0])) * (p3[0] - p4[0])) - (
        (p1[0] - p2[0]) * ((p3[0] * p4[1]) - (p3[1] * p4[0])))) / dnt

        p_y = ((((p1[0] * p2[1]) - (p1[1] * p2[0])) * (p3[1] - p4[1])) - (
        (p1[1] - p2[1]) * ((p3[0] * p4[1]) - (p3[1] * p4[0])))) / dnt

        return True, p_x, p_y, "INTERSECTING"
    else:
        return False, np.nan, np.nan, "PARALLEL"
