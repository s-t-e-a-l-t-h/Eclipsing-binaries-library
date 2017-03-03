from scipy.spatial import KDTree
import operator
import numpy as np
import objects.Projection as Projection
import objects.Sat as Sat
import objects.Function as Fn

def center_of_mass(faces=None):
    return np.array([(tri[0] + tri[1] + tri[2]) / 3.0 for tri in faces])


def distances_matrix(vertices):
    avsp = Fn.average_spacing(data=vertices, neighbours=5)
    tree = KDTree(vertices)
    dist = tree.sparse_distance_matrix(tree, max_distance=avsp * 4)
    dist_array = {}

    for i in range(0, len(vertices)):
        for j in range(0, len(vertices)):
            if i == j:
                continue
            if dist[(i, j)] == 0.0:
                # dist_array[(i, j)] = np.inf
                continue
            dist_array[(i, j)] = dist[(i, j)]

    sorted_dist_array = dict(sorted(dist_array.items(), key=operator.itemgetter(1)))
    del(dist_array, dist)

    return sorted_dist_array


def sim2tri(vertices, simplices):
    return [[vertices[sim[0]], vertices[sim[1]], vertices[sim[2]]] for sim in simplices]


def find_common_simplices(tsim1, tsim2):
    return [s for s in tsim1 if s in tsim2]


def find_not_common_simplices(simplices, common):
    return [s for s in simplices if s not in common]


def delaunay2d(vertices):
    from scipy.spatial import Delaunay as Del
    return Del(np.array(vertices, dtype="float64")).simplices


def delaunay3d(vertices):
    from scipy.spatial import Delaunay as Del
    vertices = np.array(vertices)
    convex_hull = Del(vertices).convex_hull
    return convex_hull, vertices[convex_hull]

class Tri:
    def __init__(self, vertices=None, norms=None):
        self.vertices = np.array(vertices)
        self.norms = np.array(norms)
        self._simplices = None
        self._hull = None

    def triangulate(self, srange=14, sindex=0):
        # srange (simplex range) (int)
        # sindex (start index) (int)

        # distances matrix
        d_mtx = distances_matrix(vertices=self.vertices)
        print("d_mtx done...")
        # first face
        fface = self.get_face(srange=srange, sindex=sindex, dmatrix=d_mtx)[0]

        # variables declar.
        tsimplices, tsimplices_queue, simplices_list = [fface], [fface], []

        while True:
            # break if computing queue is empty
            if len(tsimplices_queue) < 1:
                break

            # current vertices of triangle for finding next one
            current_tsimplices = tsimplices_queue[0]
            # run over all vertices in current triangle
            for simplex in current_tsimplices:
                new_tsimplices, intersc = None, None
                # trying to optim. doesn't work
                # if simplex not in simplices_list:
                #     simplices_list.append(simplex)
                # else:
                #     continue

                # new simplices
                nwsimplices = self.get_face(srange=srange, sindex=simplex, dmatrix=d_mtx)

                # iterate through every triangle created in current projection
                for simm in nwsimplices:
                    if simm in tsimplices:
                        # such simm already exists
                        continue

                    intersc = self.intersection(swhat=simm, swhere=tsimplices, srange=srange, dmatrix=d_mtx)

                    if intersc[0] is True:
                        use = False

                        # in many case happend, if triangles are in intersection, they have a two common points
                        # (see pic.) and then we should test, if we flip them, if there will be solution
                        #           D                      ________D
                        # A\\     /\                     A \\      \
                        #   \  \ /  \           ===>        \  \    \
                        #    \  / \  \                       \    \  \
                        #    B\/______\ C                    B\_______\ C

                        common = find_common_simplices(simm, intersc[1])
                        # if triangles have common simplices different to two, than continue,
                        # something like that cannot be flipped
                        if len(common) != 2:
                            continue

                        # otherwise, we have to create flipped case;
                        # in this case, there should be comnon only one point, so we are able to use [0]
                        ncmn_1 = find_not_common_simplices(simm, common)[0]
                        ncmn_2 = find_not_common_simplices(intersc[1], common)[0]
                        new_tsimplices = [sorted([ncmn_1, ncmn_2, common[0]]), sorted([ncmn_1, ncmn_2, common[1]])]

                        # we have to check, if this new_tsimplices is not in intersection with other,
                        # already existing triangles
                        for new_ts in new_tsimplices:
                            # there exist 2 triangles of course, but one of them should already exists in tsimplices
                            # array,
                            # should already exist, if there is intersection with it
                            if new_ts in tsimplices:
                                continue

                            # the other one, has to be checked "manualy";
                            # we will use the same code as above
                            intersc = self.intersection(swhat=new_ts, swhere=tsimplices, srange=srange, dmatrix=d_mtx)
                            if intersc[0] is False:
                                use, simm = True, new_ts
                                # we should break, because second one already exists
                                break

                    else:
                        use = True

                    if use is True:
                        tsimplices.append(simm)
                        tsimplices_queue.append(simm)
                        break
                del (intersc, nwsimplices, new_tsimplices)

            idx = tsimplices_queue.index(current_tsimplices)
            del tsimplices_queue[idx]

        self._simplices = tsimplices
        del (d_mtx, fface, tsimplices, tsimplices_queue, simplices_list)

    def nearest_to(self, to, srange, matrix):
        # closer_map is link between array closer and array self.vertices (in meaning of indices)
        closer, runner = [to], srange
        # get distance from self.vertices[sindex] to all other
        dist = self.distance_to(to, matrix)

        # append to closer srange closest points if dot product of normals are greather than 0
        for i in range(0, runner):
            # vertex of nearest point
            nsimplex = dist[i][0][1]
            if np.dot(self.norms[to], self.norms[nsimplex]) > 0:
                closer.append(nsimplex)
                continue
            # if dot product is negative, runner is incremented
            runner += 1

        return closer

    def intersection(self, swhat, swhere, srange, dmatrix):
        # vyhlada najblizsie trojuholniky z swhere k trojuholniku swhat a zisti ci sa medzi sebou tieto
        # najblizsie nepretinaju s swhat trojuholnikom

        # vezmeme jeden z bodov na trojuholniku swhat
        s = swhat[0]
        # nearest simplices to s point
        # nsimplices = sorted([dist[i][0][1] for i in range(0, srange)])
        closer = self.nearest_to(to=s, srange=srange, matrix=dmatrix)
        # simplices of nearest triangle from swhere to swhat
        nsimplices = []
        # iterate each close point
        for simplex in closer:
            # iterate each triangle already appended to final list (swhere)
            for face in swhere:
                if simplex in face and face not in nsimplices and face is not swhat:
                    nsimplices.append(face)

        # projection of triangles to 2d plane
        unique = np.unique(nsimplices)

        ps_projection = Projection.projection(n=self.norms[s],
                                              p=self.vertices[s],
                                              pp=-10,
                                              x=Fn.array_mask(array=self.vertices, mask=unique),
                                              dim="2d")

        # uprava aby sa to dalo podla indexov z povodneho pola precitat, ked sa pouzite funkcia sim2tri
        ps_projection_dict = {}
        for u, item in list(zip(unique, ps_projection)):
            ps_projection_dict[u] = item

        # project swhat triangle only to the same plane as other points (ps_projection)
        s_projection = Projection.projection(n=self.norms[s],
                                             p=self.vertices[s],
                                             pp=-10,
                                             x=Fn.array_mask(array=self.vertices, mask=swhat),
                                             dim="2d")

        # triangle 2d coo to comparsion;
        # nsimplices is in same order as twhere, so if twhere[1] is in intersect with twhat, than simplices for
        # twhere[1] are stored in nsimplices[1]
        twhere = sim2tri(vertices=ps_projection_dict, simplices=nsimplices)
        twhat = s_projection

        for i, t in list(zip(range(len(twhere)), twhere)):
            # tu pokracovat, treba s novym (mozno dalsim v zozname) trojuholnikom spravit intersekcny test s otanymi
            # teda s t-ckami v tomto loope a v momente, ked sa najde intersekcia vratit False, ze sa pretina a teda
            # je zlym kadnidatom na zaradenie, vlastne nemoze byt urcite zaradeny

            if Sat.intersection(twhat, t):
                return True, nsimplices[i]

        return False,

    def get_face(self, srange, sindex, dmatrix):
        closer = self.nearest_to(to=sindex, srange=srange, matrix=dmatrix)
        # projection of closer points to sindex normal plane
        ps_projection = Projection.projection(n=self.norms[sindex],
                                              p=self.vertices[sindex],
                                              pp=-10,
                                              x=Fn.array_mask(array=self.vertices, mask=closer),
                                              dim="2d")

        # first simplices (first triangle)
        fsimplices, return_simplices = delaunay2d(vertices=ps_projection), []

        for item in fsimplices:
            # zero in item, because sindex is added as firts value (zero valu in programmer term) to closer
            if 0 in item:
                return_simplices.append(sorted([closer[item[0]], closer[item[1]], closer[item[2]]]))
        return return_simplices

    def distance_to(self, to, matrix):
        dist, ml = {}, int(np.ceil(np.sqrt(len(self.vertices) ** 2 - len(self.vertices))))
        # povodne bolo toto za ml, funkcia bola staticmethod a v distance_matrix bolo odkomentovane np.inf riadok
        # int(np.ceil(np.sqrt(len(matrix))))
        for i in range(ml):
            try:
                if i == to:
                    continue
                dist[(to, i)] = matrix[(to, i)]
            except KeyError:
                continue

        sorted_dist_array = sorted(dist.items(), key=operator.itemgetter(1))

        return sorted_dist_array

    def load_csv(self, filename, separator, inplace=False):
        vertices, norms = [], []
        with open(filename, "r") as f:
            while True:
                data = f.readline()
                if not data:
                    break
                data = data.split(separator)
                try:
                    vertices.append([data[i] for i in range(0, 3)])
                    norms.append([data[i] for i in range(3, 6)])
                except IndexError:
                    break
        vertices, norms = np.array(vertices, dtype="float64"), np.array(norms, dtype="float64")

        if inplace:
            self.vertices = vertices
            self.norms = norms
        return vertices, norms

    def simplices(self):
        return self._simplices

    def hull(self):
        h = [[self.vertices[simplex[0]], self.vertices[simplex[1]], self.vertices[simplex[2]]] for simplex in
             self._simplices]
        self._hull = np.array(h[:]).tolist()
        return self._hull