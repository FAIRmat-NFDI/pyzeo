# distutils: language = c++
# distutils: sources = src/networkinfo.cc
# distutils: sources = src/networkio.cc
# distutils: sources = src/networkstorage.cc
# distutils: sources = src/graphstorage.cc
# distutils: sources = src/network.cc
"""
Cython extension declaration for pyzeo.
"""
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.set cimport set as cset

#=============================================================================
# geometry
cdef extern from "../geometry.h":
    cdef cppclass XYZ:
        XYZ(double, double, double) except +
        double x, y, z
        void scale (const double sc, XYZ*) 

    cdef cppclass CPoint "Point":
        CPoint(double, double, double) except +
        double vals[3]
        CPoint scale (const double) const
        CPoint operator-(CPoint) 
        CPoint operator+(CPoint)
        CPoint operator*(CPoint)


cdef class Xyz:
    """
    Cython wrapper declaration for Zeo++ XYZ class defined in geometry.h
    Contains a pointer to XYZ.
    """
    cdef XYZ* thisptr


cdef class Point:
    """
    Cython wrapper declaration for Zeo++ Point class defined in geometry.h
    Contains a pointer to c_point.
    """
    cdef CPoint* thisptr


#=============================================================================
# netinfo
cdef extern from "../networkinfo.h":
    cdef void zeo_initializeRadTable "initializeRadTable"()

    cdef void zeo_initializeCovRadTable "initializeCovRadTable"()

    cdef void zeo_initializeMassTable "initializeMassTable"()

    cdef void zeo_initializeAtomCharacterTable "initializeAtomCharacterTable"()

    cdef void zeo_initializeAtomicNumberTable "initializeAtomicNumberTable"()

    cdef void zeo_readRadTable "readRadTable"(char *filename)

    cdef void zeo_readMassTable "readMassTable"(char *filename)

    cdef double zeo_lookupRadius "lookupRadius"(string atomType, bint radial)

    cdef double zeo_lookupCovRadius "lookupCovRadius"(string atomType)

    cdef double zeo_lookupMass "lookupMass"(string atomType)

    cdef int zeo_lookupAtomicNumber "lookupAtomicNumber"(string atomType)

    cdef bint zeo_isMetal "isMetal"(string atomType)

#=============================================================================
# channel
cdef extern from "../channel.h":
    cdef cppclass CHANNEL:
        CHANNEL() except +

cdef extern from "../channel.h" namespace "CHANNEL":
    cdef c_findChannelsInDijkstraNet "findChannels"(DIJKSTRA_NETWORK*, 
            vector[bint] *, vector[CHANNEL] *)
    cdef c_findChannelsInVorNet "findChannels"(VORONOI_NETWORK*, double, 
            vector[bint] *, vector[CHANNEL] *)

cdef class Channel:
    cdef CHANNEL* thisptr

#=============================================================================
# psd
cdef extern from "../psd.h":
    cdef void c_calcPoreSizeDistr "calcPoreSizeDistr"(ATOM_NETWORK *atmnet, 
            ATOM_NETWORK *orgAtomnet, bint highAccuracy, double r_probe_chan, 
            double r_probe, int numSamples, bint excludePockets, 
            string histFile, string pointsFile, string nodeRadiiFile,
            string spheresDistFile, bint visualize, bint overlapsCheck)


#=============================================================================
# netio
cdef extern from '../networkio.h':
    cdef void parseFilename(const char* filename, char* name, char* extension)

    cdef bint checkInputFile(char* filename)

    cdef bint readCIFFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readARCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCUCFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readCSSRFile(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint readV1File(char *filename, ATOM_NETWORK *cell, bint radial)

    cdef bint writeToCSSR(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToCIF(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToV1(char * filename, ATOM_NETWORK *cell)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet, double minRad)

    cdef bint writeToNt2(char *filename, VORONOI_NETWORK *vornet)

    cdef bint writeToXYZ(char *filename, ATOM_NETWORK *cell, bint is_supercell,
                         bint is_duplicate_perimeter_atoms)

    cdef bint writeToVTK(char *filename, ATOM_NETWORK *cell)

    cdef bint writeToMOPAC(char *filename, ATOM_NETWORK *cell, bint is_supercell)

    cdef bint writeVornetToXYZ "writeToXYZ"(char *filename, VORONOI_NETWORK*, double)

#=============================================================================
# netstorage
cdef extern from "../networkstorage.h":
    cdef cppclass ATOM:
        ATOM() except +
        double x, y, z
        double radius
        #string type
        #int specialID
        double mass
        double charge 

    cdef cppclass ATOM_NETWORK:
        ATOM_NETWORK() except +
        void copy(ATOM_NETWORK*)
        double a, b, c      # Lattice parameters
        double alpha, beta, gamma   # lattice angles
        int no_atoms "numAtoms"
        vector[ATOM] atoms
        CPoint abc_to_xyz(double, double, double)
        CPoint abc_to_xyz(CPoint)
        CPoint xyz_to_abc(double, double, double)
        CPoint xyz_to_abc(CPoint)

    cdef cppclass VOR_NODE:
        VOR_NODE() except +
        VOR_NODE(double, double, double, double, vector[int])
        double x, y, z
        double rad_stat_sphere

    cdef cppclass VOR_EDGE:
        VOR_EDGE() except +
        VOR_EDGE(int, int, double, int, int, int, double)
        int origin "from" 
        int ending "to"
        double rad_moving_sphere
        double length
        int delta_uc_x, delta_uc_y, delta_uc_z


    cdef cppclass VORONOI_NETWORK:
        VORONOI_NETWORK() except +
        VORONOI_NETWORK prune(double)
        vector[VOR_NODE] nodes
        vector[VOR_EDGE] edges

    cdef bint c_substituteAtoms "substituteAtoms"(ATOM_NETWORK*, ATOM_NETWORK*,
            bint, int*, bint)

    cdef bint c_fracSubstituteAtoms "fracSubstituteAtoms"(ATOM_NETWORK*, 
            ATOM_NETWORK*, bint, double, 
            int, int*, double*, bint)


# At present  the return value of performVoronoiDecomp is void*
# Compile it after void* is changed to bool in the original source file
cdef extern from "../network.h":
    cdef bint performVoronoiDecomp(bint, ATOM_NETWORK*, VORONOI_NETWORK*, 
            vector[VOR_CELL]*, bint, vector[BASIC_VCELL]*)

    cdef void calculateFreeSphereParameters(VORONOI_NETWORK*, char*, bint)

    cdef void viewVoronoiDecomp(ATOM_NETWORK*, double, string)

    cdef void loadRadii(ATOM_NETWORK*)

    cdef void loadMass(bint, ATOM_NETWORK*)

cdef extern from "../area_and_volume.h":
    cdef void visVoro(char* name, double probeRad, int skel_a, int skel_b, int skel_c,
            VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet)

cdef class Atom:
    """
    Cython wrapper class for Zeo++ ATOM class.
    """
    cdef ATOM* thisptr

cdef class AtomNetwork:
    """ 
    Cython wrapper class for Zeo++ ATOM_NETWORK class.
    Contains a pointer to ATOM_NETWORK and a flag denoting whether radius
    for each atomic species is non-zero.
    """
    cdef ATOM_NETWORK* thisptr
    cdef bint rad_flag

cdef class VoronoiNode:
    """
    Cython wrapper class for Zeo++ VOR_NODE class.
    """
    cdef VOR_NODE* thisptr

cdef class VoronoiNetwork:
    """ 
    Cython wrapper class for Zeo++ VORONOI_NETWORK class.
    """
    cdef VORONOI_NETWORK* thisptr

#=============================================================================
# netstorage
cdef extern from "../graphstorage.h":
    cdef cppclass DIJKSTRA_NODE:
        int id
        double x, y, z
        bint active
        DIJKSTRA_NODE() except +
    cdef cppclass DIJKSTRA_NETWORK:
        DIJKSTRA_NETWORK() except +
cdef extern from "../graphstorage.h" namespace "DIJKSTRA_NETWORK":
    cdef void buildDijkstraNetwork(VORONOI_NETWORK*, DIJKSTRA_NETWORK*)


#cdef class DijkstraNode:
#    """
#    Cython wrapper class for Zeo++ DIJKSTRA_NODE class.
#    """
#    cdef DIJKSTRA_NODE* thisptr


cdef class DijkstraNetwork:
    """
    Cython wrapper class for Zeo++ DIJKSTRA_NETWORK class.
    """
    cdef DIJKSTRA_NETWORK* thisptr

#=============================================================================
# voronoicell
cdef extern from "../voronoicell.h":
    cdef cppclass VOR_FACE:
        VOR_FACE(vector[CPoint], ATOM_NETWORK*, VORONOI_NETWORK*) except +
        vector[CPoint] vertices "orderedVertices"
        vector[int] node_ids "nodeIDs"


    cdef cppclass VOR_CELL:
        VOR_CELL() except +
        vector[VOR_FACE] faces
        int num_vertices "numVertices"
        #cmap[CPoint, int, bint
        cmap[int,int] id_mappings "idMappings"
        cmap[int, vector[int]] reverse_id_mappings "reverseIDMappings"
        cmap[int, CPoint] vertex_coords "vertexCoords"
        vector[cset[int]] edge_connections "edgeConnections"

        void add_edge "addEdge" (CPoint, CPoint to)
        void add_face "addFace" (VOR_FACE)

    cdef cppclass BASIC_VCELL:
        BASIC_VCELL() except +


cdef class VorFace:
    cdef  VOR_FACE* thisptr

cdef class VorCell:
    cdef VOR_CELL* thisptr

cdef class BasicVCell:
    cdef BASIC_VCELL* thisptr

#=============================================================================
# cycle
cdef extern from "../cycle.h":
    cdef cppclass CYCLE:
        CYCLE() except +
        double length
        vector[DIJKSTRA_NODE] nodes

    cdef bint compute_4cycle(VORONOI_NETWORK*, vector[CYCLE]*, bint, int)

    cdef void centroid(CYCLE*, XYZ*, vector[int]*)

    cdef void face_center(ATOM_NETWORK*, vector[XYZ]*)


cdef class Cycle:
    """
    Cython wrapper class for Zeo++ CYCLE class.
    Contains a pointer to CYCLE
    """
    cdef CYCLE* thisptr

#=============================================================================
# cluster
cdef extern from "../cluster.h":
    cdef void simplify_ha_vornet(ATOM_NETWORK *atmnt)

    cdef void high_accuracy_vornodes_reduction(ATOM_NETWORK*, vector[XYZ]*)

    cdef void prune_high_accuracy_voronoi_network(VORONOI_NETWORK* vornet, 
            ATOM_NETWORK* atmnet, ATOM_NETWORK* ha_atmnet, double delta)

    cdef void nearest_largest_diameter_ha_vornet(VORONOI_NETWORK* ha_vornet,
            VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet, 
            VORONOI_NETWORK* red_vornet, float cutoff)

    cdef void simplify_high_accuracy_vornet(VORONOI_NETWORK* ha_vornet,
            ATOM_NETWORK* ha_atmnet, VORONOI_NETWORK* red_vornet)

    cdef void geometry_pruning(VORONOI_NETWORK* ha_vornet, 
            ATOM_NETWORK* ha_atmnet, float cutoff, VORONOI_NETWORK* red_vornet)

    cdef void ha_prune_within_atom(VORONOI_NETWORK* ha_vornet, 
            ATOM_NETWORK* atmnet, float cutoff, VORONOI_NETWORK* red_vornet)

#=============================================================================
# area_volume
cdef extern from "../area_and_volume.h":
    cdef string calcAV(ATOM_NETWORK*, ATOM_NETWORK*, bint, double, double, 
            int, bint, double, double)

    cdef string calcASA(ATOM_NETWORK*, ATOM_NETWORK*, bint, double, double,
            int, bint, bint)

    cdef double calcDensity(ATOM_NETWORK*) 

#=============================================================================
# high_accuracy
cdef extern from "../sphere_approx.h":
    cdef void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, 
            string AccSetting)

#=============================================================================
# feature
cdef extern from "../feature.h":
    cdef cppclass FEATURE (CHANNEL):
        FEATURE(vector[int] nodeIds, DIJKSTRA_NETWORK* dnet,  int dim,
                int basisVecs[3][3])
        # C++ code needs modification
        int createSegments(ATOM_NETWORK*, VORONOI_NETWORK*, DIJKSTRA_NETWORK,
                char* filename, int initIndex)

#=============================================================================
# holograms
cdef extern from "../holograms.h":
    cdef void analyze_accessible_voronoi_pre_segment(VORONOI_NETWORK *vornet, 
            float probeRad, vector[bint] *accessInfo, char *name, 
            char *bin_directory) 

#=============================================================================
# string_add
cdef extern from "../string_additions.h":
    cdef int strCmpList(vector[string] list, string str)
