# distutils: language = c++
# distutils: sources = src/networkinfo.cc
# distutils: sources = src/networkio.cc
# distutils: sources = src/graphstorage.cc

import sys
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool
from cython.operator cimport dereference as deref, preincrement as inc
cimport pyzeo.extension

#=============================================================================
# geometry
cdef class Xyz:
    """
    Class to store a point
    """
    
    def __cinit__(self, double x=0.0, double y=0.0, double z=0.0):
        self.thisptr = new XYZ(x,y,z)

    def __init__(self, double x=0.0, double y=0.0, double z=0.0):
        pass

    def __dealloc__(self):
        del self.thisptr

    property x:
        def __get__(self): return self.thisptr.x
        def __set__(self, x_in): self.thisptr.x = x_in 
    property y:
        def __get__(self): return self.thisptr.y
        def __set__(self, y_in): self.thisptr.y = y_in 
    property z:
        def __get__(self): return self.thisptr.z
        def __set__(self, z_in): self.thisptr.z = z_in 

    def scale(self, double factor):
        new_xyz = Xyz()
        self.thisptr.scale(factor, new_xyz.thisptr)
        return new_xyz

cdef class Point:
    """
    Class to store a point
    """
    
    def __cinit__(self, double x=0.0, double y=0.0, double z=0.0):
        self.thisptr = new CPoint(x,y,z)

    def __init__(self, double x=0.0, double y=0.0, double z=0.0):
        pass

    def __dealloc__(self):
        del self.thisptr

    def __repr__(self):
        return "("+str(self.x)+','+str(self.y)+','+str(self.y)+')'

    property x:
        def __get__(self): return self.thisptr.vals[0]
        def __set__(self, x_in): self.thisptr.vals[0] = x_in 
    property y:
        def __get__(self): return self.thisptr.vals[1]
        def __set__(self, y_in): self.thisptr.vals[1] = y_in 
    property z:
        def __get__(self): return self.thisptr.vals[2]
        def __set__(self, z_in): self.thisptr.vals[2] = z_in 

    #def scale(self, double scaling_factor):
    #    return self.thisptr.scale(scaling_factor)

#=============================================================================
# netinfo
#Python definitions for the cdefinitions in .pxd file
def initializeRadTable():
    """
    Populate the atomic radius table with Zeo++ default values
    """
    zeo_initializeRadTable()

def initializeCovRadTable():
    """
    Populate the covalent tradius table with Zeo++ default values
    """
    zeo_initializeCovRadTable()

def initializeMassTable():
    """
    Populate the atomic mass table with Zeo++ default values
    """
    zeo_initializeMassTable()

def initializeAtomCharacterTable():
    """
    Populate the Atom symbol table with Zeo++ default values
    """
    zeo_initializeAtomCharacterTable()

def initializeAtomicNumberTable():
    """
    Populate the atomic number table with Zeo++ default values
    """
    zeo_initializeAtomicNumberTable()

def readRadTable(filename):
    """
    Read atomic radii values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readRadTable(c_filename)

def readMassTable(filename):
    """
    Read atomic mass values from input file and replace the default values
    """
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    zeo_readMassTable(c_filename)

def lookupRadius(element):
    """"
    Args:
        element:
            Element name in conventional shorthand 
            Ex: Al for aluminum 
                Si for silicon 
    Returns:
        radius of the input element
    """
    radius = zeo_lookupRadius(element, True)
    return radius
    
def lookupCovRadius(element):
    return zeo_lookupCovRadius(element)

def lookupMass(element):
    return zeo_lookupMass(element)

def lookupAtomicNumber(element):
    return zeo_lookupAtomicNumber(element)

def isMetal(element):
    return zeo_isMetal(element)

#=============================================================================
# channel
cdef class Channel:
    """
    Python wrapper class to store the Zeo++ Channel objects derived from the Voronoi network.
    """
    def __cinit__(self):
        self.thisptr = new CHANNEL()
    def __dealloc__(self):
        del self.thisptr

    def find_bounding_atoms(self, AtomNetwork atmnet, VoronoiNetwork vornet):
        """
        Get the IDs of all atoms that bound this channel.
        An atom is considered to bound a channel if a Voronoi
        node of the channel is a member of the atom's Voronoi
        cell.

        Args:
        atmnet: AtomNetwork object
            The atom network structure of the associated channel
        vornet: VoronoiNetwork object
            The Voronoi network structure of the associated channel

        Returns:
            List of indices of the atoms that bound this channel
        """
        if not vornet.has_bvcells:
            raise ValueError("Voronoi cells not found. Run decomposition first.")

        cdef vector[int] atom_ids

        self.thisptr.findBoundingAtoms(atmnet.thisptr, vornet.bvcells, atom_ids)

        return [atom_ids[i] for i in range(atom_ids.size())]

#=============================================================================
# psd
def calc_pore_size_distribution(atmnet,  channel_radius, probe_radius, 
        mc_sampling_no, hist_file, high_accuracy=False, exclude_pockets=False, 
        points_file="", node_radii_file="", sphere_dist_file="", 
        vis_flag=False, overlap_check_flag=False):
    """
    Computes the pore size distribution histogram
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        hist_file:
           File to store the histogram
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        points_file (Default=None):
            File to store the points. Used in visualization
        node_radii_file (Default=None):
            File to store the node radi. Used in visualizationi
        sphere_dist_file (Default=None):
            Reserved for future use
        vis_flag (Default=False)
            Visualization Flag
        overlap_check_flag (Default=False)
            VisIT Visualization related Flag
    """
    atmnet_copy = (<AtomNetwork?>atmnet).copy()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    c_atmnetcp_ptr = (<AtomNetwork?>atmnet_copy).thisptr 
    cdef string chist_file = hist_file
    cdef string cpnt_file = points_file
    cdef string cnd_file = node_radii_file
    cdef string csph_file = sphere_dist_file
    c_calcPoreSizeDistr (c_atmnetcp_ptr, c_atmnet_ptr, high_accuracy,
              channel_radius,  probe_radius, mc_sampling_no, exclude_pockets,
              chist_file, cpnt_file, cnd_file, csph_file, vis_flag, 
              overlap_check_flag)

#=============================================================================
# netio
def readCiffile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not  readCIFFile(c_filename, atmnet.thisptr, radialflag):
        raise ValueError        # Find the appropriate error and return it
    return atmnet

def readArcfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readARCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCucfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCUCFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readCssrfile(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readCSSRFile(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def readV1file(filename, radialflag):
    atmnet = AtomNetwork()
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    if not readV1File(c_filename, atmnet.thisptr, radialflag):
        raise IOError
    return atmnet

def writeCssrfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCSSR(c_filename, c_atmnet):
        raise IOError

def writeCiffile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToCIF(c_filename, c_atmnet):
        raise IOError

def writeV1file(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToV1(c_filename, c_atmnet):
        raise IOError

def writeNt2file(filename, vornet, minRad = None):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    if minRad:
        if not writeToNt2(c_filename, c_vornet_ptr, minRad):
            raise IOError
    else:
        if not writeToNt2(c_filename, c_vornet_ptr):
            raise IOError

def writeXyzfile(filename, atmnet, supercell_flag, is_duplicate_perimeter_atoms):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToXYZ(c_filename, c_atmnet, supercell_flag, 
            is_duplicate_perimeter_atoms):
        raise IOError

def writeVtkfile(filename, atmnet):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToVTK(c_filename, c_atmnet):
        raise IOError

def writeMopacfile(filename, atmnet, supercell_flag):
    if isinstance(filename, unicode):
        filename = (<unicode>filename).encode('utf8')
    cdef char* c_filename = filename
    cdef ATOM_NETWORK* c_atmnet = (<AtomNetwork?>atmnet).thisptr
    if not writeToMOPAC(c_filename, c_atmnet, supercell_flag):
        raise IOError


#=============================================================================
# netstorage
cdef class Atom:
    """
    Class to store the information about atom (or ion) in a structure.
    """
    def __cinit__(self):
        self.thisptr = new ATOM()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property coords:
        def __get__(self):
            coords = [self.thisptr.x, self.thisptr.y, self.thisptr.z]
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]

    property radius:
        def __get__(self): return self.thisptr.radius
        def __set__(self, radius): 
            print("This value is not supposed to be modified")
            self.thisptr.radius = radius

    @property
    def type(self):
        """
        Returns the atom type.
        """
        return self.thisptr.type.decode('utf-8')


cdef class AtomNetwork:
    """
    Class to store and manipulate the input atom network.
    """
    #Cython wrapper for Zeo++ ATOM_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radius
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new ATOM_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def copy(self):
        """
        Create a copy of the AtomNetwork instance
        """
        newatmnet = AtomNetwork()
        self.thisptr.copy(newatmnet.thisptr)
        newatmnet.rad_flag = self.rad_flag
        return newatmnet

    @property  
    def no_atoms(self):  
        """
        Returns the number of atoms in the AtomNetwork.
        """  
        return self.thisptr.no_atoms

    @property  
    def atoms(self):  
        """
        Returns the Atom objects in the Atomnetwork.
        """  
        atom_list = []  
        for i in range(self.thisptr.no_atoms):  
            atom = Atom()  
            atom.thisptr[0] = self.thisptr.atoms[i]  
            atom_list.append(atom)  
        return atom_list

    #def relative_to_absolute(self, point):
    #    cdef CPoint* cpoint_ptr = (<Point?>point).thisptr
    #    cdef double x = cpoint_ptr.vals[0]
    #    cdef double y = cpoint_ptr.vals[1]
    #    cdef double z = cpoint_ptr.vals[2]
    #    cdef CPoint abs_point = self.thisptr.abc_to_xyz(x,y,z)
    #    return Point(abs_point.vals[0], abs_point.vals[1], 
    #            abs_point.vals[2])

    #def absolute_to_relative(self, point):
    #    cdef CPoint* cpoint_ptr = (<Point?>point).thisptr
    #    cdef CPoint rel_point = self.thisptr.xyz_to_abc(abs_point)
    #    return Point(rel_point.vals[0], rel_point.vals[1], 
    #            rel_point.vals[2])

    @classmethod
    def read_from_CIF(cls, filename, rad_flag=True, rad_file=None,
                      mass_flag=True, mass_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CIF file.
        Arguments:
            filename: 
                Input CIF file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, Zeo++ default values are used.
            mass_flag (optional):
                Flag denoting whether atomic masses are used.
                Default is True
            mass_file (optional):
                Input file containing atomic masses
                Works only when mass_flag is True.
                If mass_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        #Calls Zeo++ readCIFFile function defined in networkio.cc.
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(mass_file, unicode):
            mass_file = (<unicode>mass_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_rad_file = NULL
        if rad_flag:
            pyzeo.extension.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                pyzeo.extension.zeo_readRadTable(c_rad_file)

        cdef char* c_mass_file = NULL
        if mass_flag:
            pyzeo.extension.zeo_initializeMassTable()
            if mass_file:
                c_mass_file = mass_file
                pyzeo.extension.zeo_readMassTable(c_mass_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCIFFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError

        loadMass(mass_flag, atmnet.thisptr)

        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_ARC(cls, filename, rad_flag=True, rad_file=None,
                      mass_flag=True, mass_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a ARC file.
        Arguments:
            filename: 
                Input ARC file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
            mass_flag (optional):
                Flag denoting whether atomic masses are used.
                Default is True
            mass_file (optional):
                Input file containing atomic masses
                Works only when mass_flag is True.
                If mass_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(mass_file, unicode):
            mass_file = (<unicode>mass_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readARCFile function defined in networkio.cc.
        cdef char* c_rad_file = NULL
        if rad_flag:
            pyzeo.extension.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                pyzeo.extension.zeo_readRadTable(c_rad_file)

        cdef char* c_mass_file = NULL
        if mass_flag:
            pyzeo.extension.zeo_initializeMassTable()
            if mass_file:
                c_mass_file = mass_file
                pyzeo.extension.zeo_readMassTable(c_mass_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readARCFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError

        loadMass(mass_flag, atmnet.thisptr)

        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_CSSR(cls, filename, rad_flag=True, rad_file=None,
                       mass_flag=True, mass_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a CSSR file.
        Arguments:
            filename: 
                Input CSSR file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
            mass_flag (optional):
                Flag denoting whether atomic masses are used.
                Default is True
            mass_file (optional):
                Input file containing atomic masses
                Works only when mass_flag is True.
                If mass_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(mass_file, unicode):
            mass_file = (<unicode>mass_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readCSSRFile function defined in networkio.cc.
        cdef char* c_rad_file

        if rad_flag:
            pyzeo.extension.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                pyzeo.extension.zeo_readRadTable(c_rad_file)

        cdef char* c_mass_file = NULL
        if mass_flag:
            pyzeo.extension.zeo_initializeMassTable()
            if mass_file:
                c_mass_file = mass_file
                pyzeo.extension.zeo_readMassTable(c_mass_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readCSSRFile(c_filename, atmnet.thisptr, rad_flag):
            raise IOError

        loadMass(mass_flag, atmnet.thisptr)

        atmnet.rad_flag = rad_flag
        return atmnet

    @classmethod
    def read_from_V1(cls, filename, rad_flag=True, rad_file=None,
                     mass_flag=True, mass_file=None):
        """
        Static method to create and populate the AtomNetwork with 
        atom data from a V1 file.
        Arguments:
            filename: 
                Input V1 file name.
            rad_flag (optional):
                Flag denoting whether atomic radii are non-zero.
                Default is True
            rad_file (optional):
                Input file containing atomic radii
                Works only when rad_flag is True.
                If rad_file is not specified, default values are used.
            mass_flag (optional):
                Flag denoting whether atomic masses are used.
                Default is True
            mass_file (optional):
                Input file containing atomic masses
                Works only when mass_flag is True.
                If mass_file is not specified, Zeo++ default values are used.
        Returns:
            Instance of AtomNetwork
        """
        if isinstance(rad_file, unicode):
            rad_file = (<unicode>rad_file).encode('utf8')
        if isinstance(mass_file, unicode):
            mass_file = (<unicode>mass_file).encode('utf8')
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ readV1File function defined in networkio.cc.
        cdef char* c_rad_file = NULL
        if rad_flag:
            pyzeo.extension.zeo_initializeRadTable()
            if rad_file:       # rad_file is defined
                c_rad_file = rad_file
                pyzeo.extension.zeo_readRadTable(c_rad_file)

        cdef char* c_mass_file = NULL
        if mass_flag:
            pyzeo.extension.zeo_initializeMassTable()
            if mass_file:
                c_mass_file = mass_file
                pyzeo.extension.zeo_readMassTable(c_mass_file)

        atmnet = AtomNetwork()
        cdef char* c_filename = filename
        if not readV1File(c_filename, atmnet.thisptr, rad_flag):
            raise IOError

        loadMass(mass_flag, atmnet.thisptr)

        atmnet.rad_flag = rad_flag
        return atmnet

    def write_to_CSSR(self, filename):
        """
        Writes the atom data in AtomNetwork to a CSSR file.
        Arguments:
            filename: 
                Output CSSR file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCSSR function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCSSR(c_filename, self.thisptr):
            raise IOError

    def write_to_CIF(self, filename):
        """
        Writes the atom data in AtomNetwork to a CIF file.
        Arguments:
            filename: 
                Output CIF file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToCIF function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToCIF(c_filename, self.thisptr):
            raise IOError

    def write_to_V1(self, filename):
        """
        Writes the atom data in AtomNetwork to a V1 file.
        Arguments:
            filename: 
                Output V1 file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToV1 function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToV1(c_filename, self.thisptr):
            raise IOError

    def write_to_XYZ(self, filename, supercell_flag, 
                     is_duplicate_perimeter_atoms):
        """
        Writes the atom data in AtomNetwork to an XYZ file.
        Arguments:
            filename: 
                Output XYZ file name.
            supercell_flag:
                Flag denoting whether to write 2x2x2 supercell.
            is_duplicate_perimeter_atoms:
                Flag denoting whether perimeter atoms need to be replicated.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToXYZ function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToXYZ(c_filename, self.thisptr, supercell_flag, 
                is_duplicate_perimeter_atoms):
            raise IOError

    def write_to_VTK(self, filename):
        """
        Writes the boundary of unit cell within the AtomNetwork to a VTK file.
        Arguments:
            filename: 
                Output VTK file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        #Calls Zeo++ writeToVTK function defined in networkio.cc.
        cdef char* c_filename = filename
        if not writeToVTK(c_filename, self.thisptr):
            raise IOError

    def write_to_MOPAC(self, filename, supercell_flag):
        """
        Writes the atom data in AtomNetwork to a .mop file.
        Arguments:
            filename: 
                Output MOPAC file name.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeToMOPAC(c_filename, self.thisptr, supercell_flag):
             raise IOError

    def calculate_free_sphere_parameters(self, filename):
        """
        Computes the diameters of the largest included sphere, free sphere 
        and included sphere along free sphere path. 
        Arguments:
            filename:
                Name of file where the diameters are stored.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        vornet, edge_centers, face_centers = self.perform_voronoi_decomposition()
        cdef char* c_fname = filename
        vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
        calculateFreeSphereParameters(vornet_ptr, c_fname, False)

    def perform_voronoi_decomposition(self, saveVorCells=True):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef vector[VOR_CELL] vcells

        #print self.rad_flag
        if not performVoronoiDecomp(self.rad_flag, self.thisptr, 
                vornet.thisptr, &vcells, saveVorCells, &vornet.bvcells):
            raise ValueError # Change it to appropriate error

        vornet.has_bvcells = True

        cdef int N

        # Get the edge centers
        edge_centers = []
        cdef vector[VOR_EDGE] vedges = vornet.thisptr.edges
        cdef vector[VOR_NODE] vnodes = vornet.thisptr.nodes
        for i in range(vedges.size()):
            edge_orig =  vedges[i].origin
            edge_end =  vedges[i].ending
            o_vnode = vnodes[edge_orig]
            e_vnode = vnodes[edge_end]
            edge_center = (o_vnode.x + e_vnode.x, \
                           o_vnode.y + e_vnode.y, \
                           o_vnode.z + e_vnode.z)
            edge_center = tuple(x/2 for x in edge_center)
            if edge_center not in edge_centers:
                edge_centers.append(edge_center)



        # Get the vorcells and obtain the face centers
        face_centers = []
        cdef vector[VOR_FACE] vfaces
        cdef vector[CPoint] vertices
        cdef CPoint* cpoint_ptr 
        #cdef map[int, int] id_maps
        cdef vector[int] node_ids
        face_node_ids = set()
        for i in range(vcells.size()):
            vfaces = vcells[i].faces
            for j in range(vfaces.size()):
                node_ids = vfaces[j].node_ids
                node_id_list = []
                for k in range(node_ids.size()):
                    node_id_list.append(node_ids[k])
                node_id_set = frozenset(node_id_list)
                if not node_id_set in face_node_ids:
                    face_node_ids.add(node_id_set)
                    centroid = Point()
                    cpoint_ptr = (<Point?>centroid).thisptr
                    vertices = vfaces[j].vertices
                    for k in range(vertices.size()):
                        centroid.x = centroid.x + vertices[k].vals[0]
                        centroid.y = centroid.y + vertices[k].vals[1]
                        centroid.z = centroid.z + vertices[k].vals[2]
                    centroid.x = centroid.x/vertices.size()
                    centroid.y = centroid.y/vertices.size()
                    centroid.z = centroid.z/vertices.size()
                    face_centers.append(centroid)

        # Convert the Zeo++ Point objects in (x,y,z) tuple objects
        fcs = []
        for center in face_centers:
            cntr = (center.x,center.y,center.z)
            fcs.append(cntr)

        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet, edge_centers, fcs


cdef class VoronoiNode:
    """
    Class to store the voronoi nodes with coordinates and radius
    """
    def __cinit__(self):
        self.thisptr = new VOR_NODE()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    property coords:
        def __get__(self):
            coords = list(self.thisptr.x, self.thisptr.y, self.thisptr.z)
            return coords
        def __set__(self, coords):      # Don't set this
            """
            This variable is not supposed to be modified manually
            """
            print("This value is not supposed to be modified")
            self.thisptr.x = coords[0]
            self.thisptr.y = coords[1]
            self.thisptr.z = coords[2]

    property radius:
        def __get__(self): return self.thisptr.rad_stat_sphere
        def __set__(self, rad): 
            print("This value is not supposed to be modified")
            self.thisptr.rad_stat_sphere = rad

cdef class VoronoiNetwork:
    """
    Class to store the Voronoi network generated from Voronoi decomposition
    of atom network.
    """
    #Cython wrapper for Zeo++ VORONOI_NETWORK class.
    #Contains a pointer to ATOM_NETWORK and a flag denoting whether radisu
    #for each atomic species is non-zero. 
    def __cinit__(self):
        self.thisptr = new VORONOI_NETWORK()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def size(self):
        return self.thisptr.nodes.size()

    def prune(self, radius):
        """
        Removes the edges that do not allow a sphere to pass.
        Arguments:
            radius:
                Radius of the sphere
        Returns:
            Instance of VoronoiNetwork with edges pruned.
        """
        cdef VORONOI_NETWORK newcvornet = self.thisptr.prune(radius)
        newvornet = VoronoiNetwork()
        newvornet.thisptr = &newcvornet
        return newvornet

    def analyze_writeto_XYZ(self, name, double probeRad, atmnet, 
            int shift_x=0, int shift_y=0, int shift_z=0):
        """
        Create diagrams of 1) Voronoi network and 2) accessible Voronoi 
        network, and write the diagrams in VTK files and the Voronoi 
        networks in XYZ files. Useful for visualizing the Voronoi network.
        Args:
            name:
                Name to be used for output files.
            probeRad:
                Radius of the probe.
            atmnet:
                pyzeo.extension.AtomNetwork
            shift_x (default=0):
                Shift the accessible Voronoi network along x-axis
            shift_y (default=0):
                Shift the accessible Voronoi network along y-axis
            shift_z (default=0):
                Shift the accessible Voronoi network along z-axis
        """
        if isinstance(name, unicode):
            name = (<unicode>name).encode('utf8')

        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef char* cname = name
        visVoro(name, probeRad, shift_x, shift_y, shift_z, self.thisptr, 
                c_atmnetptr)

    def write_to_XYZ(self, filename, double cutoff_radius=0):
        """
        Write only voronoi node information to XYZ file.
        Args:
            filename:
                string
                Name of file to which voronoi node info is written.
            cutoff_radius:
                float
                Threshold radius (default=0)
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode('utf8')

        cdef char* c_filename = filename
        if not writeVornetToXYZ(c_filename, self.thisptr, 
                cutoff_radius):
            raise ValueError

    @classmethod
    def perform_voronoi_decomposition(cls, atmnet, saveVorCells=False):
        """
        Performs weighted voronoi decomposition of atoms in the AtomNetwork 
        to analyze void space and generate voronoi nodes, edges and faces.
        Arguments:
            saveVorCells (optional): 
                Flag to denote whether to save the VorCells.
                Reserved for future use, so ignore this.
        Returns:
            Instance of VoronoiNetwork
        """
        #Calls Zeo++ performVoronoiDecomp function defined in network.cc.
        vornet = VoronoiNetwork()  
        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
        cdef vector[VOR_CELL] vcells
        cdef vector[BASIC_VCELL] bvcells
        #print(self.rad_flag)
        if not performVoronoiDecomp(atmnet.rad_flag, c_atmnetptr, 
                vornet.thisptr, &vcells, saveVorCells, &bvcells):
            raise ValueError # Change it to appropriate error
        #cdef int N
        #vorcelllist = []
        #bvcelllist = []
        # define copy methods for BASIC_VCELL and VOR_CELL methods
        #for i in range(vcells.size()):
        #    pass
            #vorcell = VorCell()
            #vorcell.thisptr = &(vcells[i])
            #vorcelllist.append(vcells[i])
        #for i in range(bvcells.size()):
        #    pass
            #basicvcell = BasicVCell()
            #basicvcell.thisptr = &(bvcells[i])
            #bvcelllist.append(bvcells[i])
        return vornet

    def find_channels(self, double channel_radius):
        """
        Find channels in a Voronoi network.

        Identifies channels (sets of accessible Voronoi nodes) within the associated
        Voronoi network for a given probe radius.

        Args:
            channel_radius : float
                Radius of probe used to determine the accessibility of void space.

        Returns:
            1) py_channels: List of Channel objects representing accessible channels
            2) py_access_info: List of Booleans, where access_info[i] indicates
               if Voronoi node i is accessible for the probe
        """
        cdef vector[cpp_bool] access_info
        cdef vector[CHANNEL] c_channels

        c_findChannelsInVorNet(self.thisptr, channel_radius, &access_info, &c_channels)

        py_access_info = [access_info[i] for i in range(access_info.size())]

        py_channels = []
        for i in range(c_channels.size()):
            channel = Channel()
            channel.thisptr[0] = c_channels[i]
            py_channels.append(channel)

        return py_channels, py_access_info

def substitute_atoms(atmnet, substituteSeed, radialFlag):
    """
    Attempt to substitute every other Si atom with Al atom.
    AtomNetwork may only consist of Si and O atoms, where each Si atom 
    must be bonded to exactly 4 oxygen atoms and each oxygen atom must 
    be bonded to exactly 2 Si atoms. Raises exception if the substitution
    is not successful. 
    Args:
        atmnet:
            pyzeo.netstorage.AtomNetwork
        substiuteSeed:
            Boolean flag to specify whether the seeded Si atom is 
            substituted or not. Since only 2 configurations are possible 
            if the structure is consistent, changing this parameter enables 
            generation of all configurations. 
        radialFlag:
            Boolean flag to specify whether atomic sizes are to be used.
    Returns:
        If successful, returns AtomNetwork instance with Si replaced with Al
        and the number of substitutions. 
    """
    cdef int substitutionNo[1]
    atmnet_copy = AtomNetwork()
    c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    if not c_substituteAtoms(c_atmnet_ptr, atmnet_copy.thisptr, substituteSeed,
            substitutionNo, radialFlag):
        raise ValueError
    subNo = substitutionNo[0]
    return atmnet_copy, subNo

#=============================================================================
# graphstorage
cdef class DijkstraNetwork:
    """
    Python wrapper class to Zeo++ Djikstra Network
    """
    #cdef DIJKSTRA_NETWORK* thisptr
    def __cinit__(self):
        self.thisptr = new DIJKSTRA_NETWORK()
    @classmethod
    def from_VoronoiNetwork(vornet):
        """
        Build Dijkstra Net from input Voronoi Net
        """
        dijkstranet = DijkstraNetwork()
        c_vornet = (<VoronoiNetwork?>vornet).thisptr
        buildDijkstraNetwork(c_vornet, dijkstranet.thisptr)
        return dijkstranet 
    def __dealloc__(self):
        del self.thisptr

#=============================================================================
# voronoicell
#cdef class VorFace:
#    #cdef VOR_FACE *thiptr
#    def __cinit__(self, vertices,  atmnet, vornet):
#        cdef vector[CPoint] c_vertices = (<vertices
#        cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
#        cdef VORONOI_NETWORK* c_vornetptr = (<VoronoiNetwork?>vornet).thisptr
#        self.thisptr = new VOR_FACE(c_vertices, c_atmnetptr, c_vornetptr)
#
#    def __dealloc__(self):
#        del self.thisptr

cdef class VorCell:
    #cdef VOR_CELL *thiptr
    def __cinit__(self):
        self.thisptr = new VOR_CELL()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

cdef class BasicVCell:
    #cdef BASIC_VCELL *thisptr
    def __cinit__(self):
        self.thisptr = new BASIC_VCELL()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.thisptr 

#=============================================================================
# cycle
def compute_centroid_4cycles(vornet):
    """
    Computes the centroid of the 4 corners of quadrilateral voronoi face
    Args:
        vornet:
            pyzeo.storage.VoronoiNetwork
    Returns:
        List of centroids in [(x1,y1,z1),(x2,y2,z2),...] format
    """

    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef vector[CYCLE] cycles
    cdef vector[int] ids

    if not compute_4cycle(c_vornet_ptr, &cycles, False, 1):
        raise ValueError

    centroid_list = []
    cdef vector[CYCLE].iterator it = cycles.begin()
    cdef vector[int].iterator iit
    while it != cycles.end():
        new_xyz = Xyz()
        centroid(&(deref(it)), new_xyz.thisptr, &ids)
        iit = ids.begin()
        #print(ids.size())
        id_set = set()
        while iit != ids.end():
            id_set.add(deref(iit))
            inc(iit)
        
        centroid_list.append({'ids':id_set, 'coords':new_xyz})
        inc(it)

    return centroid_list

def compute_face_centers(atmnet):
    """
    Compute the face centers of the voronoi network 
    """
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef vector[XYZ] points
    face_center(c_atmnet_ptr, &points)

#=============================================================================
# cluster
def warning(*objs):
    print("WARNING", *objs)
#    print("WARNING", *objs, file=sys.stderr)

def simplify_highaccuracy_vornet(atmnet):
    """
    Generates and simplifies high accuracy voronoi network 
    """
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    simplify_ha_vornet(c_atmnetptr)


def reduced_highaccuracy_vornodes(atmnet):
    """
    Generates simplified hgh accuracy voronoi network
    """
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    cdef vector[XYZ] xyz_vect
    high_accuracy_vornodes_reduction(c_atmnetptr, &xyz_vect)
    # Conver to list of Xyz
    xyz_list = []

    cdef vector[XYZ].iterator it = xyz_vect.begin()
    while it != xyz_vect.end():
        new_xyz = Xyz((deref(it)).x, (deref(it)).y, (deref(it)).z) #Infefficient
        xyz_list.append(new_xyz)
        inc(it)

    return xyz_list


def pruned_highaccuracy_voronoi_network(atmnet, delta=0.5):
    """
    Prunes hgh accuracy voronoi network by removing voronoi
    nodes close to the center of the bigger atoms.
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atomnet(ha_atmnet, "MED")
    vornet,ecs,fcs = ha_atmnet.perform_voronoi_decomposition()
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_ha_atmnetptr = (<AtomNetwork?>ha_atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornetptr = (<VoronoiNetwork?>vornet).thisptr
    prune_high_accuracy_voronoi_network(c_vornetptr, c_atmnetptr, 
            c_ha_atmnetptr, delta)
    return vornet

def get_nearest_largest_diameter_highaccuracy_vornode( atmnet, delta=0.25):
    """
    Get the reduced high accuracy voronoi network where only nodes that 
    has the largest diameter and within the cutoff distance to the nodes
    of the low accuracy voronoi network are retained. A one-one mapping
    of high accuracy voronoi nodes and low accuracy nodes is obtained.

    Input:
        atmnet: AtomNetwork object
        delta: cutoff (default = 0.25 angstroms)
    Output:
        Reduced voronoi network
    """
    #generate_simplified_highaccuracy_voronoi_network(atmnet)
    ha_vornet = pruned_highaccuracy_voronoi_network(atmnet, delta=0.7)
    #print('')
    #print('**********ONE DECOMPOSITION.************')
    #print('')
    vornet,ecs,fcs = atmnet.perform_voronoi_decomposition()
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print('')
    #print('*********WORKED TILL HERE*********')
    #print('')
    nearest_largest_diameter_ha_vornet(c_ha_vornet_ptr, c_vornet_ptr, 
            c_atmnet_ptr, red_vornet.thisptr, delta)
    return red_vornet

def generate_simplified_highaccuracy_voronoi_network(atmnet,delta=0.6):
    """
    Generate a simplified high accuracy voronoi network. 
    Uses Zeo++ high accuracy network and simplifies it such that only voronoi 
    nodes that belong to different atoms of original atom network are 
    retained. There can be different no. of voronoi nodes when compared with
    the voronoi nodes obtained with regular tesselation. 
    Input:
        atmnet: AtomNetwork object
    Output:
        Simplified high accuracy voronoi network
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atomnet(ha_atmnet, "LOW")
    vornet,ecs,fcs = atmnet.perform_voronoi_decomposition()
    ha_vornet,ecs,fcs = ha_atmnet.perform_voronoi_decomposition()
    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    if node_size == ha_node_size:
        warning('No high accuracy')
    return ha_vornet        # The processing below is eliminated temporarily
    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    if node_size == ha_node_size:
        # No need for simplification
        return vornet
    
    #print('')
    #print('**********ONE DECOMPOSITION.************')
    #print('')
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_ha_atmnet_ptr = (<AtomNetwork?>ha_atmnet).thisptr
    cdef VORONOI_NETWORK* c_vornet_ptr = (<VoronoiNetwork?>vornet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print('')
    #print('*********WORKED TILL HERE*********')
    #print('')
    simplify_high_accuracy_vornet(c_ha_vornet_ptr, c_ha_atmnet_ptr, 
             red_vornet.thisptr)
    #print('red_vornet size', red_vornet.size())
    further_red_vornet = VoronoiNetwork()
    #print('')
    #print('*********WORKED TILL HERE*********')
    #print('')
    nearest_largest_diameter_ha_vornet(red_vornet.thisptr, c_vornet_ptr, 
            c_atmnet_ptr, further_red_vornet.thisptr, delta)
    #print ''
    #print '*********WORKED TILL HERE TOO*********'
    #print ''
    #print 'further_red_vornet size', further_red_vornet.size()
    return further_red_vornet
    #print '********SIMPLIFIED_VORNET_COMPLETE*******'

def prune_voronoi_network_close_node(atmnet,delta=0.1):
    """
    Generate a pruned high accuracy voronoi network. 
    Uses Zeo++ high accuracy network and simplifies it such that only voronoi 
    nodes that are farther than "delta" are retained. 
    Input:
        atmnet: AtomNetwork object
    Output:
        Simplified high accuracy voronoi network
    """
    ha_atmnet = atmnet.copy()
    high_accuracy_atomnet(ha_atmnet, "MED")
    vornet,ecs,fcs = atmnet.perform_voronoi_decomposition()
    ha_vornet,ecs,fcs = ha_atmnet.perform_voronoi_decomposition()

    node_size = vornet.size()
    ha_node_size = ha_vornet.size()
    print (node_size, ha_node_size)
    if node_size == ha_node_size:
        warning('No high accuracy')
        #return vornet
    
    cdef ATOM_NETWORK* c_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef VORONOI_NETWORK* c_ha_vornet_ptr = (<VoronoiNetwork?>ha_vornet).thisptr
    red_vornet = VoronoiNetwork()
    #print ''
    #print '*********WORKED TILL HERE*********'
    #print ''
    geometry_pruning(c_ha_vornet_ptr, c_atmnet_ptr, delta, 
            red_vornet.thisptr)
    print (red_vornet.size())
    pruned_vornet = VoronoiNetwork()
    ha_prune_within_atom(red_vornet.thisptr, c_atmnet_ptr,
            delta, pruned_vornet.thisptr)
    print (pruned_vornet.size())
    return pruned_vornet

#=============================================================================
# area_volume
def volume(atmnet, channel_radius, probe_radius, 
        mc_sampling_no, high_accuracy=False, high_accuracy_atmnet=None, 
        exclude_pockets=True, low_dist_range=-1, high_dist_range=-1):
    """
    Calculates the volume of channels and pockets in a given strucutre.
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        high_accuracy_atmnet (Default=None):
            pyzeo.netstorage.AtomNetwork
            Optional high accuracy AtomNetwork. If not given and high_accuracy
            flag is set to True, then it is computed and returned.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        low_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
        high_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
    Returns:
        1) string containing channel and pocket volumes
        2) if high_accuracy=True and no input high_accuracy_atmnet is given,
           returns high_accuracy_atmnet for future use.
    """
    if high_accuracy and not high_accuracy_atmnet:
        high_accuracy_atmnet = atmnet.copy()
        high_accuracy_atomnet(high_accuracy_atmnet)
        ret_high_acc_atmnet = True
    else:
        ret_high_acc_atmnet = False

    if high_accuracy_atmnet and not high_accuracy:
        high_accuracy = True

    cdef ATOM_NETWORK* c_org_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_atmnet_ptr
    if high_accuracy_atmnet:
        c_atmnet_ptr = (<AtomNetwork?>high_accuracy_atmnet).thisptr
    else:
        tmp_atmnet = atmnet.copy()
        c_atmnet_ptr = (<AtomNetwork?>tmp_atmnet).thisptr

    vol_str = calcAV(c_atmnet_ptr, c_org_atmnet_ptr, high_accuracy, 
            channel_radius, probe_radius, mc_sampling_no, exclude_pockets,
            low_dist_range, high_dist_range)
    #print vol_str
    if ret_high_acc_atmnet:
        return vol_str, high_accuracy_atmnet
    else:
        return vol_str

    #lines = vol_str.split('\n')
    #for line in lines:
    #    if "Number_of_pockets" in line
    #        fields = line.split(" ")
    #        print fields[1], fields[3]

def surface_area(atmnet, channel_radius, probe_radius, 
        mc_sampling_no, high_accuracy=False, high_accuracy_atmnet=None, 
        exclude_pockets=True, extended_output=False):

    """
    Calculates the surface area of channels and pockets in a given strucutre.
    Args:
        atmnet:
            zoe.storage.AtomNetwork
        channel_radius:
            Radius of probe used to determine the accessibility of void space.
        probe_radius:
            Radius of probe used in Monte Carlo (MC) sampling of surface.
        mc_sampling_no:
            No. of MC samples per atom
        high_accuracy (Default=False):
            Optional flag to use high accuracy.
        high_accuracy_atmnet (Default=None):
            pyzeo.netstorage.AtomNetwork
            Optional high accuracy AtomNetwork. If not given and high_accuracy
            flag is set to True, then it is computed and returned.
        exclude_pockets (Default=True):
            Optional flag to include pockets.
        low_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
        high_dist_range(Default=-1):
            Use if you know the C++ Zeo++ code.
    Returns:
        1) string containing channel and pocket surface area
        2) if high_accuracy=True and no input high_accuracy_atmnet is given,
           returns high_accuracy_atmnet for future use.
    """
    if high_accuracy and not high_accuracy_atmnet:
        high_accuracy_atmnet = atmnet.copy()
        high_accuracy_atomnet(high_accuracy_atmnet)
        ret_high_acc_atmnet = True
    else:
        ret_high_acc_atmnet = False

    if high_accuracy_atmnet and not high_accuracy:
        high_accuracy = True

    cdef ATOM_NETWORK* c_org_atmnet_ptr = (<AtomNetwork?>atmnet).thisptr
    cdef ATOM_NETWORK* c_atmnet_ptr
    if high_accuracy_atmnet:
        c_atmnet_ptr = (<AtomNetwork?>high_accuracy_atmnet).thisptr
    else:
        tmp_atmnet = atmnet.copy()
        c_atmnet_ptr = (<AtomNetwork?>tmp_atmnet).thisptr

    sa_str = calcASA(c_atmnet_ptr, c_org_atmnet_ptr, high_accuracy,
            channel_radius, probe_radius, mc_sampling_no, exclude_pockets,
            extended_output)
    if ret_high_acc_atmnet:
        return sa_str, high_accuracy_atmnet
    else:
        return sa_str

#=============================================================================
# high_accuracy
_accuracy_kw = {
        "OCC","FCC","ACC","AQC","DDH","TIH","ICH","ICC","RIH","S4","S10","S20",
        "S30","S40","S50","S100","S500","S1000","S10000","DEF","HI","MED","LOW"
        }
def high_accuracy_atomnet(atmnet, accuracy_setting="LOW"):
    """
    Increases the accuracy of voronoi decomposition by replacing big
    atoms (spheres) with a number of small spheres.
    *** Modifies atmnet argument in place ***
    Args:
        atmnet:
            pyzeo.netstorage.AtomNetwork
            Is modified in place.
        accuracy_setting: 
            String specifying the accuracy settings.
            Possible choices are "OCC","FCC","ACC","AQC","DDH",
            "TIH","ICH","ICC","RIH","S4","S10","S20","S30","S40","S50",
            "S100","S500","S1000","S10000","DEF","HI","MED","LOW".
            Default is "DEF".
    """
    if not accuracy_setting in _accuracy_kw:
        raise ValueError("Accuracy setting not understood")
    cdef ATOM_NETWORK* c_atmnetptr = (<AtomNetwork?>atmnet).thisptr
    if isinstance(accuracy_setting, unicode):
        accuracy_setting = (<unicode>accuracy_setting).encode('utf8')
    cdef string acc_set = accuracy_setting
    setupHighAccuracyAtomNetwork(c_atmnetptr, acc_set)
