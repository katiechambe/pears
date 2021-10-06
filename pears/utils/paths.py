""" 
Collection of all paths used in analysis

Dependencies:
-------------
None
"""

__author__ = "Katie Chamberlain"
__date__   = "September 2021"

class SetupPaths:
    """ 
    Defines the paths that are read into each file
    """

    def __init__(
        self, 
        basedir="/xdisk/gbesla/mig2020/rsgrps/gbeslastudents/"
        ):

        self.path_basedir = basedir
        self.path_home = self.path_basedir + "katie/"
        self.path_pears = self.path_home + "pears/"

        # directories for illustris data
        self.path_illustris = self.path_basedir + "Illustris/"
        self.path_illustristng = self.path_basedir + "IllustrisTNG/"
        # group catalog paths
        self.path_illustrisdark = self.path_illustris + "GroupCatalogsDark/"
        self.path_illustrishydro = self.path_illustris + "GroupCatalogsHydro/"
        self.path_tngdark = self.path_illustristng + "TNG100-1-Dark/"
        self.path_tnghydro = self.path_illustristng + "TNG100-1/"
        # merger trees
        self.path_illustrisdark_trees = self.path_illustris + "Illustris-1-Dark-MergerTree/"
        self.path_illustrishydro_trees = self.path_illustris + "Illustris-1-MergerTree/"
        self.path_tngdark_trees = self.path_tngdark + "postprocessing/"
        self.path_tnghydro_trees = self.path_tnghydro + "postprocessing/"

        # pears directories
        self.path_data = self.path_pears + "data/"
        self.path_groups = self.path_data + "groups/"
        self.path_subhalos = self.path_data + "subhalos/"
        self.path_maxmass = self.path_data + "max_masses/"
        self.path_am_mass = self.path_data + "am_masses/"
        self.path_snapdata = self.path_data + "snapdata/"
