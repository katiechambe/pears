class SetupPaths:
    """ 
    Defines the paths that are read into each file
    """

    def __init__(
        self, 
        basedir="/xdisk/gbesla/mig2020/rsgrps/gbeslastudents/"
        ):

        # directories for group data
        self.path_basedir = basedir
        self.path_home = self.path_basedir + "katie/"
        self.path_illustris = self.path_basedir + "Illustris/"
        self.path_illustristng = self.path_basedir + "IllustrisTNG/"

        # directories in my home 
        self.path_pears = self.path_home + "pears/"

        # group catalog paths
        self.path_illustrisdark = self.path_illustris + "GroupCatalogsDark/"
        self.path_illustrishydro = self.path_illustris + "GroupCatalogsHydro/"
        self.path_tngdark = self.path_illustristng + "TNG100-1-Dark/"
        self.path_tnghydro = self.path_illustristng + "TNG100-1/"

        # merger tree paths
        self.path_illustrisdark_trees = self.path_illustris + "Illustris-1-Dark-MergerTree/"
        self.path_illustrishydro_trees = self.path_illustris + "Illustris-1-MergerTree/"
        self.path_tngdark_trees = self.path_tngdark + "MergerTree/"
        self.path_tnghydro_trees = self.path_tnghydro + "MergerTree/"

