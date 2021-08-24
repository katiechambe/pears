class SetupPaths:
    """ 
    Defines the paths that are read into each file
    """

    def __init__(
        self, 
        basedir="/xdisk/gbesla/mig2020/rsgrps/gbeslastudents/"
        ):

        # directories for group data
        self.basedir = basedir
        self.path_home = self.basedir + "katie/"
        self.path_illustris = self.basedir + "Illustris/"
        self.path_illustristng = self.basedir + "IllustrisTNG/"

        # directories in my home 
        self.path_pears = self.home + "pears/"

        # group catalog paths
        self.path_illustrisdark = self.illustris + "GroupCatalogsDark/"
        self.path_illustrishydro = self.illustris + "GroupCatalogsHydro/"
        self.path_tngdark = self.illustristng + "TNG100-1-Dark/"
        self.path_tnghydro = self.illustristng + "TNG100-1/"

        # merger tree paths
        self.path_illustrisdark_trees = self.illustris + "Illustris-1-Dark-MergerTree/"
        self.path_illustrishydro_trees = self.illustris + "Illustris-1-MergerTree/"
        self.path_tngdark_trees = self.tngdark + "MergerTree/"
        self.path_tnghydro_trees = self.tnghydro + "MergerTree/"

