class SetupPaths:
    """ 
    Defines the paths that are read into each file
    """

    def __init__(
        self, 
        basedir="/xdisk/gbesla/mig2020/rsgrps/gbeslastudents/"
        ):

        # directories for group data
        self.basedir = self.basedir
        self.home = self.basedir + "katie/"
        self.illustris = self.basedir + "Illustris/"
        self.illustristng = self.basedir + "IllustrisTNG/"

        # directories in my home 
        self.pears = self.home + "pears/"

        # group catalog paths
        self.illustrisdark = self.illustris + "GroupCatalogsDark/"
        self.illustrishydro = self.illustris + "GroupCatalogsHydro/"
        self.tngdark = self.illustristng + "TNG100-1-Dark/"
        self.tnghydro = self.illustristng + "TNG100-1/"

        # merger tree paths
        self.illustrisdark_trees = self.illustris + "Illustris-1-Dark-MergerTree/"
        self.illustrishydro_trees = self.illustris + "Illustris-1-MergerTree/"
        self.tngdark_trees = self.tngdark + "MergerTree/"
        self.tnghydro_trees = self.tnghydro + "MergerTree/"

