"""
Savanna River F-Basin, 2012 ASCEM demo;

Script to convert the units of geochemistry parameters from Sergio
Bea's document "Savannah River F-Area Seepage Basins Geochemical
system description" (ascem_f_area_chem.pdf) into pflotran/amanzi
units.

Specific surface area: [m^2/g] ---> [cm^2/cm^3 bulk]
Site denisty: [sites/nm^2] ---> [moles/m^3 bulk]
CEC: [sites/nm^2] ---> [moles/m^3/bulk] 

"""

class Mineral(object):
    
    def __init__(self):
        self.name = None
        self._volume_fraction = None # [vol mineral / vol bulk]
        self._density = None # [g/cm^3]
        self._specific_surface_area_m2g = None # [m^2/g]
        self._specific_surface_area_cm2cm3 = None # [cm^2/cm^3]
        self._site_density_sitenm2 = None # [sites / nm^2]
        self._site_density_molm3bulk = None # [moles / m^3 bulk]
        self._cec_sitenm2 = None # [sites / nm^2]
        self._cec_molm3bulk = None # [moles / m^3 bulk]

    def __repr__(self):
        data = """  {0}:
    volume fraction = {1} [-]
    density = {2} [g/cm^3]
    specific surface area = {3} [m^2/g]
    site density = {4} [sites/nm^2]
    cec = {5} [sites/nm^2]
""".format(self.name, self._volume_fraction, self._density,
           self._specific_surface_area_m2g,
           self._site_density_sitenm2,
           self._cec_sitenm2)
        return data

    def setup(self,
              name = None,
              volume_fraction = 0.0,
              specific_surface_area_m2g = 0.0,
              density = 0.0,
              site_density_sitenm2 = 0.0,
              cec_sitesnm2 = 0.0):
        self.name = name
        self._volume_fraction = volume_fraction
        self._specific_surface_area_m2g = specific_surface_area_m2g
        self._density = density
        self._site_density_sitenm2 = site_density_sitenm2
        self._cec_sitenm2 = cec_sitesnm2

        self.convert_specific_surface_area()
        self.convert_site_density()
        self.convert_cec()

    def convert_specific_surface_area(self):
        """
        convert m^2/g to cm^2/cm^3 bulk
        """
        area_convert = 100.0**2 # (100cm/m)^2
        self._specific_surface_area_cm2cm3 = \
            (self._volume_fraction * self._specific_surface_area_m2g * 
             self._density * area_convert)

    def ssa_cm2cm3(self):
        return self._specific_surface_area_cm2cm3

    def convert_site_density(self):
        """
        convert site/nm^2 to mol / m^3
        """
        area_convert = (1.0e9)**2 # (10^9 nm/m)^2
        volume_convert = 100.0**3 # (100cm/m)**3
        mole_convert = 1.0 / 6.02e23 # 1 mole/6.02e23 sites

        self._site_density_molm3bulk = \
            (area_convert * volume_convert * mole_convert *
             self._site_density_sitenm2 * 
             self._specific_surface_area_m2g * self._density *
             self._volume_fraction)
        
    def site_density_mm3b(self):
        return self._site_density_molm3bulk

    def convert_cec(self):
        """
        convert site/nm^2 to mol / m^3

        NOTE: same conversion as site density
        """
        area_convert = (1.0e9)**2 # (10^9 nm/m)^2
        volume_convert = 100.0**3 # (100cm/m)**3
        mole_convert = 1.0 / 6.02e23 # 1 mole/6.02e23 sites

        self._cec_molm3bulk = \
            (area_convert * volume_convert * mole_convert *
             self._cec_sitenm2 * 
             self._specific_surface_area_m2g * self._density *
             self._volume_fraction)

    def cec_mm3b(self):
        return self._cec_molm3bulk

if __name__ == "__main__":

    print("""
ASCEM 2012 F-Area Phase 2 Demo : Geochemistry unit conversions from
data provided by Sergio Bea to pflotran/amanzi input
""")
    minerals = []
    quartz = Mineral()
    quartz.setup(name = "Quartz",
                 volume_fraction = 0.88,
                 specific_surface_area_m2g = 0.14,
                 density = 2.648,
                 site_density_sitenm2 = 10.0)
    minerals.append(quartz)

    kaolinite = Mineral()
    kaolinite.setup(name = "Kaolinite",
                    volume_fraction = 0.11,
                    specific_surface_area_m2g = 20.71,
                    density = 2.594,
                    site_density_sitenm2 = 2.3,
                    cec_sitesnm2 = 0.28)
    minerals.append(kaolinite)

    goethite = Mineral()
    goethite.setup(name = "Goethite",
                   volume_fraction = 0.016,
                   specific_surface_area_m2g = 16.22,
                   density = 4.268,
                   site_density_sitenm2 = 3.0,
                   cec_sitesnm2 = 0.0)
    minerals.append(goethite)

    print("Original data :\n")
    for m in minerals:
        print(m)

    print("\nConverted data :")
    print("\n-- Mineral SSA")
    for m in minerals:
        print("  {0} = {1} [cm^2/cm^3 bulk]".format(m.name, m.ssa_cm2cm3()))

    print("\n-- Surface Complexation Site Density")
    for m in minerals:
        print("  {0} = {1} [moles/m^3 bulk]".format(m.name, m.site_density_mm3b()))

    print("\n-- Cation Exchange Capacity")
    for m in minerals:
        print("  {0} = {1} [moles/m^3 bulk]".format(m.name, m.cec_mm3b()))


