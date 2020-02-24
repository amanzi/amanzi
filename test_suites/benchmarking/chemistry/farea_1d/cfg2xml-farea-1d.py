"""
This script is a one-off attempt to generate the xml for the ASCEM
2012 F-Area 1-D full chemistry benchmarking problem. It isn't intended
to be general, maintainable or even useable by anyone.... If you're
thinking of copying or using this for another purpose, don't.

Notes:

 - Convert the geochemistry constrains from batch_chem cfg format into
   a 1-D amanzi-u problem and write the xml file. We start with two
   cfg files, one for the IC and one for the upstream BC.

 - We create a bunch of string templates containing the boilerplate xml

 - Use template substitution to insert values extracted from the batch_chem files. 

 - Concatenate a bunch of templated strings in the correct order to
   get a string with a complete xml input file in it. Write it to the
   file.

"""

import sys
import ConfigParser

from string import Template

boilerplate_main_xml = Template(
"""
<ParameterList name="Main">
  <Parameter name="Amanzi Input format version" type="string" value="1.0.0"/>
  <ParameterList name="General Description">
    <Parameter name="Model ID" type="string" value="$description"/>
  </ParameterList>

  <ParameterList name="Execution Control">
    <Parameter name="Flow Model" type="string" value="Steady State Saturated"/>
    <Parameter name="Transport Model" type="string" value="On"/>
    <Parameter name="Chemistry Model" type="string" value="On"/>
    <ParameterList name="Time Integration Mode">
      <ParameterList name="Initialize To Steady">
        <Parameter name="Start" type="double" value="0.0"/>
        <Parameter name="End" type="double" value="1.5778463e9"/>
        <Parameter name="Switch" type="double" value="0.0"/>
	<Parameter name="Steady Initial Time Step" type="double" value="1.5768e+5"/>
        <Parameter name="Transient Initial Time Step" type="double" value="1.5768e+5"/>
      </ParameterList>
    </ParameterList>
    <Parameter name="Verbosity" type="string" value="High"/>
    <ParameterList name="Numerical Control Parameters">
      <ParameterList name="Unstructured Algorithm">
        <Parameter name="Transport Integration Algorithm" type="string" value="Explicit First-Order"/>
        <Parameter name="steady max iterations" type="int" value="15"/>
        <Parameter name="steady min iterations" type="int" value="10"/>
        <Parameter name="steady limit iterations" type="int" value="20"/>
        <Parameter name="steady nonlinear tolerance" type="double" value="1e-12"/>
        <Parameter name="steady max time step" type="double" value="6.0e+10"/>
        <Parameter name="steady max preconditioner lag iterations" type="int" value="4"/>
        <Parameter name="steady time step reduction factor" type="double" value="0.8"/>
        <Parameter name="steady time step increase factor" type="double" value="1.25"/>
        <Parameter name="transient max iterations" type="int" value="15"/>
        <Parameter name="transient min iterations" type="int" value="10"/>
        <Parameter name="transient limit iterations" type="int" value="20"/>
        <Parameter name="transient nonlinear tolerance" type="double" value="1e-12"/>
        <Parameter name="transient max time step" type="double" value="6.0e+10"/>
        <Parameter name="transient max preconditioner lag iterations" type="int" value="4"/>
        <Parameter name="transient time step reduction factor" type="double" value="0.8"/>
        <Parameter name="transient time step increase factor" type="double" value="1.25"/>
        <Parameter name="linear solver tolerance" type="double" value="1e-22"/>
      </ParameterList>
    </ParameterList> <!-- Numerical Control Parameters -->
  </ParameterList> <!-- Execution Control Parameters -->
""")

domain_xml = \
"""
  <ParameterList name="Domain">
    <Parameter name="Spatial Dimension" type="int" value="3"/>
  </ParameterList>
  <ParameterList name="Mesh">
    <ParameterList name="Unstructured">
      <ParameterList name="Generate Mesh">
        <ParameterList name="Uniform Structured">
          <Parameter name="Number of Cells" type="Array int" value="{100, 1, 1}"/>
          <Parameter name="Domain Low Corner" type="Array double" value="{0.0, 0.0,  0.0}" />
          <Parameter name="Domain High Corner" type="Array double" value="{100.0, 1.0, 1.0}" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList> <!-- Mesh -->

  <ParameterList name="Regions">
    <ParameterList name="Entire Domain">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{0.0,0.0,0.0}"/>
        <Parameter name="High Coordinate" type="Array double" value="{100.0,1.0,1.0}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="West">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{0.0,0.0,0.0}"/>
        <Parameter name="High Coordinate" type="Array double" value="{0.0,1.0,1.0}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="East">
      <ParameterList name="Region: Box">
        <Parameter name="Low Coordinate" type="Array double" value="{100.0,0.0,0.0}"/>
        <Parameter name="High Coordinate" type="Array double" value="{100.0,1.0,1.0}"/>
      </ParameterList>
    </ParameterList>
  </ParameterList> <!-- Regions -->
"""

phase_definitions_template = Template(
"""
  <ParameterList name="Phase Definitions">
    <ParameterList name="Aqueous">
      <ParameterList name="Phase Properties">
        <ParameterList name="Viscosity: Uniform">
          <Parameter name="Viscosity" type="double" value="1.002e-3 "/>
        </ParameterList>
        <ParameterList name="Density: Uniform">
          <Parameter name="Density" type="double" value="998.2 "/>
        </ParameterList>
      </ParameterList>
      $aqueous_phase_solutes
    </ParameterList>
    $solid_phase
  </ParameterList> <!-- Phase Definitions -->
""")

phase_components_template = Template(
"""
      <ParameterList name="Phase Components">
        <ParameterList name="Water">
          <Parameter name="Component Solutes" type="Array string" value="{$solutes}"/>
        </ParameterList>
      </ParameterList>
""")

solid_phase_template = Template(
"""
    <ParameterList name="Solid">
      <Parameter name="Minerals" type="Array string" value="{$minerals}"/>
      <Parameter name="Sorption Sites" type="Array string" value="{$sorption_sites}"/>
    </ParameterList>
""")


material_properties_template = Template(
"""
  <ParameterList name="Material Properties">
    <ParameterList name="Soil">
      <Parameter name="Assigned Regions" type="Array string" value="{Entire Domain}"/>
      <ParameterList name="Porosity: Uniform">
        <Parameter name="Value" type="double" value="0.25"/>
      </ParameterList>
      <ParameterList name="Intrinsic Permeability: Uniform">
        <Parameter name="Value" type="double" value="1.E-12"/>
      </ParameterList>
$mineralogy
$cation_exchange_capacity
$surface_complexation_sites
    </ParameterList>
  </ParameterList> <!-- Material Properties -->
""")

mineralogy_template = Template(
"""
      <ParameterList name="Mineralogy">
$mineralogy_list
      </ParameterList> <!-- Mineralogy -->
""")

mineral_template = Template(
"""
	<ParameterList name="$name">
	  <Parameter name="Volume Fraction" type="double" value="$volume_fraction"/>
	  <Parameter name="Specific Surface Area" type="double" value="$specific_surface_area"/>
	</ParameterList>
""")

cation_exchange_capacity_template = Template(
"""
      <Parameter name="Cation Exchange Capacity" type="double" value="$value"/>
""")

surface_complexation_sites_template = Template(
"""
      <ParameterList name="Surface Complexation Sites">
$site_list
      </ParameterList> <!-- Surface Complexation Sites -->
""")

surface_site_template = Template(
"""
	<ParameterList name="$name">
	  <Parameter name="Site Density" type="double" value="$value"/>
	</ParameterList>
""")

solute_ic_species_template = Template(
"""
            <ParameterList name="$name">
              <ParameterList name="IC: Uniform Concentration">
                <Parameter name="Value" type="double" value="$value"/>
                <Parameter name="Free Ion Guess" type="double" value="$free_guess"/>
              </ParameterList>
              <Parameter name="Concentration Units" type="string" value="Molarity"/>
            </ParameterList>
""")

solute_ic_template = Template(
"""
      <ParameterList name="Solute IC">
        <ParameterList name="Aqueous">
          <ParameterList name="Water">
$data
          </ParameterList>
        </ParameterList>
      </ParameterList>
""")

initial_conditions_template = Template(
"""
  <ParameterList name="Initial Conditions">
    <ParameterList name="Initial Condition">
      <Parameter name="Assigned Regions" type="Array string" value="{Entire Domain}"/>
      <ParameterList name="IC: Uniform Pressure">
        <Parameter name="Phase" type="string" value="Aqueous"/>
        <Parameter name="Value" type="double" value="201325.0"/>
      </ParameterList>
$solute_ic
    </ParameterList>
  </ParameterList> <!-- Initial Conditions -->
""")

bc_template = Template(
"""
  <ParameterList name="Boundary Conditions">
    <ParameterList name="West BC">
      <Parameter name="Assigned Regions" type="Array string" value="{West}"/>
      <ParameterList name="BC: Flux">
        <Parameter name="Times" type="Array double" value="{0.0, 1.5778463e9}"/>
        <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
        <Parameter name="Inward Mass Flux" type="Array double" value="{7.91317859e-6, 7.91317859e-6}"/>
<!--        <Parameter name="Inward Volumetric Flux" type="Array double" value="{7.927447996e-9,7.927447996e-9}"/>	 		-->
      </ParameterList>
$west_bc
    </ParameterList> <!-- West BC -->

    <ParameterList name="East BC">
      <Parameter name="Assigned Regions" type="Array string" value="{East}"/>
      <ParameterList name="BC: Uniform Pressure">
        <Parameter name="Times" type="Array double" value="{0.0, 1.5778463e9}"/>
        <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
	<Parameter name="Values" type="Array double" value="{201325.0, 201325.0}"/>	 		
      </ParameterList>
$east_bc
    </ParameterList> <!-- East BC -->
  </ParameterList> <!-- Boundary Conditions -->
""")

solute_bc_uniform_conc_template = Template(
"""
            <ParameterList name="$name">
              <ParameterList name="BC: Uniform Concentration">
                <Parameter name="Times" type="Array double" value="{0.0, 1.5778463e9}"/>
                <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
                <Parameter name="Values" type="Array double" value="{$value, $value}"/>
              </ParameterList>
              <Parameter name="Concentration Units" type="string" value="Molar Concentration"/>
            </ParameterList>
""")

solute_bc_zero_gradient_template = Template(
"""
            <ParameterList name="$name">
              <ParameterList name="BC: Zero Gradient">
                <Parameter name="Times" type="Array double" value="{0.0, 1.5778463e9}"/>
                <Parameter name="Time Functions" type="Array string" value="{Constant}"/>
                <Parameter name="Values" type="Array double" value="{$value, $value}"/>
              </ParameterList>
              <Parameter name="Concentration Units" type="string" value="Molar Concentration"/>
            </ParameterList>
""")

solute_bc_template = Template(
"""
      <ParameterList name="Solute BC">
        <ParameterList name="Aqueous">
          <ParameterList name="Water">
$data
          </ParameterList>
        </ParameterList>
      </ParameterList> <!-- Solute BC -->
""")

output_template = Template(
"""
  <ParameterList name="Output">
    <ParameterList name="Time Macros">
      <ParameterList name="Every_0.05_year">
	<Parameter name="Start_Period_Stop" type="Array double" value="{0., 1.5768e6, -1}"/>
      </ParameterList>
      <ParameterList name="Every_year">
	<Parameter name="Start_Period_Stop" type="Array double" value="{0.0, 31556926.0, -1}"/>
      </ParameterList>
      <ParameterList name="Every_0.1_year">
	<Parameter name="Start_Period_Stop" type="Array double" value="{0.0, 3155692.6, -1}"/>
      </ParameterList>
    </ParameterList>
    <ParameterList name="Visualization Data">
      <Parameter name="File Name Base" type="string" value="farea-1d"/>
      <Parameter name="File Name Digits" type="string" value="5"/>
      <Parameter name="Time Macro" type="Array string" value="{Every_year}"/>
<!--    <Parameter name="Time Macro" type="Array string" value="{Every_0.1_year}"/> -->
<!-- 	<Parameter name="Time Macro" type="string" value="Every_0.05_year"/> -->
    </ParameterList>
  </ParameterList> <!-- Output -->
""")

chemistry_template = Template(
"""
  <ParameterList name="Chemistry">
    <Parameter name="Thermodynamic Database Format" type="string" value="simple" />
    <Parameter name="Thermodynamic Database File" type="string" value="$bgd_file" />
    <Parameter name="Verbosity" type="Array string" value="{verbose}" />
    <Parameter name="Activity Model" type="string" value="debye-huckel" />
    <Parameter name="Tolerance" type="double" value="1.5e-12"/>
    <Parameter name="Maximum Newton Iterations" type="int" value="250"/>
    <Parameter name="Auxiliary Data" type="Array string" value="{pH}"/>
    <!-- 1.0 yr -->
    <Parameter name="Max Time Step (s)" type="double" value="31556926.0"/>
    <!-- 0.1 yr -->
<!--    <Parameter name="Max Time Step (s)" type="double" value="3155692.60"/> -->
    <!-- 0.05 yr -->
<!--    <Parameter name="Max Time Step (s)" type="double" value="1577846.3"/> -->
    <!-- 0.01 yr -->
<!--    <Parameter name="Max Time Step (s)" type="double" value="315569.260"/> -->
  </ParameterList> <!-- Chemistry -->
""")

finish_template = Template(
"""
</ParameterList> <!-- Main -->
""")

def initial_condition_cfg2xml(ic_cfg_file, species_names):
    """
    Extract the total list from the cfg file and write it as an amanzi
    official input spec 'Solute IC' xml section
    """
    config = ConfigParser.SafeConfigParser()
    config.read(ic_cfg_file)

    if config.has_section('total'):
        totals = config.items('total')
    else:
        raise Exception('Error: no "total" section in {0}'.format(ic_cfg_file))

    if config.has_section('free_ion'):
        free_ion = config.items('free_ion')
    else:
        raise Exception('Error: no "free_ion" section in {0}'.format(ic_cfg_file))

    initial_conditions_xml = ""
    for t,f in zip(totals, free_ion):
        if t[0] != f[0]:
            raise AssertionError('Error: total and free_ion lists are in different order!')
        species = t[0]
        value = t[1]
        free_guess = f[1]
        name = species_names[species]
        ic = solute_ic_species_template.substitute(name=name, value=value, free_guess=free_guess)
        initial_conditions_xml += ic

    solute_ic = solute_ic_template.substitute(data=initial_conditions_xml)
    ic_xml = initial_conditions_template.substitute(solute_ic=solute_ic)
    return ic_xml

def boundary_conditions(bc_cfg_file, species_names):
    uniform_conc_bc = boundary_condition_uniform_conc_cfg2xml(bc_cfg_file, species_names)
    zero_gradient_bc = boundary_condition_zero_gradient_cfg2xml(bc_cfg_file, species_names)
    bc_xml = bc_template.substitute(west_bc=uniform_conc_bc,
                                    east_bc=zero_gradient_bc)
    return bc_xml

def boundary_condition_uniform_conc_cfg2xml(bc_cfg_file, species_names):
    config = ConfigParser.SafeConfigParser()
    config.read(bc_cfg_file)

    if config.has_section('total'):
        totals = config.items('total')
    else:
        raise Exception('Error: no "total" section in {0}'.format(bc_cfg_file))

    bc_uniform_xml = ""
    for t in totals:
        species = t[0]
        value = t[1]
        name = species_names[species]
        bc_uniform_xml += solute_bc_uniform_conc_template.substitute(name=name, value=value)

    solute_bc = solute_bc_template.substitute(data=bc_uniform_xml)

    return solute_bc

def boundary_condition_zero_gradient_cfg2xml(bc_cfg_file, species_names):
    config = ConfigParser.SafeConfigParser()
    config.read(bc_cfg_file)

    if config.has_section('total'):
        totals = config.items('total')
    else:
        raise Exception('Error: no "total" section in {0}'.format(bc_cfg_file))

    bc_zero_grad_xml = ""
    for t in totals:
        species = t[0]
        value = 1.0e-40
        name = species_names[species]
        bc_zero_grad_xml += solute_bc_zero_gradient_template.substitute(name=name, value=value)

    solute_bc = solute_bc_template.substitute(data=bc_zero_grad_xml)

    return solute_bc

def phase_definitions(species_names, mineral_names, surface_site_names):
    aqueous_phase_solutes_xml = phase_components(species_names)
    solid_phase_xml = solid_phase(mineral_names, surface_site_names)
    phase_definitions_xml = phase_definitions_template.substitute(
        aqueous_phase_solutes = aqueous_phase_solutes_xml,
        solid_phase = solid_phase_xml)
    return phase_definitions_xml

def phase_components(species_names):
    solutes = ''
    i = 0
    for s in species_names:
        solutes += s
        i += 1
        if i != len(species_names):
            solutes += ', '

    aqueous_phase_solutes = phase_components_template.substitute(solutes=solutes)
    return aqueous_phase_solutes

def solid_phase(mineral_names, surface_site_names):
    minerals = ''
    i = 0
    for m in mineral_names:
        minerals += m
        i += 1
        if i != len(mineral_names):
            minerals += ', '

    site_names = ''
    i = 0
    for s in surface_site_names:
        site_names += s
        i += 1
        if i != len(surface_site_names):
            site_names += ', '
    solid_phase_section = solid_phase_template.substitute(minerals=minerals,
                                                          sorption_sites=site_names)
    return solid_phase_section

def material_properties(ic_cfg_file, bgd_file, mineral_names_dict, mineral_ssa):
    mineralogy_xml = mineralogy(ic_cfg_file, bgd_file, mineral_names_dict, mineral_ssa)
    cec_xml = cec(ic_cfg_file)
    surface_complexation_xml = surface_complexation(bgd_file)
    return material_properties_template.substitute(mineralogy=mineralogy_xml,
                                                   cation_exchange_capacity=cec_xml,
                                                   surface_complexation_sites=surface_complexation_xml)


def mineralogy(ic_cfg_file, bgd_file, mineral_names_dict, mineral_ssa):
    config = ConfigParser.SafeConfigParser()
    config.read(ic_cfg_file)

    if not config.has_section('mineral'):
        raise Exception('Error: no "mineral" section in {0}'.format(ic_cfg_file))

    mineral_vf = {}
    for m in mineral_names_dict.keys():
        name = mineral_names_dict[m]
        volume_fraction = None
        if config.has_option('mineral', m):
            volume_fraction = config.getfloat('mineral', m)
        else:
            raise Exception('Error: mineral "{0}" is not in mineral section.'.format(m))
        mineral_vf[name] = volume_fraction
    # bgd files just contains ssa = 1.0 for everything. Use the hard
    # coded values instead
    #mineral_ssa = get_mineral_ssa_from_bgd(bgd_file)
    minerals_xml = ''
    for m in mineral_vf.keys():
        name = m
        volume_fraction = mineral_vf[m]
        specific_surface_area = mineral_ssa[m]
        minerals_xml += mineral_template.substitute(name=name,
                                                  volume_fraction=volume_fraction,
                                                  specific_surface_area=specific_surface_area)
    mineralogy_xml = mineralogy_template.substitute(mineralogy_list=minerals_xml)
    return mineralogy_xml

def get_mineral_ssa_from_bgd(bgd_file):
    mineral_ssa = {}
    with open(bgd_file, 'r') as bgdfile:
        minerals_section = False
        for line in bgdfile:
            if line[0] != '#' and len(line) != 1:
                if line == '<Minerals\n':
                    minerals_section = True
                elif line == '<Mineral Kinetics\n':
                    minerals_section = False
                else:
                    if minerals_section:
                        data = line.split(';')
                        rxn = data[0].split('=')
                        name = rxn[0].strip()
                        ssa = float(data[-1])
                        mineral_ssa[name] = ssa
            
    return mineral_ssa

def cec(ic_cfg_file):
    config = ConfigParser.SafeConfigParser()
    config.read(ic_cfg_file)

    if not config.has_section('ion_exchange'):
        raise Exception('Error: no "ion_exchange" section in {0}'.format(ic_cfg_file))
    if config.has_option('ion_exchange', 'X-'):
        value = config.getfloat('ion_exchange', 'X-')
        cec_xml = cation_exchange_capacity_template.substitute(value=value)
    else:
        cec_xml = ""
    return cec_xml

def surface_complexation(bgd_file):
    surface_sites = {}
    with open(bgd_file, 'r') as bgdfile:
        site_section = False
        for line in bgdfile:
            if line[0] != '#' and len(line) != 1:
                if line == '<Surface Complex Sites\n':
                    site_section = True
                elif line == '<Surface Complexes\n':
                    site_section = False
                else:
                    if site_section:
                        data = line.split(';')
                        name = data[0].strip()
                        density = float(data[1])
                        surface_sites[name] = density

    surface_site_xml = ''
    for s in surface_sites.keys():
        surface_site_xml += surface_site_template.substitute(name=s,
                                                             value=surface_sites[s])

    surface_complexation_xml = surface_complexation_sites_template.substitute(site_list=surface_site_xml)
    return surface_complexation_xml

if __name__ == "__main__":
    bgd_file = 'farea-full.bgd'
    ic_cfg_file = "farea-background.cfg"
    bc_cfg_file = "farea-seepage.cfg"
    output_file = "amanzi-u-1d-farea-full.xml"
    output_background_bc = 'farea-full-background-bc.xml'

    print("Convert batch_chem files into a 1-D amanzi-u problem")
    print("  Input files :")
    print("    database : {0}".format(bgd_file))
    print("    initial conditions : {0}".format(ic_cfg_file))
    print("    boundary conditions : {0}".format(bc_cfg_file))

    print("  Processing batch_chem files to generate xml...")

    # config parser changes upper case to lower case and we can't do a
    # nieve conversion back, so we need a translator.
    species_names_dict = {
        'h+' : 'H+',
        'al+++' : 'Al+++',
        'ca++' : 'Ca++',
        'cl-' : 'Cl-',
        'fe+++' : 'Fe+++',
        'co2(aq)' : 'CO2(aq)',
        'k+' : 'K+',
        'mg++' : 'Mg++',
        'na+' : 'Na+',
        'sio2(aq)' : 'SiO2(aq)',
        'so4--' : 'SO4--',
        'tritium' : 'Tritium',
        'no3-' : 'NO3-',
        'uo2++' : 'UO2++'
        }
    # Dict doesn't preserve order, so we hard code it...
    species_names = [
        'H+',
        'Al+++',
        'Ca++',
        'Cl-',
        'Fe+++',
        'CO2(aq)',
        'K+',
        'Mg++',
        'Na+',
        'SiO2(aq)',
        'SO4--',
        'Tritium',
        'NO3-',
        'UO2++',
        ]

    mineral_names_dict = {
        'kaolinite' : 'Kaolinite',
        'quartz' : 'Quartz',
        'goethite' : 'Goethite',
        'schoepite' : 'Schoepite',
        'gibbsite' : 'Gibbsite',
        'jurbanite' : 'Jurbanite',
        'basaluminite' : 'Basaluminite',
        'opal' : 'Opal',
        }

    mineral_names = [
        'Quartz',
        'Goethite',
        'Kaolinite',

        'Schoepite',
        'Gibbsite',
        'Jurbanite',
        'Basaluminite',
        'Opal',
        ]

    # pflotran doesn't write the mineral specific surface area for us,
    # so we need to hard code it.
    mineral_ssa = {'Quartz' : 3262.3,
                   'Kaolinite' : 59093.9,
                   'Goethite' : 11076.3,
                   'Schoepite' : 0.1,
                   'Gibbsite' : 0.1,
                   'Basaluminite' : 0.1,
                   'Opal' : 0.1,
                   'Jurbanite' : 0.1
                   }

    surface_site_names = ['>davis_OH']
    # surface_site_names = ['>dong_kOH',
    #                       '>dong_FeOH',
    #                       '>SiOH']


    complete_xml = boilerplate_main_xml.substitute(description='1-D advection with ASCEM 2012 F-Area full chemistry')
    complete_xml += domain_xml
    phase_definitions_xml = phase_definitions(species_names, mineral_names, surface_site_names)
    complete_xml += phase_definitions_xml

    material_properties_xml = material_properties(ic_cfg_file, bgd_file,
                                                  mineral_names_dict, mineral_ssa)
    complete_xml += material_properties_xml

    ic_xml = initial_condition_cfg2xml(ic_cfg_file, species_names_dict)
    complete_xml += ic_xml

    complete_xml += boundary_conditions(bc_cfg_file, species_names_dict)
    complete_xml += output_template.substitute()
    complete_xml += chemistry_template.substitute(bgd_file=bgd_file)
    complete_xml += finish_template.substitute()

    #print complete_xml

    print('  Writing full xml to "{0}"'.format(output_file))

    with open(output_file, 'w') as output:
        output.write(complete_xml)
    print('  Finished.')
    # Now we need to dump an additional file with the background
    # chemistry used as BC instead of IC.
    print('  Creating secondary xml with the background constraint as BC...')
    secondary_xml = ''
    secondary_xml += boundary_condition_uniform_conc_cfg2xml(ic_cfg_file, species_names_dict)

    print('  Writing secondary xml to {0}'.format(output_background_bc))
    with open(output_background_bc, 'w') as output:
        output.write(secondary_xml)
    print('  Finished.')
    print('Done!')
    exit()


