# using this:
#
# 1. run this to generate a dot file
# 2. run through dot2tex:
#    


import networkx as nx
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)

g = nx.DiGraph()

_colors = dict(primary='#E8928E',
               secondary='#A3D4B8', #C5FFDE', #91BDA0',
               independent='#81C9E8E0',
               residual='#BD845D',
               )

_names = dict(pressure=r"p",
              temperature=r"T",
              effective_pressure=r"max(p,p_{atm})",
              surface_effective_pressure=r"max(p_s,p_{atm})",
              saturation=r"s",
              enthalpy=r"H",
              porosity=r"\phi",
              permeability=r"\kappa",
              density=r"\rho",
              mol_frac_gas=r"\omega",
              viscosity_liquid=r"\mu",
              unfrozen_fraction=r"\eta",
              unfrozen_fraction_relperm=r"k_r^{\eta}",
              unfrozen_effective_depth=r"h^{\eta}",
              water_content=r"\Theta",
              energy=r"E",
              internal_energy=r"u",
              air_temperature=r"T^{air}",
              base_porosity=r"\phi_0",
              capillary_pressure_gas_liq=r"p_c^{gl}",
              capillary_pressure_liq_ice=r"p_c^{li}",
              cell_volume=r"d\Omega",
              elevation=r"z",
              slope_magnitude=r"S",
              incoming_shortwave_radiation=r"Q^{sw-in}_E",
              manning_coefficient=r"m_n",
              overland_conductivity=r"\kappa_s",
              ponded_depth=r"h",
              surface_conducted_energy_source=r"Q^{cond}_E",
              surface_mass_source=r"Q",
              surface_total_energy_source=r"Q^{tot}_E",
              thermal_conductivity=r"\kappa^T",
              vapor_pressure=r"vp",
              wind_speed=r"U_{wind}",
              snow_depth=r"h_{snow}",
              relative_humidity=r"RH",
              precipitation_rain=r"P^r",
              precipitation_snow=r"P^s",
              pres_elev=r"h+z",
              relative_permeability=r"k_r",
              )
_names['3d_cell_volume_surface'] = None


reg = dict()

def get_name(key):
    try:
        return _names[key]
    except KeyError:
        pass

    if key.startswith("surface_"):
        return get_name(key[8:])+r"_s"
    if key.startswith("source_"):
        return get_name(key[7:])+r"^Q"
    if key.startswith("energy_source_"):
        return get_name(key[14:])+r"^{Q_E}"
    if key.startswith("mass_source_"):
        return get_name(key[12:])+r"^{Q_m}"

    if key.endswith("_gas"):
        return get_name(key[:-4])+r"^g"
    if key.endswith("_liquid"):
        return get_name(key[:-7])+r"^l"
    if key.endswith("_ice"):
        return get_name(key[:-4])+r"^i"
    if key.endswith("_rock"):
        return get_name(key[:-5])+r"^{rock}"
    if key.startswith("molar_density"):
        return get_name(key[6:])
    if key.startswith("mass_density"):
        return get_name(key[5:])
    if key.endswith("_bar"):
        return r"\bar{"+get_name(key[:-4])+r"}"
    return key

def texify(name):
    if name is None:
        return None
    return r"$"+name+r"$"

def get_tex_name(key):
    if key is None:
        return None
    return texify(get_name(key))

def add_node(g,key,ntype,lbl=None):
    if lbl is None:
        lbl = get_tex_name(key)
    reg[key] = lbl.replace("\\","")
    if reg[key] in g.nodes() and g.node[reg[key]].has_key('ntype'):
        pass
    else:
        print "adding node:", reg[key]
        g.add_node(reg[key],texlbl=lbl,ntype=ntype,color='black',fillcolor=_colors[ntype],shape='box',style='filled')


def add_edge(g,key1,key2):
    try:
        lbl1 = reg[key1]
    except KeyError:
        lbl1 = get_tex_name(key1)
        lbl1 = lbl1.replace("\\","")

    try:
        lbl2 = reg[key2]
    except KeyError:
        lbl2 = get_tex_name(key2)
        lbl2 = lbl2.replace("\\","")

    print "adding edge:", lbl1,lbl2
    g.add_edge(lbl1,lbl2,arrowsize='2')

with open("dependency_graph.txt",'r') as fid:
    reading = False
    for line in fid:
        if line == '\n':
            reading = False
        elif not reading:
            keylist = line.strip().strip(",")
            keys = keylist.split(",")
            reading = True
        else:
            if line.strip().startswith("Type:"):
                for key in keys:
                    typestring = " ".join(line.strip().split(" ")[1:])
                    # HACK
                    key = key.replace("mass_density","density")
                    key = key.replace("molar_density","density")
                    add_node(g,key,ntype=typestring)

            elif line.strip().startswith("Dep:"):
                for key in keys:
                    dep = " ".join(line.strip().split(" ")[1:])
                    # HACK
                    key = key.replace("mass_density","density")
                    key = key.replace("molar_density","density")
                    dep = dep.replace("mass_density","density")
                    dep = dep.replace("molar_density","density")

                    assert key is not None
                    assert dep is not None
                    add_edge(g,key, dep)

if reg.has_key("water_content"):
    #g.has_node(reg["water_content"]):
    add_node(g,"res_flow","residual",lbl=texify(r"res^{flow}"))
    add_node(g,"flux","secondary",lbl=texify(r"q"))
    add_edge(g,"res_flow", "water_content")
    add_edge(g,"res_flow", "cell_volume")
    add_edge(g,"res_flow", "flux")
    add_edge(g,"flux","permeability")
    add_edge(g,"flux","relative_permeability")
    add_edge(g,"flux","pressure")
    add_edge(g,"flux","density_liquid")
    add_edge(g,"flux","viscosity_liquid")

    if reg.has_key("surface_water_content"):
        #if g.has_node(reg["surface_water_content"]):
        add_node(g,"flux_surf_subsurf","secondary",lbl=texify(r"q_{ss}"))
        add_edge(g,"res_flow","flux_surf_subsurf")
        add_edge(g,"flux_surf_subsurf","pressure")
        add_edge(g,"flux_surf_subsurf","surface_pressure")
        add_edge(g,"flux_surf_subsurf","unfrozen_fraction_relperm")

if reg.has_key("surface_water_content"):
    #g.has_node(reg["surface_water_content"]):
    add_node(g,"res_surface_flow","residual",lbl=texify(r"res^{flow}_s"))
    add_node(g,"surface_flux","secondary",lbl=texify(r"q_s"))
    add_edge(g,"res_surface_flow", "surface_water_content")
    add_edge(g,"res_surface_flow", "surface_flux")
    add_edge(g,"res_surface_flow", "surface_cell_volume")
    add_edge(g,"res_surface_flow", "surface_mass_source")
    add_edge(g,"flux","overland_conductivity")
    add_edge(g,"flux","pres_elev")

    if reg.has_key("water_content"):
        #    if g.has_node(reg["water_content"]):
        add_edge(g,"res_surface_flow","flux_surf_subsurf")


if reg.has_key("energy"):
    #if g.has_node(reg["energy"]):
    add_node(g,"res_energy","residual",lbl=texify(r"res^{energy}"))
    add_node(g,"energy_flux","secondary",lbl=texify(r"q^E"))
    add_edge(g,"res_energy","energy")
    add_edge(g,"res_energy","flux")
    add_edge(g,"res_energy","enthalpy")
    add_edge(g,"res_energy","cell_volume")
    add_edge(g,"res_energy","energy_flux")
    add_edge(g,"energy_flux", "thermal_conductivity")
    add_edge(g,"energy_flux", "temperature")

    if reg.has_key("surface_energy"):
        #    if g.has_node(reg["surface_energy"]):
        add_node(g,"energy_flux_surf_subsurf","secondary",lbl=texify(r"q^E_{ss}"))
        add_edge(g,"res_energy","energy_flux_surf_subsurf")
        add_edge(g,"energy_flux_surf_subsurf","temperature")
        add_edge(g,"energy_flux_surf_subsurf","surface_temperature")

if reg.has_key("surface_energy"):
    #if g.has_node(reg["surface_energy"]):
    add_node(g,"res_surface_energy","residual",lbl=texify(r"res^{energy}_s"))
    add_node(g,"surface_energy_flux","secondary",lbl=texify(r"q^E_s"))
    add_edge(g,"res_surface_energy","surface_energy")
    add_edge(g,"res_surface_energy","surface_flux")
    add_edge(g,"res_surface_energy","surface_enthalpy")
    add_edge(g,"res_surface_energy","surface_cell_volume")
    add_edge(g,"res_surface_energy","surface_energy_flux")
    add_edge(g,"res_surface_energy","surface_total_energy_source")
    add_edge(g,"surface_energy_flux", "surface_thermal_conductivity")
    add_edge(g,"surface_energy_flux", "surface_temperature")

    if reg.has_key("energy"):
        #    if g.has_node(reg["energy"]):
        add_edge(g,"res_surface_energy","energy_flux_surf_subsurf")


# create the AGraph
a = nx.to_agraph(g)
a.layout('dot')

# set res rank
res_nodes = [node for node in g.nodes() if g.node[node]['ntype'] == 'residual']
res_subgraph = a.add_subgraph(res_nodes)
res_subgraph.graph_attr['rank'] = 'same'
print res_subgraph

# set p,T rank
hydro_nodes = []
if reg["temperature"] in g.nodes():
    hydro_nodes.append(reg["temperature"])
if reg["surface_temperature"] in g.nodes():
    hydro_nodes.append(reg["surface_temperature"])
if reg["pressure"] in g.nodes():
    hydro_nodes.append(reg["pressure"])
if reg["surface_pressure"] in g.nodes():
    hydro_nodes.append(reg["surface_pressure"])
hydro_subgraph = a.add_subgraph(hydro_nodes)
hydro_subgraph.graph_attr['rank'] = 'same'
print hydro_subgraph

a.write("dep.dot")

# plt.figure()
# A = nx.draw_graphviz(g, 'dot')
# plt.show()
