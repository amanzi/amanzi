import xml.etree.ElementTree as ET
import sys

#
#  Update input files from 2.1.0 to 2.1.1
#
#  Transformations:
#
#  - Update the version number from 2.1.0 to 2.1.1
#  - Materials: 
#            Update dispersion model option: set type to uniform_isotropic, burnett_frind, or lichtner_kelkar_robinson
#  - Numerical Controls: 
#            Add subelement unstructured_controls or structured_controls and move options into subelemnt
#            If unstructured: 
#                  prepend subelement names with "unstr_"
#                  move cfl from unstr_linear_solver to unstr_transport_controls
#                  consolidate preconditioner options under unstr_preconditioners
#                  rename pseudo_time_integrator to unstr_initialization, 
#                        Note: empty element means use defaults, no element means do not use
#                  update BDF1 options from attributes to subelements
#                  rename restart_tolerance_factor to restart_tolerance_relaxation_factor
#            If structured: prepend subelement names with "str_"
#  - Process Kernel: 
#            Flow: move options rel_perm_method, preconditioning_strategy, discretization_method, and atmospheric_pressure 
#                  to numerical_controls->unstructured_controls->unstr_flow_controls
#            Transport: move options algorithm and sub_cycling 
#                  to numerical_controls->unstructured_controls->unstr_transport_controls
#  - Output: 
#            Vis: update time_macro to time_macros and cycle_macro to cycle_macros
#            Observations: update time_macro to time_macros 
#            Checkpoint: update cycle_macro to cycle_macros
#            Walkabout: update cycle_macro to cycle_macros
#


def usage():

    print('')
    print('  UpdateSpec_210-211 Usage:')
    print('    python update210.py oldfile.xml (newfile.xml)')
    print('')
    print('  UpdateSpec_210-211:')
    print('    - reads a 2.1.0 amanzi input file and updates it to comform to the 2.1.1 input schema')
    print('    - if the optional argument newfile.xml is not specified, the updated xml will be written out to the oldfile.xml specified')
    print('    - it is recommended not to overwrite the oldfile.xml so that updates can be compared')
    print('    - to compare the old and new files use: diff -Eb or diff -w ')
    print('')

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def read_input_tree(inputfile):

    try:
      tree = ET.parse(inputfile)
    except:
      sys.stderr.write("Error reading inputfile - %s\n" % str(inputfile))
      usage()
      sys.exit("   exiting...")

    return tree

def report(tree, outputfile):

    root = tree.getroot()
    indent(root, 1)
    tree.write(outputfile)

def v210_update(tree):

    # try reading in inputfile
    try:
      root = tree.getroot()
    except:
      sys.stderr.write("Error getting xml tree root\n")
      usage()
      sys.exit("   exiting...")

    # check version for 2.1.0
    version = root.get('version')
    if version == '2.0.0':
      #call v200_to_210 before continuing
      sys.stderr.write("Error reading input file, can not currently update from 2.0.0 input format\n")
      sys.exit("  exiting...")
    elif (version == '1.2.3' or version == '1.2.2'):
      sys.stderr.write("Error reading input file, can not currently update from 1.2.x input format\n")
      sys.exit("  exiting...")

    # update version number
    ### EIB >> update this
    root.set('version','2.1.1')

    print("")
    print("  Beginning update from old 2.1.0 to 2.1.1")
    print("")

    # continue to updating format

    # check for hydrostatic BCs and add submodel option -> didn't exist before, can't update
    # check for dispersion_tensor and update 
    mats = root.find('./materials')
    for child in mats:
        disp = child.find('mechanical_properties/dispersion_tensor')
        if (disp is not None):
            type = disp.get("type")
            if (type is None):
                print("    Adding attribute 'type' for dispersion_tensor, material ",child.get("name"))
                type = "uniform_isotropic"
                alpha = disp.get("alpha_lh")
                if (alpha is not None):
                    type = "burnett_frind"
                alpha = disp.get("alpha_th")
                if (alpha is not None):
                    type = "burnett_frind"
                alpha = disp.get("alpha_lh")
                if (alpha is not None):
                    type = "lichtner_kelkar_robinson"
                disp.set("type",type)

    ### Process Kernel numerical controls

    # get simulation type
    type =  root.get('type')

    # check for new unstructured_controls/structured_controls
    # check subelements for unstr_/str_ 
    if (type == "unstructured"):
        numctrls = root.find('./numerical_controls')
        unstr = root.find('numerical_controls/unstructured_controls')
        if (unstr is None):
            # add unstructured_controls elements
            print("    Grouping unstructured numerical controls under 'unstructured_controls' element tag")
            unstr = ET.Element('unstructured_controls')
            # move control elements under new unstructured_controls element
            skip_names = ['comments','common_controls']
            moved_names = []
            for child in numctrls:
                if child.tag not in skip_names:
                    unstr.append(child)
                    moved_names.append(child.tag)
            for tag in moved_names:
                child = numctrls.find(tag)
                if (child is not None):
                    numctrls.remove(child)
            numctrls.append(unstr)

        # check subelements for unstr_
        for child in unstr:
            if ('unstr' not in child.tag and child.tag != 'comments'):
                print("    Adding 'ustr_' prefix to unstructured_controls subelement ",child.tag)
                child.tag = "unstr_"+child.tag
            if ("steady" in child.tag):
                for gkid in child:
                    if ("pseudo" in gkid.tag):
                        if ("unstr" not in gkid.tag):
                            gkid.tag = "unstr_"+gkid.tag

    elif (type == "structured"):
        struct = root.find('./numerical_controls/structured_controls')
        if (struct is None):
            print("    Grouping structured numerical controls under 'structured_contorls' element tag")
            # add unstructured_controls elements
            numctrls = (root.find('./numerical_controls'))
            struct = ET.SubElement(numctrls,'structured_controls')
            # move control elements under new unstructured_controls element
            for child in numctrls:
                struct.append(child)
                numctrls.remove(child)
        # check subelements for str_
        prefix_list = ['steady-state_controls','transient_controls','amr_controls']
        for child in struct:
            if ( child in prefix_list):
                print("    Adding 'str_' prefix to unstructured_controls subelement ",child.tag)
                child.tag = "str_"+child.tag

    else:
        sys.stderr.write("  ERROR: root element 'amanzi_input' has attribute 'type'\n")
        sys.stderr.write("         'type' value = %s\n" % str(type))
        sys.stderr.write("          validate values are 'unstructured' or 'structured'\n")
        sys.exit("  exiting...")

    
    # the following only apply to unstructured
    if (type == "unstructured"):

      # check for rel_perm_method, precondition_strategy, and discretization_method in old location and move
      unstr_cntl = root.find('./numerical_controls/unstructured_controls')
      flow = root.find('./process_kernels/flow')

      fields = ['rel_perm_method', 'preconditioning_strategy', 'discretization_method','atmospheric_pressure']
      moved_options = []
      for option in fields:
        if option in flow.keys():
          value = flow.get(option)
          if (unstr_cntl.find('unstr_flow_controls') == None):
            fc = ET.SubElement(unstr_cntl,'unstr_flow_controls')
          else:
            fc = unstr_cntl.find('unstr_flow_controls')
          moved_options.append(option)
          newoption = ET.SubElement(fc,option)
          newoption.text = value
      for option in moved_options:
          print("    Moved flow option ",option," to new location un 'unstr_flow_controls'")
          del flow.attrib[option]

      # check for algorithm and sub_cycling in old location and move
      trans = root.find('./process_kernels/transport')
  
      fields = ['algorithm', 'sub_cycling']
      moved_options = []
      for option in fields:
        if option in trans.keys():
          value = trans.get(option)
          if (unstr_cntl.find('unstr_transport_controls') == None):
            tc = ET.SubElement(unstr_cntl,'unstr_transport_controls')
          else:
            tc = unstr_cntl.find('unstr_transport_controls')
          moved_options.append(option)
          newoption = ET.SubElement(tc,option)
          newoption.text = value
      for option in moved_options:
          print("    Moved transport option ",option," to new location in 'unstr_transport_controls'")
          del trans.attrib[option]

      # check for cfl in old location and move
      trans = root.find('./numerical_controls/unstructured_controls/unstr_transport_controls')
      linear = root.find('./numerical_controls/unstructured_controls/unstr_linear_solver')
      if (linear is not None and trans is not None):
          cfl_old = linear.find('cfl')
          if (cfl_old is not None):
              cfl = ET.SubElement(trans,'cfl')
              cfl.text = cfl_old.text
              linear.remove(cfl_old)
              print("    Moved cfl to new location un 'unstr_transport_controls'")


      ### Preconditioners

      # check for preconditioner, check formating, saving any subelements, update formating
      precon_names = ['hypre_amg','trilinos_ml','block_ilu']

      old_precon = root.find('./numerical_controls/unstructured_controls/unstr_preconditioners')

      precon = ET.SubElement(unstr_cntl,'unstr_preconditioners')
      hypre = ET.SubElement(precon,'hypre_amg')
      trilinos = ET.SubElement(precon,'trilinos_ml')
      block = ET.SubElement(precon,'block_ilu')

      if (old_precon is not None):
          if ( old_precon.get('name') in precon_names):
              name = old_precon.get('name')
              # move any subelement options
              if (name == 'hypre_amg'):
                for op in list(old_precon):
                  hypre.append(op)
              if (name == 'trilinos_ml'):
                for op in list(old_precon):
                  trilinos.append(op)
              if (name == 'block_ilu'):
                for op in list(old_precon):
                  block.append(op)
              # remove existing preconditioner element
              unstr_cntl.remove(old_precon)
              print("    Updating preconditioner formating/definition")

      steady = root.find('./numerical_controls/unstructured_controls/unstr_steady-state_controls')
      if (steady is not None):
          pre = steady.find('preconditioner')
          if (pre is not None):
            if ( pre.get('name') in precon_names):
              print("    Updating preconditioner formating/definition in 'unstr_steady-state_controls'")
              name = pre.get('name')
              # move any subelement options
              if (name == 'hypre_amg'):
                for op in list(pre):
                  hypre.append(op)
              if (name == 'trilinos_ml'):
                for op in list(pre):
                  trilinos.append(op)
              if (name == 'block_ilu'):
                for op in list(pre):
                  block.append(op)
              # remove existing preconditioner element
              steady.remove(pre)
              # add new preconditioner element
              new_pre = ET.SubElement(steady,'preconditioner')
              new_pre.text = name

      pseudo = root.find('./numerical_controls/unstructured_controls/unstr_steady-state_controls/unstr_pseudo_time_integrator')
      # not everyone added the unstr_ prefix at this level, so let's catch it now
      if (pseudo is None):
          pseudo = root.find('./numerical_controls/unstructured_controls/unstr_steady-state_controls/pseudo_time_integrator')
          if (pseudo is not None):
              pseudo.tag = 'unstr_'+pseudo.tag
      if (pseudo is not None):
          pre = pseudo.find('preconditioner')
          if (pre is not None):
            if ( pre.get('name') in precon_names):
              print("    Updating preconditioner formating/definition in 'pseudo_time_integrator'")
              name = pre.get('name')
              # move any subelement options
              if (name == 'hypre_amg'):
                for op in list(pre):
                  hypre.append(op)
              if (name == 'trilinos_ml'):
                for op in list(pre):
                  trilinos.append(op)
              if (name == 'block_ilu'):
                for op in list(pre):
                  block.append(op)
              # remove existing preconditioner element
              pseudo.remove(pre)
              # add new preconditioner element
              new_pre = ET.SubElement(pseudo,'preconditioner')
              new_pre.text = name

      trans  = root.find('./numerical_controls/unstructured_controls/unstr_transient_controls')
      if (trans is not None):
          pre = trans.find('preconditioner')
          if (pre is not None):
            if ( pre.get('name') in precon_names):
              name = pre.get('name')
              print("    Updating preconditioner formating/definition in 'unstr_transient_controls'")
              # move any subelement options
              if (name == 'hypre_amg'):
                for op in list(pre):
                  hypre.append(op)
              if (name == 'trilinos_ml'):
                for op in list(pre):
                  trilinos.append(op)
              if (name == 'block_ilu'):
                for op in list(pre):
                  block.append(op)
              # remove existing preconditioner element
              trans.remove(pre)
              # add new preconditioner element
              new_pre = ET.SubElement(trans,'preconditioner')
              new_pre.text = name

      linear = root.find('./numerical_controls/unstructured_controls/unstr_linear_solver')
      if (linear is not None):
          pre = linear.find('preconditioner')
          if (pre is not None):
            if ( pre.get('name') in precon_names):
              name = pre.get('name')
              print("    Updating preconditioner formating/definition in 'unstr_linear_solver'")
              # move any subelement options
              if (name == 'hypre_amg'):
                for op in list(pre):
                  hypre.append(op)
              if (name == 'trilinos_ml'):
                for op in list(pre):
                  trilinos.append(op)
              if (name == 'block_ilu'):
                for op in list(pre):
                  block.append(op)
              # remove existing preconditioner element
              linear.remove(pre)
              # add new preconditioner element
              new_pre = ET.SubElement(linear,'preconditioner')
              new_pre.text = name


      ### Rename pseudo_time_integrator to initialization
      # added parameter clipping pressure, didn't exist so can't update
      pseudo = root.find('./numerical_controls/unstructured_controls/unstr_steady-state_controls/unstr_pseudo_time_integrator')
      if (pseudo is not None):
          pseudo.tag = 'unstr_initialization'
          print("    Updating 'unstr_pseudo_time_integrator' to 'unstr_initialization'")

      ### Update BDF1 attributes to elements
      bdf1 =  root.find('./numerical_controls/unstructured_controls/unstr_transient_controls/bdf1_integration_method')
      moved_list = []
      if (bdf1 is not None):
          for attr in bdf1.attrib:
              name = attr
              value = bdf1.get(attr)
              new_elem = ET.SubElement(bdf1,name)
              new_elem.text = value
              moved_list.append(name)
          if (len(moved_list) > 0):
              print("    Updating 'bdf1_integration_method' options from attributes to elements")
          for name in moved_list:
              del bdf1.attrib[name]

      ### Rename restart_tolerance_factor
      steady = root.find('./numerical_controls/unstructured_controls/unstr_steady-state_controls')
      if (steady is not None):
          for child in steady:
              if child.tag == 'restart_tolerance_factor':
                  child.tag = 'restart_tolerance_relaxation_factor'
                  print("    Renaming (steady) 'restart_tolerance_factor' to ",child.tag)
      bdf1 = root.find('./numerical_controls/unstructured_controls/unstr_transient_controls/bdf1_integration_method')
      if (bdf1 is not None):
          for child in bdf1:
              if child.tag == 'restart_tolerance_factor':
                  child.tag = 'restart_tolerance_relaxation_factor'
                  print("    Renaming (transient) 'restart_tolerance_factor' to ",child.tag)


    ### Output

    # check for observations and update time_macro to time_macros
    vis = root.find('./output/vis')

    if (vis != None):
      tm = vis.find('time_macro')
      if tm is not None:
        print("    Updating 'vis' time_macro to time_macros")
        tms = ET.SubElement(vis,'time_macros')
        tms.text = tm.text
        vis.remove(tm)
      cm = vis.find('cycle_macro')
      if cm is not None:
        print("    Updating 'vis' cycle_macro to cycle_macros")
        cms = ET.SubElement(vis,'cycle_macros')
        cms.text = cm.text
        vis.remove(cm)

    # check for observations and update time_macro to time_macros
    obs = root.find('./output/observations/liquid_phase')

    if (obs != None):
      for ob in obs:
        tm = ob.find('time_macro')
        if tm is not None:
          print("    Updating 'observations' time_macro to time_macros")
          tms = ET.SubElement(ob,'time_macros')
          tms.text = tm.text
          ob.remove(tm)

    # check for checkpoint and update cycle_macro to cycle_macros
    ckpt = root.find('./output/checkpoint')
    if (ckpt != None):
      cm = ckpt.find('cycle_macro')
      if cm is not None:
        print("    Updating 'checkpoint' cycle_macro to cycle_macros")
        cms = ET.SubElement(ckpt,'cycle_macros')
        cms.text = cm.text
        ckpt.remove(cm)

    # check for walkabout and update cycle_macro to cycle_macros
    walk = root.find('./output/walkabout')
    if (walk != None):
      cm = walk.find('cycle_macro')
      if cm is not None:
        print("    Updating 'walkabout' cycle_macro to cycle_macros")
        cms = ET.SubElement(walk,'cycle_macros')
        cms.text = cm.text
        walk.remove(cm)

    return tree


    
if __name__=="__main__":
  
    try:
      inputfile = sys.argv[1]
      try:
        outputfile = sys.argv[2]
      except:
        outputfile = inputfile
        print(">> no output filename specified, will overwrite inputfile")
    except:
      sys.stderr.write("Error reading arguments\n")
      usage()
      sys.exit("   exiting...")

    tree = read_input_tree(inputfile)
    new_tree = v210_update(tree)
    report(new_tree, outputfile)

    print("")
    print("  Use 'diff -w' or 'diff -Eb' to compare xml files")
    print("")
   
