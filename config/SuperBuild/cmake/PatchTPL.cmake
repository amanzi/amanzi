#[[
  patch_tpl(name out_dir source_dir timestamps_dir patch_files)

  Creates a sh script and cmake file that patches the source_dir for a TPL.

  Arguments:
    - name: name of the package
    - out_dir: where to put the resulting sh script
    - source_dir: source to patch
    - timestamps_dir: log from the script will be put here
    - patch_files: variable name containing a list of patch files

  Outputs:
    - ${name}_PATCH_COMMAND is set in the parent scope.

  Example (note, do not ${} patch files!):
    create_patch_sh_file(seacas
                         ${SEACAS_prefix_dir}
                         ${SEACAS_source_dir}
                         ${SEACAS_stamp_dir}
                         SEACAS_patch_files)
]]
function(patch_tpl name out_dir source_dir timestamps_dir patch_files)
  message(STATUS "Applying ${name} patches:")
  foreach(item IN LISTS ${patch_files})
    message(STATUS "  - ${item}")
  endforeach()

  # yes, this looks funny -- patch_files is a variable storing the
  # name of a list variable, so we dereference it twice to get the
  #  actual list value.
  set(local_patch_files ${${patch_files}})

  # build the patch.sh file
  configure_file(
    ${SuperBuild_TEMPLATE_FILES_DIR}/generic-patch-step.sh.in
    ${out_dir}/${name}-patch-step.sh
    @ONLY
  )

  # configure the cmake patch step
  set(sh_patch ${out_dir}/${name}-patch-step.sh)
  configure_file(
    ${SuperBuild_TEMPLATE_FILES_DIR}/generic-patch-step.cmake.in
    ${out_dir}/${name}-patch-step.cmake
    @ONLY
  )

  # set the patch command
  set(${name}_PATCH_COMMAND ${CMAKE_COMMAND} -P ${out_dir}/${name}-patch-step.cmake PARENT_SCOPE)

endfunction()