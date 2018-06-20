#
# FIND_SPACK_INSTALL_LOCATION
#
# Usage:
# FIND_SPACK_INSTALL_LOCATION(pkg_name output_name)
#
# Given a package target name, return the full path of the
# installation prefix.

function(FIND_SPACK_INSTALL_LOCATION pkg_name output_name)

#  ${SPACK_BINARY} find -ld ${pkg_name} | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g" | awk '/^[a-z]/{print $1}' | xargs -n1 ${SPACK_BINARY} find -p / | awk '/^[^-]/{print $2";"}' | tr -d '\n'

#set(${output_name} ${SPACK_BINARY} find -ld ${pkg_name} | sed -r "s/\\x1B\[\([0-9]{1,2}\(;[0-9]{1,2}\)?\)?[mGK]//g" | awk '/^[a-z]/{print $1}' | xargs -n1 ${SPACK_BINARY} find -p / | awk '/^[^-]/{print $2";"}' | tr -d '\n' PARENT_SCOPE )
  
#execute_process( COMMAND ${SPACK_BINARY} find -ld ${pkg_name} | sed -r "s/\\x1B\[\([0-9]{1,2}\(;[0-9]{1,2}\)?\)?[mGK]//g" | awk '/^[a-z]/{print $1}' | xargs -n1 ${SPACK_BINARY} find -p / | awk '/^[^-]/{print $2";"}' | tr -d '\n'
#    OUTPUT_VARIABLE temp
#)

execute_process(COMMAND ${SPACK_BINARY} find -ld ${pkg_name}
    COMMAND sed -r "s#x1B\\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]##g"
    COMMAND awk "/^[a-z,0-9]/{print $1}"
    COMMAND xargs -n1 ${SPACK_BINARY} find -p /
    COMMAND awk "NR==2,/^[^-]/{print $2\";\"}"
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE temp
)
#/^[^-]/{print $2\";\"}"
#COMMAND ${SPACK_BINARY} find -ld ${pkg_name}
#                  COMMAND sed -r "s/\\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[mGK]//g"
#		  COMMAND awk '/^[a-z]/{print $1}' 
#		  COMMAND xargs -n1 ${SPACK_BINARY} find -p / 
#		  COMMAND awk '/^[^-]/{print $2";"}' 
#        	  COMMAND tr -d '\n'
#		  OUTPUT_VARIABLE temp
#		  )

set ( ${output_name} ${temp} PARENT_SCOPE )


endfunction(FIND_SPACK_INSTALL_LOCATION)
                     
