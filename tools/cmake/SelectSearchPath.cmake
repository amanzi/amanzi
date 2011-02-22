# -*- mode: cmake -*-

# 
# Amanzi SelectSearchPath
#
include(PrintVariable)

function(select_search_path PACKNAME OUT_PATH OUT_PATH_FOUND)

    set(pack_dir          "${${PACKNAME}_DIR}")
    set(env_pack_root     "${PACKNAME}_ROOT")
    set(pack_root          $ENV{${env_pack_root}})

    #PRINT_VARIABLE(pack_dir)
    #PRINT_VARIABLE(env_pack_root)
    #PRINT_VARIABLE(pack_root)

    set(test_path "")
    if ( pack_dir )

        set(test_path "${pack_dir}")

    else(pack_dir)
        
        if (pack_root)
            set(test_path "${pack_root}")
        endif(pack_root)

    endif(pack_dir)    
    
    
    string(LENGTH "${test_path}" path_len)
    #PRINT_VARIABLE(path_len)
    if ( ${path_len} GREATER "0" )
        if ( EXISTS "${test_path}" )
            set(${OUT_PATH_FOUND} TRUE PARENT_SCOPE)
            set(${OUT_PATH} "${test_path}" PARENT_SCOPE)
        else()
            message(SEND_ERROR "The directory '${test_path}' does not exist. "
                               "Can not use this directory as a search path for ${PACKNAME}")
            set(${OUT_PATH_FOUND} FALSE PARENT_SCOPE)
        endif()    
    else()
        set(${OUT_PATH_FOUND} FALSE PARENT_SCOPE)
    endif()


endfunction()    

