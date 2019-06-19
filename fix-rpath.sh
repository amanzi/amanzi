for i in system filesystem program_options regex; do install_name_tool -change @rpath/libboost_${i}.dylib ${AMANZI_TPLS_DIR}/lib/libboost_${i}.dylib ${ATS_DIR}/bin/ats; done
