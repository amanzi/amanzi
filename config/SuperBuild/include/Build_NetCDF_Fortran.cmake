#  -*- mode: cmake -*-

#
# Build TPL: NetCDF-fortran
# 

# --- Define all the directories and common external project flags
define_external_project_args(NetCDF_Fortran TARGET netcdf-fortran
	DEPENDS ${MPI_PROJECT} NetCDF)

#set(cpp_flags_list -I${NetCDF_prefix_dir}/include)

set(cpp_flags_list -I${TPL_INSTALL_PREFIX}/include)
list(REMOVE_DUPLICATES cpp_flags_list)
build_whitespace_string(netcdf_fortran_cppflags ${cpp_flags_list})

# --- Add external project build and tie to the ZLIB build target
ExternalProject_Add(${NetCDF_Fortran_BUILD_TARGET}

	# Package dependency target
	DEPENDS ${NetCDF_Fortran_PACKAGE_DEPENDS}

	# Temporary files directory
	TMP_DIR ${NetCDF_Fortran_tmp_dir}

	# Timestamp and log directory
	STAMP_DIR ${NetCDF_Fortran_stamp_dir}

	# -- Downloads

	# Download directory
	DOWNLOAD_DIR ${TPL_DOWNLOAD_DIR}

	# URL may be a web site OR a local file
	URL ${NetCDF_Fortran_URL}

	# md5sum of the archive file
	URL_MD5 ${NetCDF_Fortran_MD5_SUM}

	# -- Configure

	# Source directory
	SOURCE_DIR ${NetCDF_Fortran_source_dir}

	CONFIGURE_COMMAND
		<SOURCE_DIR>/configure
		--prefix=<INSTALL_DIR>
		FC=${CMAKE_Fortran_COMPILER}
		FCFLAGS=${Amanzi_COMMON_FCFLAGS}
		CPPFLAGS=${netcdf_fortran_cppflags}

	# -- Build

	# Build directory 
	BINARY_DIR ${NetCDF_Fortran_build_dir}

	# $(MAKE) enables parallel builds through make
	BUILD_COMMAND $(MAKE)

	# Flag for in source builds
	BUILD_IN_SOURCE ${NetCDF_Fortran_BUILD_IN_SOURCE}

	# -- Install

	# Install directory
	INSTALL_DIR ${TPL_INSTALL_PREFIX}

	# -- Output control
	${NetCDF_Fortran_logging_args})
