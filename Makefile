HD_COMPONENT_NAME	= heacore

HD_COMPONENT_VERS	=

HD_LIBRARY_ROOT		= ${HEASP}

HD_LIBRARY_SRC_cxx	= pha.cxx phaII.cxx SPio.cxx SPutils.cxx grouping.cxx arf.cxx arfII.cxx rmf.cxx rmft.cxx table.cxx Cwrappers.cxx

HD_SHLIB_LIBS		= ${HD_LFLAGS} -l${HEAUTILS} -l${CCFITS} -l${CFITSIO} -l${HEAIO} \
			  ${SYSLIBS}

HD_CFLAGS		= ${HD_STD_CFLAGS}

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_INSTALL_HEADERS	= heasp.h Cheasp.h pha.h phaII.h SPio.h SPutils.h grouping.h arf.h arfII.h rmf.h rmft.h table.h

HD_INSTALL_LIBRARIES	= ${HD_LIBRARY_ROOT}

HD_INSTALL_HELP		= 

HD_CXXLIBS		= ${HD_LFLAGS} -${HD_STD_CXXLIBS}

HD_TEST_SUBDIRS		= 

HD_TEST_CTASK		= 

HD_TEST_CTASK_SRC_c	= 

HD_SUBDIRS		= python

include ${HD_STD_MAKEFILE}

subdir-python:
	@if [ "x${PYTHON_LIB}" != x -a "x${PYTHON_INC}" != x ]; then \
		${HD_MAKE} hd-std-subdir HD_SUBDIR=python; \
	fi
