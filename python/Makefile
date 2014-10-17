# Makefile to build the _heasp library for python.

# This Makefile does not include the SWIG command to create heasp.py and 
# heasp_wrap.cxx from heasp.i. If heasp.i is changed then these files should
# be created using
#     swig -python -c++ -classic heasp.i


HD_COMPONENT_NAME       = heacore

HD_LIBRARY_ROOT         = _heasp

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS} -I${PYTHON_INC}

HD_LIBRARY_SRC_cxx      = heasp_wrap.cxx

HD_ADD_SHLIB_LIBS       = yes

HD_SHLIB_LIBS           = ${HD_LFLAGS} -l${HEASP} -l${HEAUTILS} -l${CCFITS} -l${CFITSIO} -l${HEAIO} \
                          ${PYTHON_LIB} ${SYSLIBS}

HD_INSTALL_PYTHON_LIBS  = lib${HD_LIBRARY_ROOT}${SHLIB_SUFFIX} \
                          heasp.py

HD_PYTHON_SUBDIR        = 

HD_CLEAN                = ${HD_BLD_PYTHON_LIB}/${HD_PYTHON_SUBDIR}/${HD_LIBRARY_ROOT}.so

HD_INSTALL_EXTRA        = install-extra-libname

include ${HD_STD_MAKEFILE}

install-extra-libname:
	@if [ -f ${HD_TOP_PYTHON_LIB}/${HD_PYTHON_SUBDIR}/lib${HD_LIBRARY_ROOT}${SHLIB_SUFFIX} ]; then \
            ${HD_INSTALL} \
               ${HD_TOP_PYTHON_LIB}/${HD_PYTHON_SUBDIR}/lib${HD_LIBRARY_ROOT}${SHLIB_SUFFIX} \
               ${HD_TOP_PYTHON_LIB}/${HD_PYTHON_SUBDIR}/${HD_LIBRARY_ROOT}${PY_SHLIB_SUFFIX} ${HD_LN_S}; \
        fi
