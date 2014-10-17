HD_COMPONENT_NAME	= heacore

HD_COMPONENT_VERS	=

HD_CXXTASK              = phaIIbin
#HD_CXXTASK              = rmfexample
#HD_CXXTASK              = tableexample

HD_CXXTASK_SRC_cxx      = phaIIbin.cxx
#HD_CXXTASK_SRC_cxx      = rmfexample.cxx
#HD_CXXTASK_SRC_cxx      = tableexample.cxx

HD_CXXFLAGS             = ${HD_STD_CXXFLAGS}

HD_CXXLIBS		= ${HD_LFLAGS} -l${CCFITS} -l${HEASP} ${HD_STD_CXXLIBS}

HD_INSTALL_TASKS	= ${HD_CXXTASK}

include ${HD_STD_MAKEFILE}

