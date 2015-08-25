
# Compile a LArAnalysis dictionary and build the shared object...
include             (CetRootCint)
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
set                 (PACKAGE LArAnalysis)
cet_rootcint        (${PACKAGE})

FILE( GLOB src_files LArAnalysis.* LArValidation.* LArPID.* LArTracking.* )

art_make_library( LIBRARY_NAME ${PACKAGE}
        	  SOURCE ${src_files} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
                  LIBRARIES 
                        ${ROOT_BASIC_LIB_LIST}
                )

install_headers()
install_fhicl()
install_source()
