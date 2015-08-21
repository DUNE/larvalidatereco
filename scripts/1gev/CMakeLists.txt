set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{ROOTSYS}/lib)
art_make(	BASENAME_ONLY
		LIBRARY_NAME Content
		LIB_LIBRARIES
                        RecoBase
                        RecoBaseArt
                        Simulation
                        Geometry
                        Geometry_service
                        RecoObjects
                        LArProperties_service
                        DetectorProperties_service
                        Utilities
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
			LArPandoraInterface
			Framework
			LArAnalysis
         )

install_headers()
install_fhicl()
install_source()
