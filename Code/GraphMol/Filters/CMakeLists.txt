remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_FIRFILTERS_BUILD)
rdkit_library(FirFilters
	FirFilters.cpp
	LINK_LIBRARIES MolStandardize MolHash  SmilesParse GraphMol RDGeneral)

      rdkit_headers(FirFilters.h DEST GraphMol/FirFilters)

#rdkit_catch_test(testCondensedGraphRxn catch_main.cpp catch_tests.cpp 
#LINK_LIBRARIES CondensedGraphRxn ChemReactions SubstructMatch  SmilesParse )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()