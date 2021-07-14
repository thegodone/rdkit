remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_CONDENSEDGRAPHRXN_BUILD)
rdkit_library(CondensedGraphRxn
	CondensedGraphRxn.cpp RxnCleaning.cpp FirCGR.cpp
	LINK_LIBRARIES MolStandardize MolHash ChemReactions ChemTransforms SmilesParse SubstructMatch GraphMol RDGeneral)

      rdkit_headers(CondensedGraphRxn.h FirCGR.h DEST GraphMol/CondensedGraphRxn)

rdkit_catch_test(testCondensedGraphRxn catch_main.cpp catch_tests.cpp 
LINK_LIBRARIES CondensedGraphRxn ChemReactions SubstructMatch  SmilesParse )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
