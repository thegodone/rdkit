remove_definitions(-DRDKIT_CONDENSEDGRAPHRXN_BUILD)
rdkit_python_extension(rdCondensedGraphRxn
                       rdCondensedGraphRxn.cpp
		       DEST Chem/CondensedGraphRxn
		       LINK_LIBRARIES CondensedGraphRxn
GraphMol RDGeometryLib RDGeneral DataStructs ChemReactions) 

rdkit_python_extension(rdFirCRS
                       FirCRSPy.cpp
		       DEST Chem/CondensedGraphRxn
		       LINK_LIBRARIES CondensedGraphRxn
GraphMol RDGeometryLib RDGeneral DataStructs ChemReactions) 

add_pytest(pyCondensedGraphRxn ${CMAKE_CURRENT_SOURCE_DIR}/testCondensedGraphRxn.py)
