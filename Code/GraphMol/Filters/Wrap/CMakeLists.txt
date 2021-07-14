remove_definitions(-DRDKIT_FIRFILTERS_BUILD)
rdkit_python_extension(rdFirFilters
                       FirFiltersPy.cpp
		       DEST Chem/FirFilters
		       LINK_LIBRARIES FirFilters GraphMol RDGeneral) 

#add_pytest(pyCondensedGraphRxn ${CMAKE_CURRENT_SOURCE_DIR}/testCondensedGraphRxn.py)