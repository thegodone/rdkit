#include <GraphMol/ChemReactions/Reaction.h>
#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Filters/FirFilters.h>

namespace python = boost::python;
namespace ff = RDKit::FirFilters;

namespace RDKit {
namespace FirFilters {
namespace Aux {

    bool IsOrganicSmiles(const std::string &smiles) {
        ROMol* mol = RDKit::SmilesToMol( smiles );
        return mol ? IsOrganic(mol) : false;
    }    
}
}
}

BOOST_PYTHON_MODULE(rdFirFilters) {
    
    python::def("IsOrgMol", &ff::IsOrganic, (python::arg("mol")),
               "Method checks if a molecule is organic. Organic elements are {H,B,C,N,O,F,P,S,Cl,Br,I}.");
    
    python::def("IsOrgSmiles", &ff::Aux::IsOrganicSmiles, (python::arg("smiles")),
                "Method checks if a SMILES string defines an organic molecule. Organic elements are {H,B,C,N,O,F,P,S,Cl,Br,I}.");
        
    
}