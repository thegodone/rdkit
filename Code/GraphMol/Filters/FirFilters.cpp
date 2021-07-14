#include "FirFilters.h"
#include <boost/shared_ptr.hpp>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GraphMol.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/utils.h>

namespace RDKit {
namespace FirFilters {

bool IsOrganicAtom(Atom* atom) {
    switch(atom->getAtomicNum()) {
        case 1:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 15:
        case 16:
        case 17:
        case 35:
        case 54:
            return true;
        default:
            return false;
    }
}
    
bool IsOrganic(ROMol* mol) {
    bool isorganic = true;
    for (Atom* atom : mol->atoms()) {
        if (isorganic)
            isorganic = IsOrganicAtom(atom);
    }    
    return isorganic;
}
    
}
}