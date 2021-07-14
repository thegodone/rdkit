/*
 * (c) RUUD, Firmenich SA, 07.2020
 */

#include <RDGeneral/export.h>
#ifndef RD_FIRFILTER_H
#define RD_FIRFILTER_H
#endif

namespace RDKit {
class ROMol;
class RWMol;
class Atom;
 
namespace FirFilters {
    
RDKIT_FIRFILTERS_EXPORT bool IsOrganic(ROMol* mol);
    
} /* End of namespace FirFilters */
} /* End of namespace RDKit */