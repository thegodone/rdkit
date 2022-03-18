//
//  Copyright (C) 2020 Guillaume GODIN @  Firmenich
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_CONDENSEDGRAPHRXN_H
#define RD_CONDENSEDGRAPHRXN_H

#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <iostream>

namespace RDKit {
class ROMol;
class RWMol;
class Atom;
class ChemicalReaction;

namespace CondensedGraphRxn {

RDKIT_CONDENSEDGRAPHRXN_EXPORT std::string CRSwriter(const std::string smart, bool doRandom, unsigned int randomSeed
						     , bool aromatize, bool signature, bool charges, int radius,
						     bool addRingInfo);

RDKIT_CONDENSEDGRAPHRXN_EXPORT std::string CRSreader(RWMol *molR, const std::string crs, bool canonical, bool setAtomMap);

RDKIT_CONDENSEDGRAPHRXN_EXPORT std::string RXNCleaning(std::string rxnsmart);

RDKIT_CONDENSEDGRAPHRXN_EXPORT std::string RXNCompleteMapping(std::string sma, bool debug, bool addleavinggroups);

RDKIT_CONDENSEDGRAPHRXN_EXPORT bool getRXNComparison(std::string rxnsma, std::string rxncoresma);

RDKIT_CONDENSEDGRAPHRXN_EXPORT bool getRXNCompTotal(std::string rxnsma, std::string rxncoresma);

RDKIT_CONDENSEDGRAPHRXN_EXPORT std::string setRXNCompAtomMaps(std::string rxnsma, std::string rxncoresma);
    
RDKIT_CONDENSEDGRAPHRXN_EXPORT void transferMappedNum(ROMol *mol1, ROMol *mol2);

RDKIT_CONDENSEDGRAPHRXN_EXPORT bool getRxnDirection(const ChemicalReaction *rxn);

}  // namespace CondensedGraphRxn
}
#endif
