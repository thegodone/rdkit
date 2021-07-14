//
//  Copyright (C) 2001-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <sstream>
#include <string>
#include <iomanip>

#include "Bond.h"
#include "Atom.h"
#include "ROMol.h"
#include <RDGeneral/Invariant.h>
#include <map>
namespace RDKit {

Bond::Bond() : RDProps() { initBond(); };

Bond::Bond(BondType bT) : RDProps() {
  initBond();
  d_bondType = bT;
};

Bond::Bond(const Bond &other) : RDProps(other) {
  // NOTE: we do *not* copy ownership!
  dp_mol = nullptr;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  d_stereo = other.d_stereo;
  if (other.dp_stereoAtoms) {
    dp_stereoAtoms = new INT_VECT(*other.dp_stereoAtoms);
  } else {
    dp_stereoAtoms = nullptr;
  }
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
}

Bond::~Bond() { delete dp_stereoAtoms; }

Bond &Bond::operator=(const Bond &other) {
  if (this == &other) {
    return *this;
  }
  dp_mol = other.dp_mol;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  delete dp_stereoAtoms;
  if (other.dp_stereoAtoms) {
    dp_stereoAtoms = new INT_VECT(*other.dp_stereoAtoms);
  } else {
    dp_stereoAtoms = nullptr;
  }
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
  d_props = other.d_props;

  return *this;
}

Bond *Bond::copy() const {
  auto *res = new Bond(*this);
  return res;
}

void Bond::setOwningMol(ROMol *other) {
  // FIX: doesn't update topology
  dp_mol = other;
}

unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
  PRECONDITION(d_beginAtomIdx == thisIdx || d_endAtomIdx == thisIdx,
               "bad index");
  if (d_beginAtomIdx == thisIdx) {
    return d_endAtomIdx;
  } else if (d_endAtomIdx == thisIdx) {
    return d_beginAtomIdx;
  }
  // we cannot actually get down here
  return 0;
}

void Bond::setBeginAtomIdx(unsigned int what) {
  if (dp_mol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  d_beginAtomIdx = what;
};

void Bond::setEndAtomIdx(unsigned int what) {
  if (dp_mol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  d_endAtomIdx = what;
};

void Bond::setBeginAtom(Atom *at) {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  setBeginAtomIdx(at->getIdx());
}
void Bond::setEndAtom(Atom *at) {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  setEndAtomIdx(at->getIdx());
}

Atom *Bond::getBeginAtom() const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  return dp_mol->getAtomWithIdx(d_beginAtomIdx);
};
Atom *Bond::getEndAtom() const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  return dp_mol->getAtomWithIdx(d_endAtomIdx);
};
Atom *Bond::getOtherAtom(Atom const *what) const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");

  return dp_mol->getAtomWithIdx(getOtherAtomIdx(what->getIdx()));
};

double Bond::getBondTypeAsDouble() const {
  switch (getBondType()) {
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      return 0;
      break;
    case SINGLE:
      return 1;
      break;
    case DOUBLE:
      return 2;
      break;
    case TRIPLE:
      return 3;
      break;
    case QUADRUPLE:
      return 4;
      break;
    case QUINTUPLE:
      return 5;
      break;
    case HEXTUPLE:
      return 6;
      break;
    case ONEANDAHALF:
      return 1.5;
      break;
    case TWOANDAHALF:
      return 2.5;
      break;
    case THREEANDAHALF:
      return 3.5;
      break;
    case FOURANDAHALF:
      return 4.5;
      break;
    case FIVEANDAHALF:
      return 5.5;
      break;
    case AROMATIC:
      return 1.5;
      break;
    case CGRSD:
      return 2.0;
      break;
    case CGRST:
      return 3.0;
      break;
    case CGRSA:
      return 1.5;
      break;
    case CGRSN:
      return 0.0;
      break;
    case CGRNS:
      return 1.0;
      break;
     case CGRDS:
      return 1.0;
      break;
    case CGRDT:
      return 3.0;
      break;
    case CGRDA:
      return 1.5;
      break;
    case CGRDN:
      return 0.0;
      break;
    case CGRND:
      return 2.0;
      break;
    case CGRTS:
      return 1.0;
      break;
    case CGRTD:
      return 2.0;
      break;
    case CGRTA:
      return 1.5;
      break;
    case CGRTN:
      return 0.0;
      break;
    case CGRNT:
      return 3.0;
      break;
    case CGRAS:
      return 1.0;
      break;
    case CGRAD:
      return 2.0;
      break;
    case CGRAT:
      return 3.0;
      break;
    case CGRAN:
      return 0.0;
      break;
    case CGRNA:
      return 1.5;
      break;
    case DATIVEONE:
      return 1.0;
      break;  // FIX: this should probably be different
    case DATIVE:
      return 1.0;
      break;  // FIX: again probably wrong
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
}

double Bond::getValenceContrib(const Atom *atom) const {
  switch (getBondType()) {
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      return 0;
      break;
    case SINGLE:
      return 1;
      break;
    case DOUBLE:
      return 2;
      break;
    case TRIPLE:
      return 3;
      break;
    case QUADRUPLE:
      return 4;
      break;
    case QUINTUPLE:
      return 5;
      break;
    case HEXTUPLE:
      return 6;
      break;
    case ONEANDAHALF:
      return 1.5;
      break;
    case TWOANDAHALF:
      return 2.5;
      break;
    case THREEANDAHALF:
      return 3.5;
      break;
    case FOURANDAHALF:
      return 4.5;
      break;
    case FIVEANDAHALF:
      return 5.5;
      break;
    case AROMATIC:
      return 1.5;
      break;
    case CGRSD:
      return 2.0;
      break;
    case CGRST:
      return 3.0;
      break;
    case CGRSA:
      return 1.5;
      break;
    case CGRSN:
      return 0.0;
      break;
    case CGRNS:
      return 1.0;
      break;
    case CGRDS:
      return 1.0;
      break;
    case CGRDT:
      return 3.0;
      break;
    case CGRDA:
      return 1.5;
      break;
    case CGRDN:
      return 0.0;
      break;
    case CGRND:
      return 2.0;
      break;
    case CGRTS:
      return 1.0;
      break;
    case CGRTD:
      return 2.0;
      break;
    case CGRTA:
      return 1.5;
      break;
    case CGRTN:
      return 0.0;
      break;
    case CGRNT:
      return 3.0;
      break;
    case CGRAS:
      return 1.0;
      break;
    case CGRAD:
      return 2.0;
      break;
    case CGRAT:
      return 3.0;
      break;
    case CGRAN:
      return 0.0;
      break;
    case CGRNA:
     return 1.5;
      break;   
    case DATIVEONE:
      if (atom->getIdx() == getEndAtomIdx()) {
        return 1.0;
      } else {
        return 0.0;
      }
      break;
    case DATIVE:
      if (atom->getIdx() == getEndAtomIdx()) {
        return 1.0;
      } else {
        return 0.0;
      }
      break;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
}
    
unsigned int Bond::getNumAtomMaps() const {
    PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    unsigned int db = dp_mol->getAtomWithIdx(d_beginAtomIdx)->getAtomMapNum() > 0;
    unsigned int de = dp_mol->getAtomWithIdx(d_endAtomIdx)->getAtomMapNum() > 0;
    return db+de;
}


 
std::map<unsigned int, unsigned int > Bond::getMappedAtomsNum( bool &mappedonly) const{
  
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  std::map<unsigned int, unsigned int > res;
  unsigned int db = dp_mol->getAtomWithIdx(d_beginAtomIdx)->getAtomMapNum();
  unsigned int de = dp_mol->getAtomWithIdx(d_endAtomIdx)->getAtomMapNum();
  if (mappedonly) {
	if (db>0 ){
		res[db]  = d_beginAtomIdx ;
	}
	if (de>0 ){
                res[de]  = d_endAtomIdx ;
        }
 } 
 else {
   res[d_beginAtomIdx] = db ;
   res[d_endAtomIdx] = de;
 }


  return res;
}

std::string Bond::getBondTypeWithAtoms(bool withmapnum) const{
 PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
std::stringstream ss;
  unsigned int anb = dp_mol->getAtomWithIdx(d_beginAtomIdx)->getAtomicNum();
  unsigned int ane = dp_mol->getAtomWithIdx(d_endAtomIdx)->getAtomicNum();
  unsigned int bondtypeval;
  switch (getBondType()) { 
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      bondtypeval= 0;
      break;
    case SINGLE:
      bondtypeval= 1;
      break;
    case DOUBLE:
      bondtypeval= 2;
      break;
    case TRIPLE:
      bondtypeval = 3;
      break;
    case QUADRUPLE:
      bondtypeval= 4;
      break;
    case QUINTUPLE:
      bondtypeval= 5;
      break;
    case HEXTUPLE:
      bondtypeval= 6;
      break; 
   case AROMATIC:
      bondtypeval= 100;
      break;
   default:
      bondtypeval = 0;
   }
 ss << std::setfill('0') << std::setw(3) << bondtypeval; // bondtype number
 ss << std::setfill('0') << std::setw(3) << std::min(anb,ane); // lower atom number
 ss << std::setfill('0') << std::setw(3) << std::max(anb,ane); // higher atom number

if (withmapnum) {
 unsigned int de = dp_mol->getAtomWithIdx(d_endAtomIdx)->getAtomMapNum();
 unsigned int db = dp_mol->getAtomWithIdx(d_beginAtomIdx)->getAtomMapNum();
 ss << std::setfill('0') << std::setw(3) << std::min(db,de); // lowe atom map number
 ss << std::setfill('0') << std::setw(3) << std::max(db,de); // higher atom map number
}
 
 return ss.str();

}

void Bond::setQuery(QUERYBOND_QUERY *what) {
  //  Bonds don't have queries at the moment because I have not
  //  yet figured out what a good base query should be.
  //  It would be nice to be able to do substructure searches
  //  using molecules alone, so it'd be nice if we got this
  //  issue resolved ASAP.
  RDUNUSED_PARAM(what);
  PRECONDITION(0, "plain bonds have no Query");
}

Bond::QUERYBOND_QUERY *Bond::getQuery() const {
  PRECONDITION(0, "plain bonds have no Query");
  return nullptr;
};

bool Bond::Match(Bond const *what) const {
  bool res;
  if (getBondType() == Bond::UNSPECIFIED ||
      what->getBondType() == Bond::UNSPECIFIED) {
    res = true;
  } else {
    res = getBondType() == what->getBondType();
  }
  return res;
};

void Bond::expandQuery(Bond::QUERYBOND_QUERY *what,
                       Queries::CompositeQueryType how, bool maintainOrder) {
  RDUNUSED_PARAM(what);
  RDUNUSED_PARAM(how);
  RDUNUSED_PARAM(maintainOrder);
  PRECONDITION(0, "plain bonds have no query");
};

void Bond::initBond() {
  d_bondType = UNSPECIFIED;
  d_dirTag = NONE;
  d_stereo = STEREONONE;
  dp_mol = nullptr;
  d_beginAtomIdx = 0;
  d_endAtomIdx = 0;
  df_isAromatic = 0;
  d_index = 0;
  df_isConjugated = 0;
  dp_stereoAtoms = nullptr;
};

void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
      "bgnIdx not connected to begin atom of bond");
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
      "endIdx not connected to end atom of bond");

  INT_VECT &atoms = getStereoAtoms();
  atoms.clear();
  atoms.push_back(bgnIdx);
  atoms.push_back(endIdx);
};

};  // namespace RDKit

std::ostream &operator<<(std::ostream &target, const RDKit::Bond &bond) {
  target << bond.getIdx() << " ";
  target << bond.getBeginAtomIdx() << "->" << bond.getEndAtomIdx();
  target << " order: " << bond.getBondType();
  if (bond.getBondDir()) {
    target << " dir: " << bond.getBondDir();
  }
  if (bond.getStereo()) {
    target << " stereo: " << bond.getStereo();
    if (bond.getStereoAtoms().size() == 2) {
      const auto &ats = bond.getStereoAtoms();
      target << " stereoAts: (" << ats[0] << " " << ats[1] << ")";
    }
  }
  target << " conj?: " << bond.getIsConjugated();
  target << " aromatic?: " << bond.getIsAromatic();

  return target;
}
