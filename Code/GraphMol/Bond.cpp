//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
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
  double res;
  switch (getBondType()) {
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      res = 0;
      break;
    case SINGLE:
      res = 1;
      break;
    case DOUBLE:
      res = 2;
      break;
    case TRIPLE:
      res = 3;
      break;
    case QUADRUPLE:
      res = 4;
      break;
    case QUINTUPLE:
      res = 5;
      break;
    case HEXTUPLE:
      res = 6;
      break;
    case ONEANDAHALF:
      res = 1.5;
      break;
    case TWOANDAHALF:
      res = 2.5;
      break;
    case THREEANDAHALF:
      res = 3.5;
      break;
    case FOURANDAHALF:
      res = 4.5;
      break;
    case FIVEANDAHALF:
      res = 5.5;
      break;
    case AROMATIC:
      res = 1.5;
      break;
    case CGRSD:
      res = 2.0;
      break;
    case CGRST:
      res = 3.0;
      break;
    case CGRSA:
      res = 1.5;
      break;
    case CGRSN:
      res = 0.0;
      break;
    case CGRNS:
      res = 1.0;
      break;
     case CGRDS:
      res = 1.0;
      break;
    case CGRDT:
      res = 3.0;
      break;
    case CGRDA:
      res = 1.5;
      break;
    case CGRDN:
      res = 0.0;
      break;
    case CGRND:
      res = 2.0;
      break;
    case CGRTS:
      res = 1.0;
      break;
    case CGRTD:
      res = 2.0;
      break;
    case CGRTA:
      res = 1.5;
      break;
    case CGRTN:
      res = 0.0;
      break;
    case CGRNT:
      res = 3.0;
      break;
    case CGRAS:
      res = 1.0;
      break;
    case CGRAD:
      res = 2.0;
      break;
    case CGRAT:
      res = 3.0;
      break;
    case CGRAN:
      res = 0.0;
      break;
    case CGRNA:
      res = 1.5;
      break;
    case DATIVEONE:
      res = 1.0;
      break;  // FIX: this should probably be different
    case DATIVE:
      res = 1.0;
      break;  // FIX: again probably wrong
    case HYDROGEN:
      res = 0.0;
      break;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
  return res;
}

double Bond::getValenceContrib(const Atom *atom) const {
  double res;
  if ((getBondType() == DATIVE || getBondType() == DATIVEONE) &&
      atom->getIdx() != getEndAtomIdx()) {
    res = 0.0;
  } else {
    res = getBondTypeAsDouble();
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

    
    
uint8_t getTwiceBondType(const Bond &b) {
  switch (b.getBondType()) {
    case Bond::UNSPECIFIED:
    case Bond::IONIC:
    case Bond::ZERO:
      return 0;
      break;
    case Bond::SINGLE:
      return 2;
      break;
    case Bond::DOUBLE:
      return 4;
      break;
    case Bond::TRIPLE:
      return 6;
      break;
    case Bond::QUADRUPLE:
      return 8;
      break;
    case Bond::QUINTUPLE:
      return 10;
      break;
    case Bond::HEXTUPLE:
      return 12;
      break;
    case Bond::ONEANDAHALF:
      return 3;
      break;
    case Bond::TWOANDAHALF:
      return 5;
      break;
    case Bond::THREEANDAHALF:
      return 7;
      break;
    case Bond::FOURANDAHALF:
      return 9;
      break;
    case Bond::FIVEANDAHALF:
      return 11;
      break;
    case Bond::AROMATIC:
      return 3;
      break;
    case Bond::DATIVEONE:
      return 2;
      break;  // FIX: this should probably be different
    case Bond::DATIVE:
      return 2;
      break;  // FIX: again probably wrong
    case Bond::HYDROGEN:
      return 0;
      break;
		  
    case Bond::CGRSD:
      return 4;
      break;
    case Bond::CGRST:
      return 6;
      break;
    case Bond::CGRSA:
      return 3;
      break;
    case Bond::CGRSN:
      return 0;
      break;
    case Bond::CGRNS:
      return 2;
      break;
    case Bond::CGRDS:
      return 2;
      break;
    case Bond::CGRDT:
      return 6;
      break;
    case Bond::CGRDA:
      return 3;
      break;
    case Bond::CGRDN:
      return 0;
      break;
    case Bond::CGRND:
      return 4;
      break;
    case Bond::CGRTS:
      return 2;
      break;
    case Bond::CGRTD:
      return 4;
      break;
    case Bond::CGRTA:
      return 3;
      break;
    case Bond::CGRTN:
      return 0;
      break;
    case Bond::CGRNT:
      return 6;
      break;
    case Bond::CGRAS:
      return 2;
      break;
    case Bond::CGRAD:
      return 4;
      break;
    case Bond::CGRAT:
      return 6;
      break;
    case Bond::CGRAN:
      return 0;
      break;
    case Bond::CGRNA:
      return 3;
      break;
		  
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
}
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
