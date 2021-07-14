#include <cstdlib>
#include <iostream>
#include "EState.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace Descriptors {
namespace {

int getPeriod(int AtomicNum) {
  if (AtomicNum <= 2) {
    return 1;
  } else if (AtomicNum <= 10) {
    return 2;
  } else if (AtomicNum <= 18) {
    return 3;
  } else if (AtomicNum <= 36) {
    return 4;
  } else if (AtomicNum <= 54) {
    return 5;
  } else if (AtomicNum <= 86) {
    return 6;
  } else {
    return 7;
  }
}
    
int getGroup(int AtomicNum) {
    if (AtomicNum == 1 or AtomicNum == 3 or AtomicNum == 11 or AtomicNum == 19 or AtomicNum == 37 or AtomicNum == 55 or AtomicNum == 87) {
        return 1;
    }

    if (AtomicNum == 4 or AtomicNum == 12 or AtomicNum == 20 or AtomicNum == 38 or  AtomicNum == 56 or AtomicNum == 88) {
        return 2;
    }
    
    if (AtomicNum == 21 or AtomicNum == 39 or  AtomicNum == 71 or AtomicNum == 103) {
        return 3;
    }
    
    if (AtomicNum == 22 or AtomicNum == 40 or  AtomicNum == 72 or AtomicNum == 104) {
        return 4;
    }
    
    if (AtomicNum == 23 or AtomicNum == 41 or  AtomicNum == 73 or AtomicNum == 105) {
        return 5;
    }
    
    if (AtomicNum == 24 or AtomicNum == 42 or  AtomicNum == 74 or AtomicNum == 106) {
        return 6;
    }

    if (AtomicNum == 25 or AtomicNum == 43 or  AtomicNum == 75 or AtomicNum == 107) {
        return 7;
    }
    
    if (AtomicNum == 26 or AtomicNum == 44 or  AtomicNum == 76 or AtomicNum == 108) {
        return 8;
    }
    
    if (AtomicNum == 27 or AtomicNum == 45 or  AtomicNum == 77 or AtomicNum == 109) {
        return 9;
    }
    
    if (AtomicNum == 28 or AtomicNum == 46 or  AtomicNum == 78 or AtomicNum == 110) {
        return 10;
    }
    
    if (AtomicNum == 29 or AtomicNum == 47 or  AtomicNum == 79 or AtomicNum == 111) {
        return 11;
    }
    
    if (AtomicNum == 30 or AtomicNum == 48 or  AtomicNum == 80 or AtomicNum == 112) {
        return 12;
    }
    
    if (AtomicNum == 5 or  AtomicNum == 13 or AtomicNum == 31  or AtomicNum == 49  or AtomicNum == 81   or AtomicNum == 113) {
        return 13;
    }
    
    if (AtomicNum == 6 or  AtomicNum == 14 or AtomicNum == 32  or AtomicNum == 50  or AtomicNum == 82   or AtomicNum == 114) {
        return 14;
    }
    
    if (AtomicNum == 7 or  AtomicNum == 15 or AtomicNum == 33  or AtomicNum == 51  or AtomicNum == 83   or AtomicNum == 115) {
        return 15;
    }
    
    if (AtomicNum == 8 or  AtomicNum == 16 or AtomicNum == 34  or AtomicNum == 52  or AtomicNum == 84   or AtomicNum == 116) {
        return 16;
    }
    
    if (AtomicNum == 9 or  AtomicNum == 17 or AtomicNum == 35  or AtomicNum == 53  or AtomicNum == 85  or AtomicNum == 117) {
        return 17;
    }
    
    if (AtomicNum == 2 or  AtomicNum == 10 or AtomicNum == 18 or AtomicNum == 36  or AtomicNum == 54  or AtomicNum == 86  or AtomicNum == 118) {
        return 18;
    }

    return 0;
    
}    
    

std::vector<double> getIState(const RDKit::ROMol& mol, bool Hs_default_one) {
  int numAtoms = mol.getNumAtoms();
    
  std::vector<double>  Is(numAtoms, 0.0);

  if (Hs_default_one) {
      std::fill(Is.begin(), Is.end(), 1.0);
  }
  
  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom* atom = mol.getAtomWithIdx(i);
    int atNum = atom->getAtomicNum();
    double d = (double)atom->getDegree();
    if (d > 0 && atNum > 1) {
      int h = atom->getTotalNumHs(true); 
      int dv = RDKit::PeriodicTable::getTable()->getNouterElecs(atNum) -
               h;  
      int N = getPeriod(atNum);  
      
      if (d > 0) {
          Is[i] =  (4.0 / (N * N) * dv + 1.0) / d ;
      }
    }
  }

  return Is;
}

std::vector<double> getEState(const RDKit::ROMol& mol, int confId, bool istopolographic) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> Is = getIState(mol, false);

  double tmp, p;
  
  double* dist;
  // swish between topological or topographic distance (bond distance or euclidian distance)  
  if (istopolographic) {
      dist = MolOps::get3DDistanceMat(mol, confId, false, false, nullptr);
  }
  else {
      dist = MolOps::getDistanceMat(mol, false, false);
   }

    /*
    int length = sizeof(dist)/sizeof(dist[0]);  
    for (int k = 0; k < length; k++) {   
       dist[k]+=1.0;  
    } 
    */

    // is the 3D distance need also the + 1 ????
    
  std::vector<double> accum(numAtoms, 0.0);

  for (int i = 0; i < numAtoms; i++) {
    for (int j = i + 1; j < numAtoms; j++) {
   
      p = dist[i * numAtoms + j] + 1.0 ;     // + 1 is here for Topologic python rdkit matching!
      if (p < 1e6) {
        tmp = (Is[i] - Is[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }

  for (int i = 0; i < numAtoms; i++) {
    Is[i] += accum[i];
  }

  return Is;
}

} // end of anonymous namespace
    
void GetEStateTopographical(const ROMol &mol, std::vector<double> &res, int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  unsigned int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(numAtoms);

  res = getEState(mol, confId, true);
}
    
void GetEStateTopological(const ROMol &mol, std::vector<double> &res) {
  unsigned int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(numAtoms);

  res = getEState(mol, 0, false);
}        
    
void GetIState(const ROMol &mol, std::vector<double> &res, bool Hs_default_one) {
  unsigned int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(numAtoms);

  res = getIState(mol, Hs_default_one);
}
    
void GetAtomsGroup(const ROMol &mol,  std::vector<int> &res) { 
    unsigned int numAtoms = mol.getNumAtoms();

    res.clear();
    res.resize(numAtoms);
    
     for (int i = 0; i < numAtoms; ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        int atNum = atom->getAtomicNum();
        res[i]=  getGroup(atNum);
     }
}
    
    
void GetAtomsPeriod(const ROMol &mol, std::vector<int> &res) { 
    unsigned int numAtoms = mol.getNumAtoms();

    res.clear();
    res.resize(numAtoms);
    
    for (int i = 0; i < numAtoms; ++i) {
        const RDKit::Atom* atom = mol.getAtomWithIdx(i);
        int atNum = atom->getAtomicNum();
        res[i] =  getPeriod(atNum);
     }
        
}   
    
    
}  // namespace Descriptors
}  // namespace RDKit
