//
//  Copyright (C) 2020 Guillaume GODIN
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Adding ATOM FEATURES descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <iostream>

#include "AtomFeat.h"
#include "MolData3Ddescriptors.h"
#include "Data3Ddescriptors.h"

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <GraphMol/Atom.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MolOps.h>
#include <cmath>
#include <vector>

namespace RDKit {

namespace Descriptors {

namespace {

std::vector<Atom::ChiralType> RS{Atom::CHI_UNSPECIFIED,
	                         Atom::CHI_TETRAHEDRAL_CW,
                                 Atom::CHI_TETRAHEDRAL_CCW, 
				 Atom::CHI_OTHER};
std::vector<std::string> Symbols{"B", "C",  "N", "O", "S", "F", "Si",
                                 "P", "Cl", "Br", "I", "H"}; // [12 most representative atoms!]
std::vector<Atom::HybridizationType> HS{Atom::SP, Atom::SP2, Atom::SP3,
                                        Atom::SP3D, Atom::SP3D2};


void AtomFeatVectorori(const RDKit::Atom* atom, const ROMol* mol,
                    std::vector<double>& feats, bool addchiral) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(mol, "bad mol");

  if (addchiral) {
    feats.reserve(53);
  } else {
    feats.reserve(49);
  }

  // initiate ring info if not already done
  if (!mol->getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(*mol);
  }

  // one hot atom symbols
  std::string s = atom->getSymbol();
  bool inlist = false;
  int indx = 0;
  for (auto ind : Symbols) {
    if (ind == s) {
      feats[indx] = 1;
      inlist = true;
    } else {
      feats[indx] = 0;
    }
    ++indx;
  }

  // write UNK type if not found in the symbol list
  feats[indx] = (inlist ? 0 : 1);
  ++indx;

  // one hot degree (why 0) [0...6]
  int d = atom->getDegree();
  for (int i = 0; i < 7; i++) {
    feats[indx] = d == i;
    ++indx;
  }

  Atom::HybridizationType hs = atom->getHybridization();
  // one hot hybridization type
  for (auto hsquery : HS) {
    feats[indx] = hs == hsquery;
    ++indx;
  }

  // one hot  Implicit Valence [0...6]
  int IV = atom->getImplicitValence();
  for (int i = 0; i < 7; ++i) {
    feats[indx] = IV == i;
    ++indx;
  }

  // one hot  getFormalCharge [-1,0,1]
  int fc = atom->getFormalCharge();
  for (int i = -1; i < 2; i++) {
    feats[indx] = fc == i;
    ++indx;
  }

  // get small mid ring size [3...8]
  int atomid = atom->getIdx();
  for (unsigned int i = 3; i < 9; i++) {
    feats[indx] = mol->getRingInfo()->isAtomInRingOfSize(atomid, i);
    ++indx;
  }

  // Is in ring
  int isInRing = queryIsAtomInRing(atom);
  feats[indx++] =(double) isInRing;


  // Is aromatic
  feats[indx] = (atom->getIsAromatic());
  ++indx;

  // one hot  Total NumH [0...4]
  unsigned int toth = atom->getTotalNumHs(false);
  for (unsigned int i = 0; i < 5; ++i) {
    feats[indx] = toth == i;
    ++indx;
  }

  // add numatoms
  feats[indx] = 1. / mol->getNumAtoms();
  ++indx;

  // put if here
  if (addchiral) {
    Atom::ChiralType rs = atom->getChiralTag();
    // one hot getChiralTag type
    for (auto rsquery : RS) {
      feats[indx] = rs == rsquery;
      ++indx;
    }
  }
}


void AtomFeatEmb(const RDKit::Atom* atom, const ROMol* mol,
                    std::vector<double>& feats) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(mol, "bad mol");

  feats.reserve(10);
  
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }
  int indx = 0;  
  // AtomNumber list using Symbol
  // one hot atom symbols
  int atno = atom->getAtomicNum();
  if (atno>119) {atno = 0;}
  feats[indx++] = (double) atno;

  std::cout << "AtNo:" << atno << std::endl;

  Atom::ChiralType rs = atom->getChiralTag();
  int i = 1;
  feats[indx] =(double) 0;
  for (auto rsquery : RS) {
      if (rs == rsquery) {
         feats[indx] =(double) i;
         break;
      }
      i=i+1;
  }
  indx++;

  std::cout << "chiral:" << rs << ", i :" << i  << std::endl;

  int d = atom->getTotalDegree();
  if (d>10) {d = 11;} 
  feats[indx++] = (double)d;

  int fc = atom->getFormalCharge();
  if (fc>5) {fc=6;}
  if (fc<-5) {fc=6;}
  feats[indx++] = (double)fc;

  // TotalNumHs of the atom
  int nHs = atom->getTotalNumHs();
  if (nHs>8) {nHs=9;}
  feats[indx++] = (double)nHs;


  int nEls = atom->getNumRadicalElectrons();
  if (nEls>4) {nEls=5;}
  feats[indx++] =(double) nEls;


  Atom::HybridizationType hs = atom->getHybridization();  
  // one hot hybridization type
  i = 1;
  feats[indx] = (double)0;
  for (auto hsquery : HS) {
    if (hs == hsquery) {
        feats[indx] =(double) i;
	break;
    }
    i=i+1;
  }
  indx++;
 

  feats[indx++] = (double)atom->getIsAromatic();

  int isInRing = queryIsAtomInRing(atom); 
  feats[indx++] =(double) isInRing;

  // adding ring size
  feats[indx] = (double)0;
  if (isInRing>0) {
     int atomid = atom->getIdx();     
     for (unsigned int i = 3; i < 9; i++) {
        if (mol->getRingInfo()->isAtomInRingOfSize(atomid, i)) {
           feats[indx] = (double)i;
           break;
        }
     }
  }
  indx++;
  
}


void AtomFeatVector(const RDKit::Atom* atom, const ROMol* mol,
                    std::vector<double>& feats, bool addchiral, bool add3dfeatures) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(mol, "bad mol");

  
  if (addchiral && !add3dfeatures) {
    feats.reserve(62);
  }
  if (addchiral && add3dfeatures) {
    feats.reserve(70);
  }
  if (!addchiral && !add3dfeatures) {
    feats.reserve(58);
  }
  if (!addchiral && add3dfeatures) {
     feats.reserve(66);
  }
    
  Data3Ddescriptors data3D;

  double* relativeMw = data3D.getMW();
  double* relativePol = data3D.getPOL();
  double* relativeVdW = data3D.getVDW();
  double* rcov = data3D.getRCOV();
  double* relativeNeg = data3D.getNEG();
  double* absionpol = data3D.getIonPOL();
  
  
  // initiate ring info if not already done
  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::findSSSR(*mol);
  }

  // AtomNumber list using Symbol
  // one hot atom symbols
  std::string s = atom->getSymbol();
  bool inlist = false;
  int indx = 0; // what is the position of the first one ???
  for (auto ind : Symbols) {
    if (ind == s) {
      feats[indx++] = 1;
      inlist = true;
    } else {
      feats[indx++] = 0;
    }
  }
  // write UNK type if not found in the symbol list
  feats[indx++] = (inlist ? 0 : 1); // last symbol : ie 13th position

  // Total Degree ie neig atoms + Hs link to the atom!
  inlist= false;
  int d = atom->getTotalDegree(); // 14th to 23th positions
  for (int i = 0; i < 11; i++) {
    if (d==i) {
       feats[indx++] = 1;
       inlist = true;
     } else {
       feats[indx++] = 0;
     }
  }
  // write UNK type if not found in the degree
  feats[indx++] = (inlist ? 0 : 1); // 24th position


  // Hybridization of the atom
  inlist= false;
  Atom::HybridizationType hs = atom->getHybridization(); // 25th to 29th 
  // one hot hybridization type
  for (auto hsquery : HS) {
    if (hs == hsquery) {
        feats[indx++] = 1;
        inlist = true;
    } else {
        feats[indx++] = 0;
    }
  }
  // write UNK type if not found in the Hybridization
  feats[indx++] = (inlist ? 0 : 1); // 30th position

  // TotalNumHs of the atom
  inlist= false;
  int nHs = atom->getTotalNumHs(); // 31th to 36th position
  for (int i = 0; i < 7; ++i) {
    if (nHs == i) {
      feats[indx++] = 1;
      inlist =true;
    } else {
      feats[indx++] = 0;
    }
  }
  feats[indx++] = (inlist ? 0 : 1); // 37th position
  
  // Formal Charge of the atom
  inlist= false; 
  int fc = atom->getFormalCharge(); // 38th to 42th position
  for (int i = -2; i < 3; i++) {
    if (fc==i) {
      inlist= true;
      feats[indx++] = 1;
    } else {
      feats[indx++] = 0;
    }
  }
  feats[indx++] = (inlist ? 0 : 1); // 43th position

  // Ring size of the atom (can be in more than one ring...) 
  inlist= false;
  int nbrings = 0;
  // one hot ring size 
  int atomid = atom->getIdx(); // 44th to 49th
  for (unsigned int i = 3; i < 9; i++) {
    if ( mol->getRingInfo()->isAtomInRingOfSize(atomid, i)) {
      if (!inlist) {
      feats[indx++] = 1;
      }
      nbrings +=1;
      inlist= true;
    } else {
      feats[indx++] = 0;
    }
  } // 50th


  int isInRIng = queryIsAtomInRing(atom);  // new method to get IsInRing!
  feats[indx++] = isInRIng; // (inlist ? 1 : 0); // 57th

  // adding radicals // 51th to 53st
  inlist= false;
  int nEls = atom->getNumRadicalElectrons();
  for (int i = 0; i < 3; i++) {
    if (nEls==i) {
      inlist= true;
      feats[indx++] = 1;
    } else {
      feats[indx++] = 0;
    }
  }
  feats[indx++] = (inlist ? 0 : 1); // 54th position

  if (addchiral) {
    // 57th to 60th
  Atom::ChiralType rs = atom->getChiralTag(); // by default in the OGB code
    // one hot getChiralTag type
    for (auto rsquery : RS) {
	    if (rs == rsquery) {
              feats[indx++] = 1;
	    } else {
	      feats[indx++] = 0;
	    }
    }
  }

    // Aromaticity of the atom // 55th
   feats[indx++] = atom->getIsAromatic(); 
  //// start to be not long!
   // add numatoms  // 56nd
   feats[indx++] = 1. / mol->getNumAtoms();
 
    // put if here need to divide by max value!
    if (add3dfeatures) {
	// 61th to 68th
        int atNum = atom->getAtomicNum();

        int lookupidx =atNum - 1;
        // maybe the indx++ look slower but simpler to write
        feats[indx++] = relativeMw[lookupidx] / 22.5645;
        feats[indx++] = relativePol[lookupidx] / 33.8636;
        feats[indx++] = relativeVdW[lookupidx] / 4.2328;
        feats[indx++] = rcov[lookupidx] / 3.4211;
        feats[indx++] = relativeNeg[lookupidx] / 1.6364;
        feats[indx++] = absionpol[lookupidx] / 2.1835;

        MolData3Ddescriptors moldata3D;


        int QN = moldata3D.GetPrincipalQuantumNumber(atNum);
        feats[indx++] = QN / 7.0;
                                 
        // IState
	
        
        int degree = atom->getDegree();  // number of substituants (heavy of not?)
        if (degree > 0 && atNum > 1) {
           int h = atom->getTotalNumHs(
               true);  // caution getTotalNumHs(true) to count h !!!!
           int dv = RDKit::PeriodicTable::getTable()->getNouterElecs(atNum) -
                    h;  // number of valence (explicit with Hs)
           double d = (double)degree - h;             // degree-h
           if (d > 0) {
		// normalized the output replace 4.0 by 1.0 on numerator
                feats[indx++] =  (1.0 / (QN * QN) * dv + 1.0) / d;
            }
	   else {
	         feats[indx++] = 0.0;
	   }
        } else {

          feats[indx++] = 0.0;
	}
        

    }
}
}  // end of anonymous namespace

// entry point
void AtomFeatVect(const ROMol& mol, std::vector<double>& res, int atomid,
                  bool addchiral, bool add3dfeatures, bool emb, bool ori) {
  res.clear();
    


  if (addchiral && !add3dfeatures) {
      res.resize(62);
  }

  if (addchiral && add3dfeatures) {
      res.resize(70);
  }

  if (!addchiral && !add3dfeatures) {
      res.resize(58);
  }

  if (!addchiral && add3dfeatures) {
      res.resize(66);
  }

  if (emb) {  res.resize(10);
  }

  if (ori && !addchiral) { res.resize(49);
  }

  if (ori && addchiral) { res.resize(53);
  }


  if (!emb) {
	  if(!ori) {
              AtomFeatVector(mol.getAtomWithIdx(atomid), &mol, res, addchiral, add3dfeatures);
	  }
	  else {
	      AtomFeatVectorori(mol.getAtomWithIdx(atomid), &mol, res, addchiral);
	  }
  } else {
     AtomFeatEmb(mol.getAtomWithIdx(atomid), &mol, res);
  }
  
  }

}  // namespace Descriptors
}  // namespace RDKit
