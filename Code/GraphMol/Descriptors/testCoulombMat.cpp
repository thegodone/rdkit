//  Created by Guillaume GODIN
//  Copyright (C) 2012-2018 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <boost/tokenizer.hpp>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Descriptors/CoulombMat.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include "/Users/GVALMTGG/Github/rdkit/Code/GraphMol/Descriptors/slice.h"
#include "/Users/GVALMTGG/Github/rdkit/Code/GraphMol/Descriptors/sort.h"
#include <Eigen/Dense>

using namespace Eigen;

struct BagMatrix {
  VectorXd UCM;
  std::vector<std::string > BagTag;
};


VectorXd getUpperMat(MatrixXd mat) {

  VectorXd res(mat.rows()*(mat.cols()+1)/2);
  Index size = mat.rows();
  Index offset = 0;
  for(Index j=0; j<mat.cols(); ++j) {
      res.segment(offset,size) = mat.col(j).tail(size);
      offset += size;
      size--;
  }
  return res;
}


 BagMatrix getCoulombMatBags(Eigen::VectorXd numbers, int numatoms, Eigen::MatrixXd Distance3D,  std::vector<std::string > S, int alpha) {
  // 3D distance matrix

  Eigen::MatrixXd ProdV = numbers*numbers.transpose(); // outer products of vector ie (z[i] * z[j])

  if (alpha != 1) {
    Distance3D = Distance3D.array().pow(alpha); // power distance alpha
  }

  Eigen::MatrixXd MyVal = ProdV.cwiseQuotient(Distance3D); // ratio top / dist
  MyVal.diagonal() = 0.5 * numbers.array().pow(2.4); // set the diagonal fix values 
  // get the BagsCode
  std::vector<std::string > BagsCodes;
  std::list<std::pair<std::string, double > > Res;
  std::string s1, s2, key;
  for (unsigned int i = 0; i< numatoms; i++) {
    s1 = S[i];
    BagsCodes.push_back(s1); // adding the diagonal once!

    for (unsigned int j = i+1; j< numatoms; j++ ){
      s2 = S[j];
      if (s1 > s2) {
        key = s2+s1;
      }
      else {
        key = s1+s2;
      }
      BagsCodes.push_back(key);
    }
  }    

  Eigen::VectorXd RES = getUpperMat(MyVal);

  return BagMatrix{RES , BagsCodes}; 
}


void getBagOfBonds(const RDKit::ROMol &mol, std::vector<double> &res, double *dist3D, 
unsigned int numAtoms, std::map<std::string, unsigned int> MaxBags, int alpha) {

  // initialize the variables for the getCoulombMat
  int numatoms= mol.getNumAtoms();
  double *z = new double[numatoms];
  std::vector<std::string > S;
  for (int i=0; i< numatoms; i++){
      z[i] = mol.getAtomWithIdx(i)->getAtomicNum();
      S.push_back(mol.getAtomWithIdx(i)->getSymbol());
  }

  Eigen::VectorXd numbers = Map<VectorXd>(z, numatoms); // convert the number array to vector
  
  Eigen::MatrixXd Distance3D = Map<MatrixXd>(dist3D, numatoms, numatoms); // convert the result array to matrix (1 column)

  BagMatrix CMBags = getCoulombMatBags(numbers,  numatoms,  Distance3D, S ,alpha);
  // return the BagTag and UCM values as "upper Diagonal elements" only

  // extract the bag CM values index from BagMatrix into sorted vectors padded using MaxBags order & position index!
  std::string key;
  unsigned int sizemax;
  std::vector<double> val;
  for (std::map<std::string, unsigned int>::const_iterator MBags = MaxBags.begin();
       MBags != MaxBags.end(); ++MBags) {
         // take the Master Code Atom Pair
         key = MBags->first;
         std::vector<double> mybagarray;
         unsigned int i=0;
         // loop other the CMBags Code Atom Pair
         for (const auto& element : CMBags.BagTag){ 
            if (element == key) {
                mybagarray.push_back(CMBags.UCM[i]); 
            }
            i++;
         } 

         sizemax = MBags->second;
         //std::cout << key << ": found :" << mybagarray.size() << ", max : " << sizemax << "\n";

         std::vector<double> b;
         if (sizemax-mybagarray.size() > 0) {
           b.clear();
           b.resize(sizemax-mybagarray.size());
         }
         else {
           b.clear();
           b.resize(0);
         }
         if (mybagarray.size() > 0 )
         {
           // descending order
           std::sort(mybagarray.begin(), mybagarray.end());

           std::reverse(mybagarray.begin(), mybagarray.end());
       
           mybagarray.insert(std::end(mybagarray), std::begin(b), std::end(b));
         }
         else {
           mybagarray = b;
         }
        
        val.insert(std::end(val), std::begin(mybagarray), std::end(mybagarray));
  }

    res = val;
}



std::list<std::pair<std::string, unsigned int> > EnumeratesAtomsPair(const RDKit::ROMol &mol, 
    bool separateIsotopes, bool abbreviateHIsotopes) {
  std::map<std::string, unsigned int> counts;
  unsigned int nHs = 0;
  std::string key;
  const RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
  for (RDKit::ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    int atNum = (*atomIt)->getAtomicNum();
    key = table->getElementSymbol(atNum);
    if (separateIsotopes) {
      unsigned int isotope = (*atomIt)->getIsotope();
      if (abbreviateHIsotopes && atNum == 1 && (isotope == 2 || isotope == 3)) {
        if (isotope == 2)
          key = "D";
        else
          key = "T";
      }
    }
    if (counts.find(key) != counts.end()) {
      counts[key] += 1;
    } else {
      counts[key] = 1;
    }
    nHs += (*atomIt)->getTotalNumHs();
  }

  if (nHs) {
    key = "H";
    if (counts.find(key) != counts.end()) {
      counts[key] += nHs;
    } else {
      counts[key] = nHs;
    }
  }


  // create the Diagonal single atoms count

  std::list<std::pair<std::string, unsigned int> > ks;
  for (std::map<std::string, unsigned int>::const_iterator countIt = counts.begin();
       countIt != counts.end(); ++countIt) {
    // store the diagonal (individual atoms)
    std::pair<std::string, unsigned int > key = std::make_pair(countIt->first, countIt->second);
    ks.push_back(key);
  }


  unsigned int nbElements = counts.size();
  // enumerate Heterogenous Atom Pairs (simplely number product)
  // heterogens pairs
  for (unsigned int i = 0; i < nbElements - 1 ; i++) {
      std::map<std::string, unsigned int>::const_iterator it1 = counts.begin();
      std::advance(it1,i);
      std::pair<std::string, unsigned int> Values1 = *it1;
      for (unsigned int j = i+1; j < nbElements ; j++) {
        std::map<std::string, unsigned int>::const_iterator it2= counts.begin();
        std::advance(it2,j);
        std::pair<std::string, unsigned int> Values2 = *it2;
        std::string PairString;
        // sorting string order
        if (Values1.first > Values2.first) {
           PairString = Values2.first+Values1.first;
        }
        else {
            PairString = Values1.first+Values2.first;
        }
        unsigned int nbElementPair = Values1.second*Values2.second;
        std::pair<std::string,unsigned int > key = std::make_pair(PairString, nbElementPair);
        ks.push_back(key);
      }
  }


  for (std::map<std::string, unsigned int>::const_iterator countIt = counts.begin();
       countIt != counts.end(); ++countIt) {
    //"homologues":
    unsigned int nbElementAtom = countIt->second;
    unsigned int res =nbElementAtom*(nbElementAtom-1)/2;
    std::string PairString = countIt->first+countIt->first;
    if (res>0) {
      std::pair<std::string,unsigned int > key = std::make_pair(PairString, res);
      ks.push_back(key);
    }
  }
  

  return ks;

}



std::vector<std::string> tokenize(const std::string &s) {
        boost::char_separator<char> sep(", \n\r\t");
        boost::tokenizer<boost::char_separator<char>> tok(s, sep);
        std::vector<std::string> tokens;
        std::copy(tok.begin(), tok.end(), std::back_inserter<std::vector<std::string> >(tokens));
  return tokens;
}

bool startswith(std::string prefix, std::string text) {
    if (prefix.length() > 0 && text.length() > prefix.length()) {
        int i = 0;
        while (i < prefix.length()) {
            if (text[i] != prefix[i]) return false;
            i++;
        }
        return true;
    }

    else return false;
}

void testMOPAC2016outreader() {
  std::cout << "=>start test Coulomb\n";
  double x,y,z, atomcharge;
  std::vector<double> X,Y,Z, atomCharges;
  std::string pathName = getenv("RDBASE");

  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/benzene.out";

  std::ifstream instrm(fName.c_str());
  std::string line;
  std::vector<std::string> tokens;
  bool checkout = false;
  bool charge = false;

  RDKit::RWMol *mol = new RDKit::RWMol();
  int i=0;

  while (std::getline(instrm, line))  {
    if(startswith(" --", line)) {
          checkout = true;  
    }
    if (checkout) {
      if (std::strstr(line.c_str(),"  CARTESIAN COORDINATES") != NULL) {
        std::getline(instrm, line); // empty line
        std::getline(instrm, line); // next line

        tokens = tokenize(line); // get the list of words
        while (tokens.size() == 5) {
          // create the atoms of the molecule
          mol->addAtom(new RDKit::Atom(RDKit::PeriodicTable::getTable()->getAtomicNumber(tokens[1])));         // atom 0
          x = std::atof((char*)tokens[2].c_str());
          y = std::atof((char*)tokens[3].c_str());
          z = std::atof((char*)tokens[4].c_str());
          X.push_back(x);
          Y.push_back(y);
          Z.push_back(z);
          if (!std::getline(instrm, line))
              break;
          tokens = tokenize(line);
          i++;
        }
      }
    }
  }

  // GET THE Mulliken CHARGES FROM MOPAC
  std::ifstream instrm2(fName.c_str());
  while (std::getline(instrm2, line))  {
    if (std::strstr(line.c_str(),"  NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS") != NULL) {
      std::getline(instrm2, line); // empty line
      std::getline(instrm2, line); // next line
      if (std::strstr(line.c_str(),"CHARGE") != NULL) {
        for (int j=0;j<i;j++) {
          if (!std::getline(instrm2, line))
              break;
          tokens = tokenize(line); // get the list of words
          atomcharge = std::atof((char*)tokens[2].c_str());
          atomCharges.push_back(atomcharge);
          std::cout << atomcharge << ".";
        }  
      }   
    }
   }
  
  std::cout << "\n";
  

  // conformer definition 
  mol->addConformer(new RDKit::Conformer(mol->getNumAtoms()));
  // need a second loop for 3D points to the created conformer!
  for (int j=0;j<mol->getNumAtoms();j++) {
        mol->getConformer().setAtomPos(j, RDGeom::Point3D(X[j],Y[j],Z[j]));
  }
  // get the 3D distance from RDKit molecule!
  double* dist3D = RDKit::MolOps::get3DDistanceMat(*mol, -1, false, true);
  Eigen::MatrixXd Distance3D =  Eigen::Map< Eigen::MatrixXd>(dist3D, mol->getNumAtoms(), mol->getNumAtoms()); // convert the result array to matrix (1 column)

  std::cout << "3D distance:\n";
  std::cout << Distance3D << "\n";
  
}




void testEnumerateAtomPairs(){
    std::cout << "===================== Testing EnumerateAtomPairs =======================\n";
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCCN", "OCC(C)C(C)CN",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CCOC",
                           "c1ccccc1O",  "CCCl",     "CCBr",     "CCI",  "CCCBr",
                           "OC(=O)c1ccncc1C(=O)O",   "EOS"};
    unsigned int idx = 0;
    // caution fix size for the moment!
    std::vector< std::list < std::pair<std::string, unsigned int> > > myData(14);
    std::list<std::pair<std::string, unsigned int> > fulldata;
    while (sdata[idx] != "EOS") {
      RDKit::ROMol *mol;
      mol = RDKit::SmilesToMol(sdata[idx]);
      std::list<std::pair<std::string, unsigned int> > data = EnumeratesAtomsPair(*mol, false, false);
      myData.push_back(data);
      fulldata.merge(data);
      ++idx;
      delete mol;
    }

    std::map<std::string, unsigned int> Global;
    std::string key;

    for (auto const& i :  fulldata) {
        key= i.first;
        if (Global.find(key) != Global.end()) {
          if (i.second > Global[key]) {
            Global[key] = i.second;
          }   
        }         
        else {
            Global[key] = i.second;
        }
    }

    unsigned int resize = 0;
        std::cout << "===================== vis Master Pack ========================\n";

    for (auto const& i :  Global) {
        std::cout << i.first << ":" << i.second << ",";
        resize +=i.second;
    }
    std::cout << "\n";
    idx = 0;
    while (sdata[idx] != "EOS") {
          RDKit::ROMol *mol1;
          std::cout << sdata[idx] << "\n";
          mol1 = RDKit::SmilesToMol(sdata[idx]);
          RDKit::ROMol *mol = RDKit::MolOps::addHs(*mol1);
          int numatoms = mol->getNumAtoms();
          
          // Original distance geometry embedding
          RDKit::DGeomHelpers::EmbedMolecule(*mol, 0, 1234);
          RDKit::UFF::UFFOptimizeMolecule(*mol);

          // new Riniker and Landrum CSD-based method
          // using the parameters class
          //RDKit::DGeomHelpers::EmbedParameters params(RDKit::DGeomHelpers::ETKDG);
          //params.randomSeed = 1234;
          //RDKit::DGeomHelpers::EmbedMolecule(*mol2, params);
          
          double *dist3D = RDKit::MolOps::get3DDistanceMat(*mol, -1, false, true);

          std::vector<double>  res;

          getBagOfBonds(*mol, res, dist3D,  numatoms,  Global, 1);

          std::cout << "===================== vis bag of bond ========================\n";

          for (unsigned int i=0; i<res.size(); i++) {
              std::cout << res[i]  << ",";
          }
          std::cout << "\n";
          res.clear();
          res.resize(resize);
          ++idx;
          delete mol;
          delete mol1;

    }
} 


int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  testMOPAC2016outreader();
  testEnumerateAtomPairs();
}