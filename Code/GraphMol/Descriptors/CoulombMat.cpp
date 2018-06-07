//
//  Copyright (c) 2018, Guillaume GODIN
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "CoulombMat.h"
#include <random>
#include "slice.h"
#include "sort.h"
#include "MolData3Ddescriptors.h"
#include <Eigen/Dense>

using namespace Eigen;


struct BagMatrix {
  VectorXd UCM;
  std::vector<std::string > BagTag;
};

namespace RDKit {
namespace Descriptors {
namespace {

std::mt19937 mt;


std::unique_ptr<double[]> randnormal(const int nrolls) {

  std::normal_distribution<double> distribution(0.0,1.0);
  double *result = new double[nrolls];
  for (int i=0; i<nrolls; ++i) {
    result[i]=distribution(mt);
  }
  return std::unique_ptr<double[]>(result);
}

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

 MatrixXd getCoulombMat(VectorXd numbers, int numatoms, MatrixXd Distance3D,  int alpha, bool localdecay) {
  // 3D distance matrix strange result compare to initial data

  MatrixXd ProdV = numbers*numbers.transpose(); // outer products of vector ie (z[i] * z[j])

  if (localdecay) {
      // need to reshape two matrix (first row  Disctance3D rowwise and first column Disctance3D colwise )
      MatrixXd Rows = Distance3D.row(0);
      MatrixXd R = Rows.replicate(numatoms,1);
      MatrixXd C = Rows.transpose().replicate(1,numatoms);
      // add the distance to the first atom to the distance matrix
      Distance3D +=R;
      Distance3D +=C;
  }
   
  if (alpha != 1) {
    Distance3D = Distance3D.array().pow(alpha); // power distance alpha
  }

  MatrixXd MyVal = ProdV.cwiseQuotient(Distance3D); // ratio top / dist
  MyVal.diagonal() = 0.5 * numbers.array().pow(2.4); // set the diagonal fix values 

  return MyVal; 
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

 std::map<std::string, unsigned int> getBoBSizefromSmiles(std::vector<std::string> smiles) {
    unsigned int idx = 0;
    std::vector< std::list < std::pair<std::string, unsigned int> > > myData(smiles.size());
    std::list<std::pair<std::string, unsigned int> > fulldata;
    while (idx < smiles.size()) {
      RDKit::ROMol *mol;
      mol = RDKit::SmilesToMol(smiles[idx]);
      std::list<std::pair<std::string, unsigned int> > data = EnumeratesAtomsPair(*mol, false, false);
      myData.push_back(data);
      fulldata.merge(data); // durty need to compress not merge but... it works ;-)
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
    return Global;
}


void getBagOfBonds(const RDKit::ROMol &mol, std::vector<std::vector<double> > &res, double *dist3D, 
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

    res.push_back(val);
}





MatrixXd permuteEM(MatrixXd Mat, MatrixXi IX){
  // permutation can take an index 
  PermutationMatrix<Dynamic,Dynamic> perm(IX);
  MatrixXd PEM = perm.transpose() * Mat * perm;
  return PEM;
}


void getLocalCoulombMats(const RDKit::ROMol &mol, std::vector< std::vector<double> > &res, 
      double *dist3D,  int alpha, double rcut, bool localdecay, bool reduced, int nbmats) {

    int numatoms = mol.getNumAtoms();
    double *z = new double[numatoms];

    for (int i=0; i< numatoms; i++){
          z[i] = mol.getAtomWithIdx(i)->getAtomicNum();
    }
    MatrixXd Distance3D = Map<MatrixXd>(dist3D, numatoms, numatoms); // convert the result array to matrix (1 column)

    for (int i = 0; i < numatoms; i++) {

      // get the colum and sort it (for increase distance order)
      MatrixXd Row = Distance3D.col(i);
      MatrixXd Y;
      MatrixXi IX;
      igl::sort(Row,1,true,Y,IX); // sort and find index

      int localatomsnumber=1; 
      for (int j=1; j < numatoms; j++) {
         if (Y(j)>rcut or j>nbmats-1) { // stop if than nbmats neighbours atoms already seen!!!
           localatomsnumber = j;
           break;
         }
       } 
 
      // there are pos rcut "pos" keep atoms
      ArrayXi ri(localatomsnumber); 
      double *localz = new double[localatomsnumber];

      for (int j=0; j < localatomsnumber; j++) {
          ri(j)=IX(j);
          localz[j]=z[IX(j)];
      }
 

      // definition of submatrix size
      MatrixXd subMat(localatomsnumber,localatomsnumber);
      igl::slice(Distance3D, ri, ri, subMat);
  
      VectorXd localnumbers = Map<VectorXd>(localz, localatomsnumber); // convert the number array to vector
   
      MatrixXd localCoulomb(localatomsnumber,localatomsnumber);

      localCoulomb = getCoulombMat( localnumbers, localatomsnumber, subMat, alpha, localdecay);

      if (reduced) {

        MatrixXd MRES = MatrixXd::Zero(nbmats,nbmats);
        MRES.block(0,0,localatomsnumber,localatomsnumber) = localCoulomb; // CAUTION: need to add a check on the side of MRES lower than the FRES

        // block( startRow,  startCol,  blockRows,  blockCols)
        MatrixXd RES = MatrixXd::Zero(1,2*nbmats-1); // it's already a row vector!
        RES.block(0,0,1,nbmats) = MRES.row(0);
        RES.block(0,nbmats,1,nbmats-1) = MRES.block(1,1,nbmats,nbmats).diagonal().transpose();
        std::vector<double> lcm(RES.data(), RES.data() + RES.rows() * RES.cols()); // return the matrix to vectors
        res[i] = lcm;

      }
      if (!reduced) {
        // force the first row to be the larger norm ie Inf!
        double storelc = localCoulomb(0,0);
        MatrixXd EigenCMnorm = localCoulomb.rowwise().norm(); // compute row norms
        EigenCMnorm(0,0) /= 0; // set position one inchanged by applying infinity value
        MatrixXd Yn;
        MatrixXi IXn;

        // sort ascending norm to get permutation indexes
        igl::sort(EigenCMnorm,1,false,Yn,IXn); // caution sort descending see article
        
        MatrixXd FRES=permuteEM(localCoulomb, IXn);
        // expending the matrix using block copy operation P.block(i, j, rows, cols)
        // for extending the matrix than taking the upper vector the m*(m+1)/2 

        MatrixXd MRES = MatrixXd::Zero(nbmats,nbmats);
        MRES.block(0,0,localatomsnumber,localatomsnumber) = FRES; // CAUTION: need to add a check on the side of MRES lower than the FRES

        VectorXd RES = getUpperMat(MRES);

        std::vector<double> lcm(RES.data(), RES.data() + RES.rows() * RES.cols()); // return the matrix to vectors
        
        res[i] = lcm;
      }

    }

}

void getRandCoulombMats(const RDKit::ROMol &mol, std::vector< std::vector<double> > &res, double *dist3D, 
unsigned int  numAtoms,  int nbmats, int seed, int padding, bool localdecay, bool sorted, bool eigenval, int alpha) {

  mt.seed(seed);
  // initialize the variables for the getCoulombMat
  int numatoms= mol.getNumAtoms();
  double *z = new double[numatoms];

  for (int i=0; i< numatoms; i++){
      z[i] = mol.getAtomWithIdx(i)->getAtomicNum();
  }

  VectorXd numbers = Map<VectorXd>(z, numatoms); // convert the number array to vector
  
  MatrixXd Distance3D = Map<MatrixXd>(dist3D, numatoms, numatoms); // convert the result array to matrix (1 column)

  MatrixXd CMM = getCoulombMat(numbers, numatoms, Distance3D, alpha, localdecay); // not sure localdecay is required there!

  MatrixXd EigenCMnorm = CMM.rowwise().norm(); // compute row norms
    
  for (unsigned int j=0; j< nbmats; j++) {
    // prepare result of sort matrix using igl header
    MatrixXd Y;
    MatrixXi IX;
    if (sorted) {
    
      igl::sort(EigenCMnorm,1,false,Y,IX); // sort and find index

    }
    else {
      std::unique_ptr<double[]> result = randnormal(numAtoms); // check if this is always the same or not ;-)

      MatrixXd Eigenresult = Map<MatrixXd>(result.get(), numAtoms, 1); // convert the result array to matrix (1 column)

      // sort ascending norm+result and get the indexes
      // need to add this in the compilation ... headers
      igl::sort(EigenCMnorm+Eigenresult,1,true,Y,IX); // sort and find index
    }

    
    // last step : permute the matrix using the IX indexes (only is seed>=0 else return the orignal matrix!)
    MatrixXd RRES=permuteEM(CMM, IX);

    // we can add padding using block method too see the local example
    // caution padding >= max(numAtoms)
    MatrixXd MRES = MatrixXd::Zero(padding,padding);

    MRES.block(0,0,numAtoms,numAtoms) = RRES; // CAUTION: need to add a check on the side of MRES lower than the FRES

    VectorXd RES(padding);
    if (eigenval) {
      // get eigen vectors 
      Eigen::EigenSolver<MatrixXd> es(MRES);

      RES = es.eigenvalues().real();

    }
    // than we compress symmetric matrxi using upper triangle into a vector
    else  {
       VectorXd RES1= getUpperMat(MRES);
       RES.resize(RES1.size());
       RES = RES1;
    }
    //VectorXd RES= getUpperMat(MRES);
    
    std::vector<double> rcm(RES.data(), RES.data() + RES.rows() * RES.cols());

    
    res[j] = rcm;
  }
}



void getCoulombMats(const ROMol &mol, std::vector< std::vector<double> > &res,
              unsigned int numAtoms, int confId, unsigned int nbmats,  int alpha , int seed, int padding,
               double rcut, bool local, bool decaying, bool reduced, bool sorted, bool eigenval, bool bob, 
               std::vector<std::string> smiles) {
  // 3D distance matrix
  double *dist3D = MolOps::get3DDistanceMat(mol, confId, false, true);

  if (local) {
    res.clear();
    res.resize(numAtoms);
    getLocalCoulombMats(mol, res, dist3D, alpha, rcut, decaying, reduced, nbmats);
  } 
  else if (!bob) { 
    res.clear();
    res.resize(nbmats);
    getRandCoulombMats(mol, res, dist3D, numAtoms, nbmats, seed, padding, decaying, sorted, eigenval, alpha);

  }
  else if (bob) {
    res.clear();
    res.resize(1);
    std::cout << ".";
    std::map<std::string, unsigned int>  MaxBags = getBoBSizefromSmiles(smiles);
    getBagOfBonds(mol, res, dist3D, numAtoms,  MaxBags, alpha);

  }
}
}  // end of anonymous namespace

void CoulombMat(const ROMol &mol, std::vector<std::vector<double>>  &res, int confId,  int nbmats,
 int seed, int padding, double rcut, bool local, bool decaying, bool reduced, bool sorted, bool eigenval,
  int alpha, bool bob, std::vector<std::string>  smiles) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  unsigned int numAtoms = mol.getNumAtoms();
  // check that here or in the next call ?
  res.clear();
  res.resize(nbmats);
  bob = true;
  getCoulombMats(mol, res, numAtoms, confId, nbmats, alpha, seed, padding, rcut, local, decaying, reduced, sorted, eigenval, bob, smiles);
}
}  // namespace Descriptors
}  // namespace RDKit
