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
#include "CoulombMat.h"
#include <random>
#include "sort.h"
#include "MolData3Ddescriptors.h"
#include <Eigen/Dense>

using namespace Eigen;

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

std::unique_ptr<double[]> getCoulombMat(const RDKit::ROMol &mol, double *dist3D) {
  // 3D distance matrix strange result compare to initial data
  int numatoms= mol.getNumAtoms();
  double *cm = new double[numatoms*numatoms];
  double *z = new double[numatoms];

  for (int i=0; i< numatoms; i++){
      z[i] = mol.getAtomWithIdx(i)->getAtomicNum();
  }

  for (int i=0; i<numatoms; i++) {
      for (int j=0; j<numatoms; j++) {
          if (i == j){
              cm[i+j*numatoms] = 0.5 * pow (z[i],  2.4);
          }
          else if (i < j){
              cm[i+j*numatoms] = (z[i] * z[j]) / dist3D[i+j*numatoms];
              cm[j+i*numatoms] = cm[i+j*numatoms];
          }
      }
  }
  return std::unique_ptr<double[]>(cm);
}

Eigen::MatrixXd permuteEM(Eigen::MatrixXd Mat, Eigen::MatrixXi IX){
  // permutation can take an index 
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(IX);
  Eigen::MatrixXd PEM = perm.transpose() * Mat * perm;
  return PEM;
}

void getRandCoulombMats(const RDKit::ROMol &mol, std::vector< std::vector<double> > &res, double *dist3D, unsigned int  numAtoms,  int nbmats, int seed) {

  mt.seed(seed);

  for (unsigned int j=0; j< nbmats; j++) {

    std::unique_ptr<double[]> result = randnormal(numAtoms); // check if this is always the same or not ;-)

    Eigen::MatrixXd Eigenresult = Eigen::Map<Eigen::MatrixXd>(result.get(), numAtoms, 1); // convert the result array to matrix (1 column)

    std::unique_ptr<double[]> CM = getCoulombMat(mol, dist3D); // get the non permuted Coulomb Matrix
    Eigen::MatrixXd CMM = Eigen::Map<Eigen::MatrixXd>(CM.get(), numAtoms, numAtoms);

    Eigen::MatrixXd EigenCMnorm = CMM.rowwise().norm(); // compute row norms
    // prepare result of sort matrix using igl header
    Eigen::MatrixXd Y;
    Eigen::MatrixXi IX;

    // sort ascending norm+result and get the indexes
    // need to add this in the compilation ... headers
    igl::sort(EigenCMnorm+Eigenresult,1,true,Y,IX); // sort and find index

    // last step : permute the matrix using the IX indexes (only is seed>=0 else return the orignal matrix!)
    Eigen::MatrixXd RES=permuteEM(CMM, IX);
    std::vector<double> rcm(RES.data(), RES.data() + RES.rows() * RES.cols());
    
    res[j] = rcm;
  }
}

void getCoulombMats(const ROMol &mol, std::vector< std::vector<double> > &res,
              unsigned int numAtoms, int confId, unsigned int nbmats, int seed) {
  // 3D distance matrix
  double *dist3D = MolOps::get3DDistanceMat(mol, confId, false, true);

  res.clear();
  res.resize(nbmats);

  getRandCoulombMats(mol, res, dist3D, numAtoms, nbmats, seed);
}


}  // end of anonymous namespace

void CoulombMat(const ROMol &mol, std::vector<std::vector<double>>  &res, int confId,  int nbmats, int seed) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  unsigned int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(nbmats);

  getCoulombMats(mol, res, numAtoms, confId, nbmats, seed);
}
}  // namespace Descriptors
}  // namespace RDKit
