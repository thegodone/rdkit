//
//  Copyright (c) 2021, Guillaume GODIN
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
//

#include <GraphMol/RDKitBase.h>
#include <cmath>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "ACSF.h"
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/Core>

using namespace Eigen;
namespace RDKit {
namespace Descriptors {

void ACSF(const ROMol &mol, std::vector<std::vector<float>>  &res, std::vector<unsigned int> Zi, int confId, double rCut) {    
    PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers");

    // we only consider the G1,G2 & G4 functions (now G3 "cos" & G5 "simpler G4")
    // the g2/g4 parameters are statics ie the output dimension is static too!!!
    std::vector<std::vector<double>> g2Params{{1.,1.},{1.,2.},{1.,3.}};
    std::vector<std::vector<double>> g4Params{{1.,1.,1.},{1.,2.,1.},{1.,1.,-1.},{1.,2.,-1.}}; 
    

    int numAtoms = mol.getNumAtoms();
    // get the atomicNumbers of each atoms
    std::vector<unsigned int> atomicNumbers;
    // get the uniques atomicNumbers
    for( int i=0;i<numAtoms;++i){
         const auto at = mol.getAtomWithIdx(i);
         atomicNumbers.push_back((unsigned int)at->getAtomicNum());  
    }
    
    // define the output dimension (static)
    int nACSF = Zi.size();
    int nACSFPairs = nACSF*(nACSF+1)/2;
    
    int nG2 = g2Params.size();
    int nG4 = g4Params.size();
    
    // only compute G1,2 & 4 (don't consider cos G3 and G5 ACSF)
    int numcols = ((1+nG2)*nACSF+(nG4)*nACSFPairs);
    
    // create the zero matrix output
    MatrixXf results = MatrixXf::Zero(numAtoms, numcols);
    
    // create mapping "atomic Zi & index" in order to get correct column position for accumulation computation
    std::unordered_map<unsigned int, int> atomicNumberMap;  
    int i = 0;
    for (unsigned int Z : Zi) {
        atomicNumberMap[Z] = i;
        ++i;
    }  
    
    // get the distance matrix: rdkit is in double we need to cast after in float
    // precision will be "soso" between pure "float and double"
    double* dist3D =
       MolOps::get3DDistanceMat(mol, confId, false, true);  // 3D distance matrix
    Map<MatrixXd> dm(dist3D, numAtoms, numAtoms);

    // Calculate the symmetry function values for every specified atom    
    for( int i=0;i<numAtoms;++i){
        for( int j=0;j<numAtoms;++j){

            if (i == j) {continue;}
                float r_ij = (float)dm(i,j);
                if (r_ij<= rCut) {
                    float fc_ij = 0.5*(cos(r_ij*M_PI/rCut)+1);
                    // get proper offset using mapping to accumulate G(z,i) "sum over all atoms having same "z" value"
                    int idx_j = atomicNumberMap[atomicNumbers[j]];
                    
                    int offset = idx_j * (1+nG2);  // Skip G1, G2 types that are not the ones of atom

                    // first accumulation for G1
                    results(i,offset) +=fc_ij;   // G1
                    // go to next position
                    offset++;

                    float eta;
                    float Rs;                
                    // G2 (there is several eta/Rs pair to compute... inside)
                    for (auto params : g2Params) {
                        eta = params[0];
                        Rs = params[1];
                        // second accumulation for G2
                        results(i,offset) += exp(-eta * (r_ij - Rs)*(r_ij - Rs)) * fc_ij;   // G2
                        // go to next position
                        offset++;

                    }
                     // G4 (there is several eta/zeta/lambda triplets to compute... inside)
                     for( int k=0;k<numAtoms;++k){
                         if (k==i || k>=j) { continue;}
                         

                         float r_ik =  (float)dm(i,k);

                         if (r_ik <= rCut) {
                             float r_jk =  (float)dm(j,k);

                             float fc_ik = 0.5*(cos(r_ik*M_PI/rCut)+1);

                             float r_ij_square =  r_ij*r_ij;
                             float r_ik_square =  r_ik*r_ik;
                             float r_jk_square =  r_jk*r_jk;
                             int idx_k = atomicNumberMap[atomicNumbers[k]];
                             float costheta = 0.5/(r_ij*r_ik) * (r_ij_square+r_ik_square-r_jk_square);

                             int its;
                             if (idx_j >= idx_k) {
                                 its = (idx_j*(idx_j+1))/2 + idx_k;
                             } else  {
                                 its = (idx_k*(idx_k+1))/2 + idx_j;
                             }
                             
                             offset = nACSF * (1+nG2); // Skip this atoms G1 G2
                             offset += its * (nG4);    // skip G4 type that are not the ones of atom

                             if (r_jk> rCut) {  offset += nG4; } // don't accumulate over the rcut distance 
                             else {
                                float fc_jk = 0.5*(cos(r_jk*M_PI/rCut)+1);
                                float fc4 = fc_ij*fc_ik*fc_jk;

                                float eta;
                                float zeta;
                                float lambda;
                                float gauss;
                                for (auto params : g4Params) {
                                    eta = params[0];
                                    zeta = params[1];
                                    lambda = params[2];
                                    gauss = exp(-eta*(r_ij_square+r_ik_square+r_jk_square)) * fc4; 
                                    results(i,offset) += 2*pow(0.5*(1 + lambda*costheta), zeta) * gauss; // G4
                                    offset++;
                                    if (i==0) {
                                    }
                                }
                             }
                         }
                    }
                }
            }  
       }
    
      // return output into vector of vector of double 
      // define number of rows
      res.resize(numAtoms);
        
      for( int i=0;i<numAtoms;++i) {
            VectorXf resvecXf =  results.row(i);
            std::vector<float> resvec(resvecXf.data(), resvecXf.data() + numcols);
            // define number of cols
            res[i].resize(numcols); 
            // assign Eigen vector to results
            res[i] = resvec;
      }
  }

}  // namespace Descriptors
}  // namespace RDKit
