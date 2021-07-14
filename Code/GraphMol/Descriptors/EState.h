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

#include <RDGeneral/export.h>
#include <string>
#include <vector>

#ifndef ESTATERDKIT_H_APRIL2021
#define ESTATERDKIT_H_APRIL2021

namespace RDKit {
class ROMol;
namespace Descriptors {

const std::string EStateVersion = "1.0.0";
#ifdef RDK_BUILD_DESCRIPTORS3D

RDKIT_DESCRIPTORS_EXPORT void GetEStateTopographical(const ROMol &mol, std::vector<double> &res,
                                  int confId);    
#endif
    
RDKIT_DESCRIPTORS_EXPORT void GetEStateTopological(const ROMol &mol, std::vector<double> &res);
  
RDKIT_DESCRIPTORS_EXPORT void GetIState(const ROMol &mol, std::vector<double> &res, 
                                        bool Hs_default_one=false);
  
RDKIT_DESCRIPTORS_EXPORT void GetAtomsGroup(const ROMol &mol, std::vector<int> &res);
    
RDKIT_DESCRIPTORS_EXPORT void GetAtomsPeriod(const ROMol &mol, std::vector<int> &res);

}  // namespace Descriptors
}  // namespace RDKit
#endif


