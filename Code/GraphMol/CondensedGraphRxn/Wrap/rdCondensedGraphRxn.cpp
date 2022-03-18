//
//  Copyright (C) 2020 Guillaume GODIN @  Firmenich
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include <GraphMol/ChemReactions/Reaction.h>
#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/CondensedGraphRxn/CondensedGraphRxn.h>

namespace python = boost::python;
using namespace RDKit;
using namespace RDKit::CondensedGraphRxn;

namespace {

std::string pyObjectToString(python::object input) {
  python::extract<std::string> ex(input);
  if (ex.check()) {
    return ex();
  }
  std::wstring ws = python::extract<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

std::string CRSWriter(python::object ismart, bool doRandom , unsigned int randomSeed,
		bool aromatize, bool signature, bool charges, int radius, bool addRingInfo)  {
  const std::string smart = pyObjectToString(ismart);
  return CRSwriter( smart , doRandom, randomSeed, aromatize, signature, charges, radius, addRingInfo);
}

std::string CRSReader(python::object icrs, bool canonical, bool setAtomMap ) {
  const std::string crs = pyObjectToString(icrs);
  RDKit::RWMol *molR;
  std::string res;
  try {
    molR = SmilesToMol(crs, 0, false);
  } catch (...) {
    molR = nullptr;
  }
  if (molR) {
    res = CRSreader( molR , crs, canonical, setAtomMap );
    delete molR;
  }
  else {
    res =  "CRS SmilesParse error";
  }

  return res;
}

std::string RxnCompleteMapping(python::object  ismiles, bool debug, bool addleavinggroups){	
const std::string smiles = pyObjectToString(ismiles);
return  RXNCompleteMapping( smiles, debug, addleavinggroups);
}

bool GetRXNComp(python::object  irxnsma, python::object  irxncoresma ){
        
        const std::string rxnsma = pyObjectToString(irxnsma);
        const std::string rxncoresma = pyObjectToString(irxncoresma);
        bool res = getRXNComparison(rxnsma, rxncoresma);
        return res;
    }

bool GetRXNTotalComp(python::object  irxnsma, python::object  irxncoresma ){
        
        const std::string rxnsma = pyObjectToString(irxnsma);
        const std::string rxncoresma = pyObjectToString(irxncoresma);
        bool res = getRXNCompTotal(rxnsma, rxncoresma);
        return res;
    }
    
std::string SetRXNAtomMap(python::object  irxnsma, python::object  irxncoresma){
        const std::string rxnsma = pyObjectToString(irxnsma);
        const std::string rxncoresma = pyObjectToString(irxncoresma);
        return setRXNCompAtomMaps(rxnsma, rxncoresma);
    
    }

}  // namespace

BOOST_PYTHON_MODULE(rdCondensedGraphRxn) {
  python::scope().attr("__doc__") =
      "Module containing functions for creating a Scaffold Network";

  python::def("CRSwriter", &CRSWriter,
              (python::arg("smart"), python::arg("doRandom") = false,
	       python::arg("randomSeed")= 0, python::arg("aromatize") = true,
	       python::arg("signature") = false,
               python::arg("charges") = false,
               python::arg("radius") = 1,
	       python::arg("addRingInfo") = true),
              "convert a rxnsmarts into a crs");

   python::def("CRSreader", &CRSReader,
              (python::arg("crs"), python::arg("canonical") = true,
	       python::arg("setAtomMap") = true),
              "convert a crs into a rxnsmarts");

      python::def("RxnCompleteMapping", &RxnCompleteMapping,
              (python::arg("smiles"), python::arg("debug")=false,
               python::arg("addleavinggroups")=false),
              "Complete Mapping of RXN string");
    
    python::def("GetRXNComp", &GetRXNComp,
                (python::arg("rxnsma"), python::arg("rxnsmacore")),
                "compare reaction with core reaction site");

    python::def("GetRXNTotalComp", &GetRXNTotalComp,
                (python::arg("rxnsma"), python::arg("rxnsmacore")),
                "compare reaction with core reaction site");

    python::def("SetRXNAtomMap", &SetRXNAtomMap,
                (python::arg("rxnsma"), python::arg("rxnsmacore")),
                "fix reaction mapping using correct core reaction site mapped");

}
