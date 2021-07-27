//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"
#include "RDGeneral/test.h"
#include <sstream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/CondensedGraphRxn/CondensedGraphRxn.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/RWMol.h>

using namespace RDKit;


TEST_CASE("reactions", "[unittest][reactionclean]") {
  std::string sma = "C[C:3](=[O:4])[OH:5].[OH:1][CH:2](COC)CC>>C[C:3](=[O:4])[O:1][CH:2](COC)CC";
  SECTION("esterification") {
    std::string res = CondensedGraphRxn::RXNCompleteMapping(sma, false, false);
    CHECK(res == "[OH:1][CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11].[C:3](=[O:4])([OH:5])[CH3:6]>>[O:1]([CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11])[C:3](=[O:4])[CH3:6].[OH2:5]");
  }
  sma ="C[C:3](=[O:4])O.[OH:1][CH:2](COC)CC>>C[C:3](=[O:4])[O:1][CH:2](COC)CC";  
  SECTION("test1") {
	  std::string res = CondensedGraphRxn::RXNCompleteMapping(sma, false, false);
    CHECK(res == "[OH:1][CH:2]([CH2:6][O:7][CH3:8])[CH2:9][CH3:10].[CH:3](=[O:4])[CH3:5]>>[O:1]([CH:2]([CH2:6][O:7][CH3:8])[CH2:9][CH3:10])[C:3](=[O:4])[CH3:5]");
  }

}


TEST_CASE("cgr", "[unittest][reactionclean]") {
  std::string sma = "[OH:1][CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11].[C:3](=[O:4])([OH:5])[CH3:6]>>[O:1]([CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11])[C:3](=[O:4])[CH3:6].[OH2:5]";
  SECTION("cgr basic") {
    std::string cgr =  CondensedGraphRxn::CRSwriter(sma, false, 0, true, false, false, 0);
    CHECK(cgr == "CCC(COC)O{!-}C(C)(=O){-!}O");
  }
 SECTION("cgr signature radius 0") {
    std::string cgr =  CondensedGraphRxn::CRSwriter(sma, false, 0, true, true, false, 0);
    CHECK(cgr == "O{!-}C{-!}O");
  }
 SECTION("cgr signature radius 1") {
    std::string cgr =  CondensedGraphRxn::CRSwriter(sma, false, 0, true, true, false, 1);
    CHECK(cgr == "CO{!-}C(C)(=O){-!}O");
  }
}

TEST_CASE("crs", "[unittest][reader]") {
  std::string crs = "CC{!-}CO";
  SECTION("crs reader") {

  RDKit::RWMol *molR;
  std::string res;
  try {
    molR = SmilesToMol(crs, 0, false);
  } catch (char *excp) {
    std::cout << "Caught " << excp;
    molR = nullptr;
  }
  catch (...) {
        std::cout << "Default Exception\n";
  }

  if (molR) {
    res = CondensedGraphRxn::CRSreader( molR , crs, false, false );
    delete molR;
  }
  else {
    res =  "CRS SmilesParse error";
  }

    std::cout << "rxn crs test\n" << res;
    CHECK(res  == "CC.CO>>CCCO");
  }
}




