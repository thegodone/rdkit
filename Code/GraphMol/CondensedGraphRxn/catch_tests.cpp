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

using namespace RDKit;


TEST_CASE("reactions", "[unittest][reactionclean]") {
  std::string sma = "C[C:3](=[O:4])[OH:5].[OH:1][CH:2](COC)CC>>C[C:3](=[O:4])[O:1][CH:2](COC)CC";
  SECTION("esterification") {
    std::string res = CondensedGraphRxn::RXNCompleteMapping(sma, false, false);
    CHECK(res == "[OH:1][CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11].[C:3](=[O:4])([OH:5])[CH3:6]>>[O:1]([CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11])[C:3](=[O:4])[CH3:6].[OH2:5]");
  }
  sma ="O=[CH:1][C@@H:2](O)[C@H:3](O)[C@H:4](O)[CH2:5]O.O=[CH:6][C@H:7](O)[C@@H:8]([OH:9])[C@H:10]([OH:11])[CH2:12][OH:13].O=[CH:14][C@H:15](O)[C@@H:16](O)[C@H:17](O)[C@H:18](O)[CH2:19][OH:20].O=[CH:21][C@H:22](O)[C@@H:23](O)[C@@H:24](O)[C@H:25](O)[CH2:26][OH:27].O=[C:28]([CH2:29]O)[C@@H:30](O)[C@H:31](O)[C@H:32](O)[CH2:33][OH:34].O[CH2:35][C@H]1O[C@@H]([C@@H]([C@H]([C@@H]1O)O)O)O[C@@H:36]1[C@H:37](O[C@H:38]([C@@H:39]([C@H:40]1O)O)O)[CH2:41]O.O[CH2:42][C@:43]1(O[C@@H:44]([C@H:45]([C@@H:46]1O)O)[CH2:47]O)O[C@H:48]1O[C@@H:49]([C@H:50]([C@@H:51]([C@H:52]1O)O)O)[CH2:53]O>>[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:48][CH2:52][CH2:51]/[CH:50]=[CH:49]\[CH2:53][CH2:35][CH2:14][CH2:15][CH2:16][CH2:17][CH2:18][C:19](=[O:20])[O:13][CH2:12][CH:10]([CH2:8][O:9][C:33](=[O:34])[CH2:32][CH2:31][CH2:30][CH2:28][CH2:29][CH2:51][CH2:50]/[CH:49]=[CH:53]\[CH2:41][CH2:37][CH2:36][CH2:40][CH2:39][CH2:38][CH2:6][CH3:7])[O:11][C:26](=[O:27])[CH2:25][CH2:24][CH2:23][CH2:22][CH2:21][CH2:38][CH2:39]/[CH:40]=[CH:36]\[CH2:37][CH2:41][CH2:47][CH2:44][CH2:45][CH2:46][CH2:43][CH3:42]";  
  SECTION("bug") {
	  std::string res = CondensedGraphRxn::RXNCompleteMapping(sma, false, false);
    CHECK(res == "[CH2:35]([OH:54])[C@H:55]1[O:56][C@H:57]([O:64][C@@H:36]2[C@@H:37]([CH2:41][OH:69])[O:65][C@@H:38]([OH:68])[C@H:39]([OH:67])[C@H:40]2[OH:66])[C@H:58]([OH:63])[C@@H:59]([OH:62])[C@@H:60]1[OH:61].[CH2:42]([C@:43]1([O:75][C@@H:48]2[C@H:52]([OH:77])[C@@H:51]([OH:78])[C@H:50]([OH:79])[C@@H:49]([CH2:53][OH:80])[O:76]2)[C@@H:46]([OH:72])[C@H:45]([OH:73])[C@@H:44]([CH2:47][OH:74])[O:71]1)[OH:70].[CH:14]([C@@H:15]([C@H:16]([C@@H:17]([C@@H:18]([CH2:19][OH:20])[OH:85])[OH:84])[OH:83])[OH:82])=[O:81].[CH:21]([C@@H:22]([C@H:23]([C@H:24]([C@@H:25]([CH2:26][OH:27])[OH:90])[OH:89])[OH:88])[OH:87])=[O:86].[C:28]([CH2:29][OH:92])([C@H:30]([C@@H:31]([C@@H:32]([CH2:33][OH:34])[OH:95])[OH:94])[OH:93])=[O:91].[CH:1]([C@H:2]([C@@H:3]([C@@H:4]([CH2:5][OH:100])[OH:99])[OH:98])[OH:97])=[O:96].[CH:6]([C@@H:7]([C@@H:8]([OH:9])[C@H:10]([OH:11])[CH2:12][OH:13])[OH:102])=[O:101]>>[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:48][CH2:52][CH2:51][CH:50]=[CH:49][CH2:53][CH2:35][CH2:14][CH2:15][CH2:16][CH2:17][CH2:18][C:19]([O:13][CH2:12][CH:10]([CH2:8][O:9][C:33]([CH2:32][CH2:31][CH2:30][CH2:28][CH2:29][CH2:51][CH2:50][CH:49]=[CH:53][CH2:41][CH2:37][CH2:36][CH2:40][CH2:39][CH2:38][CH2:6][CH3:7])=[O:34])[O:11][C:26]([CH2:25][CH2:24][CH2:23][CH2:22][CH2:21][CH2:38][CH2:39][CH:40]=[CH:36][CH2:37][CH2:41][CH2:47][CH2:44][CH2:45][CH2:46][CH2:43][CH3:42])=[O:27])=[O:20].[CH2:55]1[O:56][C@H:57]([OH:64])[C@H:58]([OH:63])[C@@H:59]([OH:62])[C@@H:60]1[OH:61].[OH2:96].[OH2:92].[OH2:97].[OH2:98].[OH2:99].[OH2:100].[OH2:101].[OH2:102].[OH2:81].[OH2:82].[OH2:94].[OH2:84].[OH2:85].[OH2:86].[OH2:87].[OH2:88].[OH2:89].[OH2:90].[OH2:91].[OH2:93].[OH2:71].[OH2:95].[OH2:54].[OH2:65].[OH2:66].[OH2:67].[OH2:68].[OH2:69].[OH2:70].[OH2:83].[OH2:72].[OH2:73].[OH2:74].[OH2:75].[OH2:76].[OH2:77].[OH2:78].[OH2:79].[OH2:80]");
  }

}


TEST_CASE("cgr", "[unittest][reactionclean]") {
  std::string sma = "[OH:1][CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11].[C:3](=[O:4])([OH:5])[CH3:6]>>[O:1]([CH:2]([CH2:7][O:8][CH3:9])[CH2:10][CH3:11])[C:3](=[O:4])[CH3:6].[OH2:5]";
  SECTION("cgr basic") {
    std::string cgr =  CondensedGraphRxn::CGRwriter(sma, false, 0, true, false, false, 0);
    CHECK(cgr == "CCC(COC)O{!-}C(C)(=O){-!}O");
  }
 SECTION("cgr signature radius 0") {
    std::string cgr =  CondensedGraphRxn::CGRwriter(sma, false, 0, true, true, false, 0);
    CHECK(cgr == "O{!-}C{-!}O");
  }
 SECTION("cgr signature radius 1") {
    std::string cgr =  CondensedGraphRxn::CGRwriter(sma, false, 0, true, true, false, 1);
    CHECK(cgr == "CO{!-}C(C)(=O){-!}O");
  }
}

