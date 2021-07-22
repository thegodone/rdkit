#
# Copyright (C) 2020 Guillaume Godin @ Firmenich
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
from rdkit.Chem.CondensedGraphRxn import rdCondensedGraphRxn as CGR

rdBase.DisableLog("rdApp.info")

class TestCondensedGraphRxn(unittest.TestCase):


  def setUp(self):
    self.cgrs = ['C{-:}O' ,
    'C{:=}O' ,
    'C{:#}O' ,
    'C{#=}O' ,
    'C{!#}O' ,
    'C{!:}O' ,
    'C{#!}O' ,
    'C{#-}O' ,
    'C{=!}O' ,
    'C{=-}O' ,
    'C{:-}O' ,
    'C{-=}O' ,
    'C{!-}O' ,
    'C{-!}O' ,
    'C{-#}O' ,
    'C{:!}O' ,
    'C{#:}O' ,
    'C{=:}O' ,
    'C{!=}O',
    'C{=#}O']
    
    self.smarts = ['[CH3:1][OH:2]>>[cH2:1][o:2]',
    '[cH2:1][o:2]>>[CH2:1]=[O:2]',
    '[cH2:1][o:2]>>[CH:1]#[O:2]',
    '[CH:1]#[O:2]>>[CH2:1]=[O:2]',
    '[CH4:1].[OH2:2]>>[CH:1]#[O:2]',
    '[CH4:1].[OH2:2]>>[cH2:1][o:2]',
    '[CH:1]#[O:2]>>[CH4:1].[OH2:2]',
    '[CH:1]#[O:2]>>[CH3:1][OH:2]',
    '[CH2:1]=[O:2]>>[CH4:1].[OH2:2]',
    '[CH2:1]=[O:2]>>[CH3:1][OH:2]',
    '[cH2:1][o:2]>>[CH3:1][OH:2]',
    '[CH3:1][OH:2]>>[CH2:1]=[O:2]',
    '[CH4:1].[OH2:2]>>[CH3:1][OH:2]',
    '[CH3:1][OH:2]>>[CH4:1].[OH2:2]',
    '[CH3:1][OH:2]>>[CH:1]#[O:2]',
    '[cH2:1][o:2]>>[CH4:1].[OH2:2]',
    '[CH:1]#[O:2]>>[cH2:1][o:2]',
    '[CH2:1]=[O:2]>>[cH2:1][o:2]',
    '[CH4:1].[OH2:2]>>[CH2:1]=[O:2]',
    '[CH2:1]=[O:2]>>[CH:1]#[O:2]',]
    return self

  def test1Basicsreader(self):
    i=0
    for s,c in zip(self.smarts, self.cgrs):
        self.assertEqual( s , CGR.CGRreader(c))
        i+=1

  def test2Basicswriter(self):
    i=0
    for s,c in zip(self.smarts, self.cgrs):
        self.assertEqual( c , CGR.CGRwriter(s))
        i+=1


  def test3Acetal(self):
    sma = ('[CH3:1][C:2]([CH3:3])=[O:4].[OH:5][CH2:6][CH2:7][OH:8]>>[CH3:1][C:2]1([CH3:3])[O:5][CH2:6][CH2:7][O:8]1.[OH2:4]')
    cgr = 'CC1(C)({=!}O){!-}OCCO{!-}1' 
    self.assertEqual (cgr, CGR.CGRwriter(sma))
    self.assertEqual (sma, CGR.CGRreader(cgr))    

  def test4Radius(self):
    sma = ('[CH3:1][C:2]([CH3:3])=[O:4].[OH:5][CH2:6][CH2:7][OH:8]>>[CH3:1][C:2]1([CH3:3])[O:5][CH2:6][CH2:7][O:8]1.[OH2:4]')
    cgr = 'CC1(C)({=!}O){!-}OCCO{!-}1'
    self.assertEqual (cgr, CGR.CGRwriter(sma))

  def testCannotParse(self):
    cgr = 'C{!-}OC(=O)CC1=C(B-!}Br)C=C(F)C=C1'
    sma = 'error'
    self.assertEqual (sma, CGR.CGRreader(cgr))

  def testCanParse(self):
    cgr = 'C{!-}OC(=O)CC1=C({-!}Br)C=C(F)C=C1'
    sma = '[CH4:1].[OH:2][C:3](=[O:4])[CH2:5][c:6]1[c:7]([Br:8])[cH:9][c:10]([F:11])[cH:12][cH:13]1>>[BrH:8].[CH3:1][O:2][C:3](=[O:4])[CH2:5][c:6]1[cH:7][cH:9][c:10]([F:11])[cH:12][cH:13]1'
    self.assertEqual (sma, CGR.CGRreader(cgr))

  def testAromatizeOnBenzeneReduction(self):
    sma = '[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[CH2:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1'
    cgr = 'C1{:-}C{:-}C{:-}C{:-}C{:-}C{:-}1'
    self.assertEqual (cgr, CGR.CGRwriter(sma))
    print(CGR.CGRreader(cgr))

  def testSignature(self):
    sma = '[CH3:1][O:2][CH3:3].[CH3:4]>>[CH2:1][O:2][CH3:3][CH3:4]'
    print("Testing signature CGR for %s"%(sma))
    cgr0 = 'C{!-}C'
    cgr1 = 'C{!-}CO'
    cgr2 = 'COC{!-}C'
    self.assertEqual (cgr0, CGR.CGRwriter(sma, signature=True, radius = 0))
    self.assertEqual(cgr1, CGR.CGRwriter(sma, signature=True, radius = 1))
    self.assertEqual(cgr2, CGR.CGRwriter(sma, signature=True, radius = 2))
    print("CGR r=0: %s"%(CGR.CGRwriter(sma, signature = True,radius=0)))
    print("CGR r=1: %s"%(CGR.CGRwriter(sma,signature=True,radius=1)))
    print("CGR r=2: %s"%(CGR.CGRwriter(sma,signature=True,radius=2)))

  def testAldolCondensation(self):
      sma = "[CH3:10][CH:9]([CH3:11])[c:6]1[cH:7][cH:8][c:3]([CH:2]=[O:1])[cH:4][cH:5]1.[CH3:15][CH2:14][CH:13]=[O:12]>>[CH3:10][CH:9]([CH3:11])[c:6]1[cH:7][cH:8][c:3](\[CH:2]=[C:14](/[CH3:15])[CH:13]=[O:12])[cH:4][cH:5]1"
      cgr = "CC(C=O){!=}C({=!}O)c1ccc(C(C)C)cc1"
      cgr0 = "C{!=}C{=!}O"
      cgr1 = "CC(C){!=}C(c){=!}O"
      cgr2 = "cc(c)C({=!}O){!=}C(C)C=O"
      cgr3 = "ccc(cc)C({=!}O){!=}C(C)C=O"
      sign = [cgr0,cgr1,cgr2,cgr3]

      # Test the full molecule
      print("Testing aldol condensation with signature-CGR:")
      self.assertEqual(cgr,CGR.CGRwriter(sma))
      for r,expt in zip(list(range(4)),sign):
          self.assertEqual(expt,CGR.CGRwriter(sma,signature=True,radius=r))
          print(CGR.CGRwriter(sma,signature=True,radius=r))

    
  def testCase1(self,display=True):
      r = 0      
      print("test reaction complete mapping")
      # Check the reverse/forward case
      sma = "C[C:3](=[O:4])O.[OH:1][CH:2](COC)CC>>C[C:3](=[O:4])[O:1][CH:2](COC)CC"
      expt = "[OH:1][CH:2]([CH2:6][O:7][CH3:8])[CH2:9][CH3:10].[CH:3](=[O:4])[CH3:5]>>[O:1]([CH:2]([CH2:6][O:7][CH3:8])[CH2:9][CH3:10])[C:3](=[O:4])[CH3:5]"      
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase2(self,display=True):
      # Check the forward/reverse case
      sma = "CC[CH:2](COC)[O:1][C:3](C)=[O:4]>>C[C:3]([OH:5])=[O:4].CC[CH:2]([OH:1])COC"
      expt = "[O:1]([CH:2]([CH2:7][CH3:6])[CH2:8][O:9][CH3:10])[C:3](=[O:4])[CH3:11].[OH2:5]>>[C:3](=[O:4])([OH:5])[CH3:11].[OH:1][CH:2]([CH2:7][CH3:6])[CH2:8][O:9][CH3:10]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
  
  def testCase3(self,display=True):    
      # Check the SN2: The empty case
      sma = "[Cl:1][CH2:2]COc1ccccc1[CH2:3][Cl:4]>>[I:5][CH2:2]COc1ccccc1[CH2:3][I:6]"
      expt = "[Cl:1][CH2:2][CH2:7][O:8][c:9]1[cH:10][cH:11][cH:12][cH:13][c:14]1[CH2:3][Cl:4].[IH:6].[IH:5]>>[CH2:2]([I:5])[CH2:7][O:8][c:9]1[cH:10][cH:11][cH:12][cH:13][c:14]1[CH2:3][I:6].[ClH:4].[ClH:1]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)

  def testCase4(self,display=True):
      # Check the aldol condensation with 2x the same input molecule
      sma = "CC[CH2:3][C:4](C)=[O:5].CCC[C:8]([CH3:10])=[O:9]>>CCC[C:8]([CH3:10])=[C:3](CC)[C:4](C)=[O:5].[OH2:9]"
      expt = "[CH2:3]([C:4](=[O:5])[CH3:16])[CH2:14][CH3:15].[C:8](=[O:9])([CH3:10])[CH2:13][CH2:12][CH3:11]>>[C:3]([C:4](=[O:5])[CH3:16])(=[C:8]([CH3:10])[CH2:13][CH2:12][CH3:11])[CH2:14][CH3:15].[OH2:9]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)

  def testCase5(self,display=True):
      # Check the aldol addition => this one is tricky because we need to fix hydrogen atoms
      sma = "CC[CH2:3][C:4](C)=[O:5].CCC[C:8]([CH3:10])=[O:9]>>CCC[C:8]([CH3:10])([OH:9])[CH:3](CC)[C:4](C)=[O:5]"
      expt = "[CH2:3]([C:4](=[O:5])[CH3:16])[CH2:14][CH3:15].[C:8](=[O:9])([CH3:10])[CH2:13][CH2:12][CH3:11]>>[CH:3]([C:4](=[O:5])[CH3:16])([C:8]([OH:9])([CH3:10])[CH2:13][CH2:12][CH3:11])[CH2:14][CH3:15]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase6(self,display=True):
      # Check the Diels-Alder with symmetry on the dienophile
      sma = "C[CH:1]=[CH:2][CH:3]=[CH:4]CC.COC(=O)[CH:5]=[CH:6]C(=O)OC>>C[CH:1]1[CH:2]=[CH:3][CH:4](CC)[CH:5](C(OC)=O)[CH:6]1C(=O)OC"
      expt = "[CH:1](=[CH:2][CH:3]=[CH:4][CH2:8][CH3:9])[CH3:7].[CH:5](=[CH:6][C:14](=[O:15])[O:16][CH3:17])[C:10]([O:11][CH3:12])=[O:13]>>[CH:1]1([CH3:7])[CH:2]=[CH:3][CH:4]([CH2:8][CH3:9])[CH:5]([C:10]([O:11][CH3:12])=[O:13])[CH:6]1[C:14](=[O:15])[O:16][CH3:17]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)

  def testCase7(self,display=True):
      # Check Friedel-Crafts alkylation
      sma = "c1cc[cH:3]cc1.C[CH:2](C)[Cl:1]>>c1cc[c:3]([CH:2](C)C)cc1.[Cl-:1]"
      expt = "[cH:3]1[cH:6][cH:5][cH:4][cH:10][cH:9]1.[Cl:1][CH:2]([CH3:7])[CH3:8]>>[CH:2]([c:3]1[cH:6][cH:5][cH:4][cH:10][cH:9]1)([CH3:7])[CH3:8].[Cl-:1]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase8(self,display=True):
      # Check Reductive amination
      sma = "[NH2:2][CH2:1]C1=CC=CC=C1.C[CH:5]1CCC[CH:6](C)[C:3]1=[O:4]>>C[CH:5]1CCC[CH:6](C)[CH:3]1[NH:2][CH2:1]C1=CC=CC=C1"
      expt = "[CH2:1]([NH2:2])[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1.[C:3]1(=[O:4])[CH:5]([CH3:7])[CH2:8][CH2:9][CH2:10][CH:6]1[CH3:11]>>[CH2:1]([NH:2][CH:3]1[CH:5]([CH3:7])[CH2:8][CH2:9][CH2:10][CH:6]1[CH3:11])[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1.[OH2:4]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase9(self,display=True):
      # Williamson ether
      sma = "CC[CH2:4][OH:3].CC[CH2:1][Br:2]>>CC[CH2:1][O:3][CH2:4]CC"
      expt = "[OH:3][CH2:4][CH2:7][CH3:8].[CH2:1]([Br:2])[CH2:6][CH3:5]>>[CH2:1]([O:3][CH2:4][CH2:7][CH3:8])[CH2:6][CH3:5].[BrH:2]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase10(self,display=True):
      # Gabriel synthesis
      # alternate string output order... a potential canonical issue ?
      sma = "[NH2:4][NH2:1].[Cl:3][CH2:2]C1=CC=CC=C1>>[NH3:4].[NH2:1][CH2:2]C1=CC=CC=C1"
      expt = "[CH2:2]([Cl:3])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1.[NH2:1][NH2:4]>>[NH2:1][CH2:2][c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1.[NH3:4].[ClH:3]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase11(self,display=True):
      # 3,3-sigmatropic rearrangement
      sma = "CC(C)[CH:4]=[CH:5][CH2:6][CH2:3][CH:2]=[CH:1]C(C)C>>CC(C)[CH:1]([CH:2]=[CH2:3])[CH:4](C(C)C)[CH:5]=[CH2:6]"
      expt = "[CH:1](=[CH:2][CH2:3][CH2:6][CH:5]=[CH:4][CH:10]([CH3:11])[CH3:12])[CH:8]([CH3:7])[CH3:9]>>[CH:1]([CH:2]=[CH2:3])([CH:4]([CH:5]=[CH2:6])[CH:10]([CH3:11])[CH3:12])[CH:8]([CH3:7])[CH3:9]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase12(self,display=True):
      # Check the Claisen: Forward/reverse case
      sma = "CC[CH2:6][C:5]([O:4][CH2:3][CH:2]=[CH2:1])=[CH:7]C(C)C>>CC[CH2:6][C:5](=[O:4])[CH:7]([CH2:3][CH:2]=[CH2:1])C(C)C"
      expt = "[CH2:1]=[CH:2][CH2:3][O:4][C:5]([CH2:6][CH2:9][CH3:8])=[CH:7][CH:10]([CH3:11])[CH3:12]>>[CH2:1]=[CH:2][CH2:3][CH:7]([C:5](=[O:4])[CH2:6][CH2:9][CH3:8])[CH:10]([CH3:11])[CH3:12]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase13(self,display=True):
      # Benzene hydrogenation
      sma = "C[c:4]1[cH:3][c:2](C)[cH:1][c:6](C)[cH:5]1>>C[CH:6]1[CH2:5][CH:4](C)[CH2:3][CH:2](C)[CH2:1]1"
      expt = "[cH:1]1[c:2]([CH3:9])[cH:3][c:4]([CH3:8])[cH:5][c:6]1[CH3:7]>>[CH2:1]1[CH:2]([CH3:9])[CH2:3][CH:4]([CH3:8])[CH2:5][CH:6]1[CH3:7]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase14(self,display=True):
      # Strecker step 1
      sma = "[NH3:1].[CH:4]#[N:5].CC[CH:2]=[O:3]>>CC[CH:2]([NH2:1])[C:4]#[N:5]"
      expt = "[CH:2](=[O:3])[CH2:7][CH3:6].[NH3:1].[CH:4]#[N:5]>>[NH2:1][CH:2]([C:4]#[N:5])[CH2:7][CH3:6].[OH2:3]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
      
  def testCase15(self,display=True):
      # Strecker step 2
      sma = "CC[CH:1](N)[C:3]#[N:5]>>CC[CH:1](N)[C:3]([OH:4])=[O:2]"
      expt = "[CH:1]([C:3]#[N:5])([CH2:7][CH3:6])[NH2:8].[OH2:4].[OH2:2]>>[CH:1]([C:3](=[O:2])[OH:4])([CH2:7][CH3:6])[NH2:8].[NH3:5]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase16(self,display=True):
      # Strecker step 1+2
      sma = "[NH3:1].[CH:4]#[N:5].C1CCCCC1C[CH:2]=[O:3]>>C1CCCCC1C[CH:2]([NH2:1])[C:4]([OH:7])=[O:6]"
      expt = "[CH:2](=[O:3])[CH2:14][CH:13]1[CH2:8][CH2:9][CH2:10][CH2:11][CH2:12]1.[OH2:7].[OH2:6].[CH:4]#[N:5].[NH3:1]>>[NH2:1][CH:2]([C:4](=[O:6])[OH:7])[CH2:14][CH:13]1[CH2:8][CH2:9][CH2:10][CH2:11][CH2:12]1.[NH3:5].[OH2:3]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase17(self,display=True):
      # Complex Reaction 1 without acid reaction mapped (major product only)
      sma = "CC(C)[OH:4].CC(C)[C:1](=O)[O:2][C:3](=O)C(C)C>>CC(C)[O:4][C:1](=O)C(C)C"
      expt = "[C:1]([O:2][C:3](=[O:12])[CH:13]([CH3:14])[CH3:15])(=[O:8])[CH:9]([CH3:10])[CH3:11].[OH:4][CH:6]([CH3:5])[CH3:7]>>[C:1]([O:4][CH:6]([CH3:5])[CH3:7])(=[O:8])[CH:9]([CH3:10])[CH3:11].[CH2:13]([CH3:14])[CH3:15].[CH4:3].[OH2:2].[OH2:12]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase18(self,display=True):
      # Baeyer-Villiger (major)
      sma = "[CH3:3][C:2](=[O:1])C1=CC=CC=C1>>[CH3:3][O:4][C:2](=[O:1])C1=CC=CC=C1"
      expt ="[O:1]=[C:2]([CH3:3])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1.[OH2:4]>>[O:1]=[C:2]([O:4][CH3:3])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=display), expt)
      
  def testCase19(self,display=True):
      # Baeyer-Villiger (full)
      sma = "C[C:6](=O)[O:5][OH:4].[CH3:3][C:2](=[O:1])C1=CC=CC=C1>>C[C:6]([OH:5])=O.[CH3:3][O:4][C:2](=[O:1])C1=CC=CC=C1"
      expt = "[O:1]=[C:2]([CH3:3])[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1.[OH:4][O:5][C:6]([CH3:7])=[O:8]>>[O:1]=[C:2]([O:4][CH3:3])[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1.[OH:5][C:6]([CH3:7])=[O:8]"
      self.assertEqual( CGR.RxnCompleteMapping(sma, debug=True), expt)

if __name__ == '__main__':
  unittest.main()
