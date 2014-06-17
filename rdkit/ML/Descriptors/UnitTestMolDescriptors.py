#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#

""" unit testing code for molecular descriptor calculators

"""
import unittest,os.path
from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy
from rdkit import Chem

class TestCase(unittest.TestCase):
  def setUp(self):
    self.descs = ['MolLogP','Chi1v']
    self.vers= ('1.1.0','1.0.0')
    self.calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.descs)
    self.testD = [
      ('CCOC',    (0.6527, 1.40403)),
      ('CC=O',    (0.2052, 0.81305)),
      ('CCC(=O)O',(0.481, 1.48839))]

  def testGetNames(self):
    self.failUnlessEqual(self.calc.GetDescriptorNames(),tuple(self.descs))
    
  def _testVals(self,calc,testD):
    for smi,vals in testD:
      mol = Chem.MolFromSmiles(smi)
      ans = numpy.array(vals)
      res = numpy.array(calc.CalcDescriptors(mol))
      self.failUnless(max(abs(res-ans))<1e-4,'bad descriptor values for SMILES %s (%s)'%(smi,str(res)))
    
  def testCalcVals(self):
    self._testVals(self.calc,self.testD)

  def testSaveState(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML/Descriptors/test_data','molcalc.dsc')
    ok = 1
    try:
      inF = open(fName,'rb')
    except:
      ok = 0
    assert ok,'problems opening saved file %s'%(fName)
    try:
      calc = cPickle.load(inF)
    except:
      ok = 0
    assert ok,'problems reading saved file %s'%(fName)
      

    self.failUnlessEqual(calc.GetDescriptorNames(),tuple(self.descs))
    self.failUnlessEqual(calc.GetDescriptorVersions(),tuple(self.vers))
    self._testVals(calc,self.testD)
    
    
if __name__ == '__main__':
  unittest.main()
