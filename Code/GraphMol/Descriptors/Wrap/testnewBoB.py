from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
import numpy as np
import time
import scipy.io as sio

# convert smile to mol
def mols(smile):
	m = Chem.MolFromSmiles(smile) #,sanitize=False)
	#Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
	ps = AllChem.ETKDG()
	ps.randomSeed = 0xf00d
	AllChem.EmbedMolecule(m,ps)
	AllChem.MMFFOptimizeMolecule(m)
	return m

# generate from rdkit coulomb matrix from smiles
def generateX(m, bagdict):
	vs = rdMD.CalcBagOfBondVector(m, -1, 1, bagdict) # true, false
	return vs

# read the file and get the data
def getX(bagdict):
    X=[]
    for smile in smi:
    	m = mols(smile)
    	x= generateX(m,bagdict)
    	X.append(x)
    return X

# read the file and get the data
def loadsmi(f):
    smi=[]
    fo = open(f)
    i=0
    for line in fo:
        i+=1
        if i>10000:
            break
        if i>1:
            rl=line.split(',')
            smi.append(rl[0])
    return smi


def bagOfbondDict(smiles):
    bobDict = rdMD.CalcBagOfBondsMap(smiles)
    return bobDict

if __name__ == "__main__":
    f='../test_data/qm7smiles.csv'
    startTime = time.time()
    smi=loadsmi(f)
    print len(smi)
    bagdict = bagOfbondDict(smi)
    print bagdict
    print ('The smi concate script took {0} second !'.format(time.time() - startTime))
    startTime = time.time()
    #X,Y = loaddata(f,smi)
    X= getX(bagdict)
    print len(X)
    print ('The generation of BofB script took {0} second !'.format(time.time() - startTime))
    print len(X)
    print "---------"
    Mat = np.concatenate(X,axis=0)
    sio.savemat('BoB.mat', {'X':Mat})

