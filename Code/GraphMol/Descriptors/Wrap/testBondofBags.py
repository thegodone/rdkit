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
def generateX(m, smi):
	vs = rdMD.CalcCoulombMat(m, bob = True,  alpha = 1, smiles = smi) # true, false
	for v in vs:
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
	return arr

# read the file and get the data
def getX(smi):
    X=[]
    for smile in smi:
    	    m = mols(smile)
    	    x= generateX(m,smi)
    	    X.append(x)
    return X

# read the file and get the data
def loadsmi(f):
    smi=[]
    fo = open(f)
    i=0
    for line in fo:
        i+=1
        if i>200:
            break
        if i>1:
            rl=line.split(',')
            smi.append(rl[0])
        if i%100==0:
            print i

    return smi

if __name__ == "__main__":
    f='/Users/GVALMTGG/Downloads/qm7smiles.csv'
    startTime = time.time()
    smi=loadsmi(f)
    print smi
    print ('The smi concate script took {0} second !'.format(time.time() - startTime))
    startTime = time.time()
    #X,Y = loaddata(f,smi)
    X= getX(smi)
    for x in X:
        print x
    print ('The generation of BofB script took {0} second !'.format(time.time() - startTime))
    print len(X)
    print "---------"
    #Mat = np.concatenate(S,axis=0)
    #sio.savemat('Xoptsortedsautoanitized.mat', {'X':Mat,'Y':Y})

