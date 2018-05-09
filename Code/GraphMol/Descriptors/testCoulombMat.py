from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('NCCCCCO')
m = Chem.AddHs(m)
ps = AllChem.ETKDG()
ps.randomSeed = 0xf00d
AllChem.EmbedMolecule(m,ps)
vs = rdMD.CalcEEMcharges(m)
print vs

vs = rdMD.CalcCoulombMat(m,-1,5,1)
print vs
