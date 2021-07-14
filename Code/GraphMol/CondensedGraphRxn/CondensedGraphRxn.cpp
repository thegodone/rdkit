#include "CondensedGraphRxn.h"
#include <boost/shared_ptr.hpp>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <unordered_set>
#include <RDGeneral/utils.h>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <vector>

/////////////////// CGR format validation ///////////////////////
extern std::string validcgrs[20] = {
    "{!-}","{!=}","{!#}","{!:}",
    "{-!}","{-=}","{-#}","{-:}",
    "{=!}","{=-}","{=#}","{=:}",
    "{#!}","{#-}","{#=}","{#:}",
    "{:!}","{:-}","{:=}","{:#}"
};

// Define a method to check if known
bool IsKnown(const std::string& bond) {
    for (int i =0; i<20;i++) {
        if (bond == validcgrs[i]) return true;
    }
    return false;
}

// Define a method to check for valid CGR
bool IsValidCGR(const char* CGR) {
    // Define the test criteria
    int open = 0;  // Check close vs open
    int close = 0; // Check close vs open
    int start = 0; // We need to cache the starting position

    // Run the tests
    std::string cgr = std::string(CGR);
    int i = 0;
    while (CGR[i]!='\0') {
        // Check we start a bond
        if (CGR[i]=='{') {
            if (++open > 1) return false;
            start = i;
        }

        // Check we end a bond
        if (CGR[i]=='}') {
           // Here we need to test => can never be bigger than open
        if (++close > open) return false;
        // Here we do test to make the spacing is correct
        if (i-start != 3) return false;
        // Now we test the characters explicitly
        if (!IsKnown(cgr.substr(start,4))) return false;
        // All done, reset open/close to 0
        open = 0;
        close = 0;
        }
        i++;
    }
    // Return both values 0
    return (open == 0 and close == 0);
}

// Just a duplicate format for std::string to use the above
bool IsValidCGRString(const std::string& cgr) {
    return IsValidCGR(cgr.c_str());
}

namespace RDKit {
namespace CondensedGraphRxn {


class Row {
private:
    int* ptr;
    int l;
public:
    Row(const int& l) {
        this->ptr = (int*) std::malloc(l*sizeof(int));
        for (int i=0;i<l;i++)
            this->ptr[i] = 0;
        this->l = l;
    }
    virtual ~Row() {
        free(this->ptr);
    }
    void Set(const int idx,int value) {
        this->ptr[idx] = value;
    }
    const int Get(const int idx) {
        return this->ptr[idx];
    }
    void Display() {
        std::cout << "";
        for (int i=0;i<this->l;i++) std::cout << " " << this->ptr[i];
        std::cout << std::endl;
    }
};

class AdjMatrix {
private:
    std::vector<int*> M;
    int row,col;
public:
    AdjMatrix(const int& row, const int& col) {
        this->row = row;
        this->col = col;
        std::vector<int*> M(row);
        for (int i=0;i<row;i++) {
            int* column = (int*) std::malloc((row)*sizeof(int));
            for (int c = 0; c < row; c++)
                column[c]=0;
            M[i] = column;
        }
        this->M = M;
    }
    virtual ~AdjMatrix() {
        for (int i=0;i<this->row;i++)
            free(this->M[i]);
    }
    void SetBO(const int& from,const int& to, const int& bo) {
        int qrow = std::min(from,to);
        int qcol = std::max(from,to);
        M[qrow][qcol] = bo;
    }
    const int GetBO(const int& from,const int& to) {
        int qrow = std::min(from,to);
        int qcol = std::max(from,to);
        return M[qrow][qcol];
    }
    void Display() {
        for (int i=0;i<row;i++) {
            std::cout << "";
            for (int j=0;j<col;j++) {
                std::cout << " " << GetBO(i,j);
            }
            std::cout << std::endl;
        }
    }
};

///////////////// bond CGR type //////////////
// This may become something like
const Bond::BondType getCGRBondType(const int bor, const int bop) {
        if (bor == bop){
           switch(bor) {
              case 1:
                   return Bond::BondType::SINGLE;
              case 2:
                    return Bond::BondType::DOUBLE;
              case 3:
                    return Bond::BondType::TRIPLE;
              case 4:
                    return Bond::BondType::AROMATIC;
              default:
                    return Bond::BondType::ZERO;
           }
        }
        else {
            switch(bor) {
                case 0:
                      if (bop == 1) return  Bond::BondType::CGRNS;
                      else if (bop == 2) return  Bond::BondType::CGRND;
                      else if (bop == 3) return  Bond::BondType::CGRNT;
                      else if (bop == 4) return  Bond::BondType::CGRNA;
                case 1:
                      if (bop == 0) return  Bond::BondType::CGRSN;
                      else if (bop == 2) return  Bond::BondType::CGRSD;
                      else if (bop == 3) return  Bond::BondType::CGRST;
                      else if (bop == 4) return  Bond::BondType::CGRSA;
                case 2:
                      if (bop == 0) return  Bond::BondType::CGRDN;
                      else if (bop == 1) return  Bond::BondType::CGRDS;
                      else if (bop == 3) return  Bond::BondType::CGRDT;
                      else if (bop == 4) return  Bond::BondType::CGRDA;
                case 3:
                      if (bop == 0) return  Bond::BondType::CGRTN;
                      else if (bop == 1) return  Bond::BondType::CGRTS;
                      else if (bop == 2) return  Bond::BondType::CGRTD;
                      else if (bop == 4) return  Bond::BondType::CGRTA;
                case 4:
                      if (bop == 0) return  Bond::BondType::CGRAN;
                      else if (bop == 1) return  Bond::BondType::CGRAS;
                      else if (bop == 2) return  Bond::BondType::CGRAD;
                      else if (bop == 3) return  Bond::BondType::CGRAT;
                default:
                      return Bond::BondType::ZERO;
             }
       }
}

////////// get involved bonds
AdjMatrix* getInvolvedBonds (const int maxmapnum, const RDKit::ChemicalReaction &rxn, bool product=false) {
   AdjMatrix* adj = new AdjMatrix(maxmapnum,maxmapnum); 
   
   MOL_SPTR_VECT v = product ? rxn.getProducts() : rxn.getReactants() ;
   
   for (auto rIt = v.begin(); rIt != v.end(); ++rIt) {
      for(auto bt : (*rIt)->bonds()){
        int fromnum = bt->getBeginAtom()->getAtomMapNum();
        int tonum = bt->getEndAtom()->getAtomMapNum();
        Bond::BondType tp = bt->getBondType() ;
        
        // caution aromaticity wins on single or double bond flags
        bool fromarom = bt->getBeginAtom()->getIsAromatic();
        bool toarom = bt->getEndAtom()->getIsAromatic();
          
        if (tp == Bond::SINGLE) {
            if (fromarom && toarom) {
                adj->SetBO(fromnum, tonum, 4);
            }
            else {
                  adj->SetBO(fromnum, tonum, 1);
            }
          }
          else if  (tp == Bond::DOUBLE) {
            if (fromarom && toarom) {
                adj->SetBO(fromnum, tonum, 4);
            }
            else {
                adj->SetBO(fromnum, tonum, 2);
            }
        }
          else if (tp == Bond::TRIPLE) {
            adj->SetBO(fromnum, tonum, 3);
          }
          else if (tp == Bond::AROMATIC) {
            adj->SetBO(fromnum, tonum, 4);
        }
      }   
  }
  return adj;    
}


    
bool isCGRbond(Bond *bond){
      Bond::BondType BT = bond->getBondType();

      switch (BT){
            case Bond::CGRSN :
            case Bond::CGRSD :
            case Bond::CGRST :
            case Bond::CGRSA :
            case Bond::CGRDN :
            case Bond::CGRDS :
            case Bond::CGRDA :
            case Bond::CGRDT :
            case Bond::CGRTN :
            case Bond::CGRTS :
            case Bond::CGRTA :
            case Bond::CGRTD :
            case Bond::CGRAN :
            case Bond::CGRAS :
            case Bond::CGRAD :
            case Bond::CGRAT :
            case Bond::CGRNS :
            case Bond::CGRND :
            case Bond::CGRNA :
            case Bond::CGRNT :
               return true;
            default :
               return false;
       }
}

    

    
    
////////////////////// CGR parts ////////////////////////////
int getMaxMapNum(ChemicalReaction &rxn) {
    // Get the highest number
    int maxmapnum = -1;
    for (auto rIt = rxn.beginReactantTemplates(); rIt != rxn.endReactantTemplates(); ++rIt) {
        for (auto at : (*rIt)->atoms()) 
            maxmapnum = std::max<int>(maxmapnum,at->getAtomMapNum());
    }
    for (auto rIt = rxn.beginProductTemplates(); rIt != rxn.endProductTemplates(); ++rIt) {
        for (auto at : (*rIt)->atoms()) 
            maxmapnum = std::max<int>(maxmapnum,at->getAtomMapNum());
    }
    return maxmapnum;
}
    
std::shared_ptr<RWMol> getCGRmol(ChemicalReaction &rxn, bool charges) {
    bool aromaticity = false;
    const int maxmapnum = getMaxMapNum(rxn)+1;
    Row* atno = new Row(maxmapnum);
    Row* newindex = new Row(maxmapnum);
    Row* aromatic = new Row(maxmapnum);
    Row* charge = new Row(maxmapnum);
    Row* isotope = new Row(maxmapnum);
    
    std::vector<int> nums; 
    int atindex = 0;
    
    
    for (auto rIt = rxn.beginReactantTemplates(); rIt != rxn.endReactantTemplates();
         ++rIt) {   
      for(auto at : (*rIt)->atoms()){
         if (at->getAtomMapNum()>0) {
            int mnum = at->getAtomMapNum();
             /*
             if (!radical) {
                 if (at->getNumRadicalElectrons()) {
                     at->setNumRadicalElectrons(0);
                     at->setNoImplicit(false);
                 }
             }
            */
              
            atno->Set(mnum, at->getAtomicNum());
            aromatic->Set(mnum, at->getIsAromatic());
            newindex->Set(mnum, atindex++);
            if (charges) {
                charge->Set(mnum, at->getFormalCharge());
            }
            isotope->Set(mnum, at->getIsotope());
            nums.push_back(mnum);
         }
      }
    }

   
    /// only keep aromaticity if conserved in product ...
    for (auto rIt = rxn.beginProductTemplates(); rIt != rxn.endProductTemplates();
          ++rIt) {
      for(auto at : (*rIt)->atoms()){
         if (at->getAtomMapNum()>0) {
           int mnum = at->getAtomMapNum();

           if (aromaticity) {
               aromatic->Set(mnum , aromatic->Get(mnum));
           }
           else {
               int getaroma = at->getIsAromatic() ? 1 : 0;
               aromatic->Set(mnum , std::min(getaroma, aromatic->Get(mnum)));
           }
        }
      }
    }

   AdjMatrix* RBond = getInvolvedBonds ( maxmapnum, rxn, false);
   AdjMatrix* PBond = getInvolvedBonds ( maxmapnum, rxn, true);

   ///////////////////////////////////////////////////////////////////////
   // Constructor CGR "RDKit molecule" only works for fully map objects //
   ///////////////////////////////////////////////////////////////////////
    
   std::shared_ptr<RWMol> mol(new RWMol()); 

   for (int mnum : nums){
       RDKit::Atom ra(atno->Get(mnum));
       ra.setIsAromatic(aromatic->Get(mnum));
       ra.setIdx(newindex->Get(mnum));
       if (charges) {
           ra.setFormalCharge(charge->Get(mnum));
       }
       
       /*
       if (!radical) {
           if (ra.getNumRadicalElectrons()) {
               ra.setNumRadicalElectrons(0);
               ra.setNoImplicit(false);
           }
       }
       */
       ra.setIsotope(isotope->Get(mnum));
       mol->addAtom(&ra);
   }
    
   for (int from : nums){
       for (int to : nums){
          if (to > from){
             int boR = RBond->GetBO(from, to);
             int boP = PBond->GetBO(from, to);
             Bond::BondType bot = getCGRBondType(boR, boP);
             if (bot == Bond::BondType::ZERO) {
                continue;
             }
             else {    
                 int fromindex = newindex->Get(from);
                 int toindex = newindex->Get(to);
                 if (!mol->getBondBetweenAtoms(fromindex,toindex)){
                      mol->addBond(fromindex, toindex, bot);
                 }
             }
          } 
        }
    }
    
    delete RBond;
    delete PBond;
    delete atno;
    delete newindex;
    delete aromatic;
    delete charge;
    delete isotope;
    return mol;
}

std::vector<unsigned int> getNeighborsPtr(std::shared_ptr<RWMol> &mol, Atom* atom) {
   std::vector<unsigned int> res;
   for(const auto &nbri: make_iterator_range(mol->getAtomBonds(atom))) {
       const RDKit::Bond *bond = (*mol)[nbri];
      res.push_back(bond->getOtherAtomIdx(atom->getIdx()));
   }
   return res;
}

void BFS(std::shared_ptr<RWMol> mol, std::vector<unsigned int> atidx, bool* important, unsigned int radius, int r=0) {
   // Check if the depth has been reached or the list is empty
   // Stop if this is the case => this will never fail, because
   // eventually we will have all atoms as important.
   if (r==radius || atidx.size()==0) 
       return;        

   // Still good, let's take the next shell of neighbors on a list
   std::vector<unsigned int> newatidx;
   for (auto a: atidx) {
       Atom* atom  = mol->getAtomWithIdx(a);
       std::vector<unsigned int> nei = getNeighborsPtr(mol,atom);
       for (unsigned int ni: nei){
             // if not identified as important add it
           if (!important[ni]) {
               newatidx.push_back(ni);
               important[ni] = 1;
           }
        }
    }
    
    // Repeat for the next shell with the new set of atoms
    return BFS(mol, newatidx, important, radius, r+=1);
}


void addAtomRingCGRIdx(std::shared_ptr<RWMol> mol, bool* important) {
    
    if( !mol->getRingInfo()->isInitialized() ) {
        RDKit::MolOps::findSSSR( *mol );
    }
    
   
    VECT_INT_VECT bondRings =  mol->getRingInfo()->bondRings();
        
    for (VECT_INT_VECT_CI ringIt = bondRings.begin(); ringIt != bondRings.end();
       ++ringIt) {
          bool ringChange = false;
          for (INT_VECT_CI bondIt = ringIt->begin(); bondIt != ringIt->end();
             ++bondIt) {
              if (isCGRbond(mol->getBondWithIdx(*bondIt))){
                  ringChange = true; 
                 // break;
              }
          }

          if (ringChange){
             for (INT_VECT_CI bondIt2 = ringIt->begin(); bondIt2 != ringIt->end();
             ++bondIt2) {
                  important[mol->getBondWithIdx(*bondIt2)->getBeginAtomIdx()]=1;
                  important[mol->getBondWithIdx(*bondIt2)->getEndAtomIdx()]=1;
              }
          }
     }
}

    
    
std::shared_ptr<RWMol> getCGRsignature(ChemicalReaction &rxn, unsigned int radius = 1,
                                       bool charges = false) {

    // Construct the molecule
    std::shared_ptr<RWMol>  mol = getCGRmol(rxn, charges);
    

    // Define an array to figure out what's important
    // This array can also be used to guide the search
    // in the shells for neighbors.
    int n = mol->getNumAtoms();
    bool *important = (bool*) std::malloc(n*sizeof(bool));
    for (int i = 0; i < n; i++) {
       important[i]= 0;
    }
    
    // get the reaction atom centers as first input
    // for the search of the relevant atoms
    std::vector<unsigned int> atidx;
    for (auto b: mol->bonds()) {
          if (isCGRbond(b)){
              
          // The bond is a cgr-bond, define the atoms
          // as the first atoms for the neighbor search.
             int from = b->getBeginAtomIdx();
             int to = b->getEndAtomIdx();
             if (!important[from]) {
                atidx.push_back(from);
                important[from] = 1;
             }
             if (!important[to]) {
                atidx.push_back(to);
                important[to] = 1;
             }
          }
        
        

     }

    // Apply the BFS to get the neighbors of the reaction
    // center as specified by the desired radius.
    // This method is recursive.
    BFS(mol,atidx,important,radius);
    
    // inject the ring change two!
    addAtomRingCGRIdx(mol,important);

    // remove the atoms not identified as 'important'
    for(int i = n-1; i >=0; i--){
        if (important[i] == 0) {
            // we can check if the bond was in ring formation destruction here to avoid to delete it from the important list
              mol->removeAtom(i);
        }
    }

    // release memory for 'important' - this was tmp control array
    free(important);
    return mol;
}

std::string getCGRwriter(ChemicalReaction &rxn, bool doRandom, unsigned int randomSeed, bool aromatize = true, bool signature = false, bool charges = false, int radius=1) {

    // Define the correct molecule for the output.
    // Here we generate the full CGR molecule or the signature CGR molecule if asked for.
    // The signature CGR defines a smaller molecule around the center of reaction only
    // as specified by the radius.
    std::shared_ptr<RWMol>  mol = signature ? getCGRsignature(rxn, radius, charges) : getCGRmol(rxn, charges);
   
    //// do random need a seed
    if (randomSeed > 0) {
      getRandomGenerator(rdcast<int>(randomSeed));
    }
    
    // set moltosmiles parameters
    bool canonical = !doRandom; 
    int rootedAtAtom = -1;
    bool doKekule = false;
    bool doIsomericSmiles = true;
    bool allBondsExplicit =  false;
    bool allHsExplicit = false;

    if (aromatize) {
         unsigned int failed;
         try {
               /// this method try to sanitize the object molecule
               // we deactivate some methods that can failed with CGR bonds
                  unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^
                          MolOps::SANITIZE_CLEANUP ^
                          MolOps::SANITIZE_PROPERTIES ^
                          MolOps::SANITIZE_KEKULIZE;

                  MolOps::sanitizeMol(*mol, failed, sanitizeOps);
         }
         catch (MolSanitizeException &) {
                  BOOST_LOG(rdInfoLog) << "Issue to sanitize Mol\n";
         }
    }

    // Split the frags and write to list
    std::vector<ROMOL_SPTR> listofmol = MolOps::getMolFrags(*mol, false);
    std::stringstream ss;
    int f = 0;
    for (ROMOL_SPTR fragmol : listofmol) {
         ss << (f++==0 ? "" : "."); 
           ss << RDKit::MolToSmiles(*fragmol, doIsomericSmiles, doKekule, rootedAtAtom, canonical, allBondsExplicit, allHsExplicit, doRandom);
    }    
    // Done: to string now!
    return ss.str();
}

/////////////// cgr reader ///////////////////// 
std::string  getCGRMolecule(RWMol *molR,const std::string cgr, bool canonical, bool setAtomMap) {
        unsigned int failedR, failedP;
    // reset AtomMapNum
    if  (setAtomMap) {
           for(auto at : molR->atoms()){
              at->setAtomMapNum(at->getIdx()+1);
           }
    }
    
    /*if (!radical) {
     for(auto at : molR->atoms()){
     if (at->getNumRadicalElectrons()) {
     at->setNumRadicalElectrons(0);
     at->setNoImplicit(false);
     }
     }
     */

  std::shared_ptr<RWMol> molP(new RWMol(*molR));
    try {
       // this method try to sanitize the object molecule
       // we deactivate some methods that can failed with CGR bonds
              unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^
                    MolOps::SANITIZE_CLEANUP ^
                    MolOps::SANITIZE_PROPERTIES ^
                    MolOps::SANITIZE_KEKULIZE;
       
              MolOps::sanitizeMol(*molR, failedR, sanitizeOps);
              MolOps::sanitizeMol(*molP, failedP, sanitizeOps);
    }
    catch (MolSanitizeException &) {
              BOOST_LOG(rdInfoLog) << "Issue to sanitize Mol from this cgr:"+cgr+"\n";
    }
    
         std::set<unsigned int> keep;         
         std::vector<std::pair<int, int> > bondRrem;
         std::vector<std::pair<int, int> > bondPrem; 
         for( unsigned int i = 0 , is = molR->getNumBonds() ; i < is ; ++i ) {
              Bond *bond = molR->getBondWithIdx(i);
              Bond *bondP = molP->getBondWithIdx(i);
              Bond::BondType BT = bond->getBondType();
              switch (BT){
                  case Bond::CGRSN :  
                    bond->setBondType(Bond::SINGLE);
                    bondPrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRSD :
                    bondP->setBondType(Bond::DOUBLE);
                    bond->setBondType(Bond::SINGLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRSA :
                    bond->setBondType(Bond::SINGLE);
                    bondP->setBondType(Bond::AROMATIC);
                    bondP->getBeginAtom()->setIsAromatic(true);
                    bondP->getEndAtom()->setIsAromatic(true);
                    bond->getBeginAtom()->setIsAromatic(false);
                    bond->getEndAtom()->setIsAromatic(false);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRST :
                    bondP->setBondType(Bond::TRIPLE);
                    bond->setBondType(Bond::SINGLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRDN :
                    bond->setBondType(Bond::DOUBLE);      
                    bondPrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRDS :
                    bond->setBondType(Bond::DOUBLE);
                    bondP->setBondType(Bond::SINGLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRDA :
                    bond->setBondType(Bond::DOUBLE);
                    bondP->setBondType(Bond::AROMATIC);
                    bondP->getBeginAtom()->setIsAromatic(true);
                    bondP->getEndAtom()->setIsAromatic(true);
                    bond->getBeginAtom()->setIsAromatic(false);
                    bond->getEndAtom()->setIsAromatic(false);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRDT :
                    bondP->setBondType(Bond::TRIPLE);
                    bond->setBondType(Bond::DOUBLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRTN :
                    bond->setBondType(Bond::TRIPLE);
                    bondPrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRTS :
                    bond->setBondType(Bond::TRIPLE);
                    bondP->setBondType(Bond::SINGLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRTA :
                    bond->setBondType(Bond::TRIPLE);
                    bondP->setBondType(Bond::AROMATIC);
                    bondP->getBeginAtom()->setIsAromatic(true);
                    bondP->getEndAtom()->setIsAromatic(true);
                    bond->getBeginAtom()->setIsAromatic(false);
                    bond->getEndAtom()->setIsAromatic(false);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                    case Bond::CGRTD :
                    bond->setBondType(Bond::TRIPLE);
                    bondP->setBondType(Bond::DOUBLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;           
                  case Bond::CGRAN :
                    bond->setBondType(Bond::AROMATIC);
                    bondPrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    bondP->getBeginAtom()->setIsAromatic(false);
                    bondP->getEndAtom()->setIsAromatic(false);
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRAS :
                    bond->setBondType(Bond::AROMATIC);
                    bondP->setBondType(Bond::SINGLE);
                    bondP->getBeginAtom()->setIsAromatic(false);
                    bondP->getEndAtom()->setIsAromatic(false);
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRAD :
                    bond->setBondType(Bond::AROMATIC);
                    bondP->setBondType(Bond::DOUBLE);
                    bondP->getBeginAtom()->setIsAromatic(false);
                    bondP->getEndAtom()->setIsAromatic(false); 
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRAT :
                    bondP->setBondType(Bond::TRIPLE);
                    bond->setBondType(Bond::AROMATIC);
                    bondP->getBeginAtom()->setIsAromatic(false);
                    bondP->getEndAtom()->setIsAromatic(false);
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRNS :
                    bondRrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    bondP->setBondType(Bond::SINGLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRND :
                    bondP->setBondType(Bond::DOUBLE);
                    bondRrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRNT :
                    bondRrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    bondP->setBondType(Bond::TRIPLE);
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  case Bond::CGRNA :
                    bondP->setBondType(Bond::AROMATIC);
                    bondP->getBeginAtom()->setIsAromatic(true);
                    bondP->getEndAtom()->setIsAromatic(true);
                    bond->getBeginAtom()->setIsAromatic(false);
                    bond->getEndAtom()->setIsAromatic(false);
                    bondRrem.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                    keep.insert(bond->getBeginAtomIdx());
                    keep.insert(bond->getEndAtomIdx());
                    break;
                  default :
                    break;    
                }       
          }

    
        
        for( unsigned int i = 0 , is = molR->getNumBonds() ; i < is ; ++i ) {
              Bond *bond = molR->getBondWithIdx(i);
              
             if (bond->getBondType() == Bond::AROMATIC) {
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
             }  
        }
    
        for( unsigned int i = 0 , is = molP->getNumBonds() ; i < is ; ++i ) {
              Bond *bond = molP->getBondWithIdx(i);
              
             if (bond->getBondType() == Bond::AROMATIC) {
                    bond->getBeginAtom()->setIsAromatic(true);
                    bond->getEndAtom()->setIsAromatic(true);
             }  
        }
    
        for (auto br = std::begin(bondRrem); br != std::end(bondRrem); ++br) {
            // remove bonds broken here
            molR->removeBond(br->first, br->second);
        }
         
        for (auto br = std::begin(bondPrem); br != std::end(bondPrem); ++br) {
            // remove bonds broken here
            molP->removeBond(br->first, br->second);
        }
    

    
    
       // fixing the charges "N+,O-" => see the available scripts in rdMolStandardize to do this!
       // for the radical we need to think that we can remove initially from the reaction SMA
       
    
    
    
       // std::cout << "just before cgr from smile generator\n";
    // can remove radical at the end or at the beginning ...
    return  RDKit::MolToSmiles( *molR, true, false, -1, canonical )  + ">>" + RDKit::MolToSmiles( *molP,true, false, -1, canonical );
}
    


std::string CGRreader(RWMol *molR, const std::string cgr, bool canonical, bool setAtomMap) {
  std::cout << cgr << " is valid: "<< IsValidCGRString( cgr) << "\n";
       
    if (IsValidCGRString( cgr)){  
         return getCGRMolecule(molR, cgr, canonical, setAtomMap);
    }
    else {
         return "";
    }    
}


std::string CGRwriter(const std::string smart, bool doRandom,  unsigned int randomSeed, bool aromatize, bool signature, bool charges, int radius){
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction( smart));
    return getCGRwriter( *rxn , doRandom , randomSeed, aromatize, signature, charges, radius );
}

}
}
