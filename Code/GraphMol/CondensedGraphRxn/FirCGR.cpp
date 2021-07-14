/*
 * FirCGR.cpp defines the code to convert reaction
 * SMARTS to the new Firmenich/Bigchem CGR-smiles
 * reaction format.
 *
 * (c) RUUD, Firmenich SA, 2020
 *
 */

// Block with own imports
#include "FirCGR.h"

// Block with RDKit imports
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolHash/MolHash.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>

// Define external
#include <sstream>

namespace RDKit {
namespace Firmenich {
namespace Reactions {
 
//! Method to read a dot-separated SMILES string to a vector of molecules.
MOL_SPTR_VECT read(const std::string &smiles) {
    MOL_SPTR_VECT molvec;
    std::stringstream ss(smiles);
    for (std::string smi; std::getline(ss,smi,'.'); ) {
        ROMol* mol = RDKit::SmilesToMol(smi);
        if (mol) {
            delete mol;
            ROMOL_SPTR molptr(RDKit::SmilesToMol(smi));
            molvec.push_back(molptr);            
        }
    }
    return molvec;
}
    
CGRWriter::CGRWriter() {
}
    
CGRWriter::~CGRWriter() {
}
  
std::vector<TmpBond> CGRWriter::getBonds(RWMOL_SPTR mol, std::map<int,int> m2i) {
    // Collect the fully mapped bonds in the product.
    // TODO: We may have to expand this.
    std::vector<TmpBond> blist;
    for (Bond* bond : mol->bonds()) {
        const int from = bond->getBeginAtom()->getAtomMapNum();
        const int to = bond->getEndAtom()->getAtomMapNum();
        if (from>0 && to>0) {
            blist.push_back(TmpBond(from,to,bond->getBondType()));                        
        }
    }
    for (Atom *atom : mol->atoms()) {
        const int mnum = atom->getAtomMapNum();
        if (mnum>0)
            m2i[mnum] = atom->getIdx();
    }
    return blist;
}
    
RWMOL_SPTR CGRWriter::copy(RWMOL_SPTR in) {
    RWMOL_SPTR cp(new RWMol());
    for (Atom* atom : in->atoms()) 
        cp->addAtom(new Atom(atom->getAtomicNum()));
    for (Bond* bond : in->bonds()) {
        const unsigned int from = bond->getBeginAtomIdx();
        const unsigned int to = bond->getEndAtomIdx();
        const Bond::BondType bt = bond->getBondType();
        cp->addBond(from,to,bt);
    }
    return cp;
}
    
Bond::BondType CGRWriter::lookup(std::vector<TmpBond> list, TmpBond query) {
    auto it = std::find(list.begin(),list.end(),query);
    return it==list.end() ? Bond::ZERO : (*it).bt;
}
    
Bond::BondType CGRWriter::propose(Bond::BondType in, Bond::BondType out) {
    Bond::BondType res = Bond::ZERO;
    switch (in) {
        case Bond::ZERO:
            switch(out) {
                case Bond::SINGLE:
                    res = Bond::CGRNS;
                    break;
                case Bond::DOUBLE:
                    res = Bond::CGRND;
                    break;
                case Bond::TRIPLE:
                    res = Bond::CGRNT;
                    break;
                case Bond::AROMATIC:
                    res = Bond::CGRNA;
                    break;
                default:
                    res = in;
                    break;
            }
            break;
        case Bond::SINGLE:
           switch(out) {
                case Bond::ZERO:
                    res = Bond::CGRSN;
                    break;
                case Bond::DOUBLE:
                    res = Bond::CGRSD;
                    break;
                case Bond::TRIPLE:
                    res = Bond::CGRST;
                    break;
                case Bond::AROMATIC:
                    res = Bond::CGRSA;
                    break;
                default:
                    res = in;
                    break;
            }
            break;
        case Bond::DOUBLE:
           switch(out) {
                case Bond::ZERO:
                    res = Bond::CGRDN;
                    break;
                case Bond::SINGLE:
                    res = Bond::CGRDS;
                    break;
                case Bond::TRIPLE:
                    res = Bond::CGRDT;
                    break;
                case Bond::AROMATIC:
                    res = Bond::CGRDA;
                    break;
                default:
                    res = in;
                    break;
            }
            break;
        case Bond::TRIPLE:
           switch(out) {
                case Bond::ZERO:
                    res = Bond::CGRTN;
                    break;
                case Bond::SINGLE:
                    res = Bond::CGRTS;
                    break;
                case Bond::DOUBLE:
                    res = Bond::CGRTD;
                    break;
                case Bond::AROMATIC:
                    res = Bond::CGRTA;
                    break;
                default:
                    res = in;
                    break;
            }
            break;
        case Bond::AROMATIC:
           switch(out) {
                case Bond::ZERO:
                    res = Bond::CGRAN;
                    break;
                case Bond::SINGLE:
                    res = Bond::CGRAS;
                    break;
                case Bond::DOUBLE:
                    res = Bond::CGRAD;
                    break;
                case Bond::TRIPLE:
                    res = Bond::CGRAT;
                    break;
                default:
                    res = in;
                    break;
            }
            break;
        default:
            res=in;
            break;
    }
    return res;
}
    
Atom* lookupAtom(RWMOL_SPTR mol, const unsigned int mnum) {
    Atom* res = nullptr;
    if (mnum>0) {
        for (Atom *atom : mol->atoms()) {
            unsigned int mymnum = atom->getAtomMapNum();
            res = mymnum == mnum ? atom : nullptr;
            if(res) break;
        }
    }
    return res;
}
    
std::string CGRWriter::write(const std::string &rxnsma, bool rnd, int root) {
    std::stringstream res;
    
    // Import the vectors for reagents and products
    try {
        
        const std::string reagsmi = rxnsma.substr(0,rxnsma.find(">"));
        RWMOL_SPTR reag(new RWMol(*RDKit::SmilesToMol(reagsmi)));
        
        const std::string prodsmi = rxnsma.substr(rxnsma.find_last_of(">")+1);
        RWMOL_SPTR prod(new RWMol(*RDKit::SmilesToMol(prodsmi)));
        
        
        if (reag && prod) {
            // Collect the mapped tables
            // TODO: Collect the map nums on both site of the equation
            std::map<int,int> Rm2i,Pm2i;
            const std::vector<TmpBond> Rbonds = getBonds(reag,Rm2i);
            const std::vector<TmpBond> Pbonds = getBonds(prod,Pm2i);

            // Collect the deleted bonds => bonds no longer present in product
            //const std::vector<TmpBond> diff = difference(Rbonds,Pbonds);

            // Take the product molecules and make a modification using
            // the information in the reagent molecules. We have 3 cases:
            // 1. Bond in both or only in product => these become type CGR{N,S,D,T,A}{S,D,T,A}
            // TODO: Probably only do this if both mapnums are in the reagents (IG sometimes missing).$
            int frag=0;
            for (Bond *bond : prod->bonds()) {
                const int from = bond->getBeginAtom()->getAtomMapNum();
                const int to = bond->getEndAtom()->getAtomMapNum();
                if (from>0 && to>0) {
                    const TmpBond query = TmpBond(from,to,bond->getBondType());
                    const Bond::BondType rbt = lookup(Rbonds, query);
                    const Bond::BondType pbt = lookup(Pbonds, query);
                    if (rbt && pbt)
                        bond->setBondType(propose(rbt,pbt));
                }        
            }

            // Step 2: We now need to rebuild bonds no longer present in the products   
           for (const TmpBond &q : Rbonds) {
                auto it = std::find(Pbonds.begin(), Pbonds.end(), q);
                if (it==Pbonds.end()) {
                    // Bond is no longer present in the products 
                    // => create bond of type CGR{S,D,T,A}N
                    // Only and only if the mapped atoms are present
                    // in the product.
                    const Bond::BondType rbt = q.bt;
                    const Bond::BondType pbt = Bond::ZERO;
                    Atom* fatom = lookupAtom(prod,q.from);
                    Atom* tatom = lookupAtom(prod,q.to);
                    if (fatom && tatom)
                        prod->addBond(fatom,tatom,propose(rbt,pbt));
                }
            }

            // Step 3: Remove all mapnums
            for (Atom *atom : prod->atoms()) {
                atom->setAtomMapNum(0);
            }

            // Define the default flags.
            bool iso=false;
            bool kek=false;
            bool bondExpl = false;
            bool hExpl = false;

            // Convert to a result
            res << RDKit::MolToSmiles(*prod, iso, kek, -1, false, bondExpl, hExpl, rnd);
        } else {
            res << "invalid";
        }
        
        // Safely delete the pointers.
        
    } catch (const std::exception &e) {
        res << "error";
    }
    
    // Done    
    return res.str();
}
        
CGRReader::CGRReader() {}

CGRReader::~CGRReader() {}
    
Bond::BondType CGRReader::propose(const Bond::BondType &bt, bool &arom, bool &dearom, bool prod) {
    Bond::BondType answer = bt;
    switch(bt) {
        case Bond::CGRNS:
            answer = prod ? Bond::SINGLE : Bond::ZERO;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRND:
            answer = prod ? Bond::DOUBLE : Bond::ZERO;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRNT:
            answer = prod ? Bond::TRIPLE : Bond::ZERO;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRNA:
            answer = prod ? Bond::AROMATIC : Bond::ZERO;
            arom = prod;
            dearom = false;
            break;
        case Bond::CGRSN:
            answer = prod ? Bond::ZERO : Bond::SINGLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRSD:
            answer = prod ? Bond::DOUBLE : Bond::SINGLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRST:
            answer = prod ? Bond::TRIPLE : Bond::SINGLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRSA:
            answer = prod ? Bond::AROMATIC : Bond::SINGLE;
            arom = prod;
            dearom = false;
            break;
        case Bond::CGRDN:
            answer = prod ? Bond::ZERO : Bond::DOUBLE;
            arom = false;
            dearom = false;            
            break;
        case Bond::CGRDS:
            answer = prod ? Bond::SINGLE : Bond::DOUBLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRDT:
            answer = prod ? Bond::TRIPLE : Bond::DOUBLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRDA:
            answer = prod ? Bond::AROMATIC : Bond::DOUBLE;
            arom = prod;
            dearom = false;
            break;
        case Bond::CGRTN:
            answer = prod ? Bond::ZERO : Bond::TRIPLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRTS:
            answer = prod ? Bond::SINGLE : Bond::TRIPLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRTD:
            answer = prod ? Bond::DOUBLE : Bond::TRIPLE;
            arom = false;
            dearom = false;
            break;
        case Bond::CGRTA:
            answer = prod ? Bond::AROMATIC : Bond::TRIPLE;
            arom = prod;
            dearom = false;
            break;
        case Bond::CGRAN:
            answer = prod ? Bond::ZERO : Bond::AROMATIC;
            arom = false;
            dearom = prod;
            break;
        case Bond::CGRAS:
            answer = prod ? Bond::SINGLE : Bond::AROMATIC;
            arom = false;
            dearom = prod;
            break;
        case Bond::CGRAD:
            answer = prod ? Bond::DOUBLE : Bond::AROMATIC;
            arom = false;
            dearom = prod;
            break;
        case Bond::CGRAT:
            answer = prod ? Bond::TRIPLE : Bond::AROMATIC;
            arom = false;
            dearom = prod;
            break;
        default:
            break;
    }
    return answer;
}
    
void map(RWMOL_SPTR mol, bool mapcgratoms, bool mapallatoms) {
    // If asked to map all, map all => this is the
    // most specific solution.
    if (mapallatoms) {
        for (Atom *atom : mol->atoms()) {
            atom->setAtomMapNum(atom->getIdx()+1);
        }
    }
    
    // Check if asked to map at least the reaction site.
    else if (mapcgratoms) {
        const size_t num = mol->getNumAtoms();
        bool* cgr = (bool*) std::malloc(num*sizeof(bool));
        for (int idx=0;idx<num;idx++) {
            cgr[idx]=0;
        }
        for (Bond *bond : mol->bonds()) {
            switch(bond->getBondType()) {
                case Bond::CGRNS:
                case Bond::CGRND:
                case Bond::CGRNT:
                case Bond::CGRNA:
                case Bond::CGRSN:
                case Bond::CGRSD:
                case Bond::CGRST:
                case Bond::CGRSA:
                case Bond::CGRDN:
                case Bond::CGRDS:
                case Bond::CGRDT:
                case Bond::CGRDA:
                case Bond::CGRTN:
                case Bond::CGRTS:
                case Bond::CGRTD:
                case Bond::CGRTA:
                case Bond::CGRAN:
                case Bond::CGRAS:
                case Bond::CGRAD:
                case Bond::CGRAT:
                    cgr[bond->getBeginAtomIdx()]=1;
                    cgr[bond->getEndAtomIdx()]=1;
                    break;
                default:
                    // Ignore
                    break;
            }
        }
        // Map those atoms.        
        unsigned int mnum = 0;
        for (Atom *atom : mol->atoms()) {
            if (cgr[atom->getIdx()]) {
                atom->setAtomMapNum(++mnum);
            }
        }
        free(cgr);
    }    
}
    
std::string CGRReader::read(const std::string &cgrsmi, bool mapcgratoms, bool mapallatoms, bool can) {
    std::stringstream ss;
    ss << "";
    try {
        RWMOL_SPTR molptr1(new RWMol(*RDKit::SmilesToMol(cgrsmi)));
        RWMOL_SPTR molptr2(new RWMol(*RDKit::SmilesToMol(cgrsmi)));
        
        ss << "";
        if (molptr1 && molptr2) {
            
            // Apply the mapping as specified
            map(molptr1, mapcgratoms, mapallatoms);
            map(molptr2, mapcgratoms, mapallatoms);
            
            // Import two times the molecule: One goes to reagents one for more products
            // For all CGR-bonds, we convert them to the starting type in reagents
            // or the product type in products.
            // TODO: Think about numbering the atoms around the CGR bonds.

            // Do the reagent
            std::vector<Bond*> erase;
            for (Bond *bond : molptr1->bonds()) {
                bool arom,dearom;
                Bond::BondType bt = propose(bond->getBondType(),arom,dearom,false);
                if (bt == Bond::ZERO) {
                    // Remove
                    erase.push_back(bond);
                } else {
                    // Change
                    bond->setBondType(bt);
                }
            }
            for (Bond* bond : erase) 
                molptr1->removeBond(bond->getBeginAtomIdx(),bond->getEndAtomIdx());

            // Do the product
            erase.clear();
            for (Bond *bond : molptr2->bonds()) {
                bool arom,dearom;
                Bond::BondType bt = propose(bond->getBondType(),arom,dearom,true);
                if (bt == Bond::ZERO) {
                    // Remove
                    erase.push_back(bond);
                } else {
                    // Change
                    bond->setBondType(bt);
                }
            }
            for (Bond* bond : erase) 
                molptr2->removeBond(bond->getBeginAtomIdx(),bond->getEndAtomIdx());

            // Write the result with export, reimport and export.
            // Probably there's a better solution but it works for now.
            const std::string reagsmi = RDKit::MolToSmiles(*molptr1, true, false, -1, can, false, false, false );
            const std::string rredo = RDKit::MolToSmiles(*RDKit::SmilesToMol(reagsmi), true, false, -1, can, false, false, false );
            const std::string prodsmi = RDKit::MolToSmiles(*molptr2, true, false, -1, can, false, false, false );
            const std::string predo = RDKit::MolToSmiles(*RDKit::SmilesToMol(prodsmi), true, false, -1, can, false, false, false );
            ss << rredo << ">>" << predo;
        } else {
            ss << "invalid";
        }
    } catch (const std::exception &e) {
        ss << "error";
    }
    // Return the response
    return ss.str();
}
    
std::vector<unsigned int> getNeighborsPtr1(RWMOL_SPTR mol, Atom* atom) {
   std::vector<unsigned int> res;
   for(const auto &nbri: make_iterator_range(mol->getAtomBonds(atom))) {
       const RDKit::Bond *bond = (*mol)[nbri];
      res.push_back(bond->getOtherAtomIdx(atom->getIdx()));
   }
   return res;
}
    
bool number(Bond::BondType bt) {
    switch (bt) {
       case Bond::CGRNS:
        case Bond::CGRND:
        case Bond::CGRNT:
        case Bond::CGRNA:
        case Bond::CGRSN:
        case Bond::CGRSD:
        case Bond::CGRST:
        case Bond::CGRSA:
        case Bond::CGRDN:
        case Bond::CGRDS:
        case Bond::CGRDT:
        case Bond::CGRDA:
        case Bond::CGRTN:
        case Bond::CGRTS:
        case Bond::CGRTD:
        case Bond::CGRTA:
        case Bond::CGRAN:
        case Bond::CGRAS:
        case Bond::CGRAD:
        case Bond::CGRAT:          
            return true;
        default:
            return false;
    }
}
 
    
//! Method to convert CGRSmiles to SMARTS reaction site.
std::string CGRReader::siteSmarts(const std::string &cgrsmiles, bool siteplusone) {
        std::stringstream ss;

        //std::cout << "go:" << std::endl;

        // Define input and output and map all atoms remove sanitizedx
        RWMOL_SPTR molptr1(new RWMol(*RDKit::SmilesToMol(cgrsmiles,0, false)));
 
        if (molptr1) {
            //std::cout << "mol1 found!" << std::endl;

            /// this method try to sanitize the object molecule
            // we deactivate some methods that can failed with CGR bonds
            unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_KEKULIZE ^ MolOps::SANITIZE_PROPERTIES;
            unsigned int failed;
            
            try {
                MolOps::sanitizeMol(*molptr1, failed, sanitizeOps);
            }
            catch (MolSanitizeException &) {
                BOOST_LOG(rdInfoLog) << "Issue to sanitize Mol step 1\n";
            }
            
            // Define a bool* array to indicate the desired atoms
            const unsigned int num = molptr1->getNumAtoms();
            bool* keep = (bool*) std::malloc(num*sizeof(bool));
            for (unsigned int idx=0; idx<num; idx++) keep[idx]=0;

            // Keep all atoms with a CGR-bond
            // Keep all it's neighbors
            for (Bond* bond : molptr1->bonds()) {
                const Bond::BondType bt = bond->getBondType();
                if (number(bt)) {
                    // Keep the atoms
                    Atom* beginAtom = bond->getBeginAtom();
                    Atom* endAtom = bond->getEndAtom();
                    keep[beginAtom->getIdx()] = 1;
                    keep[endAtom->getIdx()] = 1;
                    // Keep the neighbors of the atoms
                    if (siteplusone) {
                        std::vector<unsigned int> nbrs1 = getNeighborsPtr1(molptr1,beginAtom);
                        for (const unsigned int &n : nbrs1) {
                            
                            if (!molptr1->getAtomWithIdx(n)->getIsAromatic()) {
                                 keep[n]=1;
                            }
                        }
                        std::vector<unsigned int> nbrs2 = getNeighborsPtr1(molptr1,endAtom);
                        for (const unsigned int &n : nbrs2)  {
                            if (!molptr1->getAtomWithIdx(n)->getIsAromatic()) {
                                keep[n]=1;
                            }
                        }
                    }
                }
            }

            // Delete all undesired atoms, then export to stringstream
            unsigned int mnum = 0;
            for (int idx=num-1;idx>=0;idx--) {
                if (!keep[idx]) 
                    molptr1->removeAtom(idx);
            }
            free(keep);
            //std::cout << "mol1 start aromatic check" << std::endl;

            // We not to export and reimport right here to make sure we have standardized the molecule.
            // After the step we can import and put map numbers
            
            if (!molptr1->getRingInfo()->isInitialized()) {
                MolOps::fastFindRings(*molptr1);
            }
            
            for (Atom* at: molptr1->atoms()) {
                if (at->getIsAromatic() &&  !molptr1->getRingInfo()->numAtomRings(at->getIdx())) {
                    at->setIsAromatic(false);
                    //std::cout << "remove atomatize to" << at->getIdx() << std::endl;
                }
            }
            
            unsigned int failed1;

            try {
                MolOps::sanitizeMol(*molptr1, failed1, sanitizeOps);
            }
            catch (MolSanitizeException &) {
                BOOST_LOG(rdInfoLog) << "Issue to sanitize Mol step 2\n";
            }
            
            bool can = true;
            const std::string tmpsmi = RDKit::MolToSmiles(*molptr1, true, false, -1, can, false, false, false);
            
            //std::cout << "mol1 done aromatic check =>" << tmpsmi  << std::endl;

            // Reimport and number all atoms. remove sanitized
            RWMOL_SPTR molptr2(new RWMol(*RDKit::SmilesToMol(tmpsmi, 0 , false)));
            RWMOL_SPTR molptr3(new RWMol(*RDKit::SmilesToMol(tmpsmi, 0 , false)));
            map(molptr2, false, true);
            map(molptr3, false, true);            

            // Do the reagent
            std::vector<Bond*> erase;
            for (Bond *bond : molptr2->bonds()) {
                bool arom,dearom;
                Bond::BondType bt = propose(bond->getBondType(),arom,dearom,false);
                if (bt == Bond::ZERO) {
                    // Remove
                    erase.push_back(bond);
                } else {
                    // Change
                    bond->setBondType(bt);
                }
            }
            for (Bond* bond : erase)
                molptr2->removeBond(bond->getBeginAtomIdx(),bond->getEndAtomIdx());
          
            //std::cout << "mol2 start broken ring dearomatize check" << std::endl;
            
            if (!molptr2->getRingInfo()->isInitialized()) {
                MolOps::fastFindRings(*molptr2);
            }
            
            for (Atom* at: molptr2->atoms()) {
                if (at->getIsAromatic() &&   !molptr2->getRingInfo()->numAtomRings(at->getIdx())) {
                    at->setIsAromatic(false);
                    //std::cout << "remove atomatize to" << at->getIdx() << std::endl;
                }
            }
            
            unsigned int failed2;
            try {
                MolOps::sanitizeMol(*molptr2, failed2, sanitizeOps);
            }
            catch (MolSanitizeException &) {
                BOOST_LOG(rdInfoLog) << "Issue to sanitize step 3\n";
            }
            
            //std::cout << "mol2 done aromatic check" << std::endl;
            //std::cout << "mol3 started" << std::endl;

            // Do the product
            erase.clear();
            for (Bond *bond : molptr3->bonds()) {
                bool arom,dearom;
                Bond::BondType bt = propose(bond->getBondType(),arom,dearom,true);
                if (bt == Bond::ZERO) {
                    // Remove
                    erase.push_back(bond);
                } else {
                    // Change
                    bond->setBondType(bt);
                }
            }
            for (Bond* bond : erase) 
                molptr3->removeBond(bond->getBeginAtomIdx(),bond->getEndAtomIdx());

            //std::cout << "mol3 start aromatic check" << std::endl;

            if (!molptr3->getRingInfo()->isInitialized()) {
                MolOps::fastFindRings(*molptr3);
            }
            
            for (Atom* at: molptr3->atoms()) {
                if (at->getIsAromatic() &&   !molptr3->getRingInfo()->numAtomRings(at->getIdx())) {
                    at->setIsAromatic(false);
                    //std::cout << "remove atomatize to" << at->getIdx() << std::endl;
                }
            }
            
            unsigned int failed3;
            try {
                MolOps::sanitizeMol(*molptr3, failed3, sanitizeOps);
            }
            catch (MolSanitizeException &) {
                BOOST_LOG(rdInfoLog) << "Issue to sanitize step 4\n";
            }
            
            //std::cout << "mol3 done aromatic check" << std::endl;

            //std::cout << "sanitize :"  << failed1 << " , " << failed2 << " , " << failed3 << std::endl;
            
            // Sanitize need to be false!
            // Write the result with export, reimport and export.
            // Probably there's a better solution but it works for now.
            const std::string reagsmi = RDKit::MolToSmiles(*molptr2, true, false, -1, can, false, false, false );
            
            //std::cout << "reagsmi :"  << reagsmi << std::endl;

            const std::string rredo = RDKit::MolToSmiles(*RDKit::SmilesToMol(reagsmi,0, false), true, false, -1, can, false, false, false );
            
            //std::cout << "rredo :"  << rredo << std::endl;

            const std::string prodsmi = RDKit::MolToSmiles(*molptr3, true, false, -1, can, false, false, false );

            //std::cout << "prodsmi :"  << prodsmi << std::endl;

            
            const std::string predo = RDKit::MolToSmiles(*RDKit::SmilesToMol(prodsmi,0, false), true, false, -1, can, false, false, false );

            //std::cout << "predo :"  << predo << std::endl;
            ss << rredo << ">>" << predo;
        
        }

    return ss.str();
}
    
    
}
}
}
