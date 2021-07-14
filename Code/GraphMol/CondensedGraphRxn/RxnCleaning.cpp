#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/MolHash/MolHash.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include "CondensedGraphRxn.h"
#include <sstream>
#include <algorithm>
#include <set>

template<class T>
void printSet(std::vector<T> valset) {
    std::stringstream ss;
    ss << "{";
    int v=0;
    for (const T value : valset) {
        if (++v>1) ss << ",";
        ss << value;
    }
    ss << "}";
    std::cout << ss.str() << std::endl;
}

namespace RDKit {

    namespace CondensedGraphRxn {

        // structures definitions

        struct bondid {

            unsigned int f, t, fidx, tidx;

            std::string bt;

            bondid(unsigned int fromnum, unsigned int tonum, unsigned int fromidx, unsigned int toidx, std::string bondtype) {

                f = std::min(fromnum, tonum);

                t = std::max(fromnum, tonum);

                tidx = std::max(fromidx, toidx);

                fidx = std::min(fromidx, toidx);

                bt = bondtype;

            }

            const bool operator==(const bondid &other) {

                return f == other.f && t == other.t && bt == other.bt;

            }

        };


        struct NbrType {

            unsigned int bond, atno;

            NbrType(unsigned int at, unsigned int bt) {

                bond = bt;

                atno = at;

            }

            // comparation equal
            const bool operator==(const NbrType &other) {

                return atno == other.atno &&  bond == other.bond;

            }

            // comparation atomtype smaller
            bool operator<(const NbrType &other) const {

                return atno < other.atno || bond < other.bond;
            }

        };


        struct AtomType {

            unsigned int mnum, atno, idx;
            ROMOL_SPTR source;
            bool reactant;
            std::size_t inv;
            std::vector<uint16_t> nbrs;

            AtomType(unsigned int m, unsigned int i, unsigned int a,
                     std::size_t myinv, ROMOL_SPTR mol, std::vector<uint16_t> nei, bool r) {



                inv = myinv;
                mnum = m;
                atno = a;
                idx = i;
                source = mol;
                reactant = r;
                nbrs = nei;
            }

            // empty constructor for the map object construction
            AtomType() {};

            bool operator==(const AtomType &other) const {
                return mnum == other.mnum && atno == other.atno && inv == other.inv;
            }

            // comparation atomtype smaller
            bool operator<(const AtomType &other) const {

                return idx < other.idx;
            }



        };


        struct mapBond {

            mapBond() {};

            unsigned int from, to;

            mapBond( unsigned int f, unsigned int t) {

                from = std::min(f, t);
                to = std::max(f, t);
            }

            const bool operator==(const mapBond &other) {
                return from == other.from && to == other.to;
            }

            const bool operator<(const mapBond &other) {
                return from < other.from || to < other.to;
            }

        };

        struct matchPair {

            matchPair() {};

            RDKit::MatchVectType mvec;
            unsigned int t;

            matchPair(RDKit::MatchVectType mvec, unsigned int t) {

                this->mvec = mvec;
                this->t = t;

            }

            bool operator<(const matchPair &other) const {
                return t > other.t;
            }

        };

        void clearAtomMap(std::vector<ROMOL_SPTR> &mollist) {
            //unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_KEKULIZE ^ MolOps::SANITIZE_PROPERTIES;


            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                
                
                 
                /* unsigned int failed;
                 
                 try {
                 MolOps::sanitizeMol(*it, failed, sanitizeOps);
                 }
                 catch (MolSanitizeException &) {
                 BOOST_LOG(rdInfoLog) << "Issue to sanitize\n";
                 }
                */

                
                
                for (auto atom : (*it)->atoms()) {
                    if (atom) {
                        atom->setAtomMapNum(0);
                    }
                }
            }
        }
        
        
        
        // convert to string the results
        const std::string convertMolstoSmiles(const std::vector<ROMOL_SPTR> &listofmol, bool debug = false) {
            std::stringstream ss;
            int f = 0;
            if (debug) {
                std::cout << " ...Converting reaction frags to Smiles" << std::endl;
            }

            for (ROMOL_SPTR fragmol : listofmol) {
                ss << (f++ == 0 ? "" : ".");
                ss << RDKit::MolToSmiles(*fragmol);
            }
            if (debug) {
                std::cout << "conversion done..." << std::endl;
            }
            return ss.str();

        }

        
        
        
        // convert to string the results
        const std::string VectorSmilestoSmiles(std::vector<std::string> smis, bool debug = false) {
            std::stringstream ss;
            int f = 0;
            if (debug) {
                std::cout << " ...Converting reaction frags to Smiles" << std::endl;
            }
            
            for (std::string smi : smis) {
                ss << (f++ == 0 ? "" : ".");
                ss << smi;
            }
            if (debug) {
                std::cout << "conversion done..." << std::endl;
            }
            return ss.str();
            
        }

        

        // invariant definition of atom env neighbor
        std::size_t vhash(std::vector<uint16_t> &vec) {

            std::size_t seed = vec.size();

            for(auto &i : vec) {

                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            }

            return seed;

        }

        void getMaxMapNumList(std::vector<ROMOL_SPTR> mollist, unsigned int &maxmap) {
            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                maxmap = std::max(maxmap, (*it)->getMaxMapNum());
            }
        }



        // check if an mapped fragment is full enclosed (part) of another mapped molecule
        bool isEnclosed(ROMOL_SPTR mol, std::vector<ROMOL_SPTR> mollist, bool debug ) {
            bool enclosed = false;
            std::vector< unsigned int > rset = mol->getMappedAtomSet();

            if (debug) {
            std::cout << "mol rset:" << std::endl;

            printSet<unsigned int>(rset); 
            }
            std::sort(rset.begin(), rset.end()); 

            std::vector<ROMOL_SPTR>::iterator it;
            // double test to stop the loop on first inclusion found!
            // TO DO : must check that mol is a combination of multiple reagant and not only one like here!
            // possible solution use a while to expend the list of mol combinatorial n-1 join !
            // or use not Enclose solution mean 
            // [NH2:1][c:2]1[cH:3][cH:4][c:5]([Cl:6])[cH:7][n:8]1.[Cl:9][c:10]1[cH:11][cH:12][c:13]([CH:14]=[O:15])[cH:16][cH:17]1.[CH3:18][CH2:19][CH2:20][N:21]([CH2:22][CH2:23][CH3:24])[C:25](=[O:26])[C:27]#[CH:28]>>[CH3:18][CH2:19][CH2:20][N:21]([CH2:22][CH2:23][CH3:24])[C:25](=[O:26])[CH2:27][c:28]1[c:14]([n:1][c:2]2[cH:3][cH:4][c:5]([Cl:6])[cH:7][n:8]12)-[c:13]1[cH:16][cH:17][c:10]([Cl:9])[cH:11][cH:12]1.[OH2:15]
            


            for (it = mollist.begin(); it != mollist.end() && !enclosed; ++it) {


                std::vector< unsigned int > set = (*it)->getMappedAtomSet();
                if (debug) {
                std::cout << "mol pset:" << std::endl;
                printSet<unsigned int>(set); 
                }
                std::sort(set.begin(), set.end()); 

                enclosed = std::includes(set.begin(), set.end(), rset.begin(), rset.end());
                if (debug) {
                    std::cout << "enclosed:" << enclosed << std::endl;
                }

            }


            return enclosed;
        }

        // check if an mapped fragment is covered in destination list
        bool isCovered(ROMOL_SPTR mol, std::vector<unsigned int> dstlist, bool debug ) {
            bool covered = true;
            std::vector< unsigned int > molset = mol->getMappedAtomSet();
            
            if (debug) {
                std::cout << "is covered molset:" << std::endl;

                printSet<unsigned int>(molset); 
                printSet<unsigned int>(dstlist); 
            }
            if (std::find(molset.begin(), molset.end(), 0) != molset.end()) {
                return false;
            }

            std::sort(molset.begin(), molset.end());

            std::vector< unsigned int > diff;
            // caution of the order!!!!
            std::set_difference(molset.begin(), molset.end(), dstlist.begin(), dstlist.end(), std::inserter(diff, diff.begin()));
            if (debug) {
                std::cout << "is covered diff:" << std::endl;

                printSet<unsigned int>(diff); 
            }
            if (diff.size() > 0) {
                return false;
            }

            return covered;
        }


        void mapAllAtoms(std::vector<ROMOL_SPTR> mollist, unsigned int &maxmap) {
            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                (*it)->addMapping(maxmap);
            }
        }


        void removeMonoAtomreagents(std::vector<ROMOL_SPTR> &molR, std::vector<ROMOL_SPTR> &molP) {
            std::vector<ROMOL_SPTR>::iterator itp;
            std::vector<ROMOL_SPTR>::iterator itr;

            for (itp = molP.end(); itp-- != molP.begin();) {
                unsigned int nbat = (*itp)->getNumHeavyAtoms();
                bool era = true;
                if (nbat == 1) {
                    unsigned int datommapnum = (*itp)->getAtomWithIdx(0)->getAtomMapNum();
                    unsigned int datomicid = (*itp)->getAtomWithIdx(0)->getAtomicNum();

                    for (itr = molR.end(); itr-- != molR.begin();) {
                        unsigned int nbatr = (*itr)->getNumHeavyAtoms();
                        if (nbatr == 1 && era) {
                            unsigned int atommapnum = (*itr)->getAtomWithIdx(0)->getAtomMapNum();
                            unsigned int atomicid = (*itr)->getAtomWithIdx(0)->getAtomicNum();
                            if (datommapnum == atommapnum && datomicid == atomicid) {
                                molR.erase(itr);
                                molP.erase(itp);
                                era = false;
                            }
                        }
                    }
                }
            }
        }

        bool comparedMolSized(ROMOL_SPTR mol1, ROMOL_SPTR mol2 ) {
            return mol1->getNumAtoms() > mol2->getNumAtoms();
        }


        unsigned int getMMN(ROMOL_SPTR mol) {
            unsigned int res = 0;
            for (auto atom : mol->atoms()) {
                if (atom) {
                    unsigned int atnum = atom->getAtomMapNum();
                    res = std::max(res, atnum);
                }
            }
            return res;
        }


        bool sortByMapNum(ROMOL_SPTR mol1, ROMOL_SPTR mol2 ) {
            bool res = false;
            bool orderbymaxmap = false;
            if (mol1 == nullptr ||  mol2 == nullptr) {
                return false;
            }

            bool maxmapnum = false;

            bool numat = mol1->getNumAtoms(0) > mol2->getNumAtoms(0);
            if(orderbymaxmap) {
                maxmapnum = getMMN(mol1) > getMMN(mol2);
                return numat || maxmapnum;
            }


            return numat;
        }

        void finalCleanUp(std::vector<ROMOL_SPTR> mollist,  std::vector<ROMOL_SPTR> &reslist, bool debug = false) {

            if (debug) {
                std::cout << "Final CleanUp start:" << std::endl;
            }

            std::vector<ROMOL_SPTR> remove;
            std::vector<ROMOL_SPTR>::iterator it, it1, it2;
            for (it1 = mollist.begin(); it1 != mollist.end(); ++it1) {
                std::vector< unsigned int > set1 = (*it1)->getMappedAtomSet();
                std::sort(set1.begin(), set1.end());
                if (debug) {
                    std::cout << "mol 1:" << (*it1)->getNumAtoms(0) << std::endl;
                }
                for (it2 = mollist.begin(); it2 != mollist.end(); ++it2) {
                    if (it2 < it1) {
                        if (debug) {
                            std::cout << "mol 2:" << (*it2)->getNumAtoms(0) << std::endl;
                        }
                        std::vector< unsigned int > set2 = (*it2)->getMappedAtomSet();
                        std::sort(set2.begin(), set2.end());
                        if (std::includes(set1.begin(), set1.end(), set2.begin(), set2.end())) {
                            remove.push_back((*it2));
                        } else if (std::includes(set2.begin(), set2.end(), set1.begin(), set1.end())) {
                            remove.push_back((*it1));
                        }
                    }
                }
            }

            if (remove.size() > 0) {
                if (debug) {
                    std::cout << "Remove size > 0 find" << std::endl;
                }

                for (it = mollist.begin(); it != mollist.end(); ++it) {
                    auto result = std::find(remove.begin(), remove.end(), (*it));
                    if (result == remove.end()) {
                        if ((*it) && (*it)->getNumHeavyAtoms() > 0 ) {
                            // only add none nullptr molecules
                            reslist.push_back((*it));
                        }
                    }
                }
            } else {
                if (debug) {
                    std::cout << "copy mollist without Remove" << std::endl;
                }
                for (it = mollist.begin(); it != mollist.end(); ++it) {
                    if ((*it) == nullptr) {
                        std::cout << "bingo" << std::endl;
                    }

                    else if ((*it)->getNumHeavyAtoms() > 0) {
                        // only add none nullptr molecules
                        reslist.push_back((*it));
                    }

                }
            }

            if (debug) {
                std::cout << "Get Final list: " << reslist.size()  << std::endl;
            }
            if (reslist.size() > 1) {
                std::sort(reslist.begin(), reslist.end(), sortByMapNum);
            }
            if (debug) {
                std::cout << "Sorted fragments" << std::endl;
            }
            remove.clear();
            if (debug) {
                std::cout << "Clear remove pointer" << std::endl;
            }
        }

        std::size_t atomNbrTypeInv(const Atom *atm, std::vector<uint16_t> &nbridx) {
            ROMol::ADJ_ITER begin, end;
            ROMol parent = atm->getOwningMol();
            boost::tie(begin, end) = parent.getAtomNeighbors(atm);
            unsigned int atmidx = atm->getIdx();
            std::vector<NbrType> nbrs;
            while (begin != end) {
                const Atom *at = parent.getAtomWithIdx(*begin);
                const Bond *bond = parent.getBondBetweenAtoms(atmidx, at->getIdx());
                double bt = bond->getBondTypeAsDouble();
                unsigned int ibt = bt == 1.5 ? 4 : (int) bt;
                nbridx.push_back(at->getIdx());
                nbrs.push_back(NbrType(at->getAtomicNum(), ibt));

                ++begin;
            }

            std::sort(nbrs.begin(), nbrs.end());
            std::vector<uint16_t> vec;
            for(const NbrType &n : nbrs) {
                vec.push_back(n.atno);
                vec.push_back(n.bond);
            }

            std::size_t res = vhash(vec);
            // remove temp variables
            vec.clear();
            nbrs.clear();

            return res;
        }


        void getAtomInventory(std::vector<ROMOL_SPTR> mollist,
                              std::map<unsigned int, AtomType> &mapatoms,
                              bool &nonuniqueatnum, bool reactant = true) {
            nonuniqueatnum = false;
            std::map<unsigned int, unsigned int> AtAtno;

            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                for (auto atom : (*it)->atoms()) {
                    unsigned int atnum = atom->getAtomMapNum();
                    if (atnum > 0) {

                        std::vector<uint16_t> nbrlist;
                        std::size_t inv = atomNbrTypeInv(atom, nbrlist);
                        unsigned int atidx = atom->getIdx();
                        unsigned int atno = atom->getAtomicNum();
                        AtomType at(atnum, atidx, atno, inv, *it, nbrlist, reactant);
                        mapatoms[atnum] = at;
                        if ( AtAtno.count(atnum)) {
                            nonuniqueatnum = true;
                        } else {
                            AtAtno[atnum] = atno;
                        }
                    }
                }
            }
        }

        void getBondInventory(std::vector<ROMOL_SPTR> mollist, std::vector<bondid> &full ) {
            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                for (auto bond : (*it)->bonds()) {
                    if ( bond->getNumAtomMaps() == 2) {
                        unsigned int f = bond->getBeginAtom()->getAtomMapNum();
                        unsigned int fidx = bond->getBeginAtom()->getIdx();
                        unsigned int t = bond->getEndAtom()->getAtomMapNum();
                        unsigned int tidx = bond->getEndAtom()->getIdx();
                        std::string bt = bond->getBondTypeWithAtoms(false);
                        full.push_back(bondid(f, t, fidx, tidx, bt));
                    }
                }
            }
        }

        void getDiffAtoms(std::map<unsigned int, AtomType> ratoms,
                          std::map<unsigned int, AtomType> patoms,
                          std::map<ROMOL_SPTR, std::vector<AtomType> > &rdiff,
                          std::map<ROMOL_SPTR, std::vector<AtomType> > &pdiff) {

            std::map<unsigned int, AtomType>::iterator rirat, pirat;

            for (rirat = ratoms.begin(); rirat != ratoms.end(); ++rirat) {

                AtomType rat = rirat->second;
                ROMOL_SPTR rmol = rat.source;
                unsigned int rnum = rirat->first;


                pirat = patoms.find(rnum);

                if (pirat == patoms.end()) {
                    // add to rdiff
                    rdiff[rmol].push_back(rat);

                } else {
                    AtomType pat = pirat->second;
                    if (rat.inv != pat.inv) {
                        // add to rdiff
                        rdiff[rmol].push_back(rat);

                    }

                }
            }

            for (pirat = patoms.begin(); pirat != patoms.end(); ++pirat) {


                AtomType pat = pirat->second;
                ROMOL_SPTR pmol = pat.source;

                unsigned int pnum = pirat->first;

                rirat = ratoms.find(pnum);

                if (rirat == ratoms.end()) {
                    // add to pdiff
                    pdiff[pmol].push_back(pat);

                } else {
                    AtomType rat = rirat->second;
                    if (rat.inv != pat.inv) {
                        // add to pdiff
                        pdiff[pmol].push_back(pat);

                    }

                }
            }
        }



        void getConservedBonds(std::vector<bondid> rbond, std::vector<bondid> pbond,
                               std::vector<bondid> &conserved , std::vector<bondid> &changed) {
            std::vector<bondid>::iterator rit, pit;
            // intersection only ie on both sites
            
            conserved.clear();
            changed.clear();

            for (rit = rbond.begin(); rit != rbond.end(); ++rit) {
                pit = std::find(pbond.begin(), pbond.end(), (*rit));

                if (pit != pbond.end()) {
                    conserved.push_back(*rit);
                }
                else {
                    changed.push_back(*rit);

                }
            }
        }


        void getBtab(std::vector<ROMOL_SPTR> mollist, std::vector< mapBond > &btab) {
            std::vector<ROMOL_SPTR>::iterator it;
            for (it = mollist.begin(); it != mollist.end(); ++it) {
                for (auto bond : (*it)->bonds()) {
                    unsigned int f = bond->getBeginAtom()->getAtomMapNum();
                    unsigned int t = bond->getEndAtom()->getAtomMapNum();
                    if (t > 0 && f > 0) {

                        mapBond mb = mapBond(f, t);
                        btab.push_back(mb);
                    }
                }
            }
        }


        void getAtomNbrs(Atom *atm, std::vector< Atom *> &nbrs) {
            ROMol::ADJ_ITER begin, end;
            ROMol parent = atm->getOwningMol();
            boost::tie(begin, end) = parent.getAtomNeighbors(atm);
            unsigned int atmidx = atm->getIdx();

            while (begin != end) {
                Atom *at = parent.getAtomWithIdx(*begin);

                nbrs.push_back(at);

                ++begin;
            }

        }

        const unsigned int getNumHs( Bond::BondType BT) {

            unsigned int addH;
            switch (BT) {
            case Bond::BondType::SINGLE:
                addH = 1;
                break;
            case Bond::BondType::DOUBLE:
                addH = 2;
                break;
            case Bond::BondType::TRIPLE:
                addH = 3;
                break;
            default:
                addH = 0;
            }
            return addH;
        }


        const bool isFullyMapped(ROMOL_SPTR mol) {

            bool res = true;
            for (auto atom : mol->atoms()) {
                if (res) {
                    res = atom->getAtomMapNum() > 0;
                }
            }
            return res;
        }

        const bool isValidMapping(ROMOL_SPTR mol, std::vector< mapBond > btab) {

            bool valid = true;

            for (auto bond : mol->bonds()) {
                if (valid) {
                    unsigned int fnum = bond->getBeginAtom()->getAtomMapNum();
                    unsigned int tnum = bond->getEndAtom()->getAtomMapNum();

                    if (fnum > 0 && tnum > 0) {
                        mapBond query(fnum, tnum);
                        if (std::find(btab.begin(), btab.end(), query) == btab.end()) {
                            valid = false;
                        }
                    }
                }
            }

            return valid;

        }

        bool isAllAtomsValid(ROMOL_SPTR mol) {
            bool res = true;
            for (auto atom : mol->atoms()) {
                if (res) {
                    res = (atom);
                }
            }
            return res;
        }

        void addGroup(std::map<ROMOL_SPTR, std::vector<AtomType> > diff,
                        std::vector<bondid > conserved, std::vector<ROMOL_SPTR> &dst, 
                        unsigned int prevmax, bool debug) {

            // get the full map nums in dst molecule list
            std::vector<ROMOL_SPTR>::iterator itmol;
            // double test to stop the loop on first inclusion found!
            std::vector< unsigned int > dstnums;
            for (itmol = dst.begin(); itmol != dst.end(); ++itmol) {
                std::vector< unsigned int > molset = (*itmol)->getMappedAtomSet();
                dstnums.insert(dstnums.end(), molset.begin(), molset.end());
            }
            std::sort(dstnums.begin(), dstnums.end());

            std::map<ROMOL_SPTR, std::vector<AtomType> > ::iterator it;
            for (it = diff.begin(); it != diff.end(); ++it) {
                ROMOL_SPTR mol = it->first;
                std::vector<AtomType> changes = it->second;
                

                // unchanged subtructure in product need to be keeped
                std::shared_ptr<RWMol> cpMol( new RWMol( *mol ) );
                std::vector< AtomType>::iterator cit;

                if (debug) {
                    std::cout << "try cutting :" << MolToSmiles(*mol) << std::endl;
                }

                // std::cout << "look at the changes of the input molecule\n";
                for (cit = changes.begin(); cit != changes.end(); ++cit) {
                    AtomType cat = (*cit);
                    Atom *atom = cpMol->getAtomWithIdx(cat.idx) ;
                    unsigned int atomnum = atom->getAtomMapNum();
                    unsigned int atidx = atom->getIdx();
                    std::vector<uint16_t> nbrlist = cat.nbrs;
                    // std::cout << "find a Atom Neighbors\n";
                    for (uint16_t nidx : nbrlist) {
                        unsigned int natidx = (int) nidx;
                        Atom *nat = cpMol->getAtomWithIdx(natidx);
                        unsigned int natomnum = nat->getAtomMapNum();
                        Bond *bond = cpMol->getBondBetweenAtoms(natidx, atidx) ;
                        if (!bond) {
                            // std::cout << "No bond between " << natidx << " and " << atidx  << std::endl;
                        } else {
                            if (atomnum > 0 && natomnum > 0) {
                                bondid query(atomnum, natomnum, atidx, natidx, bond->getBondTypeWithAtoms(false));
                                if (std::find(conserved.begin(), conserved.end(), query) != conserved.end()) {
                                    if (debug) {
                                        std::cout << "Conserved: " <<  atomnum << ", " << natomnum  << std::endl;
                                    }
                                    continue;
                                }
                            }

                            if (debug) {
                                      std::cout << "cutting :" << natomnum << "," << atomnum << std::endl;
                            }
                            // std::cout << "find a bond to cut\n";
                            if (atomnum <= prevmax || natomnum <= prevmax) {
                                Bond::BondType BT = bond->getBondType();

                                if (BT != Bond::BondType::AROMATIC) {
                                    // remove the bond
                                    unsigned int addH = getNumHs(BT);
                                    cpMol->removeBond(natidx, atidx);

                                    // add Hs
                                    unsigned int newidx1, newidx2;
                                    for (unsigned int i = 0; i < addH; i++) {
                                        newidx1 = cpMol->addAtom();
                                        Atom *newatom1 = cpMol->getAtomWithIdx(newidx1);
                                        newatom1->setAtomicNum(1);
                                        cpMol->addBond(natidx, newidx1, Bond::BondType::SINGLE);

                                        newidx2 = cpMol->addAtom();
                                        Atom *newatom2 = cpMol->getAtomWithIdx(newidx2);
                                        newatom2->setAtomicNum(1);
                                        cpMol->addBond(atidx, newidx2, Bond::BondType::SINGLE);
                                    }
                                }
                            }
                        }
                    }
                }
                // add new mols to dist list

                    if (debug) {
                                      std::cout << "before FragMol list"  << std::endl;
                            }


                std::vector<ROMOL_SPTR> listofmol = MolOps::getMolFrags(*cpMol, false); // in Ruud python code value is set to true ???
                if (debug) {
                                      std::cout << "after FragMol list"  << std::endl;
                            }
                for (ROMOL_SPTR fragmol : listofmol) {
                    if (debug) {
                              std::cout << "in FragMol loop"  << std::endl;
                    }
                    if (fragmol && isAllAtomsValid(fragmol) && fragmol->getNumHeavyAtoms() > 0) {
                        //try {
                            if (debug) {
                                 std::cout << "try to removeHs" << std::endl;
                            }

                            //         bool implicitOnly = false, updateExplicitCount = true , sanitize false^
                            ROMOL_SPTR molok(MolOps::removeHs(*fragmol, false, true, false));
                            /*
                            unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_KEKULIZE ^ MolOps::SANITIZE_PROPERTIES;
                            
                            unsigned int failed;
                            
                            try {
                                MolOps::sanitizeMol(*fragmol, failed, sanitizeOps);
                            }
                            catch (MolSanitizeException &) {
                                BOOST_LOG(rdInfoLog) << "Issue to sanitize fragmol after removeHs\n";
                            }
                            */
                            if (molok && isAllAtomsValid(molok) && molok->getNumHeavyAtoms() > 0) {
                                if (!isCovered(molok, dstnums, debug)) {
                                    if (debug) {
                                        std::cout << MolToSmiles(*molok)  << std::endl;
                                        std::cout << "=> is not covered" << std::endl;
                                    }
                                    if (!isEnclosed(molok, dst, debug)) {

                                        dst.push_back(molok);
                                        if (debug) {
                                             std::cout << MolToSmiles(*molok)  << std::endl;
                                             std::cout << " => mol cut added "  << std::endl;
                                        }
                                    }
                                }
                            }

                       /* } catch ( RDKit::MolSanitizeException &e ) {
                            if (debug) {
                                std::cout << "Error Remove Hs" << std::endl;
                            }
                            continue;
                        }*/
                    }
                }
            }
        }



        void getComputeSolutions(ROMOL_SPTR mol1, ROMOL_SPTR mol2, std::vector< mapBond > btab,
                                 std::vector<RDKit::MatchVectType> matches, std::vector<matchPair> &valid_hit, bool debug ) {

            //bool found = false;
            for ( size_t i = 0 ; i < matches.size() ; ++i ) {
                // need to take best mapping match in the list of submatches
                // maybe need our own submatch function with mapnums!
                ROMOL_SPTR qmol(new ROMol(*mol2));
                ROMOL_SPTR rmol(new ROMol(*mol1));
                unsigned int transferred = 0;
                bool keep = true;

                for ( size_t j = 0 ; j < matches[i].size() && keep; ++j ) {

                    unsigned int ridx = matches[i][j].second;
                    unsigned int qidx = matches[i][j].first;
                    Atom *rAt = rmol->getAtomWithIdx(ridx);
                    Atom *qAt = qmol->getAtomWithIdx(qidx);
                    unsigned int qAtnum = qAt->getAtomMapNum();
                    unsigned int rAtnum = rAt->getAtomMapNum();
                    if (qAtnum > 0 && rAtnum > 0) {
                        if (qAtnum == rAtnum) {
                            transferred += 1;
                        } else {
                            keep = false;
                        }
                    } else {
                        if (rAtnum == 0 ) {
                            rAt->setAtomMapNum(qAtnum);
                            transferred += 1;
                        } else if (qAtnum == 0 ) {
                            qAt->setAtomMapNum(rAtnum);
                            transferred += 1;
                        }
                    }
                }

                if (debug) {
                     std::cout << "Compute solution:" << std::endl;
                     std::cout << "--------------------------------------"<< std::endl;
                     std::cout << "transferred: "<< transferred << "/" << matches[i].size() <<  std::endl;
                     std::cout << "mol1:" << MolToSmiles(*rmol) << std::endl;
                     std::cout << "mol2:" << MolToSmiles(*qmol) << std::endl;
                     std::cout << "**************************************"<< std::endl;
                }

                if (keep && transferred == matches[i].size()) {
                    if (!isValidMapping(qmol, btab)) {
                        continue;
                    }
                    if (!isValidMapping(rmol, btab)) {
                        continue;
                    }


                    valid_hit.push_back(matchPair(matches[i], transferred));

                }
            }
        }

        void match(ROMOL_SPTR mol1, ROMOL_SPTR mol2, std::vector< mapBond > btab,
                     std::vector<RDKit::MatchVectType> &matches, bool debug  ) {
            // need to sort first largest to smallest molecule

            std::vector<matchPair> valid_hit;
            getComputeSolutions( mol1, mol2,  btab, matches, valid_hit, debug );

            if (valid_hit.size() > 0) {

                if (debug) {
                  std::cout << "got at least one solution of matches" << std::endl;
                }

                // get best match
                std::sort(valid_hit.begin(), valid_hit.end());

                // get the first (ie should be the best mapping solution!)
                RDKit::MatchVectType winner = valid_hit.front().mvec;

                // need to take best mapping match in the list of submatches
                // maybe need our own submatch function with mapnums!
                for ( size_t j = 0 ; j < winner.size() ; ++j ) {
                    unsigned int ridx = winner[j].second;
                    unsigned int qidx = winner[j].first;
                    Atom *rAt = mol1->getAtomWithIdx(ridx);
                    Atom *qAt = mol2->getAtomWithIdx(qidx);
                    unsigned int qAtnum = qAt->getAtomMapNum();
                    unsigned int rAtnum = rAt->getAtomMapNum();
                    if (rAtnum == 0 ) {
                        rAt->setAtomMapNum(qAtnum);
                    } else if (qAtnum == 0 ) {
                        qAt->setAtomMapNum(rAtnum);
                    }
                }
            }
        }

        void transferMapping(std::vector< mapBond > btab, std::vector<ROMOL_SPTR> &mollist, bool debug) {
            // need to sort first largest to smallest molecule
            std::sort(mollist.begin(), mollist.end(), comparedMolSized);
            std::vector<ROMOL_SPTR>::iterator it1, it2;
            for (it1 = mollist.begin(); it1 != mollist.end(); ++it1) {
                ROMOL_SPTR mol1 = (*it1);

                for (it2 = mollist.begin(); it2 != mollist.end(); ++it2) {
                    if (it2 > it1) {
                        ROMOL_SPTR mol2 = (*it2);
                          if (debug) {
                            std::cout << "=====  is mathching ======" << std::endl;
                            std::cout << RDKit::MolToSmiles(*mol1) << std::endl;
                            std::cout << RDKit::MolToSmiles(*mol2) << std::endl;
 
                          }
                        std::vector<RDKit::MatchVectType> matches;
                        // bool uniquify = false; // to get all possible rotations symetric
                        if ( RDKit::SubstructMatch( *mol1, *mol2, matches, false ) ) {
                          if (debug) {
                            std::cout << "Yes" << std::endl;
                            
                            std::cout << "----------------------" << std::endl;
 
                          }
                            match( mol1,  mol2, btab, matches, debug );
                        } else
                        {
                          if (debug) {
                            std::cout << "No" << std::endl;

                            std::cout << "----------------------" << std::endl;


                          }
                        }
                    }
                }
            }
        }


        // check only for fully map cases!
        void getActiveReagants(std::vector<ROMOL_SPTR> allmolR, std::vector<ROMOL_SPTR> allmolP, std::vector<ROMOL_SPTR> &goodmolR ) {

            std::vector<int> mapatomsP;
            std::vector<ROMOL_SPTR>::iterator itp;
            for (itp = allmolP.begin(); itp != allmolP.end(); ++itp) {
                for (auto atom : (*itp)->atoms()) {
                    mapatomsP.push_back(atom->getAtomMapNum());
                }
            }

            std::vector<ROMOL_SPTR>::iterator it;
            for (it = allmolR.begin(); it != allmolR.end(); ++it) {
                bool keep = false;
                for (auto atom : (*it)->atoms()) {
                    int atnum =  atom->getAtomMapNum();
                    if (!keep & atnum > 0) {
                        keep = std::find(mapatomsP.begin(), mapatomsP.end(), atnum) != mapatomsP.end();
                    }
                }
                if (keep) {
                    goodmolR.push_back(*it);
                }
            }

        }

        void getAcceptMols(std::vector<ROMOL_SPTR> allmol, std::vector<ROMOL_SPTR> &goodmol ) {

            std::vector<ROMOL_SPTR>::iterator it;
            for (it = allmol.begin(); it != allmol.end(); ++it) {
                bool keep = false;
                for (auto atom : (*it)->atoms()) {
                    if (!keep) {
                        keep = atom->getAtomMapNum() > 0;
                    }
                }

                if (keep) {
                    goodmol.push_back(*it);
                }
            }
        }

        void addLG(std::vector<bondid > conserved,  std::vector<ROMOL_SPTR> rst, std::vector<ROMOL_SPTR> &dst,  bool debug) {
            if (conserved.size()>0) {
                // get list of atommapnum from products molecules
                std::vector<ROMOL_SPTR>::iterator itmol;
                std::vector< unsigned int > dstnums;
                for (itmol = dst.begin(); itmol != dst.end(); ++itmol) {
                    std::vector< unsigned int > molset = (*itmol)->getMappedAtomSet();
                    dstnums.insert(dstnums.end(), molset.begin(), molset.end());
                }
                // sort them for comparison
                std::sort(dstnums.begin(), dstnums.end());
                
                
                if (debug) {
                    std::cout << "set list of mapnum: "  << std::endl;
                }
                
                // loop over reagants and check if any atommapnum is not present in product list
                std::vector<ROMOL_SPTR>::iterator rit;
                for (rit = rst.begin(); rit != rst.end(); ++rit) {
                    ROMOL_SPTR mol = *rit;

                    std::shared_ptr<RWMol> cpMol( new RWMol( *mol ) );
                    bool changed = false;

                    for(auto bt : (*rit)->bonds()){
                        int fnum = bt->getBeginAtom()->getAtomMapNum();
                        int tnum = bt->getEndAtom()->getAtomMapNum();

                        if (fnum > 0 && tnum > 0) {
                          unsigned int fidx = bt->getBeginAtom()->getIdx();
                          unsigned int tidx = bt->getEndAtom()->getIdx();

                           bondid query(fnum, tnum, fidx, tidx, bt->getBondTypeWithAtoms(false));
                           if (std::find(conserved.begin(), conserved.end(), query) != conserved.end()) {
                                if (debug) {
                                    std::cout << "Conserved: " <<  fnum << ", " << tnum  << std::endl;
                                }
                                Bond::BondType BT = bt->getBondType();

                                if (BT != Bond::BondType::AROMATIC) {
                                    // remove the bond
                                    unsigned int addH = getNumHs(BT);
                                    cpMol->removeBond(fidx, tidx);
                                    
                                    if (debug) {
                                        std::cout << "I cut: " << fidx <<  ", " << tidx  << ": ie " << fnum << ", "  << tnum << std::endl;
                                    }
                                    
                                    changed = true;
                                    // add Hs
                                    unsigned int newidx1, newidx2;
                                    for (unsigned int i = 0; i < addH; i++) {
                                        newidx1 = cpMol->addAtom();
                                        Atom *newatom1 = cpMol->getAtomWithIdx(newidx1);
                                        newatom1->setAtomicNum(1);
                                        cpMol->addBond(fidx, newidx1, Bond::BondType::SINGLE);

                                        newidx2 = cpMol->addAtom();
                                        Atom *newatom2 = cpMol->getAtomWithIdx(newidx2);
                                        newatom2->setAtomicNum(1);
                                        cpMol->addBond(tidx, newidx2, Bond::BondType::SINGLE);
                                    }
                                }
                            }
                        }
                    }

                    // if changed done please add new Fragments to product molecules list if Enclosed is valid
                    if (changed) {
                         // add new mols to dist list
                         if (debug) {
                             std::cout << "before FragMol list"  << std::endl;
                         }
                         std::vector<ROMOL_SPTR> listofmol = MolOps::getMolFrags(*cpMol, false); 
                         if (debug) {
                             std::cout << "after FragMol list"  << std::endl;
                         }
                         for (ROMOL_SPTR fragmol : listofmol) {                  
                            ROMOL_SPTR molok(MolOps::removeHs(*fragmol, false, true, false));

                            if (!isEnclosed(molok, dst, debug)) {
                                dst.push_back(molok);
                                if (debug) {
                                    std::cout << MolToSmiles(*molok)  << std::endl;
                                    std::cout << " => mol cut added "  << std::endl;
                                }
                            }
                         }
                    }
                }
            }
        }
         

        std::string RXNCompleteMapping(std::string sma, bool debug, bool addleavinggroups) {
            std::vector<std::string> results;
            boost::split(results, sma, boost::is_any_of(">"));
            sma = results[0] + ">>" + results[2] ;
            if (debug) {
                // check data process
                std::cout << "1:" << sma << std::endl;
            }

            // use useSmile = true and aromatize/sanitize = true (parameter expose in python!)
            ChemicalReaction *rxn(RxnSmartsToChemicalReaction( sma, nullptr, true, true));
            if (debug) {
                // check data process
                std::cout << "2: reaction readed." << std::endl;
            }

            std::vector<ROMOL_SPTR> rkeep0, rkeep, pkeep;
            getAcceptMols(rxn->getReactants(), rkeep0);
            getAcceptMols(rxn->getProducts(), pkeep);


            if (debug) {
                // check data process
                std::cout << "3: got R&P" << std::endl;
            }

            // extra steps in c++ to keep only active reagants and remove monoAtoms reagants
            getActiveReagants(rkeep0, pkeep, rkeep);

            removeMonoAtomreagents(rkeep, pkeep);


            // atom changes
            std::map<unsigned int, AtomType> ratoms, patoms;

            bool reactantnondupatnum;
            bool productnondupatnum;


            getAtomInventory(pkeep, patoms, productnondupatnum, false);
            getAtomInventory(rkeep, ratoms, reactantnondupatnum, true);

            if (reactantnondupatnum || productnondupatnum) {
                if (debug) {
                    std::cout << "DUPLICATE ATNUM IN THE REACTION"  << std::endl;
                }
                return "DuplicateATNUM";

            }

            if (debug) {
                // check data process
                std::cout << "4: got Atom inv!"  << std::endl;
            }

            std::map<ROMOL_SPTR, std::vector<AtomType> > rdiff, pdiff;
            getDiffAtoms(ratoms, patoms, rdiff, pdiff);
            if (debug) {
                // check data process
                std::cout << "5: got Atom Diff!" << std::endl;
            }

            // bonds conserved no need to break them
            std::vector<bondid> rbonds, pbonds;
            getBondInventory(rkeep, rbonds );
            getBondInventory(pkeep, pbonds );
            if (debug) {
                // check data process
                std::cout << "6: got Bond Inv!" << std::endl;
            }

            std::vector<bondid> conserved;
            std::vector<bondid> changed;

            getConservedBonds(rbonds, pbonds, conserved , changed);

            if (debug) {
                // check data process
                std::cout << "7: got bond Conserved!"  << std::endl;
            }

            // step 1: Do reverse
            unsigned int maxmap = 0;
            getMaxMapNumList(rkeep, maxmap);
            unsigned int newmaxmap = maxmap;
            mapAllAtoms(pkeep, newmaxmap);

            // getMaxMapNumList(pkeep, maxmap); // this  line was commented why ?
            if (debug) {
                // check data process
                std::cout << "8: max map " << maxmap << std::endl;
            }
            if (debug) {
                // check data process
                std::cout << "9: add maxmap to products " << maxmap << std::endl;
            }

            // new order in  python code inverse step 10 & 11


            addGroup(pdiff, conserved, rkeep, maxmap, debug);

            if (debug) {
                // check data process
                std::cout << "11: groups added." << std::endl;
                std::cout <<  convertMolstoSmiles( rkeep, false) <<  ">>" << convertMolstoSmiles( pkeep, false) << std::endl;
                std::cout <<  "-----------------" << std::endl;
            }

            std::vector< mapBond > btab;
            getBtab(pkeep, btab);
            getBtab(rkeep, btab);
            if (debug) {
                // check data process
                std::cout << "10: got BTAB" << std::endl;
                std::cout <<  convertMolstoSmiles( rkeep, false) <<  ">>" << convertMolstoSmiles( pkeep, false) << std::endl;
                std::cout <<  "-----------------" << std::endl;
            }


            transferMapping(btab, rkeep, debug);
            if (debug) {
                // check data process
                std::cout << "12: got transfered!" << std::endl;
                std::cout <<  convertMolstoSmiles( rkeep, false) <<  ">>" << convertMolstoSmiles( pkeep, false) << std::endl;
                std::cout <<  "-----------------" << std::endl;
            }

            // we don't nned to do the clean up on reagant cause it may be partially unmapped
            std::vector<ROMOL_SPTR> rclean;
            finalCleanUp(rkeep,  rclean, debug);
            // return the string result
            if (debug) {
                // check data process
                std::cout << "13 : reversee results" << std::endl;
                std::string ps1 = convertMolstoSmiles( pkeep, false);
                std::string rs1 = convertMolstoSmiles( rclean, false); // replace rclean by rkeep
                std::cout <<  rs1 <<  ">>" << ps1 << std::endl;
                std::cout <<  "-----------------" << std::endl;
            }


            // step 2: do forward
            unsigned int newmaxmap2 = newmaxmap;
            mapAllAtoms(rclean, newmaxmap2); // replace rclean by rkeep
            if (debug) {
                // check data process
                std::cout << "14: reageant maxmap " << maxmap << std::endl;
            }

            btab.clear();
            getBtab(rclean, btab);
            getBtab(pkeep, btab); // replace rclean by rkeep
            if (debug) {
                // check data process
                std::cout << "15: update BTAT." << std::endl;
            }

            getAtomInventory(pkeep, patoms, productnondupatnum, false);
            getAtomInventory(rclean, ratoms, reactantnondupatnum, true); // replace rclean by rkeep

            if (reactantnondupatnum || productnondupatnum) {
                if (debug) {
                    std::cout << "DUPLICATE ATNUM IN THE REACTION"  << std::endl;
                }
                return "DuplicateATNUM";
            }


            getDiffAtoms(ratoms, patoms, rdiff, pdiff);

            // extra 2 steps introduce in the python code recently
            getBondInventory(rclean, rbonds );
            getBondInventory(pkeep, pbonds );
            if (debug) {
                // check data process
                std::cout << "6: got Bond Inv!" << std::endl;
            }

            getConservedBonds(rbonds, pbonds, conserved, changed );
            if (debug) {
                std::cout << "before addgroup" << std::endl;
            }
             
            addGroup(rdiff, conserved, pkeep, newmaxmap, debug);
            if (debug) {
                // check data process
                std::cout << "16: groups added." << std::endl;
            }
            transferMapping(btab, pkeep, debug);
            if (debug) {
                // check data process

                std::cout << "17: finish transfer!" << std::endl;
                std::string ps2 = convertMolstoSmiles( pkeep, false);
                std::string rs2 = convertMolstoSmiles( rclean, false); // relace rclean by rkeep
                std::cout <<  rs2 <<  ">>" << ps2 << std::endl;
                std::cout <<  "-----------------" << std::endl;
            }

            std::vector<ROMOL_SPTR> pclean;
            finalCleanUp(pkeep,  pclean, debug);
            if (debug) {
                // check data process
                std::cout << "18: final cleanup duplicates." << std::endl;
            }

            // there ar 3 steps missing in this code that are part of the python code
            // for LG to product before mapping
            // new method setup added from Ruud python code
            std::vector<ROMOL_SPTR> pcleanlg;

            if (addleavinggroups) {
                if (debug) {
                    std::cout << "extra leaving group : 1" << std::endl;
                }

                getBondInventory(rclean, rbonds );
                getBondInventory(pclean, pbonds );
                
                getConservedBonds(rbonds, pbonds, conserved, changed );
                if (debug) {
                    std::cout << "extra leaving group : 2" << std::endl;
                }
                addLG( changed, rclean, pclean, debug);
                if (debug) {
                    std::cout << "extra leaving group : 3" << std::endl;
                }
                finalCleanUp(pclean,  pcleanlg, debug);
                if (debug) {
                    std::cout << "extra leaving group : Done" << std::endl;
                }
            }

            std::string ps = "";
            // return the string result
            if (addleavinggroups) {
                ps = convertMolstoSmiles( pcleanlg, false);
            } 
            else {
                ps = convertMolstoSmiles( pclean, false);
            }

            std::string rs = convertMolstoSmiles( rclean, false); // replacee rclean by rkeep
            sma = rs + ">>" + ps;
            if (debug) {
                // check data process
                std::cout << "19: result" << std::endl;
                std::cout <<  rs <<  ">>" << ps << std::endl;
                std::cout <<  "-----------------" << std::endl;

            }

            rkeep.clear();
            pkeep.clear();
            pclean.clear();
            rclean.clear();
            pcleanlg.clear();
            btab.clear();
            rdiff.clear();
            pdiff.clear();
            ratoms.clear();
            patoms.clear();
            rbonds.clear();
            pbonds.clear();
            conserved.clear();
            changed.clear();
            rs.clear();
            ps.clear();
            if (debug) {
                // check data process
                std::cout << "20: clearing variables" << std::endl;

            }
            delete rxn;
            if (debug) {
                // check data process
                std::cout << "21: delete rxn object" << std::endl;
            }


            return sma;
        }

        

      
        
        bool getRXNCompTotal(std::string rxnsma, std::string rxncoresma) {
            
            std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(rxnsma));
            std::unique_ptr<ChemicalReaction> core(RxnSmartsToChemicalReaction(rxncoresma));
            
            std::vector<ROMOL_SPTR>  reactants = rxn->getReactants();
            std::vector<ROMOL_SPTR>  products = rxn->getProducts();
            
            
            // get the largest product first!
            std::sort(products.begin(), products.end(), comparedMolSized);
            
            ROMOL_SPTR mainprod = products[0];
            for(auto at : (mainprod)->atoms()){
                at->setAtomMapNum(0);
            }
            
            // get the principal product only ie largest number of atoms in the product
            const std::string mainproductsmi = RDKit::MolToSmiles(*mainprod);
            
            // resets the atom map of reactants
            for (auto it = reactants.begin(); it != reactants.end(); ++it){
                for(auto at : (*it)->atoms()){
                    at->setAtomMapNum(0);
                }
            }
            core->initReactantMatchers();

            std::vector<MOL_SPTR_VECT> coreproducts = core->runReactants(reactants);

            if (coreproducts.size()==0) {
                // inverse the order of the reaction
                std::reverse(reactants.begin(),reactants.end());
                coreproducts = core->runReactants(reactants);
            }
            
            for (auto &pset: coreproducts){
                std::sort(pset.begin(), pset.end(), comparedMolSized);
                std::string coreprodsmi = RDKit::MolToSmiles(*pset[0]);

                if (coreprodsmi == mainproductsmi) {

                        return true;
                }
            }

            return false;
            
        }
      
        std::string setRXNCompAtomMaps(std::string rxnsma, std::string rxncoresma) {
            
            std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(rxnsma));
            std::unique_ptr<ChemicalReaction> core(RxnSmartsToChemicalReaction(rxncoresma));
            
            std::vector<ROMOL_SPTR>  reactants = rxn->getReactants();
            std::vector<ROMOL_SPTR>  products = rxn->getProducts();
            
            // get the largest product first!
            std::sort(products.begin(), products.end(), comparedMolSized);
            
            ROMOL_SPTR mainprod = products[0];
            for(auto at : (mainprod)->atoms()){
                at->setAtomMapNum(0);
            }
            
            // get the principal product only ie largest number of atoms in the product
            const std::string mainproductsmi = RDKit::MolToSmiles(*mainprod);
            
            // resets the atom map of reactants
            for (auto it = reactants.begin(); it != reactants.end(); ++it){
                for(auto at : (*it)->atoms()){
                    at->setAtomMapNum(0);
                }
            }
            core->initReactantMatchers();
            
            std::vector<MOL_SPTR_VECT> coreproducts = core->runReactants(reactants);
            
            if (coreproducts.size()==0) {
                // inverse the order of the reaction
                std::reverse(reactants.begin(),reactants.end());
                coreproducts = core->runReactants(reactants);
            }
            
            MOL_SPTR_VECT coreSmiProducts;
            bool found = false;

            for (auto &pset: coreproducts){
                std::sort(pset.begin(), pset.end(), comparedMolSized);
                std::string coreprodsmi = RDKit::MolToSmiles(*pset[0]);
                
                if (coreprodsmi == mainproductsmi) {
                    coreSmiProducts = pset;
                    found = true;
                    break;
                }
            }
            
            std::vector<std::string> newProds;

            if (!found) {
                
                return "";
                
            }
            else {
                
                std::map<unsigned int, unsigned int> reactMatToReactant;
                unsigned int i = 0;
                for (const auto &react: core->getReactants()) {
                    for(auto at : (*react).atoms()){
                        unsigned int atnum = at->getAtomMapNum();
                        reactMatToReactant[atnum] = i;
                    }
                    i++;
                }
                
               for (auto &prod: coreSmiProducts){
                   for(auto at : (*prod).atoms()){
                       if (at->hasProp(common_properties::reactionMapNum)) {

                           int mapnum = at->getProp<int>(common_properties::reactionMapNum);

                           if (mapnum>0) {
                               at->setAtomMapNum(mapnum);
                               unsigned int reactIdx = at->getProp<unsigned int>(common_properties::reactantAtomIdx);

                               reactants[reactMatToReactant[mapnum]]->getAtomWithIdx(reactIdx)->setAtomMapNum(mapnum);
                               
                           }
                       }
                   }
                   std:: string prodsmi = RDKit::MolToSmiles(*prod);

                   newProds.push_back(prodsmi);
               }
                
            }
            
            std::stringstream ss;
            ss << convertMolstoSmiles(reactants);
            ss << ">>" ;
            ss << VectorSmilestoSmiles(newProds);
            
            return ss.str();
            
            
        }
        
        
        bool getRXNComparison(std::string rxnsma, std::string rxncoresma) {
            

            

            
            
            
            std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(rxnsma));
            std::unique_ptr<ChemicalReaction> core(RxnSmartsToChemicalReaction(rxncoresma));
            

            
            bool submatch = hasReactionSubstructMatch(*rxn, *core);

            return submatch;
            /*

            std::vector<ROMOL_SPTR> rkeep, rcorekeep, pkeep, pcorekeep;
            getAcceptMols(rxn->getReactants(), rkeep);
            getAcceptMols(rxn->getProducts(), pkeep);
            getAcceptMols(core->getReactants(), rcorekeep);
            getAcceptMols(core->getProducts(), pcorekeep);
            
            clearAtomMap(rkeep);
            clearAtomMap(pkeep);
            clearAtomMap(rcorekeep);
            clearAtomMap(pcorekeep);
            
            core->initReactantMatchers();
            
            MOL_SPTR_VECT reag;
            std::vector<ROMOL_SPTR>::iterator itmol;
            for (itmol = rkeep.begin(); itmol != rkeep.end(); ++itmol) {
                reag.push_back(ROMOL_SPTR(*itmol));
            }
            
            
            //std::vector<ROMOL_SPTR> coreproducts;
            
            auto coreproducts = core->runReactants(reag);

   /*
            int size = rkeep.size();
            std::vector<std::vector<int> > mypermutelist;
            for (int idx=0;idx<size;idx++) {
                std::vector<int> seq;
                permutelist(idx,size,seq, mypermutelist);
            }
            
            for (int iter; (iter < size && coreproducts.size()==0) ; iter++) {
                //apply_permutation( reag,mypermutelist[iter]);
                coreproducts = core->runReactants(reag);
            }

            
            std::cout << "found exact:" << exact << std::endl;
            bool   withoutAtomMap = false;

            if (coreproducts.size() == 0) {
                std::cout << "no reaction found:" << withoutAtomMap << std::endl;
                return withoutAtomMap;
            } else {
                // test if the largest/main product for is present in the initial reaction
                // currently TODO!!!!
                withoutAtomMap = true;
                std::cout << "reaction found:" << withoutAtomMap << std::endl;
                return withoutAtomMap;
            }
                
            */

            
            
            
        }
    // end CondensedGraphRxn

        
    }
    // end RDkit
}
