/*
 * FirCRS.h defines the header file with classes
 * to convert reaction smarts to CRS-smiles and
 * vice versa.
 *
 * This new CRS-SMILES format, proposed by 
 * Firmenich/BigChem is the smallest and most
 * concise way to write reactions.
 *
 * The format is based on the SMILES string proposed
 * by Weiniger et al. SMILES are a very powerful
 * format for machine learning. Here we propose the
 * introduction of flexible bonds to indicate the
 * modifications of a reaction.
 * 
 * Apart from the good learning nature, CRSs stand
 * out for their easy reversibility, simply by putting
 * the bonds in the flexible bonds in the opposite order.
 * E.g. propyl ethionate: 
 * Hydrolysis: CCO{!-}C(=O)({-!}O)CCC
 * Esterfication: CCO{-!}C(=O)({!-}O)CCC
 *
 * (c) RUUD, Firmenich SA, 2020
 *
 */

#include <RDGeneral/export.h>
#ifndef RD_FIRCHR_H
#define RD_FIRCHR_H
#endif

#include <GraphMol/Bond.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

#include <string>
#include <memory>
#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>

namespace RDKit {
namespace Firmenich {
namespace Reactions {
    
//! Define a struct
struct TmpBond {
    int from,to;
    Bond::BondType bt;
    TmpBond() {
        this->from = -1;
        this->to = -1;
        this->bt = Bond::ZERO;
    }
    TmpBond(int from, int to, Bond::BondType bt) {
        this->from = std::min<int>(from,to);
        this->to = std::max<int>(from,to);
        this->bt = bt;
    }
    bool operator==(const TmpBond &other) const {
        return from == other.from && to == other.to;
    }
};

//! Define a class to write CRS
class CRSWriter {
private:
    std::vector<TmpBond> getBonds(RWMOL_SPTR mol,std::map<int,int> m2i);
    Bond::BondType lookup(std::vector<TmpBond> list, TmpBond query);
    Bond::BondType propose(Bond::BondType in, Bond::BondType out);
    std::vector<TmpBond> diff(std::vector<TmpBond> R, std::vector<TmpBond> P);
    RWMOL_SPTR copy(RWMOL_SPTR in);
public:
    CRSWriter();
    virtual ~CRSWriter();
    std::string write(const std::string &rxnsma, bool rnd=false, int root=-1);
};
 
//! Define a shared pointer for the writer
typedef std::shared_ptr<CRSWriter> CRSW_SPTR;
    
//! Define a class to read CRS
class CRSReader {
private:
    Bond::BondType propose(const Bond::BondType &bt, bool &arom, bool &dearom, bool prod=false);
public:
    CRSReader();
    virtual ~CRSReader();
    std::string read(const std::string &crssmiles, bool mapcrsatoms=true, bool mapallatoms=false, bool can=true);
    std::string siteSmarts(const std::string &crssmiles, bool siteplusone=true);
};
    
//! Define a shared pointer for the reader
typedef std::shared_ptr<CRSReader> CRSR_SPTR;
    
static std::string CRS2SMA(const std::string &crssmi) {
    CRSReader *r = new CRSReader();
    const std::string res = r->read(crssmi);
    delete r;
    return res;
}
    
static std::string SMA2CRS(const std::string &smarts) {
    CRSWriter *w = new CRSWriter();
    const std::string res = w->write(smarts);
    delete w;
    return res;
}
    
} /* End of namespace reactions */
} /* End of namespace Firmenich */
} /* End of namespace rdkit */
