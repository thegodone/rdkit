#include <GraphMol/ChemReactions/Reaction.h>
#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/CondensedGraphRxn/FirCRS.h>
#include <set>

namespace python = boost::python;
namespace fir = RDKit::Firmenich::Reactions;

namespace RDKit {
namespace Firmenich {
namespace Reactions {
namespace Aux {
    
    python::dict crs2smadict(CRSR_SPTR reader, const python::list lst) {
        python::dict dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &smi){ 
            dst[smi] = reader->read(smi);
        });
        return dst;
    }
    
    python::list crs2smalist(CRSR_SPTR reader, const python::list lst) {
        python::list dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &smi){ 
            dst.append(reader->read(smi));
        });
        return dst;
    }
    
    python::dict sma2crsdict(CRSW_SPTR writer, const python::list lst) {
        python::dict dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &sma){ 
            dst[sma] = writer->write(sma);
        });
        return dst;
    }
    
    python::list sma2crslist(CRSW_SPTR writer, const python::list lst) {
        python::list dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &sma){ 
            dst.append(writer->write(sma));
        });
        return dst;
    }
    
}
}
}
}


BOOST_PYTHON_MODULE(rdFirCRS) {
    
    //! Define the CRSWriter
    python::class_<fir::CRSWriter>("CRSWriter", python::init<>())
        .def("write", &fir::CRSWriter::write, 
             (python::arg("rxnsma"),python::arg("rnd")=false,python::arg("root")=-1), 
             "Method writes a CRS-type SMILES")
        .def("sma2crsdict", &fir::Aux::sma2crsdict, (python::arg("self"),python::arg("sma")),
             "Method converts a reaction SMARTS to a dictionary[sma] = crs.")
        .def("sma2crslist", &fir::Aux::sma2crslist, (python::arg("self"),python::arg("sma")),
            "Method converts a list of reaction SMARTS to CRSs.");
    python::register_ptr_to_python<fir::CRSW_SPTR>();
    
    //! Define the CRSReader
    python::class_<fir::CRSReader>("CRSReader", python::init<>())
        .def("read", &fir::CRSReader::read, 
             (python::arg("crssmi"),python::arg("mapcrs")=true,python::arg("mapall")=false, python::arg("can")=true),
            "Method reads a CRS-smiles and converts to reaction SMARTS, mappging the CRS-atoms on request.")
        .def("crs2smadict", &fir::Aux::crs2smadict, (python::arg("self"),python::arg("crss")),
            "Method converts the list of crss to a dictionary[crs] = sma.")
        .def("crs2smalist", &fir::Aux::crs2smalist, (python::arg("self"),python::arg("crss")),
             "Method converts the list of crss to a list of smarts.")
        .def("siteSmarts", &fir::CRSReader::siteSmarts, 
             (python::arg("smi"),python::arg("plusone")=true),
             "Method creates the site SMARTS. By default with the next unchanged atom.\
              This can be optionally switched off setting the option 'plusone=False'");
    python::register_ptr_to_python<fir::CRSR_SPTR>();
    
    //! Static method to read CRS
    python::def("crs2sma", &fir::CRS2SMA, (python::arg("crs")),
               "Method converts a CRS-SMILES to SMARTS.");
    
    //! Static method to write CRS from smarts
    python::def("sma2crs", &fir::SMA2CRS, (python::arg("sma")),
               "Method converts smarts to crs.");
    
}
