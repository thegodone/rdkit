#include <GraphMol/ChemReactions/Reaction.h>
#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/CondensedGraphRxn/FirCGR.h>
#include <set>

namespace python = boost::python;
namespace fir = RDKit::Firmenich::Reactions;

namespace RDKit {
namespace Firmenich {
namespace Reactions {
namespace Aux {
    
    python::dict cgr2smadict(CGRR_SPTR reader, const python::list lst) {
        python::dict dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &smi){ 
            dst[smi] = reader->read(smi);
        });
        return dst;
    }
    
    python::list cgr2smalist(CGRR_SPTR reader, const python::list lst) {
        python::list dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &smi){ 
            dst.append(reader->read(smi));
        });
        return dst;
    }
    
    python::dict sma2cgrdict(CGRW_SPTR writer, const python::list lst) {
        python::dict dst;
        python::stl_input_iterator<std::string> beg(lst), end;
        std::for_each(beg,end,[&](const std::string &sma){ 
            dst[sma] = writer->write(sma);
        });
        return dst;
    }
    
    python::list sma2cgrlist(CGRW_SPTR writer, const python::list lst) {
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


BOOST_PYTHON_MODULE(rdFirCGR) {
    
    //! Define the CGRWriter
    python::class_<fir::CGRWriter>("CGRWriter", python::init<>())
        .def("write", &fir::CGRWriter::write, 
             (python::arg("rxnsma"),python::arg("rnd")=false,python::arg("root")=-1), 
             "Method writes a CGR-type SMILES")
        .def("sma2cgrdict", &fir::Aux::sma2cgrdict, (python::arg("self"),python::arg("sma")),
             "Method converts a reaction SMARTS to a dictionary[sma] = cgr.")
        .def("sma2cgrlist", &fir::Aux::sma2cgrlist, (python::arg("self"),python::arg("sma")),
            "Method converts a list of reaction SMARTS to CGRs.");
    python::register_ptr_to_python<fir::CGRW_SPTR>();
    
    //! Define the CGRReader
    python::class_<fir::CGRReader>("CGRReader", python::init<>())
        .def("read", &fir::CGRReader::read, 
             (python::arg("cgrsmi"),python::arg("mapcgr")=true,python::arg("mapall")=false, python::arg("can")=true),
            "Method reads a CGR-smiles and converts to reaction SMARTS, mappging the CGR-atoms on request.")
        .def("cgr2smadict", &fir::Aux::cgr2smadict, (python::arg("self"),python::arg("cgrs")),
            "Method converts the list of cgrs to a dictionary[cgr] = sma.")
        .def("cgr2smalist", &fir::Aux::cgr2smalist, (python::arg("self"),python::arg("cgrs")),
             "Method converts the list of cgrs to a list of smarts.")
        .def("siteSmarts", &fir::CGRReader::siteSmarts, 
             (python::arg("smi"),python::arg("plusone")=true),
             "Method creates the site SMARTS. By default with the next unchanged atom.\
              This can be optionally switched off setting the option 'plusone=False'");
    python::register_ptr_to_python<fir::CGRR_SPTR>();
    
    //! Static method to read CGR
    python::def("cgr2sma", &fir::CGR2SMA, (python::arg("cgr")),
               "Method converts a CGR-SMILES to SMARTS.");
    
    //! Static method to write CGR from smarts
    python::def("sma2cgr", &fir::SMA2CGR, (python::arg("sma")),
               "Method converts smarts to cgr.");
    
}
