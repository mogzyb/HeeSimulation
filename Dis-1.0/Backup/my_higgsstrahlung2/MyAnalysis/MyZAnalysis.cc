// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Beam.hh"

#include <algorithm>

namespace Rivet {


  /// Analysis looking at various distributions of final state particles
  class MyZAnalysis : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MyZAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      FinalState fs;
      ZFinder zf(fs,Cuts::open(),13,91.2-5,91.2+5);
      declare(zf, "ZF");
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(zf);
      declare(jet_fs, "FS");
      declare(ChargedFinalState(jet_fs), "CFS");
  
      // Histograms   
      book(_histZMom,   "ZM",          100,    0, 100);
      book(_histZ_mass, "Z_mass",       100, 71.0, 111.0);
      book(_histZ_pT,   "Z_pT",
	   logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_histZ_pT_peak, "Z_pT_peak", 100, 40.0, 60.0);
      book(_histZ_y,       "Z_y",       100, -4.0, 4.0);
      book(_histZ_phi ,    "Z_phi",     100,  0.0, TWOPI);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Charged + neutral final state
      const ZFinder& zf     = apply<ZFinder>(event, "ZF");
      const FinalState& fs  = apply<FinalState>(event, "FS");
      if(zf.bosons().size()!=1) vetoEvent;
      Particles allp = fs.particlesByPt();
      const FourMomentum& zmom = zf.bosons()[0].momentum();

      // everything to do with the Z
      _histZ_mass->fill(zmom.mass()/GeV);
      _histZ_pT->fill(zmom.pT()/GeV);
      _histZ_pT_peak->fill(zmom.pT()/GeV);
      _histZ_y->fill(zmom.rapidity());
      _histZ_phi->fill(zmom.phi());

      // Here I initialise the four-momentum.
      FourMomentum hmom(0.,0.,0.,0.);
		if (zmom.pT()>0.0) _histZMom->fill(zmom.pT()/GeV);     		      
    }


    /// Finalize
    void finalize() {
      normalize(_histZMom);
      normalize(_histZ_mass);
      normalize(_histZ_pT);
      normalize(_histZ_pT_peak);
      normalize(_histZ_y);
      normalize(_histZ_phi);
    }

    //@}


  private:

    /// @name Histograms
    Histo1DPtr  _histZMom, _histZ_mass, _histZ_pT, _histZ_pT_peak, _histZ_y, _histZ_phi;
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MyZAnalysis);

}

