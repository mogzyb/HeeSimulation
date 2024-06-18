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
  class MyFSAnalysis : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MyFSAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      FinalState fs;
      VetoedFinalState jet_fs(fs);
      
      declare(Beam(), "Beams");
      ZFinder zf(fs,Cuts::open(),11,91.2-5,91.2+5);
      declare(zf, "ZF");
      jet_fs.addVetoOnThisFinalState(zf);
      declare(jet_fs, "FS");
      declare(ChargedFinalState(jet_fs), "CFS");
      declare(Thrust(jet_fs), "Thrust");
      declare(Sphericity(jet_fs), "Sphericity");
      declare(ParisiTensor(jet_fs), "Parisi");
      declare(Hemispheres(Thrust(jet_fs)), "Hemispheres");

      // Histograms   
	   book(_h_XS, "XS", {-0.5, 0.5});
	   //book(_h_XSUp, "XSUp", {-0.5, 0.5});
	   //book(_h_XSLow, "XSLow", {-0.5, 0.5});      
      book(_histMass,       "Mass",      100,  0, 300);      
      book(_histMult,       "Mult",      100,  0, 100);
      book(_histMultCh,     "MultCh",    100,  0, 100);
      book(_histPt,         "Pt",        100,  0, 150);
      book(_histPtCh,       "PtCh",      100,  0, 150);
      book(_histE,          "E",         100,  0, 150);
      book(_histECh,        "ECh",       100,  0, 150);
      book(_histEtaSumEt,   "EtaSumEt",   100,  0, 5);
      book(_histEta,        "Eta",        100, -5, 5);
      book(_histEtaCh,      "EtaCh",      100, -5, 5);
      book(_histRapidity,   "Rapidity",   100, -5, 5);
      book(_histRapidityCh, "RapidityCh", 100, -5, 5);
      book(_histPhi,        "Phi",        100,  0, TWOPI);
      book(_histPhiCh,      "PhiCh",      100,  0, TWOPI);

      book(_tmphistEtaPlus,    "TMP/EtaPlus",    100, 0, 5);
      book(_tmphistEtaMinus,   "TMP/EtaMinus",   100, 0, 5);
      book(_tmphistEtaChPlus,  "TMP/EtaChPlus",  100, 0, 5);
      book(_tmphistEtaChMinus, "TMP/EtaChMinus", 100, 0, 5);
      book(_tmphistRapPlus,    "TMP/RapPlus",    100, 0, 5);
      book(_tmphistRapMinus,   "TMP/RapMinus",   100, 0, 5);
      book(_tmphistRapChPlus,  "TMP/RapChPlus",  100, 0, 5);
      book(_tmphistRapChMinus, "TMP/RapChMinus", 100, 0, 5);
      book(_histEtaPMRatio,        "EtaPMRatio");
      book(_histEtaChPMRatio,      "EtaChPMRatio");
      book(_histRapidityPMRatio,   "RapidityPMRatio");
      book(_histRapidityChPMRatio, "RapidityChPMRatio");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    	
      const size_t numWeights = event.weights().size();
      const vector<pair<double,double>> xsecs = event.crossSections();
      for (size_t m = 0; m < numWeights; ++m) {
        #if defined RIVET_ENABLE_HEPMC_3 || defined HEPMC_HAS_CROSS_SECTION
        size_t idx = (xsecs.size() == numWeights)? m : 0;
        const double xs    = xsecs[idx].first;
        const double xserr = xsecs[idx].second;
        _h_XS.get()->_getPersistent(m)->point(0).setY(xs, xserr);
        # endif
      }    	
    	
      // Charged + neutral final state
      const ZFinder& zf     = apply<ZFinder>(event, "ZF");
      const FinalState& fs  = apply<FinalState>(event, "FS");
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if(zf.bosons().size()!=1) vetoEvent;
      Particles allp = fs.particlesByPt();
      const FourMomentum& fmom = allp[0].momentum();

      _histMass->fill(fmom.p());
      // this is the multiplicity of everything but the Z (I think)
      MSG_DEBUG("Total multiplicity = " << fs.size());
      _histMult->fill(fs.size());

      for (const Particle& p : fs.particles()) {
	// all the other things you fill
        _histEta->fill(p.eta());
        _histEtaSumEt->fill(p.abseta(), p.Et());
        _histRapidity->fill(p.rap());
        _histPt->fill(p.pT()/GeV);
        _histE->fill(p.E()/GeV);
        _histPhi->fill(p.phi());
        (p.eta() > 0 ? _tmphistEtaPlus : _tmphistEtaMinus)->fill(p.abseta());
        (p.rap() > 0 ? _tmphistRapPlus : _tmphistRapMinus)->fill(p.absrap());
	//if (zmom.pT()>0.0) _histZMom->fill(zmom.pT()/GeV);     		  
      } 

      // Same for the charged FS particles only      
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histMultCh->fill(cfs.size());
      for (const Particle& p : cfs.particles()) {
        _histEtaCh->fill(p.eta());
        _histRapidityCh->fill(p.rap());
        _histPtCh->fill(p.pT()/GeV);
        _histECh->fill(p.E()/GeV);
        _histPhiCh->fill(p.phi());
        (p.eta()>0? _tmphistEtaChPlus:_tmphistEtaChMinus)->fill(p.abseta());
        (p.rap()>0? _tmphistRapChPlus:_tmphistRapChMinus)->fill(p.absrap());
      }     
    }


    /// Finalize
    void finalize() {
      normalize(_histMult);
      normalize(_histEta);
      normalize(_histRapidity);
      normalize(_histMass);
      normalize(_histPt);
      normalize(_histE);
      normalize(_histPhi);
      normalize(_histMultCh);
      normalize(_histEtaCh);
      normalize(_histRapidityCh); 
      normalize(_histPtCh);
      normalize(_histECh);     
      
      divide(_tmphistEtaPlus, _tmphistEtaMinus, _histEtaPMRatio);
      divide(_tmphistEtaChPlus, _tmphistEtaChMinus, _histEtaChPMRatio);
      divide(_tmphistRapPlus, _tmphistRapMinus, _histRapidityPMRatio);
      divide(_tmphistRapChPlus, _tmphistRapChMinus, _histRapidityChPMRatio);
      
      #if !defined RIVET_ENABLE_HEPMC_3 && !defined HEPMC_HAS_CROSS_SECTION
      _h_XS->point(0).setY(xs, xserr);
      #endif
      
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histMult, _histEta, _histRapidity, _histPt, _histE, _histPhi;
    Histo1DPtr _histMultCh,  _histEtaCh, _histRapidityCh, _histPtCh, _histECh, _histPhiCh;
    Histo1DPtr _histMass;
    Profile1DPtr _histEtaSumEt;
    Scatter2DPtr _histEtaPMRatio, _histEtaChPMRatio, _histRapidityPMRatio, _histRapidityChPMRatio;
    Scatter2DPtr _h_XS;
    //@}

    /// @name Temporary histos used to calculate +/- rapidity ratio plots
    //@{
    Histo1DPtr _tmphistEtaPlus, _tmphistEtaMinus, _tmphistEtaChPlus, _tmphistEtaChMinus;
    Histo1DPtr _tmphistRapPlus, _tmphistRapMinus, _tmphistRapChPlus, _tmphistRapChMinus;
    //@}



  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MyFSAnalysis);

}
