// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

#include <algorithm>

namespace Rivet {


  /// Analysis looking at various distributions of final state particles
  class MyAnalysis : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MyAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections

      
      FinalState fs;
      
      ZFinder zf(fs,Cuts::open(),13,91.2-5,91.2+5);
      //ZFinder hf(fs,Cuts::open(),5,125.4-5,125.4+5);


      declare(zf, "ZF");
      //declare(hf, "HF");

      
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(zf);
      //jet_fs.addVetoOnThisFinalState(hf); 

     /// ZFinder zf();
      ///declare(zf, "ZF");
      ///VetoedFinalState fs(Cuts::abseta < 5 && Cuts::pT > 500*MeV);
      ///fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF"));
      declare(fs, "FS");
      declare(jet_fs, "ZFS");
      declare(ChargedFinalState(jet_fs), "CFS");

      // Histograms      /// @todo Choose E/pT ranged based on input energies... can't do anything about kin. cuts, though
      book(_histMult   , "Mult", 100, 0, 175);
      book(_histMultCh , "MultCh", 100, 0, 175);

      book(_histPt   , "Pt", 300, 0, 150);
      book(_histPtCh , "PtCh", 300, 0, 150);

      book(_histE   , "E", 100, 0, 150);
      book(_histECh , "ECh", 100, 0, 150);
     

      book(_histEtaSumEt , "EtaSumEt", 25, 0, 5);

      book(_histEta    , "Eta", 50, -5, 5);
      book(_histEtaCh  , "EtaCh", 50, -5, 5);
      book(_tmphistEtaPlus, "TMP/EtaPlus", 25, 0, 5);
      book(_tmphistEtaMinus, "TMP/EtaMinus", 25, 0, 5);
      book(_tmphistEtaChPlus, "TMP/EtaChPlus", 25, 0, 5);
      book(_tmphistEtaChMinus, "TMP/EtaChMinus", 25, 0, 5);

      book(_histRapidity    , "Rapidity", 50, -5, 5);
      book(_histRapidityCh  , "RapidityCh", 50, -5, 5);
      book(_tmphistRapPlus, "TMP/RapPlus", 25, 0, 5);
      book(_tmphistRapMinus, "TMP/RapMinus", 25, 0, 5);
      book(_tmphistRapChPlus, "TMP/RapChPlus", 25, 0, 5);
      book(_tmphistRapChMinus, "TMP/RapChMinus", 25, 0, 5);

      book(_histPhi    , "Phi", 50, 0, TWOPI);
      book(_histPhiCh  , "PhiCh", 50, 0, TWOPI);

      book(_histEtaPMRatio , "EtaPMRatio");
      book(_histEtaChPMRatio , "EtaChPMRatio");
      book(_histRapidityPMRatio , "RapidityPMRatio");
      book(_histRapidityChPMRatio , "RapidityChPMRatio");

      book(_histZMom , "ZM", 100, 0, 100);
      
      book(_histZ_mass ,"Z_mass", 50, 71.0, 111.0);
      book(_histZ_pT ,"Z_pT", logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_histZ_pT_peak ,"Z_pT_peak", 25, 40.0, 60.0);
      book(_histZ_y ,"Z_y", 40, -4.0, 4.0);
      book(_histZ_phi ,"Z_phi", 25, 0.0, TWOPI);

      
      //book(_h_H_mass ,"H_mass", 50, 119.7, 120.3);
      //book(_h_H_pT ,"H_pT", logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      //book(_h_H_pT_peak ,"H_pT_peak", 25, 0.0, 25.0);
      //book(_h_H_y ,"H_y", 40, -4, 4);
      //book(_h_H_phi ,"H_phi", 25, 0.0, TWOPI);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Charged + neutral final state
     
      const ZFinder& zf= apply<ZFinder>(event, "ZF");
      if(zf.bosons().size()!=1) vetoEvent;
     
      //const ZFinder& hf = apply<ZFinder>(event, "HF");
      //if (hf.bosons().size() != 1) vetoEvent;

      const FourMomentum& zmom = zf.bosons()[0].momentum();
      //const FourMomentum& hmom = hf.bosons()[0].momentum();
      
      const FinalState& fs = apply<FinalState>(event, "FS");
      MSG_DEBUG("Total multiplicity = " << fs.size());
      _histMult->fill(fs.size());
      for (const Particle& p : fs.particles()) {
        _histEta->fill(p.eta());
        _histEtaSumEt->fill(p.abseta(), p.Et());
        (p.eta() > 0 ? _tmphistEtaPlus : _tmphistEtaMinus)->fill(p.abseta());
        //
        _histRapidity->fill(p.rap());
        (p.rap() > 0 ? _tmphistRapPlus : _tmphistRapMinus)->fill(p.absrap());
        //
        _histPt->fill(p.pT()/GeV);
        _histE->fill(p.E()/GeV);
        _histPhi->fill(p.phi());
	
	if (zmom.pT()>0.0) _histZMom->fill(zmom.pT()/GeV);

      }

      // Same for the charged FS particles only
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histMultCh->fill(cfs.size());
      for (const Particle& p : cfs.particles()) {
        _histEtaCh->fill(p.eta());
        (p.eta() > 0 ? _tmphistEtaChPlus : _tmphistEtaChMinus)->fill(p.abseta());
        //
        _histRapidityCh->fill(p.rap());
        (p.rap() > 0 ? _tmphistRapChPlus : _tmphistRapChMinus)->fill(p.absrap());
        //
        _histPtCh->fill(p.pT()/GeV);
        _histECh->fill(p.E()/GeV);
        _histPhiCh->fill(p.phi());
      }
     
      //_h_H_mass->fill(hmom.mass()/GeV);
      //_h_H_pT->fill(hmom.pT()/GeV);
      //_h_H_pT_peak->fill(hmom.pT()/GeV);
      //_h_H_y->fill(hmom.rapidity());
      //_h_H_phi->fill(hmom.phi());

      _histZ_mass->fill(zmom.mass()/GeV);
      _histZ_pT->fill(zmom.pT()/GeV);
      _histZ_pT_peak->fill(zmom.pT()/GeV);
      _histZ_y->fill(zmom.rapidity());
      _histZ_phi->fill(zmom.phi());


    }


    /// Finalize
    void finalize() {
      normalize(_histMult); normalize(_histEta); normalize(_histRapidity); normalize(_histZMom); 
      normalize(_histZ_mass); normalize(_histZ_pT); normalize(_histZ_pT_peak); normalize(_histZ_y); normalize(_histZ_phi);
      //normalize(_h_H_mass); normalize(_h_H_pT); normalize(_h_H_pT_peak); normalize(_h_H_y); normalize(_h_H_phi);
      normalize(_histPt); normalize(_histE); normalize(_histPhi);
      normalize(_histMultCh); normalize(_histEtaCh); normalize(_histRapidityCh); 
      normalize(_histPtCh); normalize(_histECh); normalize(_histPhiCh);
      divide(_tmphistEtaPlus, _tmphistEtaMinus, _histEtaPMRatio);
      divide(_tmphistEtaChPlus, _tmphistEtaChMinus, _histEtaChPMRatio);
      divide(_tmphistRapPlus, _tmphistRapMinus, _histRapidityPMRatio);
      divide(_tmphistRapChPlus, _tmphistRapChMinus, _histRapidityChPMRatio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histMult, _histEta, _histRapidity, _histPt, _histE, _histPhi, _histZMom, _histZ_mass, _histZ_pT, _histZ_pT_peak, _histZ_y, _histZ_phi;
    Histo1DPtr _histMultCh,  _histEtaCh, _histRapidityCh, _histPtCh, _histECh, _histPhiCh; //_h_H_mass, _h_H_pT, _h_H_pT_peak, _h_H_y, _h_H_phi;
    Profile1DPtr _histEtaSumEt;
    Scatter2DPtr _histEtaPMRatio, _histEtaChPMRatio, _histRapidityPMRatio, _histRapidityChPMRatio;
    //@}

    /// @name Temporary histos used to calculate +/- rapidity ratio plots
    //@{
    Histo1DPtr _tmphistEtaPlus, _tmphistEtaMinus, _tmphistEtaChPlus, _tmphistEtaChMinus;
    Histo1DPtr _tmphistRapPlus, _tmphistRapMinus, _tmphistRapChPlus, _tmphistRapChMinus;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MyAnalysis);

}
