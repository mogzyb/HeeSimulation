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
  class MyHAnalysis : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MyHAnalysis);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      FinalState fs;
      
      
      declare(Beam(), "Beams");
      ZFinder zf(fs,Cuts::open(),11,91.2-5,91.2+5);
      declare(zf, "ZF");
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(zf);
      declare(jet_fs, "FS");
      declare(ChargedFinalState(jet_fs), "CFS");
      declare(Thrust(jet_fs), "Thrust");
      declare(Sphericity(jet_fs), "Sphericity");
      declare(ParisiTensor(jet_fs), "Parisi");
      declare(Hemispheres(Thrust(jet_fs)), "Hemispheres");
  
      // Histograms   
      book(_histH_massAdd,    "H_massAdd",    100,    0, 240);
      book(_histH_massSub,    "H_massSub",    100,    0, 240);
      book(_histH_massDif,    "H_massDif",    100,    0, 240);
      
      book(_histH_T,    "H_T",    100,    0.5, 1);
      book(_histH_TGPeak,    "H_TGPeak",    100,    0.9, 1);
      book(_histH_TM,    "H_TM",    100,    0, 1);
      book(_histH_Tm,    "H_Tm",    100,    0, 1);
      book(_histH_O,    "H_O",    100,    0, 1);
      book(_histH_S,    "H_S",    100,    0, 1);
      book(_histH_A,    "H_A",    100,    0, 1);
      //book(_histH_C,    "H_C",    100,    0, 1);
      //book(_histH_BMax,    "H_BMax",    100,    0, 1);
      //book(_histH_BMin,    "H_BMin",    100,    0, 1);
      //book(_histH_BSum,    "H_BSum",    100,    0, 1);
      //book(_histH_BDiff,    "H_BDiff",    100,    0, 1);
      
      book(_histH_LTP,    "H_LTP",    100,    -0.001, 0.001);
      
    
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Charged + neutral final state
      const ZFinder& zf     = apply<ZFinder>(event, "ZF");
      const FinalState& fs  = apply<FinalState>(event, "FS");
      if(zf.bosons().size()!=1) vetoEvent;
      Particles allp = fs.particlesByPt();
      const FourMomentum& zmom = zf.bosons()[0].momentum();
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();

      // Here I initialise the four-momentum.
      FourMomentum hmom(0.,0.,0.,0.);
      for (const Particle& p : fs.particles()) {
			// adding the four-momentum of the particle.
			hmom += p.momentum();	  
      } 
 
      _histH_massAdd->fill(hmom.mass()/GeV);
      _histH_massSub->fill((beams.first.p3().mod() + beams.second.p3().mod()-(zf.constituents()[0].momentum().p())-(zf.constituents()[1].momentum().p()))/GeV);
      _histH_massDif->fill(((beams.first.p3().mod() + beams.second.p3().mod()-(zf.constituents()[0].momentum().p())-(zf.constituents()[1].momentum().p()))/hmom.mass())/GeV);

      // Lorentz boost to CM frame
      FourMomentum CM = zmom+hmom;    
      
      LorentzTransform LT = LorentzTransform::mkFrameTransform(hmom);
		Particles higgsDecayProducts;
		for(Particle p: fs.particles()) {
			p.transformBy(LT);
			higgsDecayProducts.push_back(p);
		}

		FourMomentum LThmom(0.,0.,0.,0.);
      for (const Particle& p : higgsDecayProducts) {
			LThmom += p.momentum();
		}

		_histH_LTP->fill(sqrt(sqr(LThmom.p())));

		Thrust thrust;
		thrust.calc(higgsDecayProducts);
		_histH_T->fill(thrust.thrust());
		_histH_TGPeak->fill(thrust.thrust());
		_histH_TM->fill(thrust.thrustMajor());
		_histH_Tm->fill(thrust.thrustMinor());
		_histH_O->fill(thrust.oblateness());
      
      Sphericity sphericity;
      sphericity.calc(higgsDecayProducts);
      _histH_S->fill(sphericity.sphericity());
      _histH_A->fill(sphericity.aplanarity());
      
      const size_t numWeights = event.weights().size();
      for (size_t m = 0; m < numWeights; ++m) {
        const double weight = event.weights()[m];
      }
      
		//const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");      
      
      //ParisiTensor parisi(higgsDecayProducts);
      //parisi.calc(higgsDecayProducts);
      //_histH_C->fill(parisi.C());
      
		//const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");      
      
		//Hemispheres hemi(thrust);
		//for (const Particle& p : higgsDecayProducts) {
		//hemi.calc(p);
		//_histH_BMax->fill(hemi.Bmax());
      //_histH_BMin->fill(hemi.Bmin());
      //_histH_BSum->fill(hemi.Bsum());
      //_histH_BDiff->fill(hemi.Bdiff());
      //}
    }


    /// Finalize
    void finalize() {
		double NormVal=crossSection()*5000/femtobarn/sumOfWeights();

      normalize(_histH_massAdd,NormVal);
      normalize(_histH_massSub,NormVal);
      normalize(_histH_massDif,NormVal);
      normalize(_histH_T,NormVal);
      //scale(_histH_T, br);
      normalize(_histH_TGPeak,NormVal);
      //scale(_histH_TGPeak, br);
      normalize(_histH_TM,NormVal);
      normalize(_histH_Tm,NormVal);
      normalize(_histH_O,NormVal);  
      normalize(_histH_S,NormVal);   
      normalize(_histH_A,NormVal);
      //normalize(_histH_C);  
      //normalize(_histH_BMax);  
		//normalize(_histH_BMin);
		//normalize(_histH_BSum);
		//normalize(_histH_BDiff);
		normalize(_histH_LTP,NormVal);  
		  

      
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histH_T, _histH_TGPeak, _histH_TM, _histH_Tm, _histH_O, _histH_S, _histH_A;
    //Histo1DPtr _histH_BMax, _histH_BMin, _histH_BSum, _histH_BDiff;
    //Histo1DPtr _histH_C;
    Histo1DPtr _histH_massAdd, _histH_massSub, _histH_massDif;
    Histo1DPtr _histH_LTP;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MyHAnalysis);

}

