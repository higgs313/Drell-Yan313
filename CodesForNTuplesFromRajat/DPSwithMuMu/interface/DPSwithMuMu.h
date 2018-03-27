#include <memory>
#include <map>
#include <cmath>
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include "TH1D.h"
#include "TH2D.h"


/////////////////
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/Common/interface/View.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

typedef std::vector<edm::InputTag> vtag;
using namespace edm;
using namespace std;
using namespace reco;

class DPSwithMuMu : public edm::EDAnalyzer {
	private:  


		HLTConfigProvider hltConfig;

		//Events

		std::vector<int> T_mu_charge;
		std::vector<float>T_mu_pfreliso,T_mu_pt,T_mu_eta,T_mu_phi,T_mu_et,T_mu_pz,T_mu_mass,T_mu_energy,T_mu_dz,T_mu_d0, T_mu_beamSpot_d0,T_mu_beamSpot_dz;
		std::vector<float>T_mu_lead_vtx_pos;
		// gen particles
		std::vector<float> T_gen_part_pt,T_gen_part_eta,T_gen_part_phi,T_gen_part_energy;
		std::vector<int> T_gen_part_pdgid,T_gen_part_status,T_gen_part_mother_pdgid;
		// @@@@@@@@@@@@@@@@@@@@@@@@Generals
		std::vector<float> weights;

		std::vector<std::string>trigger_name;
		int nvtx;  int nvertices;
		std::vector<int> T_nvertices;
		//  std::vector<bool> T_Top_tag;
		int T_Event_RunNumber,T_Event_EventNumber,T_Event_LuminosityBlock;

		//  std::vector<std::string>all_trigger_name;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%% MET of the event
		double T_MET_ET;
		double T_MET_Phi;
		double T_MET_Pt;
		//  double T_MET_Eta;


		double T_MET_enup, T_MET_endown;
		double  T_MET_jecup, T_MET_jecdown, T_MET_jesup, T_MET_jesdown, T_MET_muup, T_MET_mudown;//, T_MET_elup, T_MET_eldown;
		double T_MET_phi_enup, T_MET_phi_endown;
		double  T_MET_phi_jecup, T_MET_phi_jecdown, T_MET_phi_jesup, T_MET_phi_jesdown, T_MET_phi_muup, T_MET_phi_mudown;// T_MET_phi_elup, T_MET_phi_eldown;

		//&&&&&&&&&&&&&&&&&& JETS

		std::vector<float> T_Jet_Pt,T_Jet_Phi,T_Jet_Eta,T_Jet_Area,T_Jet_Pz;
		std::vector<float> T_gen_Jet_Pt,T_gen_Jet_Phi,T_gen_Jet_Eta,T_gen_Jet_Energy,T_gen_Jet_Pz;
		std::vector<bool> flag_loose_jetid,flag_tight_jetid;
		std::vector<bool> flag_loose_chs_jetid,flag_tight_chs_jetid;
		std::vector<int> T_Jet_parton_flavor;
		int nJets;

		// root file to store histograms
		TFile*  rootFile_;

		//Tree

		TTree* mytree_;
		double deltaRpf_;
		bool printDebug_;    
		bool isMC_;//,isdataD_;
                bool isMPI_;//,isdataD_;

		vector<int> T_gen_mu_charge;
		vector<float> T_gen_Jet_Px, T_gen_Jet_Py, T_gen_mu_energy, T_gen_mu_mass, T_gen_mu_et, T_gen_mu_pz, T_gen_mu_px, T_gen_mu_py, T_gen_part_px, T_gen_part_py,T_gen_part_pz, T_mu_Glob, T_mu_Chi2, T_mu_validfraction, validfraction, T_mu_localPos, T_mu_kinkFinder, T_mu_segcomp,T_mu_trkLayers, T_mu_px, T_mu_py, T_mu_lead_vtx_d0,T_mu_lead_vtx_dz, T_Jet_JBP,T_Jet_CSV,T_Jet_CSVMVA; 
		vector<float> T_Jet_Px,T_Jet_Py, T_Jet_Energy,T_Jet_Mass, T_Jet_rawF, T_Jet_jecUnc;
		vector<float> trg1_Pt, trg1_Eta, trg1_Phi, trg1_Energy;
		vector<float> trg2_Pt, trg2_Eta, trg2_Phi, trg2_Energy;
		vector<float> jet_pu_pt, jet_sig_pt, jet_total_pt;
		vector<float> T_genjet_match_px, T_genjet_match_py, T_genjet_match_pz, T_genjet_match_energy;
		float MUF;
		int chargedMult;
		bool fired1,fired2,fired3,fired4;
		vector<float> T_Jet_puMVA;
		std::vector<const pat::Muon*> MuonVector;
		//	reco::MuonCollection IdentifiedMuons;

		//		std::vector<const pat::Muon*> GoodMuonVector;
		std::vector<const reco::GenParticle*> genMuonVector;
		std::vector<const reco::GenParticle*> GoodgenMuonVector;


	public:
		// default constructor
		explicit DPSwithMuMu(const edm::ParameterSet& );



		// default destructor
		~DPSwithMuMu();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);  
	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;
		virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);  
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void beginEvent();
		virtual void endEvent();
		virtual float deltaPhi(float , float);
		const reco::Candidate* IsFromID (const reco::Candidate* part, int targetPDGId);
		edm::Service<TFileService> tfile; 


		std::string outputFile_; // output file
		TH1F *hEvents;
		TH1D *hEvents_MCNLO;
		float aMCatNLOweight;
		std::vector<double> T_Rho;
		std::vector<TLorentzVector>*gen_Jet_P4;
		std::vector<TLorentzVector>*gen_Muon_P4;
		std::vector<TLorentzVector>*gen_Part_P4;
		std::vector<TLorentzVector>*reco_Muon_P4;
		std::vector<TLorentzVector>*reco_Jet_P4;
		std::vector<bool>*trigDecision;
		edm::EDGetTokenT<edm::View<reco::GenJet>> theGenJetTag;
		edm::EDGetTokenT<edm::View<reco::GenParticle> > src_;
		edm::EDGetTokenT<edm::View<pat::Jet>> theJetTag;

		edm::EDGetTokenT<reco::BeamSpot> beamSpotTag;
		edm::EDGetTokenT<edm::View<pat::Muon>> theMuonTag;
		edm::EDGetTokenT<pat::METCollection> theMetTag;
		edm::EDGetTokenT<vector<Vertex>> theVtxTag;
		edm::EDGetTokenT<double> theRhoTag;
		edm::EDGetTokenT<vector<PileupSummaryInfo>> thePUTag;
		edm::EDGetTokenT<GenEventInfoProduct> theGenTag;
		edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
		edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
		edm::EDGetTokenT<LHEEventProduct> theLHEPTag;


		float  _lheHt;
		int  _lheNOutPartons;
		int  _lheNOutB;
		int  _lheNOutC;
		struct ComparePt {
			bool operator()( const pat::Muon* t1, const pat::Muon* t2 ) const {
				return t1->pt()> t2->pt();
			}
		};
		ComparePt ptComparator;

		struct ComparegenPt {
			bool operator()( const reco::GenParticle* r1, const reco::GenParticle* r2 ) const {
				return r1->pt()> r2->pt();
			}
		};
		ComparegenPt genptComparator;

};

DPSwithMuMu::DPSwithMuMu(const edm::ParameterSet& iConfig ){

	//	triggerBits_            = iConfig.getParameter<edm::InputTag>("bits");
	isMC_                   = iConfig.getParameter<bool>("isMC");
        isMPI_                   = iConfig.getParameter<bool>("isMPI");
	// debug
	printDebug_             = iConfig.getParameter<bool>("printDebug");
	outputFile_           = iConfig.getParameter<std::string>("outputFile");
	rootFile_             = TFile::Open(outputFile_.c_str(),"RECREATE");
	// edm::Service<TFileService> fs;
	hEvents = new TH1F("hEvents","hEvents",100,0,100);
	hEvents_MCNLO = new TH1D("hEvents_MCNLO","hEvents_MCNLO",2,0,2);
	theGenJetTag = consumes<edm::View<reco::GenJet>> (iConfig.getParameter<edm::InputTag>("genjetCollection"));
	src_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("src"));
	theJetTag  = consumes<edm::View<pat::Jet>> (iConfig.getParameter<edm::InputTag>("jetCollection"));
	beamSpotTag = consumes<reco::BeamSpot> (iConfig.getParameter<edm::InputTag>("beamSpot"));
	theMuonTag = consumes<edm::View<pat::Muon>> (iConfig.getParameter<edm::InputTag>("muonCollection"));
	theMetTag = consumes<pat::METCollection> (iConfig.getParameter<edm::InputTag>("metCollection"));
	theVtxTag = consumes<vector<Vertex>> (iConfig.getParameter<edm::InputTag>("vtxCollection"));
	triggerBits_ = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultsLabel"));
	theRhoTag  = consumes<double> (iConfig.getParameter<edm::InputTag>("rhoCollection"));
	thePUTag = consumes<vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("puCollection"));
	theGenTag  = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genCollection"));
	triggerObjects_  = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerSet"));
	theLHEPTag          =consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>("lhepCollection"));

}

DPSwithMuMu::~DPSwithMuMu()
{
	delete rootFile_;
} 

