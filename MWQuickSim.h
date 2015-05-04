/////////////////////////////////////////////////////////////////////////////////
// Class: MWQuickSim
// 
// Author: Matthias Schott (Uni Mainz)
//
// Class for a quick simulation of the ATLAS Detector for W and Z events
// It takes as input the Truth Leptons and provides information on
// - Is event reconstruction (Taking into account muon efficiencies & acceptance)
// - kinematic cuts on lepton-PT, invariant Mass, Hadronic Recoil, M_T for the nominal
//   Z and W boson selection of the W-Mass analyses
// - Smearing of the Muon-Pt
// - Smearing of the Hadronic Recoil (parameterized approach, we smear HR_X and HR_Y)
// - Calculation of derived quantities such as transverse mass (M_T) and ETMiss
//
// Matthias Schott (Uni. Mainz)
//
// Example Implementation
//
// MWQuickSim	quickSim;
// while (getNextEvent()) {
//	...
//	TLorentzVector	truthMuon = GetTruthMuon();
//	TLorentzVector	truthNeutrino = GetTruthNeutrino();
//	...
//	if (quickSim.reconstructWEvent(truthMuon, truthNeutrino)) {
//		histMuon->Fill(quickSim.RecoMuon.Pt());
//		histETMiss->Fill(quickSim.RecoETMiss.Pt());
//		histHR->Fill(quickSim.RecoHadronicRecoil.Pt());
//		histMT->Fill(quickSim.RecoTransverseMass.Pt());
//	}
// }
//
/////////////////////////////////////////////////////////////////////////////////

#include "TRandom.h"

class MWQuickSim {
	public:
		/// Constructor
		MWQuickSim() {
			/// Muon Efficiency
			m_histMuonEfficiency = new TH1F("ATLASQUICKSIM_MuonEfficiency", "", 11, -2.5, 2.5);
			m_histMuonEfficiency->SetDirectory(0);
			m_histMuonEfficiency->SetBinContent(1, 0.95);
			m_histMuonEfficiency->SetBinContent(2, 0.95);
			m_histMuonEfficiency->SetBinContent(3, 0.90);
			m_histMuonEfficiency->SetBinContent(4, 0.98);
			m_histMuonEfficiency->SetBinContent(5, 0.98);
			m_histMuonEfficiency->SetBinContent(6, 0.80);
			m_histMuonEfficiency->SetBinContent(7, 0.98);
			m_histMuonEfficiency->SetBinContent(8, 0.98);
			m_histMuonEfficiency->SetBinContent(9, 0.90);
			m_histMuonEfficiency->SetBinContent(10, 0.95);
			m_histMuonEfficiency->SetBinContent(11, 0.95);

			/// Muon PT Resolution
			m_muonPTResolution 	= 0.0287;
			
			/// HR Bias and Resolution
			m_HRBias 		= 100;
			m_HRWidth_Pol0 		= 17230/1.41;
			m_HRWidth_Pol1 		= 0;
			m_HRWidth_Pol2 		= 0;
		}

		// Destructor
		~MWQuickSim() {
			clearRecoKinematics();
			delete m_histMuonEfficiency;
		}

		/// Function that returns true, if Z boson event is selected and fills the internal reco-Variables
		/// via Smearing Functions
		bool	reconstructZEvent(const TLorentzVector &vecMuon1, const TLorentzVector &vecMuon2) {
			/// Clear Reco Event Record
			clearRecoKinematics();
			
			/// Muon Efficiency
			if (fabs(vecMuon1.Eta())>2.5) return false;
			if (fabs(vecMuon2.Eta())>2.5) return false;

			if (gRandom->Uniform(1.0)>m_histMuonEfficiency->GetXaxis()->FindBin(vecMuon1.Eta())) return false;
			if (gRandom->Uniform(1.0)>m_histMuonEfficiency->GetXaxis()->FindBin(vecMuon2.Eta())) return false;

			/// Muon Kinematics
			RecoMuon1.SetPtEtaPhiM((1.0-gRandom->Gaus(0, m_muonPTResolution))*vecMuon1.Pt(), vecMuon1.Eta(), vecMuon1.Phi(), vecMuon1.M());
			RecoMuon2.SetPtEtaPhiM((1.0-gRandom->Gaus(0, m_muonPTResolution))*vecMuon2.Pt(), vecMuon2.Eta(), vecMuon2.Phi(), vecMuon2.M());

			/// Hadronic Recoil;
			TLorentzVector MCBoson = vecMuon1 + vecMuon2;
			double ptBGeV = MCBoson.Pt()/1000.;
			if (ptBGeV>50) ptBGeV = 50;
			double resPX 	= gRandom->Gaus(m_HRBias, m_HRWidth_Pol0 + m_HRWidth_Pol1*ptBGeV + m_HRWidth_Pol2*ptBGeV*ptBGeV);
			double resPY 	= gRandom->Gaus(m_HRBias, m_HRWidth_Pol0 + m_HRWidth_Pol1*ptBGeV + m_HRWidth_Pol2*ptBGeV*ptBGeV);
			RecoHadronicRecoil.SetPxPyPzE(resPX-MCBoson.Px(), resPY-MCBoson.Py(), 0, sqrt(pow(resPX-MCBoson.Px(),2)+pow(resPY-MCBoson.Py(),2)));

			/// Derived Quantities
			RecoETMiss.SetPxPyPzE(-(RecoMuon1.Px()+RecoMuon2.Px()-RecoHadronicRecoil.Px()), -(RecoMuon1.Py()+RecoMuon2.Py()+RecoHadronicRecoil.Py()), 0, sqrt(pow((RecoMuon1.Px()+RecoMuon2.Px()+RecoHadronicRecoil.Px()),2)+pow((RecoMuon1.Py()+RecoMuon2.Py()+RecoHadronicRecoil.Py()),2)));
			RecoTransverseMass = sqrt(2.0*RecoMuon1.Pt()*RecoETMiss.Pt()*(1.0-TMath::Cos(RecoMuon1.DeltaPhi(RecoETMiss))));

			/// Kinematic Cuts
			if (RecoMuon1.Pt()<20000) return false;
			if (RecoMuon2.Pt()<20000) return false;
			if (fabs(91000-MCBoson.M())>25000) return false;

			isReconstructedZ = true;
			return true;
		}

		/// Function that returns true, if W boson event is selected and fills the internal reco-Variables
		/// via Smearing Functions
		bool	reconstructWEvent(const TLorentzVector &vecMuon, const TLorentzVector &vecNeutrino) {
			/// Clear Reco Event Record
			clearRecoKinematics();
			
			/// Muon Efficiency
			if (fabs(vecMuon.Eta())>2.5) return false;
			if (gRandom->Uniform(1.0)>m_histMuonEfficiency->GetXaxis()->FindBin(vecMuon.Eta())) return false;

			/// Muon Kinematics
			RecoMuon.SetPtEtaPhiM((1.0-gRandom->Gaus(0, m_muonPTResolution))*vecMuon.Pt(), vecMuon.Eta(), vecMuon.Phi(), vecMuon.M());

			/// Hadronic Recoil;
			TLorentzVector MCBoson = vecMuon + vecNeutrino;
			double ptBGeV = MCBoson.Pt()/1000.;
			if (ptBGeV>50) ptBGeV = 50;
			double resPX 	= gRandom->Gaus(m_HRBias, m_HRWidth_Pol0 + m_HRWidth_Pol1*ptBGeV + m_HRWidth_Pol2*ptBGeV*ptBGeV);
			double resPY 	= gRandom->Gaus(m_HRBias, m_HRWidth_Pol0 + m_HRWidth_Pol1*ptBGeV + m_HRWidth_Pol2*ptBGeV*ptBGeV);
			RecoHadronicRecoil.SetPxPyPzE(resPX-MCBoson.Px(), resPY-MCBoson.Py(), 0, sqrt(pow(resPX-MCBoson.Px(),2)+pow(resPY-MCBoson.Py(),2)));
			/// Derived Quantities
			RecoETMiss.SetPxPyPzE(-(RecoMuon.Px()+RecoHadronicRecoil.Px()), -(RecoMuon.Py()+RecoHadronicRecoil.Py()), 0, sqrt(pow((RecoMuon.Px()+RecoHadronicRecoil.Px()),2)+pow((RecoMuon.Py()+RecoHadronicRecoil.Py()),2)));
			RecoTransverseMass = sqrt(2.0*RecoMuon.Pt()*RecoETMiss.Pt()*(1.0-TMath::Cos(RecoMuon.DeltaPhi(RecoETMiss))));

			/// Kinematic Cuts
			if (RecoMuon.Pt()<20000) return false;
			if (RecoETMiss.Pt()<20000) return false;
			if (RecoTransverseMass<20000) return false;
			if (RecoHadronicRecoil.Pt()>60000) return false; 
			
			isReconstructedW = true;
			return true;		
		}
	
	/// Private Functions
	private:
		void clearRecoKinematics() {
			RecoMuon.SetPxPyPzE(0,0,0,104);
			RecoMuon1.SetPxPyPzE(0,0,0,104);
			RecoMuon2.SetPxPyPzE(0,0,0,104);
			RecoETMiss.SetPxPyPzE(0,0,0,104);
			RecoHadronicRecoil.SetPxPyPzE(0,0,0,104);
			isReconstructedZ = false;
			isReconstructedW = false;
		}


		TH1F		*m_histMuonEfficiency;
		Double_t	m_muonPTResolution;
		Double_t	m_HRBias;
		Double_t	m_HRWidth_Pol0;
		Double_t	m_HRWidth_Pol1;
		Double_t	m_HRWidth_Pol2;

	/// Variables that contain reconstructed Information
	public:
		TLorentzVector	RecoMuon;				// W-Boson Decay Muon
		TLorentzVector	RecoMuon1;				// Z-Boson Decay Muon 1
		TLorentzVector	RecoMuon2;				// Z-Boson Decay Muon 2
		TLorentzVector	RecoETMiss;				// Reconstructed ETMiss
		TLorentzVector	RecoHadronicRecoil;		// Reconstructed Hadronic Recoil	
		Double_t	RecoTransverseMass;		// Reconstructed Transverse Mass (w.r.t to the first muon in case of Z Boson eventd)
		bool		isReconstructedZ;		// Flag if events is a reconstructed Z boson event
		bool		isReconstructedW;		// Flag if events is a reconstructed W boson event

};
