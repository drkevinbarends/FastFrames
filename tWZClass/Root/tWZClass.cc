#include "tWZClass/tWZClass.h"

#include "FastFrames/DefineHelpers.h"
#include "FastFrames/UniqueSampleID.h"

#include "Math/GenVector/Boost.h"
#include "Math/Math.h"
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>

using ROOT::VecOps::RVec;
using ROOT::VecOps::Take;
using ROOT::Math::Boost;

ROOT::RDF::RNode tWZClass::defineVariables(
  ROOT::RDF::RNode mainNode,
  const std::shared_ptr<Sample>& /*sample*/,
  const UniqueSampleID& id) 
{

  // You can also use the UniqueSampleID object to apply a custom defione
  // based on the sample and the subsample
  //   sample->name(): is the name of the sample defined in the config
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  // Define the mainNode with the campaign
  // const std::string campaign = id.campaign();
  // mainNode = MainFrame::systematicStringDefine(
  //   mainNode,
  //   "campaign_NOSYS",
  //   campaign
  // );

  // If running ML model
  // Function
  auto skipEvent = [](
      const unsigned long long& eventNumber
  ) {
    int result = 0;
    if (eventNumber % 2 != 0) { result = 1; }
  
    return result;
  };
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "skipEvent_NOSYS",
      skipEvent,
      {"eventNumber"}
  );

  /* 
    ======================================================
        Store which channel is being processed 
    ======================================================
  */
  // Function
  auto storeChannel = [](
      const char& channel
  ) {
    int result;
    if (channel) {
      // If channel is defined, store it
      result = 1;
    } else {
      // If channel is not defined, return empty string
      result = 0;
    }
    
    // Store the channel as a string
    return result;
  };
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "eeee_NOSYS",
      storeChannel,
      {"pass_eeee_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "eeemu_NOSYS",
      storeChannel,
      {"pass_eeemu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "eemumu_NOSYS",
      storeChannel,
      {"pass_eemumu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "emumumu_NOSYS",
      storeChannel,
      {"pass_emumumu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mumumumu_NOSYS",
      storeChannel,
      {"pass_mumumumu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "eee_NOSYS",
      storeChannel,
      {"pass_eee_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "eemu_NOSYS",
      storeChannel,
      {"pass_eemu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "emumu_NOSYS",
      storeChannel,
      {"pass_emumu_NOSYS"}
  );
  // Variable
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mumumu_NOSYS",
      storeChannel,
      {"pass_mumumu_NOSYS"}
  );

  /* 
    ======================================================
        Store which trigger is being processed 
    ======================================================
  */
  // Function
  auto storeTrigger = [](
      const bool& trigger
  ) {
    int result;
    if (trigger) {
      // If trigger is defined, store it
      result = 1;
    } else {
      // If trigger is not defined, return empty string
      result = 0;
    }

    // Store the trigger as a string
    return result;
  };
  // Triggers
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_mu50_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_mu50"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_mu40_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_mu40"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_mu26_ivarmedium_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_mu26_ivarmedium"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_mu20_iloose_L1MU15_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_mu20_iloose_L1MU15"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e24_lhmedium_L1EM20VH_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e24_lhmedium_L1EM20VH"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e26_lhtight_nod0_ivarloose_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e26_lhtight_nod0_ivarloose"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e60_lhmedium_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e60_lhmedium"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e60_lhmedium_nod0_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e60_lhmedium_nod0"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e120_lhloose_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e120_lhloose"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HLT_e140_lhloose_nod0_NOSYS",
      storeTrigger,
      {"trigPassed_HLT_e140_lhloose_nod0"}
  );

  /* 
    ======================================================
        Lepton triggers that passed
    ======================================================
  */
  // Function
  auto electronTriggers = [id](
    const int& e24_lhmedium,
    const int& e26_lhtight_nod0_ivarloose,
    const int& e60_lhmedium,
    const int& e60_lhmedium_nod0,
    const int& e120_lhloose,
    const int& e140_lhloose_nod0
  ) {
    int result = 0;
    if (id.campaign() == "2015") {
      if (e24_lhmedium == 1 || e60_lhmedium == 1 || e120_lhloose == 1) {
        result = 1;
      }
    }
    if (id.campaign() == "mc20a") {
      if (e24_lhmedium == 1 || e26_lhtight_nod0_ivarloose == 1 || e60_lhmedium == 1 || e60_lhmedium_nod0 == 1 || e120_lhloose == 1 || e140_lhloose_nod0 == 1) {
        result = 1;
      }
    }
    if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
      if (e26_lhtight_nod0_ivarloose == 1 || e60_lhmedium_nod0 == 1 || e140_lhloose_nod0 == 1) {
        result = 1;
      }
    }

    return result;
  };
  // Variable - Electron
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_triggers_NOSYS",
      electronTriggers,
      {
        "HLT_e24_lhmedium_L1EM20VH_NOSYS", "HLT_e26_lhtight_nod0_ivarloose_NOSYS", 
        "HLT_e60_lhmedium_NOSYS", "HLT_e60_lhmedium_nod0_NOSYS", 
        "HLT_e120_lhloose_NOSYS", "HLT_e140_lhloose_nod0_NOSYS"
      }
  );
  // Function
  auto muonTriggers = [id](
    const int& mu50,
    const int& mu40,
    const int& mu26_ivarmedium,
    const int& mu20_iloose_L1MU15
  ) {
    int result = 0;
    if (id.campaign() == "2015") {
      if (mu20_iloose_L1MU15 == 1 || mu40 == 1) {
        result = 1;
      }
    }
    if (id.campaign() == "mc20a") {
      if (mu20_iloose_L1MU15 == 1 || mu26_ivarmedium == 1 || mu40 == 1 || mu50 == 1) {
        result = 1;
      }
    }
    if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
      if (mu26_ivarmedium == 1 || mu50 == 1) {
        result = 1;
      }
    }

    return result;
  };
  // Variable - Muon
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu_triggers_NOSYS",
      muonTriggers,
      {
        "HLT_mu50_NOSYS", "HLT_mu40_NOSYS", 
        "HLT_mu26_ivarmedium_NOSYS", "HLT_mu20_iloose_L1MU15_NOSYS"
      }
  );

  /* 
    ======================================================
        Create a pdgID variable for electrons and muons 
    ======================================================
  */
  // Function
  auto electronPdgID = [] (
    const std::vector<float>& el_charge
  ) {
    const int size = el_charge.size();

    std::vector<float> result;
    for(int i = 0; i<size; i++){
      result.push_back( 11. * -1. * el_charge.at(i) );
    }

    return result;
  };
  // Variable - Electron
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_pdgID_NOSYS",
      electronPdgID,
      {"el_charge"}
  );
  // Function
  auto muonPdgID = [] (
    const std::vector<float>& mu_charge
  ) {
    const int size = mu_charge.size();

    std::vector<float> result;
    for(int i = 0; i<size; i++){
      result.push_back( 13. * -1. * mu_charge.at(i) );
    }

    return result;
  };
  // Variable - Muon
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu_pdgID_NOSYS",
      muonPdgID,
      {"mu_charge"}
  );

  /* 
    ====================================
        Convert pt and e to GeV 
    ====================================
  */
  // Function
  auto convertMeVToGeV_Vectors = [](
      const std::vector<float>& data
  ) {
      const int size = data.size();

      std::vector<float> result;
      for (int i = 0; i < size; ++i) {
        result.push_back( data.at(i) / 1.e3 );
      }

      return result;
  };
  // Variables
  // Jet
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_pt_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"jet_pt_NOSYS"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_e_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"jet_e_NOSYS"}
  );
  // Electron
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_pt_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"el_pt_NOSYS"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_e_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"el_e_NOSYS"}
  );
  // Muon
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu_pt_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"mu_pt_NOSYS"}
  );
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu_e_GeV_NOSYS",
      convertMeVToGeV_Vectors,
      {"mu_e_NOSYS"}
  );

  // Function
  auto convertMeVToGeV_Scalar = [](
      float& data
  ) {
      float result = data/ 1.e3;

      return result;
  };
  // Variable - MET
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "met_met_GeV_NOSYS",
      convertMeVToGeV_Scalar,
      {"met_met_NOSYS"}
  );

  /* 
    ===========================================
        Jet pt and eta cuts
    ===========================================
  */
  // Function - apply jet pt and eta cuts
  auto applyJetPtEtaCut = [](
      const std::vector<float>& jet_pt,
      const std::vector<float>& jet_eta
  ) {
      const int size = jet_pt.size();

      int result = 1;
      for (int i = 0; i < size; i++) {
        if (std::abs(jet_eta.at(i)) > 2.5 || jet_pt.at(i) < 20.) {
          result = 0;
          break;
        }
      }
      
      return result;
  };
  // Variable - Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "passedJetPtEtaCut_NOSYS",
      applyJetPtEtaCut,
      {"jet_pt_GeV_NOSYS", "jet_eta"}
  );

  /* 
    ===========================================
        Electron pt and eta cuts
    ===========================================
  */
  // Function - apply electron pt and eta cuts
  auto applyElectronPtEtaCut = [](
      const std::vector<float>& electron_pt,
      const std::vector<float>& electron_eta
  ) {
      const int size = electron_pt.size();

      int result = 1;
      for (int i = 0; i < size; i++) {
        if (
            std::abs(electron_eta.at(i)) > 2.47
            || electron_pt.at(i) < 10.
        ) {
            result = 0;
            break;
          }
      }
      
      return result;
  };
  // Variable - Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "passedElectronPtEtaCut_NOSYS",
      applyElectronPtEtaCut,
      {"el_pt_GeV_NOSYS","el_eta"}
  );
  /* 
    ===========================================
        Muon pt and eta cuts
    ===========================================
  */
  // Function - apply muon pt and eta cuts
  auto applyMuonPtEtaCut = [](
      const std::vector<float>& muon_pt,
      const std::vector<float>& muon_eta
  ) {
      const int size = muon_pt.size();

      int result = 1;
      for (int i = 0; i < size; i++) {
        if (
            std::abs(muon_eta.at(i)) > 2.5
            || muon_pt.at(i) < 10.
        ) {
            result = 0;
            break;
        }
      }

      return result;
  };
  // Variable - Muons 
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "passedMuonPtEtaCut_NOSYS",
      applyMuonPtEtaCut,
      {"mu_pt_GeV_NOSYS","mu_eta"}
  );

  /*
    ======================================================
        Calculate deltaR between electrons and jets 
    ======================================================
  */
  // Function - calculate deltaR
  auto calculateDeltaR = [](
      const std::vector<float>& jet_eta,
      const std::vector<float>& jet_phi,
      const std::vector<float>& el_eta,
      const std::vector<float>& el_phi
  ) {
      const int nJets = jet_eta.size();
      const int nElectrons = el_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nJets; ++i) {
          for (int j = 0; j < nElectrons; ++j) {
              float dEta = jet_eta[i] - el_eta[j];
              float dPhi = jet_phi[i] - el_phi[j];
              deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_el_deltaR_NOSYS",
      calculateDeltaR,
      {"jet_eta", "jet_phi", "el_eta", "el_phi"}
  );
  /*
    ======================================================
        Calculate deltaR between muons and jets 
    ======================================================
  */ 
  // Function - calculate deltaR
  auto calculateDeltaR_muon = [](
      const std::vector<float>& jet_eta,
      const std::vector<float>& jet_phi,
      const std::vector<float>& mu_eta,
      const std::vector<float>& mu_phi
  ) {
      const int nJets = jet_eta.size();
      const int nMuons = mu_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nJets; ++i) {
          for (int j = 0; j < nMuons; ++j) {
              float dEta = jet_eta[i] - mu_eta[j];
              float dPhi = jet_phi[i] - mu_phi[j];
              deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_mu_deltaR_NOSYS",
      calculateDeltaR_muon,
      {"jet_eta", "jet_phi", "mu_eta", "mu_phi"}
  );
  /* 
    ======================================================
        Calculate deltaR between electrons and muons 
    ======================================================
  */
  // Function - calculate deltaR
  auto calculateDeltaR_electron_muon = [](
      const std::vector<float>& el_eta,
      const std::vector<float>& el_phi,
      const std::vector<float>& mu_eta,
      const std::vector<float>& mu_phi
  ) {
      const int nElectrons = el_eta.size();
      const int nMuons = mu_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nElectrons; ++i) {
          for (int j = 0; j < nMuons; ++j) {
              float dEta = el_eta[i] - mu_eta[j];
              float dPhi = el_phi[i] - mu_phi[j];
              deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_mu_deltaR_NOSYS",
      calculateDeltaR_electron_muon,
      {"el_eta", "el_phi", "mu_eta", "mu_phi"}
  );

  /* 
    ======================================================
        Calculate deltaR between electrons and electrons 
    ======================================================
  */
  // Function - calculate deltaR
  auto calculateDeltaR_electron_electron = [](
      const std::vector<float>& el_eta,
      const std::vector<float>& el_phi
  ) {
      const int nElectrons = el_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nElectrons; ++i) {
          for (int j = 0; j < nElectrons; ++j) {
              if (i != j) { // Avoid self-comparison
                  float dEta = el_eta[i] - el_eta[j];
                  float dPhi = el_phi[i] - el_phi[j];
                  deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
              }
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el_el_deltaR_NOSYS",
      calculateDeltaR_electron_electron,
      {"el_eta", "el_phi"}
  );

  /* 
    ======================================================
        Calculate deltaR between muons and muons 
    ======================================================
  */
  // Function - calculate deltaR
  auto calculateDeltaR_muon_muon = [](
      const std::vector<float>& mu_eta,
      const std::vector<float>& mu_phi
  ) {
      const int nMuons = mu_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nMuons; ++i) {
          for (int j = 0; j < nMuons; ++j) {
              if (i != j) { // Avoid self-comparison
                  float dEta = mu_eta[i] - mu_eta[j];
                  float dPhi = mu_phi[i] - mu_phi[j];
                  deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
              }
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu_mu_deltaR_NOSYS",
      calculateDeltaR_muon_muon,
      {"mu_eta", "mu_phi"}
  );

  /* 
    ======================================================
        Calculate deltaR between jets and jets 
    ======================================================
  */
  // Function - calculate deltaR
  auto calculateDeltaR_jet_jet = [](
      const std::vector<float>& jet_eta,
      const std::vector<float>& jet_phi
  ) {
      const int nJets = jet_eta.size();

      std::vector<float> deltaR;
      for (int i = 0; i < nJets; ++i) {
          for (int j = 0; j < nJets; ++j) {
              if (i != j) { // Avoid self-comparison
                  float dEta = jet_eta[i] - jet_eta[j];
                  float dPhi = jet_phi[i] - jet_phi[j];
                  deltaR.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
              }
          }
      }

      return deltaR;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_jet_deltaR_NOSYS",
      calculateDeltaR_jet_jet,
      {"jet_eta", "jet_phi"}
  );

  /* 
    ======================================================
       Determine which electron triggers passed and which electron fired the trigger
    ======================================================
  */
  // Function - Determine which electron triggers passed and which electron fired the trigger
  auto electronFiredTriggers = [id](
    const int& e24_lhmedium,
    const int& e26_lhtight_nod0_ivarloose,
    const int& e60_lhmedium,
    const int& e60_lhmedium_nod0,
    const int& e120_lhloose,
    const int& e140_lhloose_nod0,
    const std::vector<float>& el_pt,
    const std::vector<char>& el_tight,
    const std::vector<char>& el_trigMatched_HLT_e24_lhmedium_L1EM20VH,
    const std::vector<char>& el_trigMatched_HLT_e26_lhtight_nod0_ivarloose,
    const std::vector<char>& el_trigMatched_HLT_e60_lhmedium,
    const std::vector<char>& el_trigMatched_HLT_e60_lhmedium_nod0,
    const std::vector<char>& el_trigMatched_HLT_e120_lhloose,
    const std::vector<char>& el_trigMatched_HLT_e140_lhloose_nod0
  ) {
    int result = 0;
    if (id.campaign() == "2015") {
      if (e24_lhmedium == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e24_lhmedium_L1EM20VH[i] && el_pt[i] > 26. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e60_lhmedium == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e60_lhmedium[i] && el_pt[i] > 62. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e120_lhloose == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e120_lhloose[i] && el_pt[i] > 122. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }
    if (id.campaign() == "mc20a") {
      if (e24_lhmedium == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e24_lhmedium_L1EM20VH[i] && el_pt[i] > 26. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e26_lhtight_nod0_ivarloose == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e26_lhtight_nod0_ivarloose[i] && el_pt[i] > 28. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e60_lhmedium == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e60_lhmedium[i] && el_pt[i] > 62. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e60_lhmedium_nod0 == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e60_lhmedium_nod0[i] && el_pt[i] > 62. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e120_lhloose == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e120_lhloose[i] && el_pt[i] > 122. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e140_lhloose_nod0 == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e140_lhloose_nod0[i] && el_pt[i] > 142. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }
    if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
      if (e26_lhtight_nod0_ivarloose == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e26_lhtight_nod0_ivarloose[i] && el_pt[i] > 28. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e60_lhmedium_nod0 == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e60_lhmedium_nod0[i] && el_pt[i] > 62. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (e140_lhloose_nod0 == 1) {
        for (size_t i = 0; i < el_pt.size(); ++i) {
          if (el_trigMatched_HLT_e140_lhloose_nod0[i] && el_pt[i] > 142. && el_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }

    return result;
  };
  // Variable - Electron
  if (id.campaign() == "2015") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_firedTriggers_NOSYS",
        electronFiredTriggers,
        {
          "HLT_e24_lhmedium_L1EM20VH_NOSYS", "HLT_e24_lhmedium_L1EM20VH_NOSYS", 
          "HLT_e60_lhmedium_NOSYS", "HLT_e60_lhmedium_NOSYS", 
          "HLT_e120_lhloose_NOSYS", "HLT_e120_lhloose_NOSYS",
          "el_pt_GeV_NOSYS", "el_select_tight_NOSYS",
          "el_trigMatched_HLT_e24_lhmedium_L1EM20VH", 
          "el_trigMatched_HLT_e24_lhmedium_L1EM20VH", 
          "el_trigMatched_HLT_e60_lhmedium", 
          "el_trigMatched_HLT_e60_lhmedium", 
          "el_trigMatched_HLT_e120_lhloose", 
          "el_trigMatched_HLT_e120_lhloose"
        }
    );
  }
  if (id.campaign() == "mc20a") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_firedTriggers_NOSYS",
        electronFiredTriggers,
        {
          "HLT_e24_lhmedium_L1EM20VH_NOSYS", "HLT_e26_lhtight_nod0_ivarloose_NOSYS", 
          "HLT_e60_lhmedium_NOSYS", "HLT_e60_lhmedium_nod0_NOSYS", 
          "HLT_e120_lhloose_NOSYS", "HLT_e140_lhloose_nod0_NOSYS",
          "el_pt_GeV_NOSYS", "el_select_tight_NOSYS",
          "el_trigMatched_HLT_e24_lhmedium_L1EM20VH", 
          "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", 
          "el_trigMatched_HLT_e60_lhmedium", 
          "el_trigMatched_HLT_e60_lhmedium_nod0", 
          "el_trigMatched_HLT_e120_lhloose", 
          "el_trigMatched_HLT_e140_lhloose_nod0"
        }
    );
  }
  if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_firedTriggers_NOSYS",
        electronFiredTriggers,
        {
          "HLT_e26_lhtight_nod0_ivarloose_NOSYS", "HLT_e26_lhtight_nod0_ivarloose_NOSYS", 
          "HLT_e60_lhmedium_nod0_NOSYS", "HLT_e60_lhmedium_nod0_NOSYS", 
          "HLT_e140_lhloose_nod0_NOSYS", "HLT_e140_lhloose_nod0_NOSYS",
          "el_pt_GeV_NOSYS", "el_select_tight_NOSYS",
          "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", 
          "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", 
          "el_trigMatched_HLT_e60_lhmedium_nod0", 
          "el_trigMatched_HLT_e60_lhmedium_nod0", 
          "el_trigMatched_HLT_e140_lhloose_nod0", 
          "el_trigMatched_HLT_e140_lhloose_nod0"
        }
    );
  }
  /* 
    ======================================================
        Determine which muon triggers passed and which muon fired the trigger
    ======================================================
  */
  // Function - Determine which muon triggers passed and which muon fired the trigger
  auto muonFiredTriggers = [id](
    const int& mu50,
    const int& mu40,
    const int& mu26_ivarmedium,
    const int& mu20_iloose_L1MU15,
    const std::vector<float>& mu_pt,
    const std::vector<char>& mu_tight,
    const std::vector<char>& mu_trigMatched_HLT_mu50,
    const std::vector<char>& mu_trigMatched_HLT_mu40,
    const std::vector<char>& mu_trigMatched_HLT_mu26_ivarmedium,
    const std::vector<char>& mu_trigMatched_HLT_mu20_iloose_L1MU15
  ) {
    int result = 0;
    if (id.campaign() == "2015") {
      if (mu20_iloose_L1MU15 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu20_iloose_L1MU15[i] && mu_pt[i] > 22. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (mu40 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu40[i] && mu_pt[i] > 42. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }
    if (id.campaign() == "mc20a") {
      if (mu20_iloose_L1MU15 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu20_iloose_L1MU15[i] && mu_pt[i] > 22. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (mu26_ivarmedium == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu26_ivarmedium[i] && mu_pt[i] > 28. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (mu40 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu40[i] && mu_pt[i] > 42. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (mu50 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu50[i] && mu_pt[i] > 52. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }
    if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
      if (mu26_ivarmedium == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu26_ivarmedium[i] && mu_pt[i] > 28. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
      if (mu50 == 1) {
        for (size_t i = 0; i < mu_pt.size(); ++i) {
          if (mu_trigMatched_HLT_mu50[i] && mu_pt[i] > 52. && mu_tight[i]) {
            result = 1;
            break;
          }
        }
      }
    }
    return result;
  };
  // Variable - Muon
  if (id.campaign() == "2015") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_firedTriggers_NOSYS",
        muonFiredTriggers,
        {
          "HLT_mu40_NOSYS", "HLT_mu40_NOSYS", 
          "HLT_mu20_iloose_L1MU15_NOSYS", "HLT_mu20_iloose_L1MU15_NOSYS",
          "mu_pt_GeV_NOSYS", "mu_select_tight_NOSYS",
          "mu_trigMatched_HLT_mu40", 
          "mu_trigMatched_HLT_mu40", 
          "mu_trigMatched_HLT_mu20_iloose_L1MU15", 
          "mu_trigMatched_HLT_mu20_iloose_L1MU15"
        }
    );
  }
  if (id.campaign() == "mc20a") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_firedTriggers_NOSYS",
        muonFiredTriggers,
        {
          "HLT_mu50_NOSYS", "HLT_mu40_NOSYS", 
          "HLT_mu26_ivarmedium_NOSYS", "HLT_mu20_iloose_L1MU15_NOSYS",
          "mu_pt_GeV_NOSYS", "mu_select_tight_NOSYS",
          "mu_trigMatched_HLT_mu50", 
          "mu_trigMatched_HLT_mu40", 
          "mu_trigMatched_HLT_mu26_ivarmedium", 
          "mu_trigMatched_HLT_mu20_iloose_L1MU15"
        }
    );
  }
  if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_firedTriggers_NOSYS",
        muonFiredTriggers,
        {
          "HLT_mu50_NOSYS", "HLT_mu50_NOSYS", 
          "HLT_mu26_ivarmedium_NOSYS", "HLT_mu26_ivarmedium_NOSYS",
          "mu_pt_GeV_NOSYS", "mu_select_tight_NOSYS",
          "mu_trigMatched_HLT_mu50", 
          "mu_trigMatched_HLT_mu50", 
          "mu_trigMatched_HLT_mu26_ivarmedium", 
          "mu_trigMatched_HLT_mu26_ivarmedium"
        }
    );
  }

  /*
    ======================================================
        Combine electron variables into one vector 
    ======================================================
  */
  // Function - Data
  auto combineElectronVariablesData = [id](
      const std::vector<float>& pT,
      const std::vector<float>& e,
      const std::vector<float>& eta,
      const std::vector<float>& phi,
      const std::vector<float>& charge,
      const std::vector<char>& medium_noniso,
      const std::vector<char>& medium,
      const std::vector<char>& tight_noniso,
      const std::vector<char>& tight,
      const std::vector<char>& HLT_e24_e26,
      const std::vector<char>& HLT_e60,
      const std::vector<char>& HLT_e120_e140
  ) {
      // Convert all variables to float
      std::vector<float> e_float(e.begin(), e.end());
      std::vector<float> eta_float(eta.begin(), eta.end());
      std::vector<float> phi_float(phi.begin(), phi.end());
      std::vector<float> charge_float(charge.begin(), charge.end());
      std::vector<float> medium_noniso_float(medium_noniso.begin(), medium_noniso.end());
      std::vector<float> medium_float(medium.begin(), medium.end());
      std::vector<float> tight_noniso_float(tight_noniso.begin(), tight_noniso.end());
      std::vector<float> tight_float(tight.begin(), tight.end());
      std::vector<float> HLT_e24_e26_float(HLT_e24_e26.begin(), HLT_e24_e26.end());
      std::vector<float> HLT_e60_float(HLT_e60.begin(), HLT_e60.end());
      std::vector<float> HLT_e120_e140_float(HLT_e120_e140.begin(), HLT_e120_e140.end());

      const int size = pT.size();

      // Create a vector of indices to sort the other variables
      std::vector<int> indices(size);
      for (int i = 0; i < size; ++i) {
          indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), [&pT](int a, int b) {
          return pT[a] > pT[b];
      });

      // Create one large vector to store all the variables in order of electron
      std::vector<float> all_variables;
      for (int i = 0; i < size; ++i) {
          all_variables.push_back(pT[indices[i]]);
          all_variables.push_back(e_float[indices[i]]);
          all_variables.push_back(eta_float[indices[i]]);
          all_variables.push_back(phi_float[indices[i]]);
          all_variables.push_back(11.);  // pdgFlavour for electrons
          all_variables.push_back(charge_float[indices[i]]);
          all_variables.push_back(0.);
          if (tight_float[indices[i]] == 1) {
            all_variables.push_back(4.);
          } else if (tight_noniso_float[indices[i]] == 1) {
            all_variables.push_back(3.);
          } else if (medium_float[indices[i]] == 1) {
            all_variables.push_back(2.);
          } else if (medium_noniso_float[indices[i]] == 1) {
            all_variables.push_back(1.);
          } else {
            all_variables.push_back(0.);
          }
          if (id.campaign() == "2015") {
            if (HLT_e24_e26_float.size() == 1. && pT[indices[i]] > 24.)
                all_variables.push_back(1.);
            else if (HLT_e60_float.size() == 1. && pT[indices[i]] > 62.)
              all_variables.push_back(1.);
            else if (HLT_e120_e140_float.size() == 1. && pT[indices[i]] > 122.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }
          if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018") {
            if (HLT_e24_e26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_e60_float.size() == 1. && pT[indices[i]] > 62.)
              all_variables.push_back(1.);
            else if (HLT_e120_e140_float.size() == 1. && pT[indices[i]] > 142.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }
      }

      return all_variables;
  };
  // Variable - Electron
  if (id.campaign() == "2015") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_variables_NOSYS",
        combineElectronVariablesData,
        {
            "el_pt_GeV_NOSYS", "el_e_GeV_NOSYS", "el_eta", "el_phi",
            "el_charge", "el_select_medium_noniso_NOSYS", "el_select_medium_NOSYS", "el_select_tight_noniso_NOSYS", "el_select_tight_NOSYS",
            "el_trigMatched_HLT_e24_lhmedium_L1EM20VH", "el_trigMatched_HLT_e60_lhmedium", "el_trigMatched_HLT_e120_lhloose"
        }
    );
  }
  if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_variables_NOSYS",
        combineElectronVariablesData,
        {
            "el_pt_GeV_NOSYS", "el_e_GeV_NOSYS", "el_eta", "el_phi",
            "el_charge", "el_select_medium_noniso_NOSYS", "el_select_medium_NOSYS", "el_select_tight_noniso_NOSYS", "el_select_tight_NOSYS",
            "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", "el_trigMatched_HLT_e60_lhmedium_nod0", "el_trigMatched_HLT_e140_lhloose_nod0"
        }
    );
  }
  // Function - MC
  auto combineElectronVariablesMC = [id](
      const std::vector<float>& pT,
      const std::vector<float>& e,
      const std::vector<float>& eta,
      const std::vector<float>& phi,
      const std::vector<float>& charge,
      const std::vector<int>& IFF_Class,
      const std::vector<char>& medium_noniso,
      const std::vector<char>& medium,
      const std::vector<char>& tight_noniso,
      const std::vector<char>& tight,
      const std::vector<char>& HLT_e24,
      const std::vector<char>& HLT_e26,
      const std::vector<char>& HLT_e60,
      const std::vector<char>& HLT_e60_nod0,
      const std::vector<char>& HLT_e120,
      const std::vector<char>& HLT_e140
  ) {
      // Convert all variables to float
      std::vector<float> e_float(e.begin(), e.end());
      std::vector<float> eta_float(eta.begin(), eta.end());
      std::vector<float> phi_float(phi.begin(), phi.end());
      std::vector<float> charge_float(charge.begin(), charge.end());
      std::vector<float> IFF_Class_float(IFF_Class.begin(), IFF_Class.end());
      std::vector<float> medium_noniso_float(medium_noniso.begin(), medium_noniso.end());
      std::vector<float> medium_float(medium.begin(), medium.end());
      std::vector<float> tight_noniso_float(tight_noniso.begin(), tight_noniso.end());
      std::vector<float> tight_float(tight.begin(), tight.end());
      std::vector<float> HLT_e24_float(HLT_e24.begin(), HLT_e24.end());
      std::vector<float> HLT_e26_float(HLT_e26.begin(), HLT_e26.end());
      std::vector<float> HLT_e60_float(HLT_e60.begin(), HLT_e60.end());
      std::vector<float> HLT_e60_nod0_float(HLT_e60_nod0.begin(), HLT_e60_nod0.end());
      std::vector<float> HLT_e120_float(HLT_e120.begin(), HLT_e120.end());
      std::vector<float> HLT_e140_float(HLT_e140.begin(), HLT_e140.end());

      const int size = pT.size();

      // Create a vector of indices to sort the other variables
      std::vector<int> indices(size);
      for (int i = 0; i < size; ++i) {
          indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), [&pT](int a, int b) {
          return pT[a] > pT[b];
      });

      // Create one large vector to store all the variables in order of electron
      std::vector<float> all_variables;
      for (int i = 0; i < size; ++i) {
          all_variables.push_back(pT[indices[i]]);
          all_variables.push_back(e_float[indices[i]]);
          all_variables.push_back(eta_float[indices[i]]);
          all_variables.push_back(phi_float[indices[i]]);
          all_variables.push_back(11.);  // pdgFlavour for electrons
          all_variables.push_back(charge_float[indices[i]]);
          all_variables.push_back(IFF_Class_float[indices[i]]);
          if (tight_float[indices[i]] == 1) {
            all_variables.push_back(4.);
          } else if (tight_noniso_float[indices[i]] == 1) {
            all_variables.push_back(3.);
          } else if (medium_float[indices[i]] == 1) {
            all_variables.push_back(2.);
          } else if (medium_noniso_float[indices[i]] == 1) {
            all_variables.push_back(1.);
          } else {
            all_variables.push_back(0.);
          }
          if (id.campaign() == "mc20a") {
            if (HLT_e24_float.size() == 1. && pT[indices[i]] > 26.)
                all_variables.push_back(1.);
            else if (HLT_e26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_e60_float.size() == 1. && pT[indices[i]] > 62.)
              all_variables.push_back(1.);
            else if (HLT_e60_nod0_float.size() == 1. && pT[indices[i]] > 62.)
              all_variables.push_back(1.);
            else if (HLT_e120_float.size() == 1. && pT[indices[i]] > 122.)
              all_variables.push_back(1.);
            else if (HLT_e140_float.size() == 1. && pT[indices[i]] > 142.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }
          if (id.campaign() == "mc20d" || id.campaign() == "mc20e") {
            if (HLT_e26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_e60_nod0_float.size() == 1. && pT[indices[i]] > 62.)
              all_variables.push_back(1.);
            else if (HLT_e140_float.size() == 1. && pT[indices[i]] > 142.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }

      }

      return all_variables;
  };
  // Variable - Electron
  if (id.campaign() == "mc20a") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_variables_NOSYS",
        combineElectronVariablesMC,
        {
            "el_pt_GeV_NOSYS", "el_e_GeV_NOSYS", "el_eta", "el_phi",
            "el_charge", "el_IFFClass", "el_select_medium_noniso_NOSYS", "el_select_medium_NOSYS", "el_select_tight_noniso_NOSYS", "el_select_tight_NOSYS",
            "el_trigMatched_HLT_e24_lhmedium_L1EM20VH", "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", "el_trigMatched_HLT_e60_lhmedium",
            "el_trigMatched_HLT_e60_lhmedium_nod0", "el_trigMatched_HLT_e120_lhloose", "el_trigMatched_HLT_e140_lhloose_nod0"
        }
    );
  }
  if (id.campaign() == "mc20d" || id.campaign() == "mc20e") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "el_variables_NOSYS",
        combineElectronVariablesMC,
        {
            "el_pt_GeV_NOSYS", "el_e_GeV_NOSYS", "el_eta", "el_phi",
            "el_charge", "el_IFFClass", "el_select_medium_noniso_NOSYS", "el_select_medium_NOSYS", "el_select_tight_noniso_NOSYS", "el_select_tight_NOSYS",
            "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", "el_trigMatched_HLT_e60_lhmedium_nod0",
            "el_trigMatched_HLT_e60_lhmedium_nod0", "el_trigMatched_HLT_e140_lhloose_nod0", "el_trigMatched_HLT_e140_lhloose_nod0"
        }
    );
  }

  /*
    ======================================================
      Total loose electrons from el_variables_NOSYS
    ======================================================
  */
  // Function - Count electrons
  auto countLooseElectrons = [](
      const std::vector<float>& el_variables
  ) {
      const int size = el_variables.size()/9;

      return size;
  };
  // Variable - Count Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nLooseElectrons_NOSYS",
      countLooseElectrons,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Total medium_noniso electrons from el_variables_NOSYS
    ======================================================
  */
  // Function - Count electrons
  auto countMediumNonIsoElectrons = [](
      const std::vector<float>& el_variables
  ) {
      const int size = el_variables.size();

      int count = 0;
      for (int i = 0; i < size/9; ++i) {
          if (el_variables[i * 9 + 7] >= 1.) { // Check if medium_noniso
              count++;
          }
      }
      return count;
  };
  // Variable - Count Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nMediumNonIsoElectrons_NOSYS",
      countMediumNonIsoElectrons,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Total medium electrons from el_variables_NOSYS
    ======================================================
  */
  // Function - Count electrons
  auto countMediumElectrons = [](
      const std::vector<float>& el_variables
  ) {
      const int size = el_variables.size();

      int count = 0;
      for (int i = 0; i < size/9; ++i) {
          if (el_variables[i * 9 + 7] == 2. || el_variables[i * 9 + 7] == 4.) { // Check if medium
              count++;
          }
      }
      return count;
  };
  // Variable - Count Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nMediumElectrons_NOSYS",
      countMediumElectrons,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Total tight_noniso electrons from el_variables_NOSYS
    ======================================================
  */
  // Function - Count electrons
  auto countTightNonIsoElectrons = [](
      const std::vector<float>& el_variables
  ) {
      const int size = el_variables.size();

      int count = 0;
      for (int i = 0; i < size/9; ++i) {
          if (el_variables[i * 9 + 7] >= 3.) { // Check if tight_noniso
              count++;
          }
      }
      return count;
  };
  // Variable - Count Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nTightNonIsoElectrons_NOSYS",
      countTightNonIsoElectrons,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Total tight electrons from el_variables_NOSYS
    ======================================================
  */
  // Function - Count electrons
  auto countTightElectrons = [](
      const std::vector<float>& el_variables
  ) {
      const int size = el_variables.size();

      int count = 0;
      for (int i = 0; i < size/9; ++i) {
          if (el_variables[i * 9 + 7] == 4.) { // Check if tight
              count++;
          }
      }
      return count;
  };
  // Variable - Count Electrons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nTightElectrons_NOSYS",
      countTightElectrons,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Extract 1st electron pt (i.e., leading electron)
    ======================================================
  */
  // Function - Extract leading electron pt
  auto leadingElectronPt = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.empty()) return -999.0f; // Return -999 if no electrons

      // The first element in the vector is the leading electron's pt
      return el_variables[0];
  };
  // Variable - Leading Electron Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el0_pt_NOSYS",
      leadingElectronPt,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Extract 1st electron eta (i.e., leading electron)
    ======================================================
  */
  // Function - Extract leading electron eta
  auto leadingElectronEta = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.empty()) return -999.0f; // Return -999 if no electrons

      // The third element in the vector is the leading electron's eta
      return el_variables[2];
  };
  // Variable - Leading Electron Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el0_eta_NOSYS",
      leadingElectronEta,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 1st electron phi (i.e., leading electron)
    ======================================================
  */
  // Function - Extract leading electron phi
  auto leadingElectronPhi = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.empty()) return -999.0f; // Return -999 if no electrons

      // The fourth element in the vector is the leading electron's phi
      return el_variables[3];
  };
  // Variable - Leading Electron Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el0_phi_NOSYS",
      leadingElectronPhi,
      {"el_variables_NOSYS"}
  );

  /*
    ======================================================
      Extract 2nd electron pt (i.e., subleading electron)
    ======================================================
  */
  // Function - Extract subleading electron pt
  auto subleadingElectronPt = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 2) return -999.0f; // Return -999 if no second electron

      // The first element in the second electron's vector is the subleading electron's pt
      return el_variables[9];
  };
  // Variable - Subleading Electron Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el1_pt_NOSYS",
      subleadingElectronPt,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 2nd electron eta (i.e., subleading electron)
    ======================================================
  */
  // Function - Extract subleading electron eta
  auto subleadingElectronEta = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 2) return -999.0f; // Return -999 if no second electron

      // The third element in the second electron's vector is the subleading electron's eta
      return el_variables[11];
  };
  // Variable - Subleading Electron Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el1_eta_NOSYS",
      subleadingElectronEta,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 2nd electron phi (i.e., subleading electron)
    ======================================================
  */
  // Function - Extract subleading electron phi
  auto subleadingElectronPhi = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 2) return -999.0f; // Return -999 if no second electron

      // The fourth element in the second electron's vector is the subleading electron's phi
      return el_variables[12];
  };
  // Variable - Subleading Electron Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el1_phi_NOSYS",
      subleadingElectronPhi,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd electron pt (i.e., third electron)
    ======================================================
  */
  // Function - Extract third electron pt
  auto thirdElectronPt = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 3) return -999.0f; // Return -999 if no third electron

      // The first element in the third electron's vector is the third electron's pt
      return el_variables[18];
  };
  // Variable - Third Electron Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el2_pt_NOSYS",
      thirdElectronPt,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd electron eta (i.e., third electron)
    ======================================================
  */
  // Function - Extract third electron eta
  auto thirdElectronEta = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 3) return -999.0f; // Return -999 if no third electron

      // The third element in the third electron's vector is the third electron's eta
      return el_variables[20];
  };
  // Variable - Third Electron Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el2_eta_NOSYS",
      thirdElectronEta,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd electron phi (i.e., third electron)
    ======================================================
  */
  // Function - Extract third electron phi
  auto thirdElectronPhi = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 3) return -999.0f; // Return -999 if no third electron

      // The fourth element in the third electron's vector is the third electron's phi
      return el_variables[21];
  };
  // Variable - Third Electron Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el2_phi_NOSYS",
      thirdElectronPhi,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th electron pt (i.e., fourth electron)
    ======================================================
  */
  // Function - Extract fourth electron pt
  auto fourthElectronPt = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth electron

      // The first element in the fourth electron's vector is the fourth electron's pt
      return el_variables[27];
  };
  // Variable - Fourth Electron Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el3_pt_NOSYS",
      fourthElectronPt,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th electron eta (i.e., fourth electron)
    ======================================================
  */
  // Function - Extract fourth electron eta
  auto fourthElectronEta = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth electron

      // The third element in the fourth electron's vector is the fourth electron's eta
      return el_variables[29];
  };
  // Variable - Fourth Electron Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el3_eta_NOSYS",
      fourthElectronEta,
      {"el_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th electron phi (i.e., fourth electron)
    ======================================================
  */
  // Function - Extract fourth electron phi
  auto fourthElectronPhi = [](
      const std::vector<float>& el_variables
  ) {
      if (el_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth electron

      // The fourth element in the fourth electron's vector is the fourth electron's phi
      return el_variables[30];
  };
  // Variable - Fourth Electron Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "el3_phi_NOSYS",
      fourthElectronPhi,
      {"el_variables_NOSYS"}
  );

  /* 
    ======================================================
        Combine muon variables into one vector 
    ======================================================
  */
  // Function - Data
  auto combineMuonVariablesData = [id](
      const std::vector<float>& pT,
      const std::vector<float>& e,
      const std::vector<float>& eta,
      const std::vector<float>& phi,
      const std::vector<float>& charge,
      const std::vector<char>& medium,
      const std::vector<char>& tight,
      const std::vector<char>& HLT_mu20_mu26,
      const std::vector<char>& HLT_mu40_mu50
  ) {
      // Convert all variables to float
      std::vector<float> e_float(e.begin(), e.end());
      std::vector<float> eta_float(eta.begin(), eta.end());
      std::vector<float> phi_float(phi.begin(), phi.end());
      std::vector<float> charge_float(charge.begin(), charge.end());
      std::vector<float> medium_float(medium.begin(), medium.end());
      std::vector<float> tight_float(tight.begin(), tight.end());
      std::vector<float> HLT_mu20_mu26_float(HLT_mu20_mu26.begin(), HLT_mu20_mu26.end());
      std::vector<float> HLT_mu40_mu50_float(HLT_mu40_mu50.begin(), HLT_mu40_mu50.end());  

      const int size = pT.size();

      // Create a vector of indices to sort the other variables
      std::vector<int> indices(size);
      for (int i = 0; i < size; ++i) {
          indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), [&pT](int a, int b) {
          return pT[a] > pT[b];
      });

      // Create one large vector to store all the variables
      std::vector<float> all_variables;
      for (int i = 0; i < size; ++i) {
          all_variables.push_back(pT[indices[i]]);
          all_variables.push_back(e_float[indices[i]]);
          all_variables.push_back(eta_float[indices[i]]);
          all_variables.push_back(phi_float[indices[i]]);
          all_variables.push_back(13.);  // pdgFlavour for muons
          all_variables.push_back(charge_float[indices[i]]);
          all_variables.push_back(0.);
          if (tight_float[indices[i]] == 1) {
            all_variables.push_back(2.);
          } else if (medium_float[indices[i]] == 1) {
            all_variables.push_back(1.);
          } else {
            all_variables.push_back(0.);
          }
          if (id.campaign() == "2015") {
            if (HLT_mu20_mu26_float.size() == 1. && pT[indices[i]] > 22.)
                all_variables.push_back(1.);
            else if (HLT_mu40_mu50_float.size() == 1. && pT[indices[i]] > 42.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }
          if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018") {
            if (HLT_mu20_mu26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_mu40_mu50_float.size() == 1. && pT[indices[i]] > 52.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger, set to 0
          }
      }

      return all_variables;
  };
  // Variable - Muon
  if (id.campaign() == "2015") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_variables_NOSYS",
        combineMuonVariablesData,
        {
            "mu_pt_GeV_NOSYS", "mu_e_GeV_NOSYS", "mu_eta", "mu_phi",
            "mu_charge", "mu_select_medium_NOSYS", "mu_select_tight_NOSYS", 
            "mu_trigMatched_HLT_mu20_iloose_L1MU15", "mu_trigMatched_HLT_mu40"
        }
    );
  }
  if (id.campaign() == "2016" || id.campaign() == "2017" || id.campaign() == "2018") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_variables_NOSYS",
        combineMuonVariablesData,
        {
            "mu_pt_GeV_NOSYS", "mu_e_GeV_NOSYS", "mu_eta", "mu_phi",
            "mu_charge", "mu_select_medium_NOSYS", "mu_select_tight_NOSYS",
            "mu_trigMatched_HLT_mu26_ivarmedium", "mu_trigMatched_HLT_mu50"
        }
    );
  }
  // Function - MC
  auto combineMuonVariablesMC = [id](
      const std::vector<float>& pT,
      const std::vector<float>& e,
      const std::vector<float>& eta,
      const std::vector<float>& phi,
      const std::vector<float>& charge,
      const std::vector<int>& IFF_Class,
      const std::vector<char>& medium,
      const std::vector<char>& tight,
      const std::vector<char>& HLT_mu20,
      const std::vector<char>& HLT_mu26,
      const std::vector<char>& HLT_mu40,
      const std::vector<char>& HLT_mu50
  ) {
      // Convert all variables to float
      std::vector<float> e_float(e.begin(), e.end());
      std::vector<float> eta_float(eta.begin(), eta.end());
      std::vector<float> phi_float(phi.begin(), phi.end());
      std::vector<float> charge_float(charge.begin(), charge.end());
      std::vector<float> IFF_Class_float(IFF_Class.begin(), IFF_Class.end());
      std::vector<float> medium_float(medium.begin(), medium.end());
      std::vector<float> tight_float(tight.begin(), tight.end());
      std::vector<float> HLT_mu20_float(HLT_mu20.begin(), HLT_mu20.end());
      std::vector<float> HLT_mu26_float(HLT_mu26.begin(), HLT_mu26.end());
      std::vector<float> HLT_mu40_float(HLT_mu40.begin(), HLT_mu40.end());
      std::vector<float> HLT_mu50_float(HLT_mu50.begin(), HLT_mu50.end());

      const int size = pT.size();

      // Create a vector of indices to sort the other variables
      std::vector<int> indices(size);
      for (int i = 0; i < size; ++i) {
          indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), [&pT](int a, int b) {
          return pT[a] > pT[b];
      });


      // Create one large vector to store all the variables
      std::vector<float> all_variables;
      for (int i = 0; i < size; ++i) {
          all_variables.push_back(pT[indices[i]]);
          all_variables.push_back(e_float[indices[i]]);
          all_variables.push_back(eta_float[indices[i]]);
          all_variables.push_back(phi_float[indices[i]]);
          all_variables.push_back(13.);  // pdgFlavour for muons
          all_variables.push_back(charge_float[indices[i]]);
          all_variables.push_back(IFF_Class_float[indices[i]]);
          if (tight_float[indices[i]] == 1) {
            all_variables.push_back(2.);
          } else if (medium_float[indices[i]] == 1) {
            all_variables.push_back(1.);
          } else {
            all_variables.push_back(0.);
          }
          if (id.campaign() == "mc20a") {
            if (HLT_mu20_float.size() == 1. && pT[indices[i]] > 22.)
                all_variables.push_back(1.);
            else if (HLT_mu26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_mu40_float.size() == 1. && pT[indices[i]] > 42.)
              all_variables.push_back(1.);
            else if (HLT_mu50_float.size() == 1. && pT[indices[i]] > 52.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger match, set to 0
          }
          if (id.campaign() == "mc20d" || id.campaign() == "mc20e") {
            if (HLT_mu26_float.size() == 1. && pT[indices[i]] > 28.)
              all_variables.push_back(1.);
            else if (HLT_mu50_float.size() == 1. && pT[indices[i]] > 52.)
              all_variables.push_back(1.);
            else
              all_variables.push_back(0.); // If no trigger match, set to 0
          }
      }

      return all_variables;
  };
  // Variable - Muon
  if (id.campaign() == "mc20a") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_variables_NOSYS",
        combineMuonVariablesMC,
        {
            "mu_pt_GeV_NOSYS", "mu_e_GeV_NOSYS", "mu_eta", "mu_phi",
            "mu_charge", "mu_IFFClass", "mu_select_medium_NOSYS", "mu_select_tight_NOSYS",
            "mu_trigMatched_HLT_mu20_iloose_L1MU15", "mu_trigMatched_HLT_mu26_ivarmedium", 
            "mu_trigMatched_HLT_mu40", "mu_trigMatched_HLT_mu50"
        }
    );
  }
  if (id.campaign() == "mc20d" || id.campaign() == "mc20e") {
    mainNode = MainFrame::systematicDefine(
        mainNode,
        "mu_variables_NOSYS",
        combineMuonVariablesMC,
        {
            "mu_pt_GeV_NOSYS", "mu_e_GeV_NOSYS", "mu_eta", "mu_phi",
            "mu_charge", "mu_IFFClass", "mu_select_medium_NOSYS", "mu_select_tight_NOSYS",
            "mu_trigMatched_HLT_mu26_ivarmedium", "mu_trigMatched_HLT_mu26_ivarmedium", 
            "mu_trigMatched_HLT_mu50", "mu_trigMatched_HLT_mu50"
        }
    );
  }

  /*
    ======================================================
      Total loose muons from mu_variables_NOSYS
    ======================================================
  */
  // Function - Count muons
  auto countLooseMuons = [](
      const std::vector<float>& mu_variables
  ) {
      const int size = mu_variables.size()/9;

      return size;
  };
  // Variable - Count Muons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nLooseMuons_NOSYS",
      countLooseMuons,
      {"mu_variables_NOSYS"}
  );

  /*
    ======================================================
      Total medium muons from mu_variables_NOSYS
    ======================================================
  */
  // Function - Count muons
  auto countMediumMuons = [](
      const std::vector<float>& mu_variables
  ) {
      const int size = mu_variables.size()/9;

      int count = 0;
      for (int i = 0; i < size; ++i) {
          if (mu_variables[i * 9 + 7] >= 1.) { // Check if medium
              count++;
          }
      }
      return count;
  };
  // Variable - Count Muons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nMediumMuons_NOSYS",
      countMediumMuons,
      {"mu_variables_NOSYS"}
  );

  /*
    ======================================================
      Total tight muons from mu_variables_NOSYS
    ======================================================
  */
  // Function - Count muons
  auto countTightMuons = [](
      const std::vector<float>& mu_variables
  ) {
      const int size = mu_variables.size();

      int count = 0;
      for (int i = 0; i < size/9; ++i) {
          if (mu_variables[i * 9 + 7] == 2.) { // Check if tight
              count++;
          }
      }
      return count;
  };
  // Variable - Count Muons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nTightMuons_NOSYS",
      countTightMuons,
      {"mu_variables_NOSYS"}
  );

  /*
    ======================================================
      Extract 1st muon pt (i.e., leading muon)
    ======================================================
  */
  // Function - Extract leading muon pt
  auto leadingMuonPt = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.empty()) return -999.0f; // Return -999 if no muons

      // The first element in the vector is the leading muon's pt
      return mu_variables[0];
  };
  // Variable - Leading Muon Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu0_pt_NOSYS",
      leadingMuonPt,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 1st muon eta (i.e., leading muon)
    ======================================================
  */
  // Function - Extract leading muon eta
  auto leadingMuonEta = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.empty()) return -999.0f; // Return -999 if no muons

      // The third element in the vector is the leading muon's eta
      return mu_variables[2];
  };
  // Variable - Leading Muon Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu0_eta_NOSYS",
      leadingMuonEta,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 1st muon phi (i.e., leading muon)
    ======================================================
  */
  // Function - Extract leading muon phi
  auto leadingMuonPhi = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.empty()) return -999.0f; // Return -999 if no muons

      // The fourth element in the vector is the leading muon's phi
      return mu_variables[3];
  };
  // Variable - Leading Muon Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu0_phi_NOSYS",
      leadingMuonPhi,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 2nd muon pt (i.e., subleading muon)
    ======================================================
  */
  // Function - Extract subleading muon pt
  auto subleadingMuonPt = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 2) return -999.0f; // Return -999 if no second muon

      // The first element in the second muon's vector is the subleading muon's pt
      return mu_variables[9];
  };
  // Variable - Subleading Muon Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu1_pt_NOSYS",
      subleadingMuonPt,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 2nd muon eta (i.e., subleading muon)
    ======================================================
  */
  // Function - Extract subleading muon eta
  auto subleadingMuonEta = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 2) return -999.0f; // Return -999 if no second muon

      // The third element in the second muon's vector is the subleading muon's eta
      return mu_variables[11];
  };
  // Variable - Subleading Muon Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu1_eta_NOSYS",
      subleadingMuonEta,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 2nd muon phi (i.e., subleading muon)
    ======================================================
  */
  // Function - Extract subleading muon phi
  auto subleadingMuonPhi = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 2) return -999.0f; // Return -999 if no second muon

      // The fourth element in the second muon's vector is the subleading muon's phi
      return mu_variables[12];
  };
  // Variable - Subleading Muon Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu1_phi_NOSYS",
      subleadingMuonPhi,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd muon pt (i.e., third muon)
    ======================================================
  */
  // Function - Extract third muon pt
  auto thirdMuonPt = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 3) return -999.0f; // Return -999 if no third muon

      // The first element in the third muon's vector is the third muon's pt
      return mu_variables[18];
  };
  // Variable - Third Muon Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu2_pt_NOSYS",
      thirdMuonPt,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd muon eta (i.e., third muon)
    ======================================================
  */
  // Function - Extract third muon eta
  auto thirdMuonEta = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 3) return -999.0f; // Return -999 if no third muon

      // The third element in the third muon's vector is the third muon's eta
      return mu_variables[20];
  };
  // Variable - Third Muon Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu2_eta_NOSYS",
      thirdMuonEta,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 3rd muon phi (i.e., third muon)
    ======================================================
  */
  // Function - Extract third muon phi
  auto thirdMuonPhi = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 3) return -999.0f; // Return -999 if no third muon

      // The fourth element in the third muon's vector is the third muon's phi
      return mu_variables[21];
  };
  // Variable - Third Muon Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu2_phi_NOSYS",
      thirdMuonPhi,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th muon pt (i.e., fourth muon)
    ======================================================
  */
  // Function - Extract fourth muon pt
  auto fourthMuonPt = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth muon

      // The first element in the fourth muon's vector is the fourth muon's pt
      return mu_variables[27];
  };
  // Variable - Fourth Muon Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu3_pt_NOSYS",
      fourthMuonPt,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th muon eta (i.e., fourth muon)
    ======================================================
  */
  // Function - Extract fourth muon eta
  auto fourthMuonEta = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth muon

      // The third element in the fourth muon's vector is the fourth muon's eta
      return mu_variables[29];
  };
  // Variable - Fourth Muon Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu3_eta_NOSYS",
      fourthMuonEta,
      {"mu_variables_NOSYS"}
  );
  /*
    ======================================================
      Extract 4th muon phi (i.e., fourth muon)
    ======================================================
  */
  // Function - Extract fourth muon phi
  auto fourthMuonPhi = [](
      const std::vector<float>& mu_variables
  ) {
      if (mu_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth muon

      // The fourth element in the fourth muon's vector is the fourth muon's phi
      return mu_variables[30];
  };
  // Variable - Fourth Muon Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "mu3_phi_NOSYS",
      fourthMuonPhi,
      {"mu_variables_NOSYS"}
  );

  /*
    ====================================================
        Count leptons from electrons and muons
    ====================================================
  */
  // Function - Count leptons
  auto countLeptons = [](
      const int& nElectrons,
      const int& nMuons
  ) {

      // Return the total count of leptons
      int nLeptons = nElectrons + nMuons;

      return nLeptons;
  };
  // Variable - Count Loose Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nLooseLeptons_NOSYS",
      countLeptons,
      {"nLooseElectrons_NOSYS", "nLooseMuons_NOSYS"}
  );
  // Variable - Count Tight Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nTightLeptons_NOSYS",
      countLeptons,
      {"nTightElectrons_NOSYS", "nTightMuons_NOSYS"}
  );

  /* 
    ======================================================
        Combine jet variables into one vector
    ======================================================
  */
  // Function
  auto combineJetVariables = [](
      const std::vector<float>& pT,
      const std::vector<float>& e,
      const std::vector<float>& eta,
      const std::vector<float>& phi,
      const std::vector<char>& bTag
  ) {
      // Convert all variables to float
      std::vector<float> e_float(e.begin(), e.end());
      std::vector<float> eta_float(eta.begin(), eta.end());
      std::vector<float> phi_float(phi.begin(), phi.end());
      std::vector<float> bTag_float(bTag.begin(), bTag.end());

      const int size = pT.size();

      // Create a vector of indices to sort the other variables
      std::vector<int> indices(size);
      for (int i = 0; i < size; ++i) {
          indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), [&pT](int a, int b) {
          return pT[a] > pT[b];
      });

      // Create one large vector to store all the variables
      std::vector<float> all_variables;
      for (int i = 0; i < size; ++i) {
        if (abs(eta_float[indices[i]]) > 2.5) continue; // Skip jets with |eta| > 2.5
        all_variables.push_back(pT[indices[i]]);
        all_variables.push_back(e_float[indices[i]]);
        all_variables.push_back(eta_float[indices[i]]);
        all_variables.push_back(phi_float[indices[i]]);
        all_variables.push_back(bTag_float[indices[i]]);
      }

      return all_variables;
  };
  // Variable - Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_variables_NOSYS",
      combineJetVariables,
      {
        "jet_pt_GeV_NOSYS", "jet_e_GeV_NOSYS", "jet_eta", "jet_phi", "jet_GN2v01_FixedCutBEff_77_select"
      }
  );

  /* 
    ======================================================
        Count the number of jets from jet_variables_NOSYS
    ======================================================
  */
  // Function - Count jets
  auto countJets = [](
      const std::vector<float>& jet_variables
  ) {
      const int size = jet_variables.size()/5; // Each jet has 5 variables

      return size;
  };
  // Variable - Count Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nJets_NOSYS",
      countJets,
      {"jet_variables_NOSYS"}
  );

  /* 
    ======================================================
        Extract leading jet pt from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract leading jet pt
  auto leadingJetPt = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.empty()) return -999.0f; // Return -999 if no jets

      // The first element in the vector is the leading jet's pt
      return jet_variables[0];
  };
  // Variable - Leading Jet Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet0_pt_NOSYS",
      leadingJetPt,
      {"jet_variables_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading jet eta from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract leading jet eta
  auto leadingJetEta = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.empty()) return -999.0f; // Return -999 if no jets

      // The second element in the vector is the leading jet's eta
      return jet_variables[2];
  };
  // Variable - Leading Jet Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet0_eta_NOSYS",
      leadingJetEta,
      {"jet_variables_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading jet phi from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract leading jet phi
  auto leadingJetPhi = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.empty()) return -999.0f; // Return -999 if no jets

      // The third element in the vector is the leading jet's phi
      return jet_variables[3];
  };
  // Variable - Leading Jet Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet0_phi_NOSYS",
      leadingJetPhi,
      {"jet_variables_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading jet pt from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract subleading jet pt
  auto subleadingJetPt = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.size()/5 < 2) return -999.0f; // Return -999 if no second jet

      // The first element in the second jet's vector is the subleading jet's pt
      return jet_variables[5];
  };
  // Variable - Subleading Jet Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet1_pt_NOSYS",
      subleadingJetPt,
      {"jet_variables_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading jet eta from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract subleading jet eta
  auto subleadingJetEta = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.size()/5 < 2) return -999.0f; // Return -999 if no second jet

      // The second element in the second jet's vector is the subleading jet's eta
      return jet_variables[7];
  };
  // Variable - Subleading Jet Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet1_eta_NOSYS",
      subleadingJetEta,
      {"jet_variables_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading jet phi from jet_variables_NOSYS
    ======================================================
  */
  // Function - Extract subleading jet phi
  auto subleadingJetPhi = [](
      const std::vector<float>& jet_variables
  ) {
      if (jet_variables.size()/5 < 2) return -999.0f; // Return -999 if no second jet

      // The third element in the second jet's vector is the subleading jet's phi
      return jet_variables[8];
  };
  // Variable - Subleading Jet Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet1_phi_NOSYS",
      subleadingJetPhi,
      {"jet_variables_NOSYS"}
  );

  /* 
    ======================================================
        Number of b-tagged jets 
    ======================================================
  */ 
  // Function
  auto countBTaggedJets = [](
      const std::vector<float>& jet_variables
  ) {
      const int size = jet_variables.size(); // Each jet has 5 variables

      int count = 0;
      for (int i = 0; i < size; i += 5) { // Each jet has 5 variables
          if (jet_variables[i + 4] == 1.) { // Check b-tagging variable
              count++;
          }
      }
      return count;
  };
  // Variable - Number of b-tagged jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "nBjets_NOSYS",
      countBTaggedJets,
      {"jet_variables_NOSYS"}
  );

  /* 
    ======================================================
        Combine electron and muon variables into one vector 
    ======================================================
  */
  // Function
  auto combineElectronsAndMuons = [](
      const std::vector<float>& electron_variables,
      const std::vector<float>& muon_variables
  ) {
      // Construct a new vector with elements from 'a', then append 'b'
      std::vector<float> result(electron_variables.begin(), electron_variables.end());
      result.insert(result.end(), muon_variables.begin(), muon_variables.end());

      return result;
  };
  // Variable - Lepton
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lepton_variables_NOSYS",
      combineElectronsAndMuons,
      {"el_variables_NOSYS", "mu_variables_NOSYS"}
  );

  /* 
    ======================================================
        Sort by Pt lepton variables 
    ======================================================
  */
  // Function
  auto sortLeptonVariablesByPt = [](
      const std::vector<float>& lepton_variables
  ) {
      // Create a vector of pairs (pt, index)
      std::vector<std::pair<float, int>> pt_index_pairs;
      for (size_t i = 0; i < lepton_variables.size(); i += 9) {
          pt_index_pairs.emplace_back(lepton_variables[i], i / 9);
      }

      // Sort the pairs by pt in descending order
      std::sort(pt_index_pairs.begin(), pt_index_pairs.end(),
                std::greater<std::pair<float, int>>());

      // Create a new vector to hold the sorted variables
      std::vector<float> sorted_lepton_variables;
      for (const auto& pair : pt_index_pairs) {
          int index = pair.second * 9; // Each lepton has 9 variables
          for (int j = 0; j < 9; ++j) {
              sorted_lepton_variables.push_back(lepton_variables[index + j]);
          }
      }

      return sorted_lepton_variables;
  };
  // Variable - Lepton variables sorted by pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lepton_variables_sortedByPt_NOSYS",
      sortLeptonVariablesByPt,
      {"lepton_variables_NOSYS"}
  );

  /* 
    ======================================================
        Extract leading lepton pt
    ======================================================
  */
  // Function
  auto extractLeadingLeptonPt = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.empty()) return -999.0f; // Return -999 if no leptons

      // The first element in the sorted lepton variables is the leading lepton's pt
      return lepton_variables[0];
  };
  // Variable - Leading lepton pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep0_pt_NOSYS",
      extractLeadingLeptonPt,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading lepton eta
    ======================================================
  */
  // Function
  auto extractLeadingLeptonEta = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.empty()) return -999.0f; // Return -999 if no leptons

      // The third element in the sorted lepton variables is the leading lepton's eta
      return lepton_variables[2];
  };
  // Variable - Leading lepton eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep0_eta_NOSYS",
      extractLeadingLeptonEta,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading lepton phi
    ======================================================
  */
  // Function
  auto extractLeadingLeptonPhi = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.empty()) return -999.0f; // Return -999 if no leptons

      // The fourth element in the sorted lepton variables is the leading lepton's phi
      return lepton_variables[3];
  };
  // Variable - Leading lepton phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep0_phi_NOSYS",
      extractLeadingLeptonPhi,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading lepton pt
    ======================================================
  */
  // Function
  auto extractSubleadingLeptonPt = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 2) return -999.0f; // Return -999 if no second lepton

      // The first element in the second lepton's vector is the subleading lepton's pt
      return lepton_variables[9];
  };
  // Variable - Subleading lepton pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep1_pt_NOSYS",
      extractSubleadingLeptonPt,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading lepton eta
    ======================================================
  */
  // Function
  auto extractSubleadingLeptonEta = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 2) return -999.0f; // Return -999 if no second lepton

      // The third element in the second lepton's vector is the subleading lepton's eta
      return lepton_variables[11];
  };
  // Variable - Subleading lepton eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep1_eta_NOSYS",
      extractSubleadingLeptonEta,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading lepton phi
    ======================================================
  */
  // Function
  auto extractSubleadingLeptonPhi = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 2) return -999.0f; // Return -999 if no second lepton

      // The fourth element in the second lepton's vector is the subleading lepton's phi
      return lepton_variables[12];
  };
  // Variable - Subleading lepton phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep1_phi_NOSYS",
      extractSubleadingLeptonPhi,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubleading lepton pt
    ======================================================
  */
  // Function
  auto extractSubsubleadingLeptonPt = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 3) return -999.0f; // Return -999 if no third lepton

      // The first element in the third lepton's vector is the subsubleading lepton's pt
      return lepton_variables[18];
  };
  // Variable - Subsubleading lepton pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep2_pt_NOSYS",
      extractSubsubleadingLeptonPt,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubleading lepton eta
    ======================================================
  */
  // Function
  auto extractSubsubleadingLeptonEta = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 3) return -999.0f; // Return -999 if no third lepton

      // The third element in the third lepton's vector is the subsubleading lepton's eta
      return lepton_variables[20];
  };
  // Variable - Subsubleading lepton eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep2_eta_NOSYS",
      extractSubsubleadingLeptonEta,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubleading lepton phi
    ======================================================
  */
  // Function
  auto extractSubsubleadingLeptonPhi = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 3) return -999.0f; // Return -999 if no third lepton

      // The fourth element in the third lepton's vector is the subsubleading lepton's phi
      return lepton_variables[21];
  };
  // Variable - Subsubleading lepton phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep2_phi_NOSYS",
      extractSubsubleadingLeptonPhi,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubsubleading lepton pt
    ======================================================
  */
  // Function
  auto extractSubsubsubleadingLeptonPt = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth lepton

      // The first element in the fourth lepton's vector is the subsubsubleading lepton's pt
      return lepton_variables[27];
  };
  // Variable - Subsubsubleading lepton pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep3_pt_NOSYS",
      extractSubsubsubleadingLeptonPt,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubsubleading lepton eta
    ======================================================
  */
  // Function
  auto extractSubsubsubleadingLeptonEta = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth lepton

      // The third element in the fourth lepton's vector is the subsubsubleading lepton's eta
      return lepton_variables[29];
  };
  // Variable - Subsubsubleading lepton eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep3_eta_NOSYS",
      extractSubsubsubleadingLeptonEta,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  /* 
    ======================================================
        Extract subsubsubleading lepton phi
    ======================================================
  */
  // Function
  auto extractSubsubsubleadingLeptonPhi = [](
      const std::vector<float>& lepton_variables
  ) {
      if (lepton_variables.size()/9 < 4) return -999.0f; // Return -999 if no fourth lepton

      // The fourth element in the fourth lepton's vector is the subsubsubleading lepton's phi
      return lepton_variables[30];
  };
  // Variable - Subsubsubleading lepton phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep3_phi_NOSYS",
      extractSubsubsubleadingLeptonPhi,
      {"lepton_variables_sortedByPt_NOSYS"}
  );

  /* 
    ======================================================
        Calculate the sum of lepton charges 
    ======================================================
  */
  // Function
  auto calculateSumOfLeptonCharges = [](
      const std::vector<float>& lepton_variables
  ) {
      const int size = lepton_variables.size();

      int sum_of_charges = 0;

      for (int i = 0; i < size; i += 9) {
          // Cast float as int to get the charge
          int charge = static_cast<int>(lepton_variables[i + 5]); // Charge is the 6th element in each lepton's variables
          // Add the charge to the sum
          sum_of_charges += charge;
      }
      return sum_of_charges;
  };
  // Variable - Sum of lepton charges
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lep_sum_charge_NOSYS",
      calculateSumOfLeptonCharges,
      {"lepton_variables_sortedByPt_NOSYS"}
  );

  /*
    ======================================================
        Store the 4 vectors of leptons as TLVs from the variables
    ======================================================
  */
  // Function
  auto storeLeptonTLVs = [](
      const std::vector<float>& lepton_variables
  ) {
      const int size = lepton_variables.size();

      std::vector<TLV> lepton_TLVs;
      for (int i = 0; i < size; i += 9) {
          lepton_TLVs.push_back(TLV{lepton_variables[i], // pT
                                    lepton_variables[i + 2], // eta
                                    lepton_variables[i + 3], // phi
                                    lepton_variables[i + 1]}); // energy
      }
      return lepton_TLVs;
  };
  // Variable - Lepton TLVs
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "lepton_TLV_NOSYS",
      storeLeptonTLVs,
      {"lepton_variables_sortedByPt_NOSYS"}
  );

  /* 
    ======================================================
        Extract lepton pdgID
    ======================================================
  */
  // Function
  auto extractLeptonPdgID = [](
      const std::vector<float>& lepton_variables
  ) {
      const int size = lepton_variables.size();

      std::vector<float> lepton_pdgID;
      for (int i = 0; i < size; i += 9) {
          lepton_pdgID.push_back(lepton_variables[i + 4]*lepton_variables[i + 5]); // pdgFlavour is the 5th element in each lepton's variables
      }
      return lepton_pdgID;
  };
  // Variable - Lepton pdgID
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sorted_lep_pdgID_NOSYS",
      extractLeptonPdgID,
      {"lepton_variables_sortedByPt_NOSYS"}
  );
  
  /*
    ======================================================
        Store the 4 vectors of jets as TLVs from the variables
    ======================================================
  */
  // Function
  auto storeJetTLVs = [](
      const std::vector<float>& jet_variables
  ) {
      const int size = jet_variables.size();

      std::vector<TLV> jet_TLVs;
      for (int i = 0; i < size; i += 5) {
          jet_TLVs.push_back(TLV{jet_variables[i], // pT
                                 jet_variables[i + 2], // eta
                                 jet_variables[i + 3], // phi
                                 jet_variables[i + 1]}); // energy
      }
      return jet_TLVs;
  };
  // Variable - Jet TLVs
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "jet_TLV_NOSYS",
      storeJetTLVs,
      {"jet_variables_NOSYS"}
  );

  /* 
    ======================================================
        Store the 4 vectors of btagged jets as TLVs
    ======================================================
  */
  // Function
  auto storeBtaggedJetTLVs = [](
      const std::vector<float>& jet_variables
  ) {
      const int size = jet_variables.size();

      std::vector<TLV> btagged_jet_TLVs;
      for (int i = 0; i < size; i += 5) {
          if (jet_variables[i + 4] == 1.) { // Check if the jet is b-tagged
              btagged_jet_TLVs.push_back(TLV{jet_variables[i], // pT
                                              jet_variables[i + 2], // eta
                                              jet_variables[i + 3], // phi
                                              jet_variables[i + 1]}); // energy
          }
      }
      return btagged_jet_TLVs;
  };
  // Variable - Btagged Jet TLVs
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet_TLV_NOSYS",
      storeBtaggedJetTLVs,
      {"jet_variables_NOSYS"}
  );

  /* 
    ======================================================
        Extract leading b-tagged jet pt from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract leading b-tagged jet pt
  auto leadingBtaggedJetTLVPt = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.empty()) return -999.0f; // Return -999 if no b-tagged jets

      float pt = bjet_TLVs[0].Pt(); // Get the pt of the leading b-tagged jet

      // Return the pt of the leading b-tagged jet
      return pt;
  };
  // Variable - Leading Btagged Jet TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet0_pt_NOSYS",
      leadingBtaggedJetTLVPt,
      {"bjet_TLV_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading b-tagged jet eta from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract leading b-tagged jet eta
  auto leadingBtaggedJetTLVEta = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.empty()) return -999.0f; // Return -999 if no b-tagged jets

      float eta = bjet_TLVs[0].Eta(); // Get the eta of the leading b-tagged jet

      // Return the eta of the leading b-tagged jet
      return eta;
  };
  // Variable - Leading Btagged Jet TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet0_eta_NOSYS",
      leadingBtaggedJetTLVEta,
      {"bjet_TLV_NOSYS"}
  );
  /* 
    ======================================================
        Extract leading b-tagged jet phi from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract leading b-tagged jet phi
  auto leadingBtaggedJetTLVPhi = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.empty()) return -999.0f; // Return -999 if no b-tagged jets

      float phi = bjet_TLVs[0].Phi(); // Get the phi of the leading b-tagged jet

      // Return the phi of the leading b-tagged jet
      return phi;
  };
  // Variable - Leading Btagged Jet TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet0_phi_NOSYS",
      leadingBtaggedJetTLVPhi,
      {"bjet_TLV_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading b-tagged jet pt from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract subleading b-tagged jet pt
  auto subleadingBtaggedJetTLVPt = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second b-tagged jet

      float pt = bjet_TLVs[1].Pt(); // Get the pt of the subleading b-tagged jet

      // Return the pt of the subleading b-tagged jet
      return pt;
  };
  // Variable - Subleading Btagged Jet TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet1_pt_NOSYS",
      subleadingBtaggedJetTLVPt,
      {"bjet_TLV_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading b-tagged jet eta from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract subleading b-tagged jet eta
  auto subleadingBtaggedJetTLVEta = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second b-tagged jet

      float eta = bjet_TLVs[1].Eta(); // Get the eta of the subleading b-tagged jet

      // Return the eta of the subleading b-tagged jet
      return eta;
  };
  // Variable - Subleading Btagged Jet TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet1_eta_NOSYS",
      subleadingBtaggedJetTLVEta,
      {"bjet_TLV_NOSYS"}
  );
  /* 
    ======================================================
        Extract subleading b-tagged jet phi from bjet_TLV_NOSYS
    ======================================================
  */
  // Function - Extract subleading b-tagged jet phi
  auto subleadingBtaggedJetTLVPhi = [](
      const std::vector<TLV>& bjet_TLVs
  ) {
      if (bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second b-tagged jet

      float phi = bjet_TLVs[1].Phi(); // Get the phi of the subleading b-tagged jet

      // Return the phi of the subleading b-tagged jet
      return phi;
  };
  // Variable - Subleading Btagged Jet TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "bjet1_phi_NOSYS",
      subleadingBtaggedJetTLVPhi,
      {"bjet_TLV_NOSYS"}
  );

  /* 
    ======================================================
        Store the 4 vectors of non-btagged jets as TLVs
    ======================================================
  */
  // Function
  auto storeNonBtaggedJetTLVs = [](
      const std::vector<float>& jet_variables
  ) {
      const int size = jet_variables.size();

      std::vector<TLV> non_btagged_jet_TLVs;
      for (int i = 0; i < size; i += 5) {
          if (jet_variables[i + 4] == 0.) { // Check if the jet is non-b-tagged
              non_btagged_jet_TLVs.push_back(TLV{jet_variables[i], // pT
                                                jet_variables[i + 2], // eta
                                                jet_variables[i + 3], // phi
                                                jet_variables[i + 1]}); // energy
          }
      }
      return non_btagged_jet_TLVs;
  };
  // Variable - Non-Btagged Jet TLVs
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet_TLV_NOSYS",
      storeNonBtaggedJetTLVs,
      {"jet_variables_NOSYS"}
  );

  /* 
    ============================================================
        Extract leading non-b-tagged jet pt from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract leading non-b-tagged jet pt
  auto leadingNonBtaggedJetTLVPt = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.empty()) return -999.0f; // Return -999 if no non-b-tagged jets

      float pt = non_bjet_TLVs[0].Pt(); // Get the pT of the leading non-b-tagged jet

      // Return the pT of the leading non-b-tagged jet
      return pt;
  };
  // Variable - Leading Non-Btagged Jet TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet0_pt_NOSYS",
      leadingNonBtaggedJetTLVPt,
      {"non_bjet_TLV_NOSYS"}
  );

  /* 
    ============================================================
        Extract leading non-b-tagged jet eta from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract leading non-b-tagged jet eta
  auto leadingNonBtaggedJetTLVEta = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.empty()) return -999.0f; // Return -999 if no non-b-tagged jets

      float eta = non_bjet_TLVs[0].Eta(); // Get the eta of the leading non-b-tagged jet

      // Return the eta of the leading non-b-tagged jet
      return eta;
  };
  // Variable - Leading Non-Btagged Jet TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet0_eta_NOSYS",
      leadingNonBtaggedJetTLVEta,
      {"non_bjet_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract leading non-b-tagged jet phi from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract leading non-b-tagged jet phi
  auto leadingNonBtaggedJetTLVPhi = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.empty()) return -999.0f; // Return -999 if no non-b-tagged jets

      float phi = non_bjet_TLVs[0].Phi(); // Get the phi of the leading non-b-tagged jet

      // Return the phi of the leading non-b-tagged jet
      return phi;
  };
  // Variable - Leading Non-Btagged Jet TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet0_phi_NOSYS",
      leadingNonBtaggedJetTLVPhi,
      {"non_bjet_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract subleading non-b-tagged jet pt from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract subleading non-b-tagged jet pt
  auto subleadingNonBtaggedJetTLVPt = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second non-b-tagged jet

      float pt = non_bjet_TLVs[1].Pt(); // Get the pT of the subleading non-b-tagged jet

      // Return the pT of the subleading non-b-tagged jet
      return pt;
  };
  // Variable - Subleading Non-Btagged Jet TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet1_pt_NOSYS",
      subleadingNonBtaggedJetTLVPt,
      {"non_bjet_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract subleading non-b-tagged jet eta from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract subleading non-b-tagged jet eta
  auto subleadingNonBtaggedJetTLVEta = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second non-b-tagged jet

      float eta = non_bjet_TLVs[1].Eta(); // Get the eta of the subleading non-b-tagged jet

      // Return the eta of the subleading non-b-tagged jet
      return eta;
  };
  // Variable - Subleading Non-Btagged Jet TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet1_eta_NOSYS",
      subleadingNonBtaggedJetTLVEta,
      {"non_bjet_TLV_NOSYS"}
  );
  /*
    ============================================================
        Extract subleading non-b-tagged jet phi from non_bjet_TLV_NOSYS
    ============================================================
  */
  // Function - Extract subleading non-b-tagged jet phi
  auto subleadingNonBtaggedJetTLVPhi = [](
      const std::vector<TLV>& non_bjet_TLVs
  ) {
      if (non_bjet_TLVs.size() < 2) return -999.0f; // Return -999 if no second non-b-tagged jet

      float phi = non_bjet_TLVs[1].Phi(); // Get the phi of the subleading non-b-tagged jet

      // Return the phi of the subleading non-b-tagged jet
      return phi;
  };
  // Variable - Subleading Non-Btagged Jet TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "non_bjet1_phi_NOSYS",
      subleadingNonBtaggedJetTLVPhi,
      {"non_bjet_TLV_NOSYS"}
  );

  /*
    ===================================================================================
        Calculate the mass of all possible opposite sign same flavour leptons
    ===================================================================================
  */
  // Function
  auto calculateOSSFMass = [](
      const std::vector<TLV>& leptons,
      const std::vector<float>& lep_pdgID
  ) {
      const int nLeptons = leptons.size();

      std::vector<float> ossf_mass_flat;
      if (nLeptons == 0) { return ossf_mass_flat; }
      else {
        for(int lep_i = 0; lep_i<nLeptons; lep_i++) {
          for(int lep_j = 0; lep_j<nLeptons; lep_j++) {
            if (lep_i >= lep_j) { continue; }
            if (lep_pdgID.at(lep_i) + lep_pdgID.at(lep_j) != 0.) { continue; } // OSSF

            TLV pair = leptons.at(lep_i) + leptons.at(lep_j);
            float mass = pair.M();

            ossf_mass_flat.push_back(mass); // Store the mass of the OSSF
            // Store the indices of the leptons involved in the OSSF pair
            ossf_mass_flat.push_back(static_cast<float>(lep_i)); // Index of first le
            ossf_mass_flat.push_back(static_cast<float>(lep_j)); // Index of second lepton
          }
        }
      }

      return ossf_mass_flat;
  };
  // Variable - OSSF Mass
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ossf_mass_NOSYS",
      calculateOSSFMass,
      {"lepton_TLV_NOSYS", "sorted_lep_pdgID_NOSYS"}
  );
  // Variable - Number of OSSF pairs
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "nOSSF_NOSYS",
      "static_cast<int>(ossf_mass_NOSYS.size() / 3)" //
  );

  /*
    ===================================================================================
        Create OSSF Flag for selection cut purposes
    ===================================================================================
  */
  // Function
  auto ossfFlag = [](
      const std::vector<float>& ossf_mass_flat
  ) {
      const int size = ossf_mass_flat.size();

      int result = 1; // Default to 1 (true)
      for (int i = 0; i < size; i+= 3) {
          if (ossf_mass_flat.at(i) < 10.) { // If any mass is less than 10 GeV, set result to 0 (false)
              result = 0;
              break;
          }
      }
      return result;
  };
  // Variable - OSSF Flag
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ossfFlag_NOSYS",
      ossfFlag,
      {"ossf_mass_NOSYS"}
  );

  /*
    ===================================================================================
        Determine if the massOSSF sits within th Z window
    ===================================================================================
  */
  // Function
  auto massOSSFInZWindow = [](
      const std::vector<float>& ossf_mass_flat
  ) {
      const double z_mass = 91.1876;
      const double minMass = z_mass - 10.;
      const double maxMass = z_mass + 10.;

      const int size = ossf_mass_flat.size();
      std::vector<float> result; // Initialize result vector 
      for (int i = 0; i < size; i+= 3) {
          if (ossf_mass_flat[i] >= minMass && ossf_mass_flat[i] <= maxMass) {
              result.push_back(ossf_mass_flat[i]); // Store the mass if it is within the Z window
          }
      }
      return result; // Return the vector of masses within the Z window
  };
  // Variable - Mass OSSF in Z Window
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "massOSSF_inZWindow_NOSYS",
      massOSSFInZWindow,
      {"ossf_mass_NOSYS"}
  );
  // Variable - Number of OSSF pairs in Z window
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "nOSSF_inZWindow_NOSYS",
      "static_cast<int>(massOSSF_inZWindow_NOSYS.size())"
  );

  /* 
    ============================================================
        Calculate difference mass of opposite sign same flavour leptons
    ============================================================
  */
  // Function
  auto calculateDifferenceMassOSSF = [](
      const std::vector<TLV>& leptons, 
      const std::vector<float>& lep_pdgID
  )  {
      const float z_mass = 91.1876;
      const int nLeptons = leptons.size();
  
      std::vector<float> difference_flat;
      if (nLeptons == 0) { return difference_flat; }
      else {
        for(int lep_i = 0; lep_i<nLeptons; lep_i++) {
          for(int lep_j = 0; lep_j<nLeptons; lep_j++) {
            if (lep_i >= lep_j) { continue; }
            if (lep_pdgID.at(lep_i) + lep_pdgID.at(lep_j) != 0.) { continue; } // OSSF

            TLV pair = leptons.at(lep_i) + leptons.at(lep_j);

            float mass = pair.M();
            float diff_mass = std::abs(pair.M() - z_mass);
            float index_i = static_cast<float>(lep_i);
            float index_j = static_cast<float>(lep_j);

            // Store the difference, indices of leptons
            int already_stored_flag = 0;
            int size_difference = difference_flat.size() / 3;
            std::vector<int> erase_indices;
            for (int i = 0; i < size_difference; i++) {
              float prev_mass = difference_flat[i*3 + 0];
              float prev_diff_mass = std::abs(prev_mass - z_mass);
              float prev_index_i = difference_flat[i*3 + 1];
              float prev_index_j = difference_flat[i*3 + 2];
              if (prev_index_i == index_i || prev_index_i == index_j || prev_index_j == index_i || prev_index_j == index_j) {
                if (prev_diff_mass > diff_mass) {
                  erase_indices.push_back(i*3+0);
                  erase_indices.push_back(i*3+1);
                  erase_indices.push_back(i*3+2);
                }
                else {
                  already_stored_flag = 1;
                }
              }
            }
            // Erase indices
            for (int i = erase_indices.size() - 1; i >= 0; i--) {
              difference_flat.erase(difference_flat.begin() + erase_indices[i]);
            }
            if (already_stored_flag == 1) { continue; }
            else {
              difference_flat.push_back(mass);
              difference_flat.push_back(index_i);
              difference_flat.push_back(index_j);
            }
          }
        }
      }
    
      return difference_flat;
  };
  // Variable  - Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "differenceMassOSSF_NOSYS",
      calculateDifferenceMassOSSF,
      {"lepton_TLV_NOSYS","sorted_lep_pdgID_NOSYS"}
  );

  /* 
    ============================================================
        Create differenceMassOSSF Flag for selection cut purposes
    ============================================================
  */
  // Function
  auto differenceMassOSSFFlag = [](
      const std::vector<float>& mass_flat
  ) {
      const int size = mass_flat.size() / 3;

      int result = 1;   
      
      for (int mass_i = 0; mass_i<size; mass_i++) {
        if (mass_flat.at(mass_i*3 + 0) <10.)  {
          result = 0;
          break;
        }
      }

      return result;
  };
  // Variable - OSSF Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "differenceMassOSSFFlag_NOSYS",
      differenceMassOSSFFlag,
      {"differenceMassOSSF_NOSYS"}
  );

  /* 
    ============================================================
        Tie massOSSF with leptons
    ============================================================
  */
  // Function
  auto tieMassOSSFWithLeptons = [](
      const std::vector<float>& mass_flat,
      const std::vector<TLV>& leptons
  ) {
      const double z_mass = 91.1876;
      const double minMass = z_mass - 10.;
      const double maxMass = z_mass + 10.;

      const int size = mass_flat.size() / 3;
      const int nLeptons = leptons.size();

      std::vector<float> result(nLeptons);
      for (int i = 0; i < size; i++) {
        if (mass_flat[i*3 + 0] < minMass || mass_flat[i*3 + 0] > maxMass) { continue; }

        int index_i = static_cast<int>(mass_flat[i*3 + 1]);
        int index_j = static_cast<int>(mass_flat[i*3 + 2]);

        //if (index_i < 0 || index_j < 0) { continue; } // Skip if indices are invalid
        //if (index_i >= 4 || index_j >= 4) { continue; } // Skip if indices are out of bounds

        result.at(index_i) = mass_flat[i*3 + 0];
        result.at(index_j) = mass_flat[i*3 + 0];
      }

      return result;
  };
  // Variables - OSSF Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "massesOSSF_NOSYS",
      tieMassOSSFWithLeptons,
      {"differenceMassOSSF_NOSYS", "lepton_TLV_NOSYS"}
  );

  /* 
    ============================================================
        Calculate Z candidates
    ============================================================
  */
  // Function
  auto extractMassOSSF = [](
    const std::vector<float>& mass_flat
  ) {
      std::vector<float> massOSSF;
      const int size = mass_flat.size() / 3;
      for (int i = 0; i < size; i++) {
        massOSSF.push_back(mass_flat[i*3 + 0]);
      }
      return massOSSF;
  };
  // Variable - OSSF Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "extractMassOSSF_NOSYS",
      extractMassOSSF,
      {"differenceMassOSSF_NOSYS"}
  );

  /* 
    ============================================================
        Calculate Z candidates
    ============================================================
  */
  // Function
  auto calculateZCandidates = [](
    const std::vector<float>& mass
  ) {
      const double z_mass = 91.1876;
      const double minMass = z_mass - 10.;
      const double maxMass = z_mass + 10.;

      const int size = mass.size();

      std::vector<float> zCandidates;
      for (int i = 0; i < size; i++) {
        if (minMass > mass.at(i) || maxMass < mass.at(i)) { continue; }
        zCandidates.push_back(mass.at(i));
      }

      return zCandidates;
  };
  // Variable - Z Candidates
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidates_NOSYS",
      calculateZCandidates,
      {"extractMassOSSF_NOSYS"}
  );

  /* 
    ============================================================
        Calculate the Z candidate mutliplicity
    ============================================================
  */
  mainNode = MainFrame::systematicStringDefine(
    mainNode,
    "nZCandidates_NOSYS",
    "static_cast<int>(ZCandidates_NOSYS.size())"
  );

  /* 
    ============================================================
        Create zCandidates off the back of zMassOSSF
    ============================================================
  */
  // Function
  auto zCandidatesFlag = [](
      const std::vector<float>& mass,
      const std::vector<TLV>& leptons
  ) {
      const int size = mass.size();
      const int nLeptons = leptons.size();

      std::vector<int> result(nLeptons);
      if (size == 0) {
        return result; // Return empty vector if no mass values
      }
      for (int mass_i = 0; mass_i<size; mass_i++) {
        if (result.at(mass_i) == 1) { continue; } // Already set as a candidate
        for (int mass_j = 0; mass_j<size; mass_j++) {
          if (mass_i >= mass_j) { continue; }
          if (mass.at(mass_i) == 0. || mass.at(mass_j) == 0.) { continue; } // Skip if mass is zero
          if (result.at(mass_j) == 1) { continue; } // Already set as a candidate

          if (mass.at(mass_i) == mass.at(mass_j)) { 
            result.at(mass_i) = 1;
            result.at(mass_j) = 1;
            break;
          }
        }
      }

      return result;
  };
  // Variable - Z Candidates
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidateFlag_NOSYS",
      zCandidatesFlag,
      {"massesOSSF_NOSYS", "lepton_TLV_NOSYS"}
  );

  /*
    ===============================================================
      Create Z candidate 4 vector
    ===============================================================
  */
  // Function
  auto createZCandidate4Vector = [](
      const std::vector<TLV>& lepton,
      const std::vector<float>& mass,
      const std::vector<int>& zCandidateFlag
  ) {
      const int size = lepton.size();

      std::vector<TLV> result;
      for (int i=0; i<size; i++){
        if (zCandidateFlag.at(i) == 0) { continue; } //Not a Z candidate
        
        for (int j=0; j<size; j++) {
          if (i>=j) {continue; }
          if (zCandidateFlag.at(j) == 0) { continue; }
          if (mass.at(i) != mass.at(j)) { continue; }
           
          result.push_back(lepton.at(i) + lepton.at(j));
        }
      }

      return result;
  };
  // Variables - Z Candidates
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate_TLV_NOSYS",
      createZCandidate4Vector,
      {"lepton_TLV_NOSYS","massesOSSF_NOSYS","ZCandidateFlag_NOSYS"}
  );
  // Function 
  auto sortByPt = [](
      const std::vector<TLV>& tlv
  ) {
      const int size = tlv.size();

      std::vector<TLV> result;
      for (int i=0; i<size; i++) {
        result.push_back(tlv.at(i));
      }
      // sort them based on pT
      std::sort(result.begin(), result.end(), [](const TLV& v1, const TLV& v2) { return v1.pt() > v2.pt(); });

      return result;
  };
  // Variable - Z Candidates
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sorted_ZCandidate_TLV_NOSYS",
      sortByPt,
      {"ZCandidate_TLV_NOSYS"}
  );

  /*
    ============================================================
        Extract leading Z candidate pt
    ============================================================
  */
  // Function - Extract leading Z candidate pt
  auto leadingZCandidateTLVPt = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.empty()) return -999.0f; // Return -999 if no Z candidates

      float pt = zCandidate_TLVs[0].Pt(); // Get the pt of the leading Z candidate

      return pt;
  };
  // Variable - Leading Z Candidate TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate0_pt_NOSYS",
      leadingZCandidateTLVPt,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract leading Z candidate eta
    ============================================================
  */
  // Function - Extract leading Z candidate eta
  auto leadingZCandidateTLVEta = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.empty()) return -999.0f; // Return -999 if no Z candidates

      float eta = zCandidate_TLVs[0].Eta(); // Get the eta of the leading Z candidate

      return eta;
  };
  // Variable - Leading Z Candidate TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate0_eta_NOSYS",
      leadingZCandidateTLVEta,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract leading Z candidate phi
    ============================================================
  */
  // Function - Extract leading Z candidate phi
  auto leadingZCandidateTLVPhi = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.empty()) return -999.0f; // Return -999 if no Z candidates

      float phi = zCandidate_TLVs[0].Phi(); // Get the phi of the leading Z candidate

      return phi;
  };
  // Variable - Leading Z Candidate TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate0_phi_NOSYS",
      leadingZCandidateTLVPhi,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract leading Z candidate mass
    ============================================================
  */
  // Function - Extract leading Z candidate mass
  auto leadingZCandidateTLVMass = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.empty()) return -999.0f; // Return -999 if no Z candidates

      float mass = zCandidate_TLVs[0].M(); // Get the mass of the leading Z candidate

      return mass;
  };
  // Variable - Leading Z Candidate TLV Mass
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate0_mass_NOSYS",
      leadingZCandidateTLVMass,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );

  /* 
    ============================================================
        Extract subleading Z candidate pt
    ============================================================
  */
  // Function - Extract subleading Z candidate pt
  auto subleadingZCandidateTLVPt = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.size() < 2) return -999.0f; // Return -999 if no second Z candidate

      float pt = zCandidate_TLVs[1].Pt(); // Get the pt of the subleading Z candidate

      return pt;
  };
  // Variable - Subleading Z Candidate TLV Pt
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate1_pt_NOSYS",
      subleadingZCandidateTLVPt,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract subleading Z candidate eta
    ============================================================
  */
  // Function - Extract subleading Z candidate eta
  auto subleadingZCandidateTLVEta = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.size() < 2) return -999.0f; // Return -999 if no second Z candidate

      float eta = zCandidate_TLVs[1].Eta(); // Get the eta of the subleading Z candidate

      return eta;
  };
  // Variable - Subleading Z Candidate TLV Eta
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate1_eta_NOSYS",
      subleadingZCandidateTLVEta,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract subleading Z candidate phi
    ============================================================
  */
  // Function - Extract subleading Z candidate phi
  auto subleadingZCandidateTLVPhi = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.size() < 2) return -999.0f; // Return -999 if no second Z candidate

      float phi = zCandidate_TLVs[1].Phi(); // Get the phi of the subleading Z candidate

      return phi;
  };
  // Variable - Subleading Z Candidate TLV Phi
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate1_phi_NOSYS",
      subleadingZCandidateTLVPhi,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );
  /* 
    ============================================================
        Extract subleading Z candidate mass
    ============================================================
  */
  // Function - Extract subleading Z candidate mass
  auto subleadingZCandidateTLVMass = [](
      const std::vector<TLV>& zCandidate_TLVs
  ) {
      if (zCandidate_TLVs.size() < 2) return -999.0f; // Return -999 if no second Z candidate

      float mass = zCandidate_TLVs[1].M(); // Get the mass of the subleading Z candidate

      return mass;
  };
  // Variable - Subleading Z Candidate TLV Mass
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ZCandidate1_mass_NOSYS",
      subleadingZCandidateTLVMass,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );

  /* 
    ============================================================
        Specific additional identifier for tWZOFSR | tWZSFSR
    ============================================================
  */
  // Function
  auto signalRegionType = [](
      const std::vector<int>& candidatesFlag,
      const std::vector<float>& leptonPdgID
  ) {
      const int size = candidatesFlag.size();

      int result = 5;    
      for (int cand_i = 0; cand_i<size; cand_i++) {
        if (result != 5) { break; }
        if (candidatesFlag.at(cand_i) == 1) { continue; }

        for (int cand_j = 0; cand_j<size; cand_j++) {
          if (candidatesFlag.at(cand_j) == 1) { continue; }
          if (cand_i >= cand_j) { continue; }

          if (leptonPdgID.at(cand_i) + leptonPdgID.at(cand_j) == 0) {
            result = 1; // Same Flavour (SF)
            break;
          }
          else {
            result = 0; // Opposite Flavour (OF)
            break;
          }
        }       
      }

      return result;
  };
  // Variable - OF or SF
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "signalRegionType_NOSYS",
      signalRegionType,
      {"ZCandidateFlag_NOSYS","sorted_lep_pdgID_NOSYS"}
  );

  /*
    =========================================
      Calculate the Scalar sum of pt
    =========================================
  */
  // Function
  auto calculateScalarSumOfPt = [](
      const std::vector<TLV>& tlv
  ) {
      const std::size_t size = tlv.size();

      float result = 0.;    
      for (std::size_t tlv_i = 0; tlv_i<size; tlv_i++) {
        result += tlv.at(tlv_i).Pt();
      }

      return result;
  };
  // Variable
  // Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "HT_NOSYS",
      calculateScalarSumOfPt,
      {"jet_TLV_NOSYS"}
  );
  // Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "LT_NOSYS",
      calculateScalarSumOfPt,
      {"lepton_TLV_NOSYS"}
  );
  // Bjets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sumBjetPt_NOSYS",
      calculateScalarSumOfPt,
      {"bjet_TLV_NOSYS"}
  );
  // Z Candidates
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "sumZCandidatePt_NOSYS",
      calculateScalarSumOfPt,
      {"sorted_ZCandidate_TLV_NOSYS"}
  );

  /*
    ===================================================
      Calculate the scalar sum of pt from two object
    ===================================================
  */
  // Function
  auto calculateScalarSumOfPtTwoObjects = [](
      const std::vector<TLV>& tlv1,
      const std::vector<TLV>& tlv2
  ) {
      const std::size_t size1 = tlv1.size();
      const std::size_t size2 = tlv2.size();

      float result = 0.;    
      for (std::size_t tlv_i = 0; tlv_i<size1; tlv_i++) {
        result += tlv1.at(tlv_i).Pt();
      }

      for (std::size_t tlv_j = 0; tlv_j<size2; tlv_j++) {
        result += tlv2.at(tlv_j).Pt();
      }

      return result;
  };
  // Variables - Jets & Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "ST_NOSYS",
      calculateScalarSumOfPtTwoObjects,
      {"jet_TLV_NOSYS","lepton_TLV_NOSYS"}
  );

  /*
    ===============================================================
      Calculate the Scalar sum of pt from Jets, Leptons and MET
    ===============================================================
  */
  // Function
  auto calculateSMT = [](
      const std::vector<TLV>& tlv1,
      const std::vector<TLV>& tlv2,
      const float& met
  ) {
      const std::size_t size1 = tlv1.size();
      const std::size_t size2 = tlv2.size();

      float result = 0.;    
      for (std::size_t tlv_i = 0; tlv_i<size1; tlv_i++) {
        result += tlv1.at(tlv_i).Pt();
      }

      for (std::size_t tlv_j = 0; tlv_j<size2; tlv_j++) {
        result += tlv2.at(tlv_j).Pt();
      }

      result += met;

      return result;
  };
  // Variables - Jets & Leptons & MET
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "SMT_NOSYS",
      calculateSMT,
      {"jet_TLV_NOSYS","lepton_TLV_NOSYS","met_met_GeV_NOSYS"}
  );

  /*
    ===============================================================
      Calculate 4 vector and components of the 4 lepton system
    ===============================================================
  */
  // Function
  auto combineLeptonTLV = [](
      const std::vector<TLV>& tlv
  ) {
      const int size = tlv.size();

      TLV result;
      for (int i = 0; i < size; i++) {
          result += tlv.at(i);
      }

      return result;
  };
  // Variables - Leptons
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "llll_TLV_NOSYS",
      combineLeptonTLV,
      {"lepton_TLV_NOSYS"}
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_pt_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_eta_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_phi_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.Phi())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_mass_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.M())"
  );

  /* 
  ======================================================
    Event truth flavour
  ======================================================
  */
  // Function
  auto eventTruthFlavour = [](
    const std::vector<int>& truthFlavour
  ){
    const int size = truthFlavour.size();

    int result = 0;
    for (int i=0; i<size; i++) {
      if (truthFlavour.at(i) == 5) { result = 5; break; }
    }
    if (result == 5) { return result; }

    for (int i=0; i<size; i++){ 
      if (truthFlavour.at(i) == 4) { result = 4; break; }
    }

    return result;
  };
  // Variable
  if (id.dsid() != 0) {
    mainNode = MainFrame::systematicDefine(
      mainNode,
      "eventTruthFlavour_NOSYS",
      eventTruthFlavour,
      {"jet_HadronConeExclTruthLabelID"}
    );
  }

  /*
    ===============================================================
      DeltaPhi betweeen different leading Leptons and MET
    ===============================================================
  */
  // Variable - Leptons and MET
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPhi_Lep0_MET_NOSYS",
      "static_cast<float>(std::abs(lep0_phi_NOSYS - met_phi_NOSYS))"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPhi_Lep1_MET_NOSYS",
      "static_cast<float>(std::abs(lep1_phi_NOSYS - met_phi_NOSYS))"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPhi_Lep2_MET_NOSYS",
      "static_cast<float>(std::abs(lep2_phi_NOSYS - met_phi_NOSYS))"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPhi_Lep3_MET_NOSYS",
      "static_cast<float>(std::abs(lep3_phi_NOSYS - met_phi_NOSYS))"
  );

  /*
    ===============================================================
      Combine the 4 lepton system with other particles
    ===============================================================
  */
  // 4Lepton + Leading order Bjet
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_Bjet0_TLV_NOSYS",
      "static_cast<ROOT::Math::PtEtaPhiEVector>((llll_TLV_NOSYS + bjet_TLV_NOSYS.at(0)))"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_Bjet0_pt_NOSYS",
      "static_cast<float>((llll_TLV_NOSYS + bjet_TLV_NOSYS.at(0)).Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_Bjet0_eta_NOSYS",
      "static_cast<float>((llll_TLV_NOSYS + bjet_TLV_NOSYS.at(0)).Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_Bjet0_mass_NOSYS",
      "static_cast<float>((llll_TLV_NOSYS + bjet_TLV_NOSYS.at(0)).M())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPt_llll_Bjet0_NOSYS",
      "static_cast<float>(std::abs(llll_TLV_NOSYS.Pt() - bjet_TLV_NOSYS.at(0).Pt()))"
  );
  // Leading Z Candidate + 4Lepton & Leading Bjet system
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "deltaPhi_ZCandidate0_llllBjet0_NOSYS",
      "static_cast<float>(std::abs(llll_Bjet0_TLV_NOSYS.Phi() - ZCandidate0_phi_NOSYS))"
  );
  // MET + 4Lepton & Leading Bjet system
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llllBjet0_MET_pt_NOSYS",
      "static_cast<float>(llll_Bjet0_TLV_NOSYS.Pt() + met_met_GeV_NOSYS)"
  );
  // MET + 4Lepton system
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_MET_pt_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.Pt() + met_met_GeV_NOSYS)"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "llll_MET_mass_NOSYS",
      "static_cast<float>(llll_TLV_NOSYS.M() + met_met_GeV_NOSYS)"
  );
  // Leptons + MET
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "l_l_l_l_MET_pt_NOSYS",
      "static_cast<float>(LT_NOSYS + met_met_GeV_NOSYS)"
  );
  // Leptons + MET + Leading Bjet
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "l_l_l_l_MET_Bjet0_pt_NOSYS",
      "static_cast<float>(LT_NOSYS + met_met_GeV_NOSYS + bjet_TLV_NOSYS.at(0).Pt())"
  );
  // MET + Leading Bjet
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "MET_Bjet0_pt_NOSYS",
      "static_cast<float>(met_met_GeV_NOSYS + bjet_TLV_NOSYS.at(0).Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "MET_Bjet0_mass_NOSYS",
      "static_cast<float>(met_met_GeV_NOSYS + bjet_TLV_NOSYS.at(0).M())"
  );

  /*
    ===============================================================
      Combine all leptons and all jets
    ===============================================================
  */
  // Function
  auto combineTwoTLVs = [](
      const std::vector<TLV>& lep_tlv,
      const std::vector<TLV>& jet_tlv
  ) {
      const int lep_size = lep_tlv.size();
      const int jet_size = jet_tlv.size();

      TLV result;
      for (int i=0; i<lep_size; i++) {
        result += lep_tlv.at(i);
      }
      for (int j=0; j<jet_size; j++) {
        result += jet_tlv.at(j);
      }

      return result;
  };
  // Variables - Leptons & Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "allLep_allJet_Sys_TLV_NOSYS",
      combineTwoTLVs,
      {"lepton_TLV_NOSYS","jet_TLV_NOSYS"}
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allJet_Sys_pt_NOSYS",
      "static_cast<float>(allLep_allJet_Sys_TLV_NOSYS.Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allJet_Sys_eta_NOSYS",
      "static_cast<float>(allLep_allJet_Sys_TLV_NOSYS.Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allJet_Sys_phi_NOSYS",
      "static_cast<float>(allLep_allJet_Sys_TLV_NOSYS.Phi())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allJet_Sys_mass_NOSYS",
      "static_cast<float>(allLep_allJet_Sys_TLV_NOSYS.M())"
  );

  /*
    ===============================================================
      Combine all jets
    ===============================================================
  */
  // Function
  auto combineTLV = [](
      const std::vector<TLV>& tlv
  ) {
      const int size = tlv.size();

      TLV result;
      for (int i=0; i<size; i++) {
        result += tlv.at(i);
      }

      return result;
  };
  //Variable - Jets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "allJets_Sys_TLV_NOSYS",
      combineTLV,
      {"jet_TLV_NOSYS"}
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allJets_Sys_pt_NOSYS",
      "static_cast<float>(allJets_Sys_TLV_NOSYS.Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allJets_Sys_eta_NOSYS",
      "static_cast<float>(allJets_Sys_TLV_NOSYS.Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allJets_Sys_phi_NOSYS",
      "static_cast<float>(allJets_Sys_TLV_NOSYS.Phi())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allJets_Sys_mass_NOSYS",
      "static_cast<float>(allJets_Sys_TLV_NOSYS.M())"
  );

  /*
    ===============================================================
      Combine all Bjets
    ===============================================================
  */
  //Variable - Bjets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "allBjets_Sys_TLV_NOSYS",
      combineTLV,
      {"bjet_TLV_NOSYS"}
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allBjets_Sys_pt_NOSYS",
      "static_cast<float>(allBjets_Sys_TLV_NOSYS.Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allBjets_Sys_eta_NOSYS",
      "static_cast<float>(allBjets_Sys_TLV_NOSYS.Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allBjets_Sys_phi_NOSYS",
      "static_cast<float>(allBjets_Sys_TLV_NOSYS.Phi())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allBjets_Sys_mass_NOSYS",
      "static_cast<float>(allBjets_Sys_TLV_NOSYS.M())"
  );

  /*
    ===============================================================
      Combine all leptons and all bjets
    ===============================================================
  */
  //Variable - Leptons & Bjets
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "allLep_allBjets_Sys_TLV_NOSYS",
      combineTwoTLVs,
      {"lepton_TLV_NOSYS","bjet_TLV_NOSYS"}
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allBjets_Sys_pt_NOSYS",
      "static_cast<float>(allLep_allBjets_Sys_TLV_NOSYS.Pt())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allBjets_Sys_eta_NOSYS",
      "static_cast<float>(allLep_allBjets_Sys_TLV_NOSYS.Eta())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allBjets_Sys_phi_NOSYS",
      "static_cast<float>(allLep_allBjets_Sys_TLV_NOSYS.Phi())"
  );
  mainNode = MainFrame::systematicStringDefine(
      mainNode,
      "allLep_allBjets_Sys_mass_NOSYS",
      "static_cast<float>(allLep_allBjets_Sys_TLV_NOSYS.M())"
  );

  /*
    ===============================================================
      DeltaR between loBjet and loZCandidate
    ===============================================================
  */
  // Function
  auto deltaR_leading = [](
      const std::vector<TLV>& tlv1,
      const std::vector<TLV>& tlv2
  ) {
      float delEta = 0.;
      float delPhi = 0.;
      float delRSqrd = 0.;

      float result = 0.;
      delEta = tlv1.at(0).Eta() - tlv2.at(0).Eta();
      delPhi = tlv1.at(0).Phi() - tlv2.at(0).Phi();
      delRSqrd = pow(delEta,2) + pow(delPhi,2);
      result = sqrt(delRSqrd);

      return result;
  };
  // Variable - DeltaR
  mainNode = MainFrame::systematicDefine(
      mainNode,
      "deltaR_Bjet0_ZCandidate0_NOSYS",
      deltaR_leading,
      {"bjet_TLV_NOSYS","sorted_ZCandidate_TLV_NOSYS"}
  );
  
  /*
    ===============================================================
      DeltaPhi between Bjet0 and ZCandidate0
    ===============================================================
  */
  // Variable - DeltaPhi
  mainNode = MainFrame::systematicStringDefine(
    mainNode,
    "deltaPhi_Bjet0_ZCandidate0_NOSYS",
    "static_cast<float>(std::abs(bjet_TLV_NOSYS.at(0).Phi() - sorted_ZCandidate_TLV_NOSYS.at(0).Phi()))"
  );

  /* 
  ================================
    End of defineVariables
  ================================
  */

  return mainNode;
}

ROOT::RDF::RNode tWZClass::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                      const std::shared_ptr<Sample>& sample,
                                                      const UniqueSampleID& id) {

  mainNode = defineVariables(mainNode,sample,id);

  return mainNode;
}

ROOT::RDF::RNode tWZClass::defineVariablesTruth(ROOT::RDF::RNode node,
                                                     const std::string& /*sample*/,
                                                     const std::shared_ptr<Sample>& /*sample*/,
                                                     const UniqueSampleID& /*sampleID*/) {
  return node;
}

ROOT::RDF::RNode tWZClass::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}

