import os
import ROOT

base_dir = "/eos/user/k/kebarend/tWZ/FastFrames/input/"

for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file.endswith(".root"):
            file_path = os.path.join(root, file)
            f = ROOT.TFile.Open(file_path, "READ")
            if not f or f.IsZombie():
                print(f"Could not open {file_path}, skipping.")
                continue
            has_reco = f.GetListOfKeys().Contains("reco")
            f.Close()
            if not has_reco:
                print(f"Deleting {file_path} (no 'reco' TTree found)")
                os.remove(file_path)
            else:
                # Add histogram if not already present
                f = ROOT.TFile.Open(file_path, "UPDATE")
                if not f.GetListOfKeys().Contains("listOfSystematics"):
                    hist = ROOT.TH1F("listOfSystematics", "listOfSystematics", 1, 0, 1)
                    hist.GetXaxis().SetBinLabel(1, "NOSYS")
                    hist.Write()
                    print(f"Added 'listOfSystematics' histogram to {file_path}")
                f.Close()