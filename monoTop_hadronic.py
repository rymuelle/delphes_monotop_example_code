#!/usr/bin/env python

#test sublime sftp
import sys

import ROOT
import math

if len(sys.argv) < 2:
  print " Usage: monoTop_hadronic.py input_file"
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

inputFile = sys.argv[1]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")
branchMET = treeReader.UseBranch("MissingET")

# Book histograms
canv = ROOT.TCanvas()
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0)
histMass = ROOT.TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0)
hist_JJB = ROOT.TH1F("massJJB", "M_{j,j,b}", 100, 0, 1500)
hist_MET_pT = ROOT.TH1F("hist_MET_pT", "MET p_{T} hadronic #chi = 1 TeV; p_{T}; counts", 100, 0, 1000)
hist_bottom_pT = ROOT.TH1F("hist_bottom_pT", "bottom p_{T} hadronic #chi = 1 TeV; p_{T}; counts", 100, 0, 1000)
hist_bottom_eta = ROOT.TH1F("hist_bottom_eta", "bottom #eta hadronic #chi = 1 TeV; #eta; counts", 100, -2.5,2.5)
hist_b_over_t = ROOT.TH1F("hist_b_over_t", "#frac{E(b)}{E(t)} hadronic #chi = 1 TeV; #frac{e(b)}{E(t)}; counts", 100, 0,1)

count = 0
count_jet = 0

total = 0.0
btagCut = 0.0
jetCut = 0.0
leptonCut = 0.0
missingEt = 0.0
mJJB = 0.0
# Loop over all events
for entry in range(0, numberOfEntries):

  
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  btagIndex=-1
  nonBtagIndex=-1
  nBtag = 0
  total = total + 1
 
  print entry 

  hist_MET_pT.Fill(branchMET.At(0).MET)

  #jet loop
  if branchJet.GetEntries() > 0:
    for i in range(branchJet.GetEntries()):
    	jet = branchJet.At(i)
	if jet.BTag ==1:
		nBtag = nBtag +1
	if jet.BTag == 1 and btagIndex < 0:
		btagIndex=i
		bJet = jet
	if jet.BTag == 0 and nonBtagIndex == 1:
		nonBtagIndex = 2
		jet2 = jet
	if jet.BTag == 0 and nonBtagIndex < 0:
		nonBtagIndex = 1
		jet1 = jet
	
  if btagIndex >-1 : hist_bottom_pT.Fill( branchJet.At(btagIndex).PT )
  if btagIndex >-1: hist_bottom_eta.Fill( branchJet.At(btagIndex).Eta )

  #btag cut
  if not (nBtag == 1 and btagIndex >-1 and branchJet.At(btagIndex).PT >70 and abs(branchJet.At(btagIndex).Eta) < 2.5): continue
 
  btagCut = btagCut + 1

# Loop over all events
  #other jet cut
  if nonBtagIndex == 2:
	  if not(jet1.PT > 30 and abs(jet1.Eta) < 2.5): continue
	  if not(jet2.PT > 30 and abs(jet2.Eta) < 2.5): continue
  
  jetCut = jetCut + 1

  nLepton_over_30 = 0
  electron_index = 0
  muon_index = 0

  #electron
#  if branchElectron.GetEntries() > 0:
#	  for i in range(branchElectron.GetEntries()):
#		  electron = branchElectron.At(i)
#		  if(electron.PT > 30 and abs(electron.Eta) < 2.1):
#			 nLepton_over_30 = nLepton_over_30 + 1 
#			 electron_index = i
#			 lepton = electron
#  #muon
#  if branchMuon.GetEntries() > 0:
#	  for i in range(branchMuon.GetEntries()):
#		  muon = branchMuon.At(i)
#		  if(muon.PT > 30 and abs(muon.Eta) < 2.1):
#			 nLepton_over_30 = nLepton_over_30 + 1 
#			 muon_index = i
#			 lepton = muon
#
#  #lepton cuts
#  if not nLepton_over_30==1: continue
#
#  # W pt cut
#  w_pt = 0
##  if btagIndex > -1 and  electron_index > -1:
##	w_pt = (lepton.P4()+bJet.P4()).Pt()
##  if btagIndex > -1 and  muon_index > -1:
#  w_pt = (lepton.P4()+bJet.P4()).Pt()
#  #w_pt = lepton.PT
#  if not w_pt > 50: continue
#
#  #phit cut:
#
#  if not abs(lepton.Phi - bJet.Phi) <1.7: continue
#
  leptonCut = leptonCut + 1
#  #Met cut
#  if not branchMET.At(0).MET > 350: continue
  missingEt = missingEt + 1 
  #mass j j b
  mass_j_j_b = (jet1.P4()+jet2.P4()+bJet.P4()).M()
  #print mass_j_j_b
  hist_JJB.Fill(mass_j_j_b)
  hist_b_over_t.Fill(bJet.P4().E()/(jet1.P4()+jet2.P4()+bJet.P4()).E())
  if mass_j_j_b > 450: continue
  
  mJJB = mJJB + 1

#
#  #mass T cut:
#  transverse_mass = 2*(lepton.P4().E())*(branchMET.At(0).P4().E())*(1-math.cos(lepton.Phi-branchMET.At(0).Phi))
#  transverse_mass= math.sqrt(transverse_mass)
#  #print transverse_mass, ((lepton.Phi-branchMET.At(0).Phi))
#
#  if not transverse_mass > 400: continue
#
#  #if not leadJetBTag: continue
#  #if not leadJetPT > 70: continue
 
  count = count + 1
print count, (count+0.0)/numberOfEntries, numberOfEntries

hist_JJB.Draw()
canv.SaveAs("hist_JJB.png")
hist_b_over_t.Draw()
canv.SaveAs("hist_b_over_t.png")
hist_bottom_eta.Draw()
canv.SaveAs("hist_bottom_eta.png")
hist_bottom_pT.Draw()
canv.SaveAs("hist_bottom_pT.png")
hist_MET_pT.Draw()
canv.SaveAs("hist_MET_pT.png")

print "total","\t", total,"\t", total/total,"\t", total/total
print "btagCut","\t", btagCut,"\t", btagCut/total,"\t", btagCut/total
print "jetCut","\t", jetCut,"\t", jetCut/total,"\t", jetCut/btagCut
print "leptonCut","\t", leptonCut,"\t", leptonCut/total,"\t", leptonCut/btagCut
print "missingEt","\t", missingEt,"\t", missingEt/total,"\t", missingEt/leptonCut
print "mJJB","\t", mJJB,"\t", mJJB/total,"\t", mJJB/missingEt
 


#    # Plot jet transverse momentum
#    histJetPT.Fill(jet.PT)
#
#    # Print jet transverse momentum
#   # print jet.PT
#
#  # If event contains at least 2 electrons
#  if branchElectron.GetEntries() > 1:
#    # Take first two electrons
#    elec1 = branchElectron.At(0)
#    elec2 = branchElectron.At(1)
#
#    # Plot their invariant mass
#    histMass.Fill(((elec1.P4()) + (elec2.P4())).M())
#
## Show resulting histograms
#histJetPT.Draw()
##histMass.Draw()
#

raw_input("Press Enter to continue...")
