#!/usr/bin/env python

import sys

import ROOT
import math

if len(sys.argv) < 2:
  print " Usage: Example1.py input_file"
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
hist_JJB = ROOT.TH1F("massJJB", "M_{j,j,b}", 100, 0, 1400)
hist_MET_pT = ROOT.TH1F("hist_MET_pT", "MET p_{T} leptonic #chi = 1 TeV; p_{T}; counts", 100, 0, 1000)
hist_bottom_pT = ROOT.TH1F("hist_bottom_pT", "bottom p_{T} leptonic #chi = 1 TeV; p_{T}; counts", 100, 0, 1000)
hist_bottom_eta = ROOT.TH1F("hist_bottom_eta", "bottom #eta leptonic #chi = 1 TeV; #eta; counts", 100, -2.5,2.5)
hist_lep_pT = ROOT.TH1F("hist_lep_pT", "leading lepton p_{T} leptonic #chi = 1 TeV; p_{T}; counts", 100, 0, 1000)
hist_lep_eta = ROOT.TH1F("hist_lep_eta", "leading lepton #eta leptonic #chi = 1 TeV; #eta; counts", 100, -2.5,2.5)
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0)
histMass = ROOT.TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0)
b_ratio = ROOT.TH1F("b_ratio", "R E(b)/E(t)", 0, 1, 100)
MTW = ROOT.TH1F("MTW", "MTW", 100, 0, 1000)

count = 0
count_jet = 0

total = 0.0
bJetCut = 0.0
nonBJet = 0.0
leptonCut = 0.0
Tau = 0.0
pTW = 0.0
deltaPhi = 0.0
deltaR = 0.0
MET = 0.0
mT = 0.0

# Loop over all events
for entry in range(0, numberOfEntries):

  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  btagIndex=-1
  nonBtagIndex=-1
  nNonBJet = 0
  nBtag = 0
  print entry
  total = total + 1
  #jet loop
  if branchJet.GetEntries() > 0:
    for i in range(branchJet.GetEntries()):
    	jet = branchJet.At(i)
        if jet.BTag ==1 and jet.PT > 70 and abs(jet.Eta) < 2.5:
		nBtag = nBtag + 1
	if jet.BTag == 1 and btagIndex < 0:
		btagIndex=i
		bJet = jet
	if jet.BTag == 0 and nonBtagIndex < 0:
		nonBtagIndex=i
	if jet.BTag == 0 and jet.PT > 30 and (jet.Eta) < 2.1:
		nNonBJet = nNonBJet + 1

  nLepton_over_30 = 0
  electron_index = 0
  muon_index = 0
  if branchElectron.GetEntries() > 0:
    for i in range(branchElectron.GetEntries()):
      electron = branchElectron.At(i)
      if(electron.PT > 30 and abs(electron.Eta) < 2.1):
       nLepton_over_30 = nLepton_over_30 + 1 
       electron_index = i
       lepton = electron
  #muon
  if branchMuon.GetEntries() > 0:
    for i in range(branchMuon.GetEntries()):
      muon = branchMuon.At(i)
      if(muon.PT > 30 and abs(muon.Eta) < 2.1):
       nLepton_over_30 = nLepton_over_30 + 1 
       muon_index = i
       lepton = muon
  #lepton cuts

  hist_MET_pT.Fill(branchMET.At(0).MET)
  hist_lep_pT.Fill(lepton.PT)
  hist_lep_eta.Fill(lepton.Eta)

  #btag cut
  #if not (nBtag == 1  and branchJet.At(btagIndex).PT >70 and abs(branchJet.At(btagIndex).Eta) < 2.5): continue
  if not (nBtag == 1  and bJet.PT >  70 and abs(bJet.Eta) < 2.5): continue
  bJetCut = bJetCut + 1

  #other jet cut
#  if  (nonBtagIndex >-1):
#    if branchJet.At(nonBtagIndex).PT > 30 or abs(branchJet.At(nonBtagIndex).Eta) > 2.5: continue
  
  if nNonBJet > 1: continue
  nonBJet = nonBJet + 1
 

  #electron
#  if branchElectron.GetEntries() > 0:
#	  for i in range(branchElectron.GetEntries()):
#		  electron = branchElectron.At(i)
#		  if(electron.PT > 30 and abs(electron.Eta) < 2.1):
#			 nLepton_over_30 = nLepton_over_30 + 1 
#			 electron_index = i
#			 lepton = electron
 # #muon
  #if branchMuon.GetEntries() > 0:
#	  for i in range(branchMuon.GetEntries()):
#		  muon = branchMuon.At(i)
#		  if(muon.PT > 30 and abs(muon.Eta) < 2.1):
#			 nLepton_over_30 = nLepton_over_30 + 1 
#			 muon_index = i
#			 lepton = muon
  #lepton cuts
  if not nLepton_over_30==1: continue
  leptonCut = leptonCut + 1

  Tau = Tau +1
  # W pt cut
  w_pt = 0
#  if btagIndex > -1 and  electron_index > -1:
#	w_pt = (lepton.P4()+bJet.P4()).Pt()
#  if btagIndex > -1 and  muon_index > -1:
  #w_pt = (lepton.P4()+bJet.P4()).Pt()
  w_pt = (lepton.P4()+bJet.P4()).Pt()

  #w_pt = lepton.PT
  if not w_pt > 50: continue
  pTW = pTW + 1

  #phit cut:

  if not abs(lepton.Phi - bJet.Phi) <1.7: continue
  deltaPhi = deltaPhi + 1

  deltaR = deltaR + 1
  #Met cut
  if not branchMET.At(0).MET > 100: continue
  MET = MET + 1

  #mass T cut:
  transverse_mass = 2*(lepton.P4().E())*(branchMET.At(0).P4().E())*(1-math.cos(lepton.Phi-branchMET.At(0).Phi))
  transverse_mass= math.sqrt(transverse_mass)
  #print transverse_mass, ((lepton.Phi-branchMET.At(0).Phi))


  MTW.Fill(transverse_mass)
  if not transverse_mass > 400: continue
  mT = mT + 1
  #if not leadJetBTag: continue
  #if not leadJetPT > 70: continue
 
  count = count + 1

MTW.Draw()
print count, (count+0.0)/numberOfEntries, numberOfEntries

print "total", "\t", total, "\t", total/total, "\t", total/total
print "bJetCut", "\t", bJetCut, "\t", bJetCut/total, "\t", bJetCut/total
print "nonBJet", "\t", nonBJet, "\t", nonBJet/total, "\t", nonBJet/bJetCut
print "leptonCut", "\t", leptonCut, "\t", leptonCut/total, "\t", leptonCut/nonBJet
print "Tau", "\t", Tau, "\t", Tau/total, "\t", Tau/leptonCut
print "pTW", "\t", pTW, "\t", pTW/total, "\t", pTW/Tau
print "deltaPhi", "\t", deltaPhi, "\t", deltaPhi/total, "\t", deltaPhi/pTW
print "deltaR", "\t", deltaR, "\t", deltaR/total, "\t", deltaR/deltaPhi
print "MET", "\t", MET, "\t", MET/total, "\t", MET/deltaR
print "mT", "\t", mT, "\t", mT/total, "\t", mT/MET


hist_MET_pT.Draw()
canv.SaveAs("hist_MET_pT.png")
hist_lep_pT.Draw()
canv.SaveAs("hist_lep_pT.png")
hist_lep_eta.Draw()
canv.SaveAs("hist_lep_eta.png")





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
