#!/usr/bin/env python

import argparse
import ROOT
import math
import os

def getTree( fNames, tName ):

  t = ROOT.TChain( tName )
  for fName in fNames:
    t.AddFile( fName )

  t.SetDirectory(0)
  return t
# getTree()


def bookAngularHistograms( Vars, Sample ):

  hList = {}
  nBins = 40
  xMin  = -1.
  xMax  = 1.
  xTitle = 'cos(#theta)'
  yTitle = ''

  for var in Vars:
    
    hName = '%s_%s' % ( Sample, var )

    hList[var] = ROOT.TH1D( hName, "%s; %s; %s" % ( hName, xTitle, yTitle ), nBins, xMin, xMax )

  return hList

# def bookAngularHistograms()

def fillHistograms( tree, hAngularList ):

  n = tree.GetEntries()
  nPassedEvents = {}
  rPassedEvents = {}
  
  for var in hAngularList.keys():
    nPassedEvents[var] = 0
    rPassedEvents[var] = 0
  
  for i in tree:
        
    for var in hAngularList.keys():
      leaf = eval('i.%s[0]' % var)
      costheta = math.cos( leaf )
      if ( i.isIn10kton == 1 ):
        hAngularList[var].Fill( costheta )
        nPassedEvents[var] += 1

  for var in hAngularList.keys():
    rPassedEvents[var] = float(nPassedEvents[var])/float(n)
    print '%s: %d/%d' %( var, nPassedEvents[var], n )

  return hAngularList
# def fillHistograms()

def makeAngularPlot( hADict, var, oDir ):

  ROOT.gStyle.SetOptStat(0)
  cName = var
  c = ROOT.TCanvas( cName , cName, 1600, 1200 )
  c.SetBottomMargin(0.15)
  
  colors = { 'nominal': 906, 'noosc_E001G_010G': 593, 'noosc_E010G_100G': 870, 'maxnutau_E001G_010G': 417, 'maxnutau_E010G_100G': 797 }
  isFirst = True

  l = ROOT.TLegend( 0.6, 0.65, 0.9, 0.85 )
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  
  # Find the Ymax
  # h = hADict['noosc_E010G_100G'][var]
  # h.Scale( 1./h.Integral() )
  # Ymax = h.GetBinContent( h.GetMaximumBin() )
  
  for sample in hADict.keys():
    h = hADict[sample][var]
    h.Scale( 1./h.Integral() )
    h.SetLineColor( colors[sample] )
    h.SetLineWidth(3)
    if isFirst:
      h.SetTitle("")
      h.GetXaxis().SetTitle("cos(#theta)")
      h.GetXaxis().SetTitleSize(0.07)
      h.GetXaxis().SetLabelSize(0.06)
      h.GetYaxis().SetLabelSize(0.06)
      h.GetYaxis().SetRangeUser(0.01, 0.05)
      h.Draw("HIST")
      isFirst = False
    else:
      h.Draw("same HIST")

    leg = ''
    if sample.find( 'noosc' ) != -1:
      leg = 'Honda'
    elif sample.find( 'maxnutau' ) != -1:
      leg = 'Honda #nu_{#tau}'
    if sample.find( 'E001G_010G' ) != -1:
      leg += ' medium E'
    elif sample.find( 'E010G_100G' ) != -1:
      leg += ' high E'
    if sample == 'nominal':
      leg = 'Nominal'
      
    l.AddEntry( h, leg  )
  
  
  l.Draw()
  c.Draw()
  plotName = '%s/%s.png' %( oDir, var )
  c.SaveAs( plotName )
  
  
# def makeAngularPlot()



if __name__ == "__main__":
  
  Flavors = [ 'nominal', 'noosc', 'maxnutau' ]
  ERanges = [ 'E001G_010G', 'E010G_100G' ]
  Fluxes  = [ 'max', 'min' ]
  
  Vars = [ 'VisibleAngle', 'VisibleNoNAngle', 'LeadingParticleAngle', 'LeadingParticleNoNAngle',
           'SmearedVisibleAngle', 'SmearedVisibleNonAngle', 'SmearedReconstructableAngle', 
           'SmearedReconstructableNoNAngle', 'LeadingSmearedAngle', 'LeadingSmearedNoNAngle',
           'LeadingSmearedReconstructableAngle', 'LeadingSmearedReconstructableNoNAngle' ]
  
  parser = argparse.ArgumentParser( description = 'Make angular distributions for the atmospheric neutrino samples.')
  parser.add_argument( '-i', dest = 'iDir', type = str, help = 'Specify the directory of the input ROOT files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'Specify the directory for the output plots.' )
  
  args = parser.parse_args()
  
  if not os.path.exists( args.oDir ):
    os.makedirs( args.oDir )  
  
  hADict = {}
  Files  = {}
  
  for flavor in Flavors:
    
    fil = 'NCFilter'
    if flavor == 'maxnutau':
      fil = 'CCTauHadronicFilter'

    for eRange in ERanges:
      
      if flavor == 'nominal' and eRange == 'E010G_100G': continue
      fKey = '%s_%s' % ( flavor, eRange )
      if flavor == 'nominal':
        fKey = flavor
      filelist = []

      for flux in Fluxes:
        
        if flavor == 'nominal':
          f = '%s/prodgenie_atmnu_%s_dune10kt_gen_g4_NCFilter_reco_ana.root' % ( args.iDir, flux )
        else:
          f = '%s/prodgenie_atmnu_%s_dune10kt_honda_%s_%s_gen_g4_%s_reco_ana.root' % ( args.iDir, flux, flavor, eRange, fil )
        
        filelist.append( f )
        
      Files[fKey] = filelist
      
  for sample in Files.keys():
    
    print 'Filling the histograms for the sample %s...' % sample
    hADict[sample] = bookAngularHistograms( Vars, sample )
    tree = getTree( Files[sample], 'MCParticles' )
    hADict[sample] = fillHistograms( tree, hADict[sample] )

  for var in Vars:
    makeAngularPlot( hADict, var, args.oDir )
