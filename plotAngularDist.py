#!/usr/bin/env python

import argparse
import ROOT
import math

def getTree( fName, tName ):

  t = ROOT.TChain( tName )
  t.AddFile( fName )

  t.SetDirectory(0)
  return t
# getTree()


def bookHistograms( Vars, Mass, Gamma ):

  hList = {}
  nBins = 50
  xMin  = -1.
  xMax  = 1.
  xTitle = 'cos(#theta)'
  yTitle = ''

  for var in Vars:
    
    hName = "%s_e%s_m%s" % ( var, str(Gamma*Mass), str(Mass) )
    
    nBins = 50
    xMin  = -1.
    xMax  = 1.
    xTitle = 'cos(#theta)'
    yTitle = ''
    
    if Mass == 0.:
      if var == 'InDMP':
        hName = 'InNuP_atmosNu'
        xTitle = 'Momentum [GeV]'
        nBins = 50
        xMin = 0.
        xMax = 20.
      elif var == 'OutDMP':
        hName = 'OutNuP_atmosNu'
        xTitle = 'Momentum [GeV]'
        nBins = 50
        xMin = 0.
        xMax = 20.

      else:
        hName = '%s_atmosNu' % var

    if var in [ 'InDMP', 'OutDMP' ]:
      xTitle = 'Momentum [GeV]'
      nBins = 50
      xMin = 0.
      xMax = 20.

      
    hList[var] = ROOT.TH1D( hName, "%s; %s; %s" % ( hName, xTitle, yTitle ), nBins, xMin, xMax )

  return hList

# def bookHistograms()

def selectEvents( tree, hList, AngularVars, isDM ):

  n = tree.GetEntries()
  nPassedEvents = {}
  rPassedEvents = {}
  
  for var in AngularVars:
    nPassedEvents[var] = 0
    rPassedEvents[var] = 0
  
  for i in tree:
    
    if isDM:
      for var in [ 'InDMP', 'OutDMP' ]:
        leaf = eval( 'i.%s[0]' % var )
        hList[var].Fill( leaf )

      for var in [ 'InDMAngle', 'OutDMAngle' ]:
        leaf = eval( 'i.%s[0]' % var )
        costheta = math.cos( leaf )
        hList[var].Fill( costheta )
    else:
      for var in [ 'InNuP', 'OutNuP' ]:
        leaf = eval( 'i.%s[0]' % var )
        hList[var.replace('Nu', 'DM')].Fill( leaf )

      for var in [ 'InNuAngle', 'OutNuAngle' ]:
        leaf = eval( 'i.%s[0]' % var )
        costheta = math.cos( leaf )
        hList[var.replace('Nu', 'DM')].Fill( costheta )

    
    for var in AngularVars:
      leaf = eval('i.%s[0]' % var)
      costheta = math.cos( leaf )
      hList[var].Fill( costheta )
      if costheta > 0.6:
        if isDM:
          nPassedEvents[var] += 1
        elif i.isIn10kton == 1:
          nPassedEvents[var] += 1


  for var in AngularVars:
    rPassedEvents[var] = float(nPassedEvents[var])/float(n)

  return hList, n, nPassedEvents, rPassedEvents
# def selectEvents()

def makePlot( hDict, Mass, Gammas, var, oDir ):

  ROOT.gStyle.SetOptStat(0)
  cName = 'm%s_%s' %( Mass, var )
  c = ROOT.TCanvas( cName , cName, 800, 600 )
  c.SetBottomMargin(0.15)
  
  colors = { '1.25': 593, '2': 906, '10': 417, 'atmos': 797 }
  isFirst = True

  l = ROOT.TLegend( 0.15, 0.6, 0.55, 0.8 )
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  
  # Find the Ymax
  gamma = Gammas[-1]
  e = gamma*Mass
  sKey = 'e%s_m%s' %( e, Mass )
  h = hDict[sKey][var]
  Ymax = h.GetBinContent( h.GetMaximumBin() )
  
  for gamma in Gammas:
    e = gamma*Mass
    if e == 25. or e == 50.: eround = int(e)
    else: eround = e
    sKey = 'e%s_m%s' %( eround, Mass )
    h = hDict[sKey][var]
    if gamma == 1.25:
      h.Scale( 5. )
    h.SetLineColor( colors[str(gamma)] )
    h.SetLineWidth(3)
    if isFirst:
      h.SetTitle("")
      # h.GetXaxis().SetTitle("cos(#theta)")
      h.GetXaxis().SetTitleSize(0.07)
      h.GetXaxis().SetLabelSize(0.06)
      h.GetYaxis().SetLabelSize(0.06)
      h.GetYaxis().SetRangeUser(0, 1.1*Ymax)
      h.Draw("HIST")
      isFirst = False
    else:
      h.Draw("same HIST")
    if gamma == 1.25:
      l.AddEntry( h, 'DM, E = %s GeV #times 5' % str(e) )
    else:
      l.AddEntry( h, 'DM, E = %s GeV' % str(e) )
  
  h = hDict['atmos'][var]
  if Mass == 5:
    h.Scale( 100./14.239 )
  h.SetLineColor( colors['atmos'] )
  h.SetLineWidth(3)
  h.Draw("same HIST")
  l.AddEntry( h, 'Atmospheric Nus #times 100' )
  
  l.Draw()
  c.Draw()
  # h.Scale( 0.014239 )
  plotName = '%s/m%s_%s.png' %( oDir, str(Mass), var )
  c.SaveAs( plotName )
  
  
# def makePlot

if __name__ == "__main__":

  Vars = [ 'InDMP', 'InDMAngle', 'OutDMP', 'OutDMAngle',
           'VisibleAngle', 'VisibleNoNAngle', 'LeadingParticleAngle', 'LeadingParticleNoNAngle',
           'SmearedVisibleAngle', 'SmearedVisibleNonAngle', 'SmearedReconstructableAngle', 'SmearedReconstructableNoNAngle',
           'LeadingSmearedAngle', 'LeadingSmearedNoNAngle', 'LeadingSmearedReconstructableAngle', 'LeadingSmearedReconstructableNoNAngle' ]
 
  AngularVars = [ 'VisibleAngle', 'VisibleNoNAngle', 'LeadingParticleAngle', 'LeadingParticleNoNAngle',
                  'SmearedVisibleAngle', 'SmearedVisibleNonAngle', 'SmearedReconstructableAngle', 'SmearedReconstructableNoNAngle',
                  'LeadingSmearedAngle', 'LeadingSmearedNoNAngle', 'LeadingSmearedReconstructableAngle', 'LeadingSmearedReconstructableNoNAngle' ]

  parser = argparse.ArgumentParser( description = 'Make angular distributions.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  
  args = parser.parse_args()



  Masses   = [ 5, 10, 20, 40 ]
  Gammas   = [ 1.25, 2, 10 ]
  hDict    = {}
  nTotal   = {}
  nPassed  = {}
  passRate = {}

  # Background
  bFile = '%s/BkgdSmearedMCParticles.root' % args.bDir
  hDict['atmos'] = bookHistograms( Vars, 0., 0. )
  bTree = getTree( bFile, 'MCParticles' )
  hDict['atmos'], nTotal['atmos'], nPassed['atmos'], passRate['atmos'] = selectEvents( bTree, hDict['atmos'], AngularVars, False )

  for Mass in Masses:
    for Gamma in Gammas:
      E = Mass * Gamma
      if E == 25. or E == 50.: Eround = int(E)
      else: Eround = E
      sFile = '%s/dune_fermion_e%s_m%s_g1_z1_fixedflux.0_Gen_g4_RecoSmear_ana.root' %( args.sDir, str(Eround), str(Mass) )
      sKey = 'e%s_m%s' %( Eround, Mass )
      hDict[sKey] = bookHistograms( Vars, Mass, Gamma )
      sTree = getTree( sFile, 'MCParticles' )
      hDict[sKey], nTotal[sKey], nPassed[sKey], passRate[sKey] = selectEvents( sTree, hDict[sKey], AngularVars, True )

  for Mass in Masses:
    for var in AngularVars:
      makePlot( hDict, Mass, Gammas, var, args.oDir )


  # Output the selection efficiency
  txtName = '%s/Efficiency.txt' % args.oDir
  txtFile = open( txtName, 'w' )
  for sample in nTotal.keys():
    n = nTotal[sample]
    if sample in [ 'atmos' ]:
      n = float( nTotal[sample] )/ 14.239
    txtFile.write( 'Sample: %s, total events: %f\n' % ( sample, n ) )
    
    for var in AngularVars:
      p = nPassed[sample][var]
      if sample in [ 'atmos' ]:
        p = float( nPassed[sample][var] )/ 14.239
      txtFile.write( '   %s: %f (%f)\n' %( var, p, passRate[sample][var] ) )
