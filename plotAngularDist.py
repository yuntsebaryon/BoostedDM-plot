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


def bookAngularHistograms( Vars, Mass, Gamma ):

  hList = {}
  nBins = 100
  xMin  = -1.
  xMax  = 1.
  xTitle = 'cos(#theta)'
  yTitle = ''

  for var in Vars:
    
    hName = "e%s_m%s_%s" % ( str(Gamma*Mass), str(Mass), var )
        
    if Mass == 0.:
      hName = 'atmosNu_%s' % var
      
    hList[var] = ROOT.TH1D( hName, "%s; %s; %s" % ( hName, xTitle, yTitle ), nBins, xMin, xMax )

  return hList

# def bookAngularHistograms()

def bookMomentumHistograms( Vars, Mass, Gamma ):

  hList = {}
  nBins = 100
  xMin  = 0.
  xMax  = 20.
  xTitle = 'Momentum [GeV]'
  yTitle = ''

  for var in Vars:
    
    hName = "e%s_m%s_%s" % ( str(Gamma*Mass), str(Mass), var )
    
    if Mass == 0.:
      hName = 'atmosNu_%s' % var
      
    hList[var] = ROOT.TH1D( hName, "%s; %s; %s" % ( hName, xTitle, yTitle ), nBins, xMin, xMax )

  return hList

# def bookMomentumHistograms()

def selectEvents( tree, hMomentumList, hAngularList, costhetaCut ):

  n = tree.GetEntries()
  nPassedEvents = {}
  rPassedEvents = {}
  
  for var in hAngularList.keys():
    if var in [ 'InParticleAngle', 'OutParticleAngle' ]: continue
    nPassedEvents[var] = 0
    rPassedEvents[var] = 0
  
  for i in tree:
    
    for var in [ 'InParticleP', 'OutParticleP' ]:
      leaf = eval( 'i.%s[0]' % var )
      hMomentumList[var].Fill( leaf )
    
    for var in hAngularList.keys():
      leaf = eval('i.%s[0]' % var)
      costheta = math.cos( leaf )
      hAngularList[var].Fill( costheta )
      if var in [ 'InParticleAngle', 'OutParticleAngle' ]: continue
      if ( i.isIn10kton == 1 ):
        nPassedEvents[var] += 1

  for var in hAngularList.keys():
    if var in [ 'InParticleAngle', 'OutParticleAngle' ]: continue
    rPassedEvents[var] = float(nPassedEvents[var])/float(n)
    print '%s: %d/%d' %( var, nPassedEvents[var], n )

  return hMomentumList, hAngularList, n, nPassedEvents, rPassedEvents
# def selectEvents()

def makeAngularPlot( hADict, Mass, Gammas, var, oDir ):

  ROOT.gStyle.SetOptStat(0)
  cName = 'm%s_%s' %( Mass, var )
  c = ROOT.TCanvas( cName , cName, 1600, 1200 )
  c.SetBottomMargin(0.15)
  
  colors = { '1.1': 906, '1.25': 593, '2': 870, '10': 417, 'atmos': 797 }
  isFirst = True

  l = ROOT.TLegend( 0.15, 0.6, 0.55, 0.8 )
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  
  # Find the Ymax
  gamma = Gammas[-1]
  e = gamma*Mass
  sKey = 'e%s_m%s' %( e, Mass )
  h = hADict[sKey][var]
  Ymax = h.GetBinContent( h.GetMaximumBin() )
  
  for gamma in Gammas:
    if ( gamma == 2 ) and Mass in [ 5, 20, 40 ]:
      continue
    e = gamma*Mass
    if e in [ 11., 22., 25., 44., 50. ]: eround = int(e)
    else: eround = e
    sKey = 'e%s_m%s' %( eround, Mass )
    h = hADict[sKey][var]
    if gamma in [ 1.1, 1.25 ]:
      h.Scale( 5. )
    h.SetLineColor( colors[str(gamma)] )
    h.SetLineWidth(3)
    if isFirst:
      h.SetTitle("")
      h.GetXaxis().SetTitle("cos(#theta)")
      h.GetXaxis().SetTitleSize(0.07)
      h.GetXaxis().SetLabelSize(0.06)
      h.GetYaxis().SetLabelSize(0.06)
      h.GetYaxis().SetRangeUser(0, 1.1*Ymax)
      h.Draw("HIST")
      isFirst = False
    else:
      h.Draw("same HIST")
    if gamma in [ 1.1, 1.25 ]:
      l.AddEntry( h, 'DM, E = %s GeV #times 5' % str(eround) )
    else:
      l.AddEntry( h, 'DM, E = %s GeV' % str(eround) )
    # print 'M = %d, E = %f, total events = %f' %( Mass, e, h.Integral() )
  
  h = hADict['atmos'][var]
  h.SetLineColor( colors['atmos'] )
  h.SetLineWidth(3)
  print 'Atmospheric neutrino, total events = %f' % h.Integral()
  h.Scale( 100./28.705 )
  l.AddEntry( h, 'Atmospheric Nus 1 Mton-year' )
  h.Draw("same HIST")
  
  l.Draw()
  c.Draw()
  plotName = '%s/m%s_%s.png' %( oDir, str(Mass), var )
  c.SaveAs( plotName )
  h.Scale( 28.705/100. )
  
  
# def makeAngularPlot()

if __name__ == "__main__":

  MomentumVars = [ 'InParticleP', 'OutParticleP' ]
 
  AngularVars = [ 'InParticleAngle', 'OutParticleAngle',
                  'VisibleAngle', 'VisibleNoNAngle', 'LeadingParticleAngle', 'LeadingParticleNoNAngle',
                  'SmearedVisibleAngle', 'SmearedVisibleNonAngle', 'SmearedReconstructableAngle', 'SmearedReconstructableNoNAngle',
                  'LeadingSmearedAngle', 'LeadingSmearedNoNAngle', 'LeadingSmearedReconstructableAngle', 'LeadingSmearedReconstructableNoNAngle' ]

  parser = argparse.ArgumentParser( description = 'Make angular distributions.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  
  args = parser.parse_args()


  Masses   = [ 5, 10, 20, 40 ]
  Gammas   = [ 1.1, 1.25, 2, 10 ]
  hPDict   = {}
  hADict   = {}
  nTotal   = {}
  nPassed  = {}
  passRate = {}


  if not os.path.exists( args.oDir ):
    os.makedirs( args.oDir )

  # Background
  bFile = [ '%s/prodgenie_atmnu_max_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir, '%s/prodgenie_atmnu_min_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir ]
  hPDict['atmos'] = bookMomentumHistograms( MomentumVars, 0., 0. )
  hADict['atmos'] = bookAngularHistograms( AngularVars, 0., 0. )
  bTree = getTree( bFile, 'MCParticles' )
  print 'Atmospheric neutrino...'
  hPDict['atmos'], hADict['atmos'], nTotal['atmos'], nPassed['atmos'], passRate['atmos'] = selectEvents( bTree, hPDict['atmos'], hADict['atmos'], 0.6 )

  for Mass in Masses:
    for Gamma in Gammas:
      
      if Mass in [ 5, 20, 40 ] and ( Gamma == 2 ):
        continue

      E = Mass * Gamma
      if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
      else: Eround = E
      if E in [ 5.5, 6.25, 12.5 ]:
        Estr = str(E)
        Estr.replace( '.', 'p' )
      else: Estr = str(Eround)
      sFile = ['%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_RecoSmear_ana.root' %( args.sDir, Estr, str(Mass) ) ]
      sKey = 'e%s_m%s' %( Eround, Mass )
      hPDict[sKey] = bookMomentumHistograms( MomentumVars, Mass, Gamma )
      hADict[sKey] = bookAngularHistograms( AngularVars, Mass, Gamma )
      sTree = getTree( sFile, 'MCParticles' )
      print 'M = %d, E = %f' %( Mass, E )
      hPDict[sKey], hADict[sKey], nTotal[sKey], nPassed[sKey], passRate[sKey] = selectEvents( sTree, hPDict[sKey], hADict[sKey], 0.6 )

  for Mass in Masses:
    for var in AngularVars:
      makeAngularPlot( hADict, Mass, Gammas, var, args.oDir )


  # Output the selection efficiency
  # txtName = '%s/Efficiency.txt' % args.oDir
  # txtFile = open( txtName, 'w' )
  # for sample in nTotal.keys():
  #   n = nTotal[sample]
  #   if sample in [ 'atmos' ]:
  #     n = float( nTotal[sample] )/ 14.239
  #   txtFile.write( 'Sample: %s, total events: %f\n' % ( sample, n ) )
    
  #   for var in AngularVars:
  #     if var in ['InParticleAngle', 'OutParticleAngle']: continue
  #     p = nPassed[sample][var]
  #     if sample in [ 'atmos' ]:
  #       p = float( nPassed[sample][var] )/ 14.239
  #     txtFile.write( '   %s: %f (%f)\n' %( var, p, passRate[sample][var] ) )
