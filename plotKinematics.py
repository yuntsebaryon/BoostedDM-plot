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


def bookHistograms( Vars, Mass, Gamma ):

  hList  = {}
  nXBins = 50
  xMin   = 0.
  xMax   = 1.
  xTitle = 'cos(#theta)'
  nYBins = 20
  yMin   = 0.
  
  if Gamma == 1.25:
    yMax = 2.
  elif Gamma in [ 0, 2 ]:
    yMax = 4.
  elif Gamma == 10:
    yMax = 80.
  else:
    yMax   = 1.
    
  yTitle = '3-Momentum [GeV]'

  for var in Vars:
    
    if Mass == 0.:
      hName = 'atmosNu_%s' % var
    else:
      hName = "e%s_m%s_%s" % ( str(Gamma*Mass), str(Mass), var )
      
    hList[var] = ROOT.TH2D( hName, "%s; %s; %s" % ( hName, xTitle, yTitle ), nXBins, xMin, xMax, nYBins, yMin, yMax )

  return hList

# def bookHistograms()

def fillHistograms( t, h, Vars ):
  
  for i in t:
    for var in Vars:
      if var == 'TrueKin':
        visibleCosth = math.cos( i.VisibleAngle[0] )
        h[var].Fill( visibleCosth, i.VisibleP[0] )
      elif var == 'RecoKin':
        recoCosth = math.cos( i.SmearedReconstructableAngle[0] )
        h[var].Fill( recoCosth, i.SmearedReconstructableP[0] )
      elif var == 'TrueNoNKin':
        visibleNoNCosth = math.cos( i.VisibleNoNAngle[0] )
        h[var].Fill ( visibleNoNCosth, i.VisibleNoNP[0] )
      elif var == 'RecoNoNKin':
        recoNoNCosth = math.cos( i.SmearedReconstructableNoNAngle[0] )
        h[var].Fill( recoNoNCosth, i.SmearedReconstructableNoNP[0] )
        
  return h

# def fillHistograms()

def makePlots( sample, hDict, Vars, oDir ):
  
  if not os.path.exists( outdir ):
    os.makedirs( outdir )
  
  for var in Vars:
    
    ROOT.gStyle.SetOptStat(0)
    cName = '%s%s' % ( sample, var )
    cTitle = '%s %s Kinematics' % ( sample, var )
    c = ROOT.TCanvas( cName, cTitle, 800, 600 )
    c.SetLeftMargin( 0.15 )
    c.SetBottomMargin( 0.15 )

    hDict[sample][var].GetXaxis().SetLabelSize(0.06)
    hDict[sample][var].GetYaxis().SetLabelSize(0.06)
    hDict[sample][var].GetXaxis().SetTitleSize(0.07)
    hDict[sample][var].GetYaxis().SetTitleSize(0.07)
    hDict[sample][var].SetTitle( '')
    hDict[sample][var].Draw('COLZ')
    c.Draw()
    plotName = '%s/%s_%s.png' % ( outdir, sample, var )
    c.SaveAs( plotName )

# def makePlots()


if __name__ == "__main__":
  
  Vars = [ 'TrueKin', 'RecoKin', 'TrueNoNKin', 'RecoNoNKin' ]
  
  parser = argparse.ArgumentParser( description = 'Make kinematic correlation plots.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  
  args = parser.parse_args()


  Masses   = [ 5, 10, 20, 40 ]
  Gammas   = [ 1.1, 1.25, 2, 10 ]
  hDict    = {}
  

  # Background
  bFile = [ '%s/prodgenie_atmnu_max_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir, '%s/prodgenie_atmnu_min_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir ]
  hDict['atmos'] = bookHistograms( Vars, 0, 0 )
  bTree = getTree( bFile, 'MCParticles' )
  print 'Atmospheric neutrino sample: Fill histograms...'
  hDict['atmos'] = fillHistograms( bTree, hDict['atmos'], Vars )
  print 'Atmospheric neutrino sample: Make plots...'
  outdir = '%s/atmos_kinematics' % args.oDir
  makePlots( 'atmos', hDict, Vars, outdir )

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
      hDict[sKey] = bookHistograms( Vars, Mass, Gamma )
      sTree = getTree( sFile, 'MCParticles' )
      print 'BDM M = %d, E = %f sample: Fill histograms...' %( Mass, E )
      hDict[sKey] = fillHistograms( sTree, hDict[sKey], Vars )
      print 'BDM M = %d, E = %f sample: Make plots...' %( Mass, E )
      outdir = '%s/%s_kinematics' % ( args.oDir, sKey )
      makePlots( sKey, hDict, Vars, outdir )
