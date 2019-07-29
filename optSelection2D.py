#!/usr/bin/env python

import os
import argparse
import ROOT
import math
import numpy

def createDir( odir ):
  if not os.path.exists( odir ):
    os.mkdir( odir )
    print "   Creating the output directory %s " % odir
  else:
    print "   The output directory already exists"
# createDir()

def getTree( fNames, tName ):

  t = ROOT.TChain( tName )
  for fName in fNames:
    t.AddFile( fName )

  t.SetDirectory(0)
  return t
# getTree()

def selectEvents( tree, Vars, pl, ph, np, costhll, costhlh, ncosthl, costhhl, costhhh, ncosthh ):

  n = tree.GetEntries()
  nFiducial = 0
  nPassed = {}
  rPassed = {}
  largeP  = {}
  smallP  = {}
  
  pBins      = numpy.linspace( pl, ph, np )
  costhlCuts = numpy.linspace( costhll, costhlh, ncosthl )
  costhhCuts = numpy.linspace( costhhl, costhhh, ncosthh )
  
  for var in Vars:
    nPassed[var] = {}
    rPassed[var] = {}
    largeP[var]  = {}
    smallP[var]  = {}
    for pBin in pBins:
      largeP[var][pBin] = {}
      smallP[var][pBin] = {}
      for costhlCut in costhlCuts:
        for costhhCut in costhhCuts:
          if costhhCut <= costhlCut: continue
          cut = '%0.2f_%0.2f' % ( costhlCut, costhhCut )
          largeP[var][pBin][cut] = 0
          smallP[var][pBin][cut] = 0
  
  for i in tree:
    if i.isIn10kton == 0: continue
    nFiducial += 1
    
    for var in Vars:
      
      p = eval( 'i.%sP[0]' % var )
      leaf = eval( 'i.%sAngle[0]' % var )
      costheta = math.cos( leaf )
      if p == 0.: continue
      
      for costhlCut in costhlCuts:
        for costhhCut in costhhCuts:
          if costhhCut <= costhlCut: continue
          if costheta < costhlCut or costheta > costhhCut: continue
          for pBin in pBins:
            cut = '%0.2f_%0.2f' % ( costhlCut, costhhCut )
            # print 'var: %s, pBin: %f, cut: %s' % ( var, pBin, cut )
            if p < pBin:
              smallP[var][pBin][cut] += 1
            else:
              largeP[var][pBin][cut] += 1

  for var in Vars:
    for pBin in pBins:
      for smallKey in smallP[var][pBin].keys():
        for largeKey in largeP[var][pBin].keys():
          cut = '%f_%s_%s' % ( pBin, smallKey, largeKey )
          nPassed[var][cut] = smallP[var][pBin][smallKey] + largeP[var][pBin][largeKey]


  for var in Vars:
    for Cut in nPassed[var].keys():
      rPassed[var][Cut] = float(nPassed[var][Cut])/float(nFiducial)
      # print nPassed[var][costhetaCut], nFiducial

  return n, nFiducial, nPassed, rPassed
# def selectEvents()

def optimizeSelection( mass, gamma, Vars, nPassed, rPassed ):

  E = mass * gamma
  if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
  else: Eround = E
  sKey = 'e%s_m%s' %( Eround, mass )
  # Scale the background events to 40kton*10 year exposure
  backgroundScale = 40./28.705
  bestCut = {}
  bestEff = {}
  bestBkg = {}
  
  for var in Vars:
    SprimeMin = 10000.
    for Cut in nPassed[sKey][var].keys():
      sEff = rPassed[sKey][var][Cut]
      bEvents = float(nPassed['atmos'][var][Cut])*backgroundScale
      Sprime = 25./ (2.*sEff) + math.sqrt( 25.*bEvents/ (sEff*sEff) + 625./(4.*sEff*sEff) )
      if Sprime < SprimeMin:
        SprimeMin = Sprime
        bestCut[var] = Cut
        bestEff[var] = sEff
        bestBkg[var] = bEvents
        bestBkgErr[var] = math.sqrt( float(nPassed['atmos'][var][Cut]) )*backgroundScale
  
  return bestCut, bestEff, bestBkg, bestBkgErr
  
# def optimizeSelection()

if __name__ == "__main__":
 
  Vars = [ 'Visible', 'VisibleNoN', 'SmearedReconstructable', 'SmearedReconstructableNoN' ]

  parser = argparse.ArgumentParser( description = 'Optimize the selection.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  
  args = parser.parse_args()

  Masses     = [ 5, 10, 20, 40 ]
  Gammas     = [ 1.1, 1.25 ]
  hDict      = {}
  nTotal     = {}
  nFiducial  = {}
  nPassed    = {}
  passRate   = {}
  costhll    = 0.2
  costhlh    = 1.
  ncosthl    = 17
  costhhl    = 0.5
  costhhh    = 1.
  ncosthh    = 11
  pl         = 0.1
  ph         = 1.
  np         = 19
  bestCut    = {}
  bestEff    = {}
  bestBkg    = {}
  bestBkgErr = {}
  
  # Create the output directory
  createDir( args.oDir )
  
  # Create the output text file
  txtName = '%s/Efficiency.txt' % ( args.oDir )
  txtFile = open( txtName, 'w' )
    

  # Background
  bFile = [ '%s/prodgenie_atmnu_max_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir, '%s/prodgenie_atmnu_min_dune10kt_gen_g4_NCFilter_RecoSmear_ana.root' % args.bDir ]
  bTree = getTree( bFile, 'MCParticles' )
  nTotal['atmos'], nFiducial['atmos'], nPassed['atmos'], passRate['atmos'] = selectEvents( bTree, Vars, pl, ph, np, costhll, costhlh, ncosthl, costhhl, costhhh, ncosthh )

  # Signal
  for Mass in Masses:
    for Gamma in Gammas:
      
      if Mass in [ 5, 20, 40 ] and ( Gamma == 2 ):
        continue
      
      E = Mass * Gamma
      if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
      else: Eround = E
      sFile = [ '%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_RecoSmear_ana.root' %( args.sDir, str(Eround), str(Mass) ) ]
      sKey = 'e%s_m%s' %( Eround, Mass )
      sTree = getTree( sFile, 'MCParticles' )
      nTotal[sKey], nFiducial[sKey], nPassed[sKey], passRate[sKey] = selectEvents( sTree, Vars, pl, ph, np, costhll, costhlh, ncosthl, costhhl, costhhh, ncosthh )
      # optimize the selection
      bestCut[sKey], bestEff[sKey], bestBkg[sKey], bestBkgErr[sKey] = optimizeSelection( Mass, Gamma, Vars, nPassed, passRate )
      txtFile.write('Mass = %s GeV, E = %s GeV\n' % ( str(Mass), str(Eround) ) )
      
      for var in Vars:
        txtFile.write( '  %s: optimal cut at %s, signal efficieny = %f, background counts = %f, background uncertainty =  %f\n' %( var, bestCut[sKey][var], bestEff[sKey][var], bestBkg[sKey][var], bestBkgErr[sKey][var] ) )
