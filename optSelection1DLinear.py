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

def selectEvents( tree, Vars, costhl, costhh, ncosth, ph ):

  n = tree.GetEntries()
  nFiducial = 0
  nPassed = {}
  rPassed = {}
  
  costhetaCuts = numpy.linspace( costhl, costhh, ncosth )
  
  for var in Vars:
    nPassed[var] = {}
    rPassed[var] = {}
    for costhetaCut in costhetaCuts:
      nPassed[var][costhetaCut] = 0
      rPassed[var][costhetaCut] = 0
  
  for i in tree:
    if i.isIn10kton == 0: continue
    nFiducial += 1
    
    for var in Vars:
      
      p = eval( 'i.%sP[0]' % var )
      leaf = eval( 'i.%sAngle[0]' % var )
      costheta = math.cos( leaf )
      
      for costhetaCut in costhetaCuts:
        # pc = ( 1. - costhetaCut )*p + ph*costheta
        # print 'p: %f, theta: %f, costheta: %f, pc: %f, ph: %f, costhl: %f' % ( p, leaf, costheta, pc, ph, costhetaCut )
        if p == 0. or p*( 1. - costhetaCut ) > ph*( costheta - costhetaCut ): continue
        nPassed[var][costhetaCut] += 1

  for var in Vars:
    for costhetaCut in costhetaCuts:
      rPassed[var][costhetaCut] = float(nPassed[var][costhetaCut])/float(nFiducial)
      # print nPassed[var][costhetaCut], nFiducial


  return n, nFiducial, nPassed, rPassed
# def selectEvents()

def optimizeSelection( mass, gamma, Vars, bgScale, nPassed, rPassed ):

  E = mass * gamma
  if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
  else: Eround = E
  sKey = 'e%s_m%s' %( Eround, mass )
  bKey = 'atmos_%s' % gamma
  # Scale the background events to 40kton*10 year exposure
  backgroundScale = bgScale* 40./28.705
  bestCut = {}
  bestEff = {}
  bestBkg = {}
  
  for var in Vars:
    SprimeMin = 10000.
    for costhetaCut in nPassed[sKey][var].keys():
      sEff = rPassed[sKey][var][costhetaCut]
      bEvents = float(nPassed[bKey][var][costhetaCut])*backgroundScale
      Sprime = 25./ (2.*sEff) + math.sqrt( 25.*bEvents/ (sEff*sEff) + 625./(4.*sEff*sEff) )
      if Sprime < SprimeMin:
        SprimeMin = Sprime
        bestCut[var] = costhetaCut
        bestEff[var] = sEff
        bestBkg[var] = bEvents
        bestBkgErr[var] = math.sqrt( float(nPassed[bKey][var][costhetaCut]) )*backgroundScale
  
  return bestCut, bestEff, bestBkg, bestBkgErr
  
# def optimizeSelection()

if __name__ == "__main__":
 
  Vars = [ 'Visible', 'VisibleNoN', 'LeadingParticle', 'LeadingParticleNoN',
           'SmearedReconstructable', 'SmearedReconstructableNoN',
           'LeadingSmearedReconstructable', 'LeadingSmearedReconstructableNoN',
           'SmearedVisible', 'SmearedVisibleNon', 'LeadingSmeared', 'LeadingSmearedNoN' ]

  parser = argparse.ArgumentParser( description = 'Optimize the selection.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  parser.add_argument( '-m', dest = 'bgScale', type = float, default = 1., 
                      help = 'The scale factor on the background events to account for additional background source.' )
  
  args = parser.parse_args()

  Masses     = [ 5, 10, 20, 40 ]
  Gammas     = [ 1.1, 1.25, 2, 10 ]
  hDict      = {}
  nTotal     = {}
  nFiducial  = {}
  nPassed    = {}
  passRate   = {}
  costhl     = 0.1
  costhh     = 0.95
  ncosth     = 18
  ph         = { 1.1: 1., 1.25: 1.6, 2: 3., 10: 80. }
  bestCut    = {}
  bestEff    = {}
  bestBkg    = {}
  bestBkgErr = {}

  # Create the output directory
  createDir( args.oDir )
  
  # Create the output text file
  txtName = '%s/1DLinear_Efficiency_scalar_bgScale%f.txt' % ( args.oDir, args.bgScale )
  txtFile = open( txtName, 'w' )
    

  # Background
  bFile = [ '%s/prodgenie_atmnu_max_dune10kt_gen_g4_NCFilter_reco_ana.root' % args.bDir, '%s/prodgenie_atmnu_min_dune10kt_gen_g4_NCFilter_reco_ana.root' % args.bDir ]
  bTree = getTree( bFile, 'MCParticles' )
  
  # Signal
  for Mass in Masses:
    for Gamma in Gammas:
      
      bKey = 'atmos_%s' % Gamma
      nTotal[bKey], nFiducial[bKey], nPassed[bKey], passRate[bKey] = selectEvents( bTree, Vars, costhl, costhh, ncosth, ph[Gamma] )
    
      
      if Mass in [ 5, 20, 40 ] and ( Gamma == 2 ):
        continue
      
      E = Mass * Gamma
      if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
      else: Eround = E
      sFile = [ '%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_reco_ana.root' %( args.sDir, str(Eround), str(Mass) ) ]
      sKey = 'e%s_m%s' %( Eround, Mass )
      sTree = getTree( sFile, 'MCParticles' )
      nTotal[sKey], nFiducial[sKey], nPassed[sKey], passRate[sKey] = selectEvents( sTree, Vars, costhl, costhh, ncosth, ph[Gamma] )
      # optimize the selection
      bestCut[sKey], bestEff[sKey], bestBkg[sKey], bestBkgErr[sKey] = optimizeSelection( Mass, Gamma, Vars, args.bgScale, nPassed, passRate )
      txtFile.write('Mass = %s GeV, E = %s GeV\n' % ( str(Mass), str(Eround) ) )
      
      for var in Vars:
        txtFile.write( '  %s: optimal cut at %f, signal efficieny = %f, background counts = %f, background uncertainty =  %f\n' %( var, bestCut[sKey][var], bestEff[sKey][var], bestBkg[sKey][var], bestBkgErr[sKey][var] ) )
