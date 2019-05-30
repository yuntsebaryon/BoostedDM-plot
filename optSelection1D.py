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

def selectEvents( tree, AngularVars, isSignal, costhl, costhh, ncosth ):

  n = tree.GetEntries()
  nFiducial = 0
  nPassed = {}
  rPassed = {}
  cutValues = numpy.arange( 0.2, 1, 0.05 )
  backgroundScale = 40./28.705
  
  costhetaCuts = numpy.linspace( costhl, costhh, ncosth )
  
  for var in AngularVars:
    nPassed[var] = {}
    rPassed[var] = {}
    for costhetaCut in costhetaCuts:
      nPassed[var][costhetaCut] = 0
      rPassed[var][costhetaCut] = 0
  
  for i in tree:
    if i.isIn10kton == 0: continue
    nFiducial += 1
    
    for var in AngularVars:
      
      leaf = eval('i.%s[0]' % var)
      costheta = math.cos( leaf )
      
      for costhetaCut in costhetaCuts:
        if costheta < costhetaCut: continue
        nPassed[var][costhetaCut] += 1

  for var in AngularVars:
    for costhetaCut in costhetaCuts:
      rPassed[var][costhetaCut] = float(nPassed[var][costhetaCut])/float(nFiducial)
      # print nPassed[var][costhetaCut], nFiducial

  # if not isSignal:
  #   nFiducial *= backgroundScale
  #   for costheta in nPassed.values():
  #     for passed in costheta.values():
  #       passed *= backgroundScale

  return n, nFiducial, nPassed, rPassed
# def selectEvents()

def optimizeSelection( mass, gamma, AngularVars, nPassed, rPassed ):

  E = mass * gamma
  if E == 25. or E == 50.: Eround = int(E)
  else: Eround = E
  sKey = 'e%s_m%s' %( Eround, mass )
  # Scale the background events to 40kton*10 year exposure
  backgroundScale = 40./28.705
  bestCut = {}
  bestEff = {}
  bestBkg = {}
  
  for var in AngularVars:
    SprimeMin = 10000.
    for costhetaCut in nPassed[sKey][var].keys():
      sEff = rPassed[sKey][var][costhetaCut]
      bEvents = float(nPassed['atmos'][var][costhetaCut])*backgroundScale
      Sprime = 25./ (2.*sEff) + math.sqrt( 25.*bEvents/ (sEff*sEff) + 625./(4.*sEff*sEff) )
      if Sprime < SprimeMin:
        SprimeMin = Sprime
        bestCut[var] = costhetaCut
        bestEff[var] = sEff
        bestBkg[var] = bEvents
        bestBkgErr[var] = math.sqrt( float(nPassed['atmos'][var][costhetaCut]) )*backgroundScale
  
  return bestCut, bestEff, bestBkg, bestBkgErr
  
# def optimizeSelection()

if __name__ == "__main__":
 
  AngularVars = [ 'VisibleAngle', 'VisibleNoNAngle', 'LeadingParticleAngle', 'LeadingParticleNoNAngle',
                  'SmearedReconstructableAngle', 'SmearedReconstructableNoNAngle',
                  'LeadingSmearedReconstructableAngle', 'LeadingSmearedReconstructableNoNAngle' ]

  parser = argparse.ArgumentParser( description = 'Optimize the selection.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )
  parser.add_argument( '-o', dest = 'oDir', type = str, help = 'The directory of the output plots.' )
  
  args = parser.parse_args()

  Masses     = [ 5, 10, 20, 40 ]
  Gammas     = [ 1.25, 2, 10 ]
  hDict      = {}
  nTotal     = {}
  nFiducial  = {}
  nPassed    = {}
  passRate   = {}
  costhl     = 0.2
  costhh     = 1.
  ncosth     = 17
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
  nTotal['atmos'], nFiducial['atmos'], nPassed['atmos'], passRate['atmos'] = selectEvents( bTree, AngularVars, False, costhl, costhh, ncosth )

  # Signal
  for Mass in Masses:
    for Gamma in Gammas:
      
      if Mass in [ 5, 20, 40 ] and ( Gamma == 2 ):
        continue
      
      E = Mass * Gamma
      if E == 25. or E == 50.: Eround = int(E)
      else: Eround = E
      sFile = [ '%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_RecoSmear_ana.root' %( args.sDir, str(Eround), str(Mass) ) ]
      sKey = 'e%s_m%s' %( Eround, Mass )
      sTree = getTree( sFile, 'MCParticles' )
      nTotal[sKey], nFiducial[sKey], nPassed[sKey], passRate[sKey] = selectEvents( sTree, AngularVars, True, costhl, costhh, ncosth )
      # optimize the selection
      bestCut[sKey], bestEff[sKey], bestBkg[sKey], bestBkgErr[sKey] = optimizeSelection( Mass, Gamma, AngularVars, nPassed, passRate )
      txtFile.write('Mass = %s GeV, E = %s GeV\n' % ( str(Mass), str(Eround) ) )
      
      for var in AngularVars:
        txtFile.write( '  %s: optimal cut at %f, signal efficieny = %f, background counts = %f, background uncertainty =  %f\n' %( var, bestCut[sKey][var], bestEff[sKey][var], bestBkg[sKey][var], bestBkgErr[sKey][var] ) )
