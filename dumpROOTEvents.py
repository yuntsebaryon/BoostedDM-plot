#!/usr/bin/env python

import argparse
import ROOT
import math
import numpy

def getTree( fNames, tName ):

  t = ROOT.TChain( tName )
  for fName in fNames:
    t.AddFile( fName )

  t.SetDirectory(0)
  return t
# getTree()

def printEvents( tree, outFilename, Vars ):

  f = open( outFilename, 'w' )

  n = tree.GetEntries()

  f.write( '#' )
  for var in Vars:
    f.write( '%s    ' % var )
  f.write ( '\n' )

  for i in tree:
    
    doChangeLine = True
    for var in Vars:
      if ( i.isIn10kton == 0 ):
        doChangeLine = False
        continue
      leaf = eval('i.%s[0]' % var)
      if ( numpy.isnan( leaf ) ):
        costheta = -1.
      else:
        costheta = math.cos( leaf )
      f.write( '  %f                     ' % costheta )
    if doChangeLine: f.write( '\n' )

  f.close()
# def printEvents()

if __name__ == "__main__":


  Vars = [ 'SmearedReconstructableAngle', 'SmearedReconstructableNoNAngle' ]
           
  parser = argparse.ArgumentParser( description = 'Make angular distributions.')
  parser.add_argument( '-s', dest = 'sDir', type = str, help = 'The directory of the input SIGNAL files.' )
  parser.add_argument( '-b', dest = 'bDir', type = str, help = 'The directory of the input BACKGROUND files.' )

  args = parser.parse_args()


  Masses   = [ 5, 10, 20, 40 ]
  Gammas   = [ 1.1, 1.25, 10 ]

  # Background
  bFile = [ '%s/prodgenie_atmnu_max_dune10kt_gen_g4_NCFilter_reco_ana.root' % args.bDir, '%s/prodgenie_atmnu_min_dune10kt_gen_g4_NCFilter_reco_ana.root' % args.bDir ]
  bTree = getTree( bFile, 'MCParticles' )
  bOut = '%s/prodgenie_atmnu_maxmin_dune10kt_gen_g4_NCFilter_reco_ana.dat' % args.bDir
  print 'Atmospheric neutrino...'
  printEvents( bTree, bOut, Vars )

  for Mass in Masses:
    for Gamma in Gammas:

      E = Mass * Gamma
      if E in [ 11., 22., 25., 44., 50. ]: Eround = int(E)
      else: Eround = E
      if E in [ 5.5, 6.25, 12.5 ]:
        Estr = str(E)
        Estr.replace( '.', 'p' )
      else: Estr = str(Eround)
      sFile = ['%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_reco_ana.root' %( args.sDir, Estr, str(Mass) ) ]
      sKey = 'e%s_m%s' %( Eround, Mass )
      sTree = getTree( sFile, 'MCParticles' )
      sOut = '%s/dune_scalar_e%s_m%s_g1_z1.0_Gen_g4_reco_ana.dat' % ( args.sDir, Estr, str(Mass) )
      print 'M = %d, E = %f' %( Mass, E )
      printEvents( sTree, sOut, Vars )
