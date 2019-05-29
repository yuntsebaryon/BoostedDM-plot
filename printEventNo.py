#!/usr/bin/env python

import argparse
import os
import ROOT

def getTree( fNames, tName ):

  t = ROOT.TChain( tName )
  for fName in fNames:
    t.AddFile( fName )

  t.SetDirectory(0)
  return t
# getTree()

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser( description = 'Print the event numbers')
  parser.add_argument( '-i', dest = 'iFile', type = str, help = 'The input files.' )
  parser.add_argument( '-o', dest = 'oFile', type = str, help = 'The output text file.' )
  
  args = parser.parse_args()
  
  if os.path.isfile( args.oFile ):
    os.remove( args.oFile )
  
  iFiles = [ args.iFile ]
  tree = getTree( iFiles, 'MCParticles' )
  
  oFile = open( args.oFile, 'w' )
  
  for i in tree:
    if i.nVisible == 0:
      oFile.write( '%d\n' % i.Event )
