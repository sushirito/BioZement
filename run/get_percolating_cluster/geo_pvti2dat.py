#!/usr/bin/env python

from sys import stdout, argv, exit
#import numpy as np
import myvtk

if len(argv)<5:
    print
    print '  Usage: '
    print '          geo_pvti2dat.py filename resolution in-out-nodes wall-nodes'
    print
    exit();
    

infile = argv[1]
outfile = infile.split('pvti')[0] + 'dat'
resolution = float(argv[2])
in_out = int(argv[3])
wall = int(argv[4])

geo = myvtk.IO.read_pvti(infile)['mode']
if in_out>0 and wall>0:
    geo = geo[in_out:-in_out, wall:-wall, wall:-wall]
geo[geo>0.5] = 1.0
myvtk.IO.write_geo_dat(geo.astype('uint'), outfile, resolution)
print '   Wrote ' + outfile





