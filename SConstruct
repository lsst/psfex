# -*- python -*-
from lsst.sconsUtils import scripts
env = scripts.BasicSConstruct("psfex", versionModuleName="python/astromatic/%s/version.py")

import os,sys
env.Append(CCFLAGS = ['-DHAVE_CONFIG_H=1','-std=c99'])
if len(os.environ['PSFEX_DIR']) == 0:
	psfexdir = os.path.normpath(os.path.join(os.path.join(os.getcwd(),'../'),'lapack_functions'))
else:
	psfexdir = os.path.join(os.environ['PSFEX_DIR'],'lapack_functions')
env.Append(CPPPATH = [Dir('.').abspath,os.path.join(psfexdir,'include')])
