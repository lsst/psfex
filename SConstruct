# -*- python -*-
from lsst.sconsUtils import scripts
env = scripts.BasicSConstruct("psfex", versionModuleName="python/astromatic/%s/version.py")

env.Append(CCFLAGS = ['-DHAVE_CONFIG_H=1','-std=c99'])
env.Append(CPPPATH = [Dir('#')])
