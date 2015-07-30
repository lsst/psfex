# -*- python -*-
from lsst.sconsUtils import scripts
env = scripts.BasicSConstruct("psfex", versionModuleName="python/astromatic/%s/version.py")

import os
if os.environ.has_key("PSFEX_DIR"):
    env.Append(CPPPATH = [Dir('.').abspath,os.path.join(os.path.realpath(os.getcwd()),'lapack_functions/include')]) # needed for config.h.  N.b. doesn't work with _wrap.cc
env.Append(CCFLAGS = ['-DHAVE_CONFIG_H=1','-std=c99'])
