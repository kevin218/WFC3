
from numpy import get_include
import os
import re
from distutils.core import setup, Extension

files = os.listdir('c_code/')
#this will filter the results for just the c files
files = list(filter(lambda x: not re.search('.+[.]c[~]$',x),files))
files = list(filter(lambda x: not re.search('[.#].+[.]c$',x),files))
files.remove(".svn")

ext_mod = []
inc = [get_include()]

for i in range(len(files)):
    exec('mod'+str(i)+'=Extension("'+files[i].rstrip('.c')+'",sources=["c_code/'+files[i]+'"],include_dirs=inc,extra_compile_args=["-fopenmp"],extra_link_args=["-lgomp"])')

    exec('ext_mod.append(mod'+str(i)+')')


setup (name='models_c',version='1.0',description='Models in c for mcmc',
       ext_modules = ext_mod)
