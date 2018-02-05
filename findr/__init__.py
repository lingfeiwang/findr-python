# Copyright 2016-2018 Lingfei Wang
# 
# This file is part of Findr.
# 
# Findr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Findr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with Findr.  If not, see <http://www.gnu.org/licenses/>.
# 
"""Python interface.
For interface class, see: findr.lib.
For examples, see: findr.examples.
"""

import ctypes
__all__=["auto","common","pij","netr","osdepend","types"]
from . import pij, netr
try: from exceptions import ValueError,OSError
except ImportError: pass

class lib:
	@staticmethod	
	def default_libpaths():
		from .common import lpaths,libfname
		from os.path import join as pjoin
		from os import environ
		if 'WINDIR' in environ:
			lp=[pjoin(environ['WINDIR'],'System32')]+lpaths
		else:
			lp=lpaths
		return [pjoin(x,libfname) for x in lp]
	def __init__(self,path=None,loglv=6,rs=0,nth=0):
		"""Links and initializes shared library.
		path:	Extra exact file location for shared library
		loglv:	Level of logging output. 1-3: Errors, 4-6: Warnings, 7-9: Infos, 10-12: Debug, 0: Default(6).
		rs:		Initial random seed. Default (0) indicates to use current time.
		nth:	Maximum number of parallel threads. Default (0) indicates to use automatically determined number of cores (not always correct).
		"""
		self.lib=None
		import logging
		from .auto import pkgname,version
		from .osdepend import fdll,typesizet
		if type(rs) is not int or type(loglv) is not int or type(nth) is not int:
			raise ValueError('Wrong input type')
		if loglv<0 or loglv>12:
			raise ValueError('Wrong log level')
		if nth<0:
			raise ValueError('Wrong number of threads')
		paths=self.default_libpaths()
		if path:
			paths=[path]+paths
		for p in paths:
			try:
				lib=fdll(p)
				ln=ctypes.CFUNCTYPE(ctypes.c_char_p)(('lib_name',lib))().decode()
				lv=ctypes.CFUNCTYPE(ctypes.c_char_p)(('lib_version',lib))().decode()
				lv1=ctypes.CFUNCTYPE(typesizet)(('lib_version1',lib))()
				lv2=ctypes.CFUNCTYPE(typesizet)(('lib_version2',lib))()
				lv3=ctypes.CFUNCTYPE(typesizet)(('lib_version3',lib))()
				pv='.'.join(map(str,version))
				if((ln!=pkgname) or (lv1!=version[0]) or (lv2!=version[1])):
					logging.warning('Located library {} {} different from python interface for {} {} at {}. Skipped.'.format(ln,lv,pkgname,pv,p))
					lib=None
					continue
				elif lv3!=version[2]:
					logging.info('Located library {} {} different from python interface for {} {} at {}, but only at patch version. Loaded.'.format(ln,lv,pkgname,pv,p))
				lib.lib_init(ctypes.c_ubyte(loglv),ctypes.c_ulong(rs),typesizet(nth))
			except:
				lib=None
				continue
			break
		if lib is None:
			raise OSError("Library not found at default path. Please install/update "+pkgname+' library and python interface, or set library path manually.')
		self.lib=lib
	def cfunc(self,*a,**ka):
		if self.lib is None:
			raise ValueError("Not initialized.")
		from .types import cfunc
		return cfunc(self.lib,*a,**ka)
	pij_gassist=pij.gassist
	pij_gassist_trad=pij.gassist_trad
	pijs_gassist=pij.gassists
	pijs_gassist_pv=pij.gassists_pv
	pij_cassist=pij.cassist
	pij_cassist_trad=pij.cassist_trad
	pijs_cassist=pij.cassists
	pijs_cassist_pv=pij.cassists_pv
	pij_rank=pij.rank
	pij_rank_pv=pij.rank_pv
	netr_one_greedy=netr.one_greedy
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

