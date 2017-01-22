# Copyright 2016, 2017 Lingfei Wang
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
__all__=["auto","common","pij","osdepend","types"]
from . import pij

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
		return map(lambda x:pjoin(x,libfname),lp)
	def __init__(self,path=None,loglv=6,rs=0,nth=0):
		"""Links and initializes shared library.
		path:	Extra exact file location for shared library
		loglv:	Level of logging output. 1-3: Errors, 4-6: Warnings, 7-9: Infos, 10-12: Debug, 0: Default(6).
		rs:		Initial random seed. Default (0) indicates to use current time.
		nth:	Maximum number of parallel threads. Default (0) indicates to use automatically determined number of cores (not always correct).
		"""
		self.lib=None
		from exceptions import ValueError, OSError
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
				libname=ctypes.CFUNCTYPE(ctypes.c_char_p)(('lib_name',lib))
				libversion=ctypes.CFUNCTYPE(ctypes.c_char_p)(('lib_version',lib))
				ln=libname()
				lv=libversion()
				pv='.'.join(map(str,version))
				if((ln!=pkgname) or (lv!=pv)):
					print 'Located library '+ln+' '+lv+' different from python interface for '+pkgname+' '+pv+' at '+p+'. Skipped.'
					lib=None
					continue
				lib.lib_init(ctypes.c_ubyte(loglv),ctypes.c_ulong(rs),typesizet(nth))
			except:
				lib=None
				continue
			break
		if lib is None:
			raise OSError("Library not found at default path. Please install/update "+pkgname+' library and python interface, or set library path manually.')
		self.lib=lib
	def cfunc(self,*a,**ka):
		from exceptions import ValueError
		if self.lib is None:
			raise ValueError("Not initialized.")
		from .types import cfunc
		return cfunc(self.lib,*a,**ka)
	pij_gassist=pij.gassist
	pij_gassist_trad=pij.gassist_trad
	pijs_gassist=pij.gassists
	pij_rank=pij.rank
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

