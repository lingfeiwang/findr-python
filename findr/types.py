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
"""Python interface for this library."""

import ctypes as c
import numpy as np
from .common import *
from .osdepend import typesizet,npulong

try: from exceptions import ValueError,NotImplementedError,NameError
except ImportError: pass

try:
	long
	isint=lambda v: type(v) is int or type(v) is long
except NameError:
	isint=lambda v: type(v) is int

class blockuc(c.Structure):
	_fields_=[("size",typesizet),("data",c_ubyte_p)]

class vectoruc(c.Structure):
	_fields_=[("size",typesizet),("stride",typesizet),("data",c_ubyte_p),("block",c.POINTER(blockuc)),("owner",c.c_int)]

class matrixuc(c.Structure):
	_fields_=[("size1",typesizet),("size2",typesizet),("tda",typesizet),("data",c_ubyte_p),("block",c.POINTER(blockuc)),("owner",c.c_int)]

class blockf(c.Structure):
	_fields_=[("size",typesizet),("data",ftype_p)]

class vectorf(c.Structure):
	_fields_=[("size",typesizet),("stride",typesizet),("data",ftype_p),("block",c.POINTER(blockf)),("owner",c.c_int)]

class matrixf(c.Structure):
	_fields_=[("size1",typesizet),("size2",typesizet),("tda",typesizet),("data",ftype_p),("block",c.POINTER(blockf)),("owner",c.c_int)]

class blockg(c.Structure):
	_fields_=[("size",typesizet),("data",gtype_p)]

class vectorg(c.Structure):
	_fields_=[("size",typesizet),("stride",typesizet),("data",gtype_p),("block",c.POINTER(blockg)),("owner",c.c_int)]

class matrixg(c.Structure):
	_fields_=[("size1",typesizet),("size2",typesizet),("tda",typesizet),("data",gtype_p),("block",c.POINTER(blockg)),("owner",c.c_int)]

def vectoruc_pin(d0,req=['A','C','W']):
	if len(d0.shape)!=1:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!='B':
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockuc(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(c_ubyte_p))
	sv=vectoruc(d.ctypes.shape[0],
		d.ctypes.strides[0]//d.dtype.type().nbytes,
		d.ctypes.data_as(c_ubyte_p),
		c.pointer(sb),
		0)
	return sv
def vectoruc_pout(d):
	raise NotImplementedError

def matrixuc_pin(d0,req=['A','C','W']):
	import numpy as np
	if len(d0.shape)!=2:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!='B':
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockuc(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(c_ubyte_p))
	sv=matrixuc(d.ctypes.shape[0],
		d.ctypes.shape[1],
		d.ctypes.strides[0]//nb,
		d.ctypes.data_as(c_ubyte_p),
		c.pointer(sb),
		0)
	return sv
def matrixuc_pout(d):
	raise NotImplementedError
	

def vectorf_pin(d0,req=['A','C','W']):
	if len(d0.shape)!=1:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockf(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(ftype_p))
	sv=vectorf(d.ctypes.shape[0],
		d.ctypes.strides[0]//d.dtype.type().nbytes,
		d.ctypes.data_as(ftype_p),
		c.pointer(sb),
		0)
	return sv
def vectorf_pout(d):
	raise NotImplementedError

def matrixf_pin(d0,req=['A','C','W']):
	import numpy as np
	if len(d0.shape)!=2:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockf(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(ftype_p))
	sv=matrixf(d.ctypes.shape[0],
		d.ctypes.shape[1],
		d.ctypes.strides[0]//nb,
		d.ctypes.data_as(ftype_p),
		c.pointer(sb),
		0)
	return sv
def matrixf_pout(d):
	raise NotImplementedError
	
def vectorg_pin(d0,req=['A','C','W']):
	if len(d0.shape)!=1:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockg(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(gtype_p))
	sv=vectorg(d.ctypes.shape,
		d.ctypes.strides[0]//d.dtype.type().nbytes,
		d.ctypes.data_as(gtype_p),
		c.pointer(sb),
		0)
	return sv
def vectorg_pout(d):
	raise NotImplementedError

def matrixg_pin(d0,req=['A','C','W']):
	import numpy as np
	if len(d0.shape)!=2:
		raise ValueError('Wrong input shape')
	if d0.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype')
	d=np.require(d0,requirements=req)
	nb=d.dtype.type().nbytes
	sb=blockg(typesizet(d.ctypes.strides[0]*d.shape[0]//nb),d.ctypes.data_as(gtype_p))
	sv=matrixg()
	sv=matrixg(d.ctypes.shape[0],
		d.ctypes.shape[1],
		d.ctypes.strides[0]//nb,
		d.ctypes.data_as(gtype_p),
		c.pointer(sb),
		0)
	return sv
def matrixg_pout(d):
	raise NotImplementedError
	

class varpipe:
	"""Variable pipe to process input and return variables for C functions"""
	def __init__(self,fin,fout):
		self.fin=fin
		self.fout=fout
	def pin(self,v):
		"""Process function input variables from python to C format"""
		return self.fin(v)
	def pout(self,v):
		"""Process function output variables from C to python format"""
		return self.fout(v)

vp_byte=varpipe(c.c_byte,lambda x:x)
vp_int=varpipe(c.c_int,lambda x:x)
vp_ulong=varpipe(c.c_ulong,lambda x:x)
vp_ulonglong=varpipe(typesizet,lambda x:x)
vp_sizet=varpipe(typesizet,lambda x:x)
vp_void=varpipe(lambda x:None,lambda x:None)
vp_char_p=varpipe(c.c_char_p,lambda x:x.value)
vp_vectorf=varpipe(vectorf_pin,vectorf_pout)
vp_matrixf=varpipe(matrixf_pin,matrixf_pout)
vp_vectorg=varpipe(vectorg_pin,vectorg_pout)
vp_matrixg=varpipe(matrixg_pin,matrixg_pout)
vp_vectoruc=varpipe(vectoruc_pin,vectoruc_pout)
vp_matrixuc=varpipe(matrixuc_pin,matrixuc_pout)
vp_vectorcf=varpipe(lambda x:vectorf_pin(x,req=['A','C']),vectorf_pout)
vp_matrixcf=varpipe(lambda x:matrixf_pin(x,req=['A','C']),matrixf_pout)
vp_vectorcg=varpipe(lambda x:vectorg_pin(x,req=['A','C']),vectorg_pout)
vp_matrixcg=varpipe(lambda x:matrixg_pin(x,req=['A','C']),matrixg_pout)

vmap={
	'byte':('i1',c.c_byte,vp_byte),
	'int':('i8',c.c_int,vp_int),
	'unsigned long':(npulong,c.c_ulong,vp_ulong),
	'size_t':('u8',typesizet,vp_sizet),
	'void':('u0',None,vp_void),
	'char*':(None,c.c_char_p,vp_char_p),
	'VECTORF*':(ftype_np,c.POINTER(vectorf),vp_vectorf),
	'MATRIXF*':(ftype_np,c.POINTER(matrixf),vp_matrixf),
	'VECTORUC*':('B',c.POINTER(vectoruc),vp_vectoruc),
	'MATRIXUC*':('B',c.POINTER(matrixuc),vp_matrixuc),
	'VECTORG*':(gtype_np,c.POINTER(vectorg),vp_vectorg),
	'MATRIXG*':(gtype_np,c.POINTER(matrixg),vp_matrixg),
	'const VECTORF*':(ftype_np,c.POINTER(vectorf),vp_vectorcf),
	'const MATRIXF*':(ftype_np,c.POINTER(matrixf),vp_matrixcf),
	'const VECTORG*':(gtype_np,c.POINTER(vectorg),vp_vectorcg),
	'const MATRIXG*':(gtype_np,c.POINTER(matrixg),vp_matrixcg)
	}
		
class cfunc:
	def __init__(self,lib,funcname,rettype='void',argtypes=[]):
		self.rt=rettype
		self.at=argtypes
		self.f=c.CFUNCTYPE(vmap[rettype][1],*[vmap[x][1] for x in argtypes])((funcname,lib))
	def __call__(self,*arg):
		a=[vmap[self.at[x]][2].pin(arg[x]) for x in range(len(arg))]
		return vmap[self.rt][2].pout(self.f(*a))







































