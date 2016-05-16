# Copyright 2016 Lingfei Wang
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
"""Python interface"""


def gassist_a(self,dg,dt,dt2,na=None,nodiag=False):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with the model E(A)->A->B, by converting log likelihoods into probabilities per A for all B.
	dg:	numpy.ndarray(ng,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(ng,ns,dtype=ftype(='<f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt,ns,dtype=ftype(='<f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:ng,:]=dt, and
		set parameter nodiag = 1.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(ng,dtype=ftype(='<f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		Since the program already expects significant eQTLs as input, p1 is
		set to constant 1 for all A.
	p2b:numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 2 bold.
		Test 2 bold calculates E(A)->A->B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p2c:numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 2 conservative.
		Test 2 conservative calculates E(A)->A->B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dg.dtype.char!=gtype_np or dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	ng=dg.shape[0]
	nt=dt2.shape[0]
	ns=dg.shape[1]
	nvx=na+1 if na else dg.max()+1
	nd=1 if nodiag else 0
	
	if nvx<2:
		raise ValueError('Invalid genotype values')
	if dt.shape!=dg.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')

	func=self.cfunc('pijs_gassist_a',rettype='int',argtypes=['const MATRIXG*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','size_t','byte'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2b=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d2c=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	ret=func(dg,dt,dt2,d1,d2b,d2c,d3,nvx,nd)
	ans={'ret':ret,'p1':d1,'p2b':d2b,'p2c':d2c,'p3':d3}
	return ans

def gassist_tot(self,dg,dt,dt2,na=None,nodiag=False):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with the model E(A)->A->B, by converting log likelihoods into probabilities altogether.
	dg:	numpy.ndarray(ng,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,nv-1.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(ng,ns,dtype=ftype(='<f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt,ns,dtype=ftype(='<f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:ng,:]=dt, and
		set parameter nodiag = 1.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(ng,dtype=ftype(='<f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		Since the program already expects significant eQTLs as input, p1 is
		set to constant 1 for all A.
	p2b:numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 2 bold.
		Test 2 bold calculates E(A)->A->B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p2c:numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 2 conservative.
		Test 2 conservative calculates E(A)->A->B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((ng,nt),dtype=ftype(='<f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dg.dtype.char!=gtype_np or dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	ng=dg.shape[0]
	nt=dt2.shape[0]
	ns=dg.shape[1]
	nvx=na+1 if na else dg.max()+1
	nd=1 if nodiag else 0
	if nvx<2:
		raise ValueError('Invalid genotype values')
	if dt.shape!=dg.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')
		
	func=self.cfunc('pijs_gassist_tot',rettype='int',argtypes=['const MATRIXG*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','size_t','byte'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2b=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d2c=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	ret=func(dg,dt,dt2,d1,d2b,d2c,d3,nvx,nd)
	ans={'ret':ret,'p1':d1,'p2b':d2b,'p2c':d2c,'p3':d3}
	return ans

def rank_a(self,dt,dt2,nodiag=False):
	"""Calculates probability of gene i correlating with gene j by converting log likelihoods into probabilities per A for all B.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='<f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='<f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, a subset of, or a superset of dt. When dt2 is a superset of (or identical with) dt, dt2 must be arranged to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and set parameter nodiag = 1. Similarly if dt2 is a subset of dt.
	nodiag:	skip diagonal regulations, i.e. regulation A--B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='<f4' by default)). Probability for A--B.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype')
	if len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	ng=dt.shape[0]
	nt=dt2.shape[0]
	ns=dt.shape[1]
	nd=1 if nodiag else 0
	
	if dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')

	func=self.cfunc('pij_rank_a',rettype='int',argtypes=['const MATRIXF*','const MATRIXF*','MATRIXF*','byte'])
	dp=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	ret=func(dt,dt2,dp,nd)
	ans={'ret':ret,'p':dp}
	return ans



















