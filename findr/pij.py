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
"""Python interface"""

def _gassists_any(self,dg,dt,dt2,name,na=None,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with multiple tests, by converting log likelihoods into probabilities per A for all B.
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	name:	actual C function name to call 
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='=f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		For nodiag=False, because the function expects significant eQTLs, p1 always return 1.
		For nodiag=True, uses diagonal elements of p2.
		Consider replacing p1 with your own (1-FDR) from eQTL discovery.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A. The earlier one is preferred.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dg.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype for genotype data: dg.dtype.char is '+dtype.dtype.char+'!='+gtype_np)
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not (type(memlimit) is int or type(memlimit) is long):
		raise ValueError('Wrong memlimit type')
	if not (na is None or type(na) is int or type(na) is long):
		raise ValueError('Wrong na type')
	if na is not None and na<=0:
		raise ValueError('Input requires na>0.')
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
	
	arglist=['const MATRIXG*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','MATRIXF*','size_t','byte','size_t']
	names=name
	dgr=np.require(dg,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d4=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d5=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	args=[dgr,dtr,dt2r,d1,d2,d3,d4,d5,nvx,nd,memlimit]
	func=self.cfunc(names,rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p1':d1,'p2':d2,'p3':d3,'p4':d4,'p5':d5}
	return ans

def gassists(self,dg,dt,dt2,na=None,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with multiple tests, by converting log likelihoods into probabilities per A for all B.
	Probabilities are converted from likelihood ratios separately for each A. This gives better
	predictions when the number of secondary targets (dt2) is large. (Check program warnings.)
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='=f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		For nodiag=False, because the function expects significant eQTLs, p1 always return 1.
		For nodiag=True, uses diagonal elements of p2.
		Consider replacing p1 with your own (1-FDR) from eQTL discovery.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A. The earlier one is preferred.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	return _gassists_any(self,dg,dt,dt2,"pijs_gassist",na=na,nodiag=nodiag,memlimit=memlimit)

def _gassist_any(self,dg,dt,dt2,name,na=None,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with the recommended combination of multiple tests.
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	name:	actual C function name to call 
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dg.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype for genotype data: dg.dtype.char is '+dtype.dtype.char+'!='+gtype_np)
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not (type(memlimit) is int or type(memlimit) is long):
		raise ValueError('Wrong memlimit type')
	if not (na is None or type(na) is int or type(na) is long):
		raise ValueError('Wrong na type')
	if na is not None and na<=0:
		raise ValueError('Input requires na>0.')
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

	func=self.cfunc(name,rettype='int',argtypes=['const MATRIXG*','const MATRIXF*','const MATRIXF*','MATRIXF*','size_t','byte','size_t'])
	d=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	dgr=np.require(dg,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	ret=func(dgr,dtr,dt2r,d,nvx,nd,memlimit)
	ans={'ret':ret,'p':d}
	return ans

def gassist(self,dg,dt,dt2,na=None,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method,
	with the recommended combination of multiple tests.
	Probabilities are converted from likelihood ratios separately for each A. This gives better
	predictions when the number of secondary targets (dt2) is large. (Check program warnings.)
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	return _gassist_any(self,dg,dt,dt2,"pij_gassist",na=na,nodiag=nodiag,memlimit=memlimit)

def gassist_trad(self,dg,dt,dt2,na=None,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with genotype data assisted method with traditional test.
	WARNING: This is not and is not intended as a loyal reimplementation of Trigger. This test does not include p1.
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	"""
	return _gassist_any(self,dg,dt,dt2,"pij_gassist_trad",na=na,nodiag=nodiag,memlimit=memlimit)

def rank(self,dt,dt2,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i correlating with gene j by converting log likelihoods into probabilities per A for all B.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, a subset of, or a superset of dt. When dt2 is a superset of (or identical with) dt, dt2 must be arranged to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and set parameter nodiag = 1. Similarly if dt2 is a subset of dt.
	nodiag:	skip diagonal regulations, i.e. regulation A--B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for A--B.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not (type(memlimit) is int or type(memlimit) is long):
		raise ValueError('Wrong memlimit type')
	ng=dt.shape[0]
	nt=dt2.shape[0]
	ns=dt.shape[1]
	nd=1 if nodiag else 0
	
	if dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')

	dp=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	arglist=['const MATRIXF*','const MATRIXF*','MATRIXF*','byte','size_t']
	args=[dtr,dt2r,dp,nd,memlimit]
	names='pij_rank'
	func=self.cfunc(names,rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p':dp}
	return ans






























