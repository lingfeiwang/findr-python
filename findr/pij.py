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
"""Python interface"""

try: from exceptions import ValueError
except ImportError: pass

def gassists_pv(self,dg,dt,dt2,na=None,memlimit=-1):
	"""Calculates p-values of gene i regulating gene j with genotype data assisted method with multiple tests.
	dg:	numpy.ndarray(nt,ns,dtype=gtype(='u1' by default)) Genotype data.
		Entry dg[i,j] is genotype i's value for sample j.
		Each value must be among 0,1,...,na.
		Genotype i must be best (and significant) eQTL of gene i (in dt).
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
		Genotype i (in dg) must be best (and significant) eQTL of gene i.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
	na:	Number of alleles the species have. It determintes the maximum number of values each genotype can take. When unspecified, it is automatically
		determined as the maximum of dg.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='f4' by default)). P-values for LLR of test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A.
	For more information on tests, see paper.
	ftype and gtype can be found in auto.py.
	
	Example: see findr.examples.geuvadis6
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	from .types import isint
	if dg.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype for genotype data: dg.dtype.char is '+dg.dtype.char+'!='+gtype_np)
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	if not (na is None or isint(na)):
		raise ValueError('Wrong na type')
	if na is not None and na<=0:
		raise ValueError('Input requires na>0.')
	ng=dg.shape[0]
	nt=dt2.shape[0]
	ns=dg.shape[1]
	nvx=na+1 if na else dg.max()+1
	
	if nvx<2:
		raise ValueError('Invalid genotype values')
	if dt.shape!=dg.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')
	
	arglist=['const MATRIXG*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','MATRIXF*','size_t','size_t']
	dgr=np.require(dg,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d4=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d5=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	args=[dgr,dtr,dt2r,d1,d2,d3,d4,d5,nvx,memlimit]
	func=self.cfunc('pijs_gassist_pv',rettype='int',argtypes=arglist)
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
	
	Example: see findr.examples.geuvadis4
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	from .types import isint
	if dg.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype for genotype data: dg.dtype.char is '+dg.dtype.char+'!='+gtype_np)
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	if not (na is None or isint(na)):
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
	dgr=np.require(dg,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d4=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d5=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	args=[dgr,dtr,dt2r,d1,d2,d3,d4,d5,nvx,nd,memlimit]
	func=self.cfunc("pijs_gassist",rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p1':d1,'p2':d2,'p3':d3,'p4':d4,'p5':d5}
	return ans

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
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	from .types import isint
	if dg.dtype.char!=gtype_np:
		raise ValueError('Wrong input dtype for genotype data: dg.dtype.char is '+dg.dtype.char+'!='+gtype_np)
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dg.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	if not (na is None or isint(na)):
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
	
	Example: see findr.examples.geuvadis2, findr.examples.geuvadis3
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
	
	Example: see findr.examples.geuvadis2, findr.examples.geuvadis3 (same format)
	"""
	return _gassist_any(self,dg,dt,dt2,"pij_gassist_trad",na=na,nodiag=nodiag,memlimit=memlimit)

def cassists_pv(self,dc,dt,dt2,memlimit=-1):
	"""Calculates p-values of gene i regulating gene j with continuous anchor data assisted method with multiple tests.
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
	memlimit:	The approximate memory usage limit in bytes for the library. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='f4' by default)). P-values for LLR of test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='f4' by default)). P-values for LLR of test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	
	Example: see findr.examples.geuvadis6 (similar format)
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np
	from .types import isint
	if dc.dtype.char!=ftype_np or dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dc.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	ng=dc.shape[0]
	nt=dt2.shape[0]
	ns=dc.shape[1]
	
	if dt.shape!=dc.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dc).sum()+np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')
	
	arglist=['const MATRIXF*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','MATRIXF*','size_t']

	dcr=np.require(dc,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d4=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d5=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	args=[dcr,dtr,dt2r,d1,d2,d3,d4,d5,memlimit]
	func=self.cfunc('pijs_cassist_pv',rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p1':d1,'p2':d2,'p3':d3,'p4':d4,'p5':d5}
	return ans
	
def _cassists_any(self,dc,dt,dt2,name,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with continuous anchor data assisted method,
	with multiple tests, by converting log likelihoods into probabilities per A for all B.
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	name:	actual C function name to call 
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='=f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		For nodiag=False, because the function expects significant anchors, p1 always return 1.
		For nodiag=True, uses diagonal elements of p2.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A. The earlier one is preferred.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np
	from .types import isint
	if dc.dtype.char!=ftype_np or dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dc.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	ng=dc.shape[0]
	nt=dt2.shape[0]
	ns=dc.shape[1]
	nd=1 if nodiag else 0
	
	if dt.shape!=dc.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dc).sum()+np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')
	
	arglist=['const MATRIXF*','const MATRIXF*','const MATRIXF*','VECTORF*','MATRIXF*','MATRIXF*','MATRIXF*','MATRIXF*','byte','size_t']
	names=name
	dcr=np.require(dc,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	d1=np.require(np.zeros(ng,dtype=dt.dtype),requirements=['A','C','O','W'])
	d2=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d3=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d4=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	d5=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	args=[dcr,dtr,dt2r,d1,d2,d3,d4,d5,nd,memlimit]
	func=self.cfunc(names,rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p1':d1,'p2':d2,'p3':d3,'p4':d4,'p5':d5}
	return ans

def cassists(self,dc,dt,dt2,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with continuous data assisted method,
	with multiple tests, by converting log likelihoods into probabilities per A for all B.
	Probabilities are converted from likelihood ratios separately for each A. This gives better
	predictions when the number of secondary targets (dt2) is large. (Check program warnings.)
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p1:	numpy.ndarray(nt,dtype=ftype(='=f4' by default)). Probability for test 1.
		Test 1 calculates E(A)->A v.s. E(A)  A. The earlier one is preferred.
		For nodiag=False, because the function expects significant anchors, p1 always return 1.
		For nodiag=True, uses diagonal elements of p2.
	p2:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 2.
		Test 2 calculates E(A)->A--B with E(A)->B v.s. E(A)->A<-B. The earlier one is preferred.
	p3:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 3.
		Test 3 calculates E(A)->A--B with E(A)->B v.s. E(A)->A->B. The latter one is preferred.
	p4:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 4.
		Test 4 calculates E(A)->A--B with E(A)->B v.s. E(A)->A  B. The earlier one is preferred.
	p5:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for test 5.
		Test 5 calculates E(A)->A--B with E(A)->B v.s. B<-E(A)->A. The earlier one is preferred.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	
	Example: see findr.examples.geuvadis4 (similar format)
	"""
	return _cassists_any(self,dc,dt,dt2,"pijs_cassist",nodiag=nodiag,memlimit=memlimit)

def _cassist_any(self,dc,dt,dt2,name,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with continuous data assisted method,
	with the recommended combination of multiple tests.
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	name:	actual C function name to call 
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np
	from .types import isint
	if dc.dtype.char!=ftype_np or dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dc.shape)!=2 or len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	ng=dc.shape[0]
	nt=dt2.shape[0]
	ns=dc.shape[1]
	nd=1 if nodiag else 0
	
	if dt.shape!=dc.shape or dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dc).sum()+np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')

	func=self.cfunc(name,rettype='int',argtypes=['const MATRIXF*','const MATRIXF*','const MATRIXF*','MATRIXF*','byte','size_t'])
	d=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	dcr=np.require(dc,requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	ret=func(dcr,dtr,dt2r,d,nd,memlimit)
	ans={'ret':ret,'p':d}
	return ans

def cassist(self,dc,dt,dt2,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with continuous data assisted method,
	with the recommended combination of multiple tests.
	Probabilities are converted from likelihood ratios separately for each A. This gives better
	predictions when the number of secondary targets (dt2) is large. (Check program warnings.)
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	
	Example: see findr.examples.geuvadis5
	"""
	return _cassist_any(self,dc,dt,dt2,"pij_cassist",nodiag=nodiag,memlimit=memlimit)

def cassist_trad(self,dc,dt,dt2,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i regulating gene j with continuous data assisted method with traditional test.
	dc:	numpy.ndarray(nt,ns,dtype=ftype(='f4' by default)) Continuous anchor data.
		Entry dc[i,j] is anchor i's value for sample j.
		Anchor i is used to infer the probability of gene i -> any other gene.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	nodiag:	skip diagonal regulations, i.e. regulation A->B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will be split into smaller chunks. If the memory limit is smaller than minimum required, calculation can fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)).
		Probability function from for recommended combination of multiple tests.
	For more information on tests, see paper.
	ftype can be found in auto.py.
	
	Example: see findr.examples.geuvadis5 (same format)
	"""
	return _cassist_any(self,dc,dt,dt2,"pij_cassist_trad",nodiag=nodiag,memlimit=memlimit)

def rank_pv(self,dt,dt2,memlimit=-1):
	"""Calculates p-values of gene i correlating with gene j by converting log likelihoods into probabilities per A for all B.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, a subset of, or a superset of dt.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). P-values for A--B.
	ftype and gtype can be found in auto.py.
	
	Example: see findr.examples.geuvadis1 (similar format)
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	from .types import isint
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if not isint(memlimit):
		raise ValueError('Wrong memlimit type')
	ng=dt.shape[0]
	nt=dt2.shape[0]
	ns=dt.shape[1]
	
	if dt2.shape[1]!=ns:
		raise ValueError('Wrong input shape')
	if np.isnan(dt).sum()+np.isnan(dt2).sum()>0:
		raise ValueError('NaN found.')

	dp=np.require(np.zeros((ng,nt),dtype=dt.dtype),requirements=['A','C','O','W'])
	dtr=np.require(dt,requirements=['A','C','O','W'])
	dt2r=np.require(dt2,requirements=['A','C','O','W'])
	arglist=['const MATRIXF*','const MATRIXF*','MATRIXF*','size_t']
	args=[dtr,dt2r,dp,memlimit]
	func=self.cfunc('pij_rank_pv',rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p':dp}
	return ans

def rank(self,dt,dt2,nodiag=False,memlimit=-1):
	"""Calculates probability of gene i correlating with gene j by converting log likelihoods into probabilities per A for all B.
	dt:	numpy.ndarray(nt,ns,dtype=ftype(='=f4' by default)) Gene expression data for A
		Entry dt[i,j] is gene i's expression level for sample j.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, a subset of, or a superset of dt. When dt2 is a superset of (or identical with) dt, dt2 must be arranged to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and set parameter nodiag = 1. Similarly if dt2 is a subset of dt.
	nodiag:	skip diagonal regulations, i.e. regulation A--B for A=B.
		This should be set to True when A is a subset of B and aligned correspondingly.
	memlimit:	The approximate memory usage limit in bytes for the library.  For datasets require a larger memory, calculation will fail with an error message. memlimit=0 defaults to unlimited memory usage.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	p:	numpy.ndarray((nt,nt2),dtype=ftype(='=f4' by default)). Probability for A--B.
	ftype and gtype can be found in auto.py.
	
	Example: see findr.examples.geuvadis1
	"""
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np,gtype_np
	from .types import isint
	if dt.dtype.char!=ftype_np or dt2.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for gene expression data')
	if len(dt.shape)!=2 or len(dt2.shape)!=2:
		raise ValueError('Wrong input shape')
	if type(nodiag) is not bool:
		raise ValueError('Wrong nodiag type')
	if not isint(memlimit):
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
	func=self.cfunc('pij_rank',rettype='int',argtypes=arglist)
	ret=func(*args)
	ans={'ret':ret,'p':dp}
	return ans


































