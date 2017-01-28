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

def one_greedy(self,dp,namax=None,nimax=None,nomax=None):
	"""Reconstructs a directed acyclic graph according to prior information of edge significance.
	This function first ranks all edges and introduce the most significant one by one, avoiding
	those that would create a loop. Optional constraints on the maximum total number of edges,
	the number of incoming or outgoing edges for every gene can be specified.
	dp:	numpy.ndarray(nt,nt,dtype=ftype(='f4' by default))
		Prior information of edge significance levels. Entry dp[i,j] is significance of edge i to j. 
		A larger values indicates the edge's presence is more probable.
		One option to obtain the prior information is to use pairwise inference methods in findr.
	dt2:numpy.ndarray(nt2,ns,dtype=ftype(='=f4' by default)) Gene expression data for B.
		dt2 has the same format as dt, and can be identical with, different from, or a superset of dt.
		When dt2 is a superset of (or identical with) dt, dt2 must be arranged
		to be identical with dt at its upper submatrix, i.e. dt2[:nt,:]=dt, and
		set parameter nodiag = 1.
	namax:	Constraint on the maximum total number of edges in the reconstructed network.
	nimax:	Constraint on the maximum number of incoming edges for each node in the reconstructed
		network.
	nomax:	Constraint on the maximum number of outgoing edges for each node in the reconstructed
		network.
	Return:	dictionary with following keys:
	ret:0 iff execution succeeded.
	net:	numpy.ndarray((nt,nt),dtype=bool). The reconstructed direct acyclic graph or network
		net[i,j]=True if an edge from i to j exists in the reconstructed network, and False otherwise.
	ftype and gtype can be found in auto.py.
	"""
	from exceptions import ValueError
	if self.lib is None:
		raise ValueError("Not initialized.")
	import numpy as np
	from .auto import ftype_np
	if dp.dtype.char!=ftype_np:
		raise ValueError('Wrong input dtype for prior matrix')
	if len(dp.shape)!=2:
		raise ValueError('Wrong input shape')
	if not (namax is None or type(namax) is int or type(namax) is long):
		raise ValueError('Wrong namax type')
	if namax is not None and namax<=0:
		raise ValueError('Input requires namax>0.')
	if namax is None:
		namax=-1
	if not (nimax is None or type(nimax) is int or type(nimax) is long):
		raise ValueError('Wrong nimax type')
	if nimax is not None and nimax<=0:
		raise ValueError('Input requires nimax>0.')
	if nimax is None:
		nimax=-1
	if not (nomax is None or type(nomax) is int or type(nomax) is long):
		raise ValueError('Wrong nomax type')
	if nomax is not None and nomax<=0:
		raise ValueError('Input requires nomax>0.')
	if nomax is None:
		nomax=-1
	
	nt=dp.shape[0]
	if nt==0:
		raise ValueError('Invalid prior dimension')
	if dp.shape[1]!=nt:
		raise ValueError('Wrong input shape')
	if np.isnan(dp).sum()>0:
		raise ValueError('NaN found.')
	func=self.cfunc('netr_one_greedy',rettype='size_t',argtypes=['const MATRIXF*','MATRIXUC*','size_t','size_t','size_t'])
	d=np.require(np.zeros((nt,nt),dtype='u1'),requirements=['A','C','O','W'])
	dpr=np.require(dp,requirements=['A','C','O','W'])
	ret=func(dpr,d,namax,nimax,nomax)
	d=d.astype(bool)
	ret=(ret==0)
	ans={'ret':ret,'net':d}
	return ans






























