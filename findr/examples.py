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
"""Examples of findr. For geuvadis dataset, see findr.examples.geuvadis1, findr.examples.geuvadis2, ... , findr.examples.geuvadis4"""

def load_geuvadis_data():
	"""This function loads downsampled data files from the Geuvadis study (Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506-511 (2013)), including expression levels of 10 miRNAs and 3000 genes for 360 European individuals. Among them, all miRNAs and 1000 genes have significant cis-eQTLs, whose haplotypes are also included. File data formats follow Findr's binary interface input/output requirement. A description of each file is available below:
	dmi.dat:	Expression levels of 10 miRNAs
	dgmi.dat:	Haplotypes of cis-eQTLs of 10 miRNAs
	dc.dat:		Continuous causal anchors for demonstration purposes, simulated from adding continuous noise to dgmi.dat
	dt.dat:		Expression levels of 1000 genes that have cis-eQTLs
	dt2.dat:	Expression levels of 3000 genes
	dgt.dat:	Haplotypes of cis-eQTLs of 1000 genes
	namest.txt:	3000 gene names"""
	from os.path import dirname,join
	from .auto import gtype_np,ftype_np
	import numpy as np
	def getdata(name,dtype,shape):
		d=join(dirname(__file__),'data','geuvadis',name)
		d=np.fromfile(d,dtype=dtype)
		d=d.reshape(*shape)
		return d
	
	ans={'dc':getdata('dc.dat',ftype_np,(10,360)),
	'dgmi':getdata('dgmi.dat',gtype_np,(10,360)),
	'dmi':getdata('dmi.dat',ftype_np,(10,360)),
	'dgt':getdata('dgt.dat',gtype_np,(1000,360)),
	'dt':getdata('dt.dat',ftype_np,(1000,360)),
	'dt2':getdata('dt2.dat',ftype_np,(3000,360))}
	
	f=open(join(dirname(__file__),'data','geuvadis','namest.txt'),'r')
	namest=[x.strip('\r\n') for x in f.readlines()]
	f.close()
	ans['namest']=namest
	return ans

def runcode(code):
	"""Run the given code line by line with printing, as list of lines, and return variable 'ans'."""
	for line in code:
		print('# '+line)
		exec(line,globals())
	print('# return ans')
	return ans
		
def geuvadis1():
	"""Perform the correlation test, without genotype data, from 10 miRNAs to 1000 genes"""
	lines=['import findr,findr.examples',
		'l=findr.lib() # or verbose: l=findr.lib(loglv=12)',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_rank(d["dmi"],d["dt"])']
	return runcode(lines)
	
def geuvadis2():
	"""Perform the novel causal inference test from 10 miRNAs to 1000 genes, using cis-eQTLs as causal anchors"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_gassist(d["dgmi"],d["dmi"],d["dt"])']
	return runcode(lines)

def geuvadis3():
	"""Perform the novel causal inference test from 1000 genes with cis-eQTLs to all 3000 genes"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_gassist(d["dgt"],d["dt"],d["dt2"],nodiag=True)']
	return runcode(lines)

def geuvadis4():
	"""Perform 5 subtests for causal inference test from 10 miRNAs to 1000 genes, using cis-eQTLs as causal anchors"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pijs_gassist(d["dgmi"],d["dmi"],d["dt"])']
	return runcode(lines)

def geuvadis5():
	"""Perform the novel causal inference test from 10 miRNAs to 1000 genes, using continuous causal anchors"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_cassist(d["dc"],d["dmi"],d["dt"])']
	return runcode(lines)

def geuvadis6():
	"""Perform 5 subtests for causal inference test from 10 miRNAs to 1000 genes, using cis-eQTLs as causal anchors and only obtaining p-values"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pijs_gassist_pv(d["dgmi"],d["dmi"],d["dt"])']
	return runcode(lines)

def geuvadis7():
	"""Constructs maximal direct acyclic graph from pairwise causal inference test probability among 1000 genes with cis-eQTLs"""
	lines=['import findr,findr.examples',
		'l=findr.lib()',
		'd=findr.examples.load_geuvadis_data()',
		'ans1=l.pij_gassist(d["dgt"],d["dt"],d["dt"],nodiag=True)',
		'ans=l.netr_one_greedy(ans1["p"])']
	return runcode(lines)





















