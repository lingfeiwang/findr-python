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
"""Examples of findr. For geuvadis dataset, see findr.examples.geuvadis1 and findr.examples.geuvadis2"""

def load_geuvadis_data():
	"""Return: dict of 'dg','dt','dt2' as part of geuvadis data."""
	from os.path import dirname,join
	from .auto import gtype_np,ftype_np
	import numpy as np
	d=join(dirname(__file__),'data','geuvadis')
	fg=join(d,'dg.dat')
	ft=join(d,'dt.dat')
	ft2=join(d,'dt2.dat')
	dg=np.fromfile(fg,dtype=gtype_np).reshape(10,360)
	dt=np.fromfile(ft,dtype=ftype_np).reshape(10,360)
	dt2=np.fromfile(ft2,dtype=ftype_np).reshape(3000,360)
	return {'dg':dg,'dt':dt,'dt2':dt2}
	
def geuvadis1():
	"""Example with geuvadis data without genotype information. Every line of execution is printed."""
	lines=['import findr,findr.examples',
		'l=findr.lib(loglv=12)',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_rank(d["dt"],d["dt2"])']
	for line in lines:
		print '# '+line
		exec line
	print '# return ans'
	return ans
	
def geuvadis2():
	"""Example with geuvadis data with genotype information. Every line of execution is printed."""
	lines=['import findr,findr.examples',
		'l=findr.lib(loglv=12)',
		'd=findr.examples.load_geuvadis_data()',
		'ans=l.pij_gassist(d["dg"],d["dt"],d["dt2"])']
	for line in lines:
		print '# '+line
		exec line
	print '# return ans'
	return ans
	
















