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

from .auto import *
import ctypes as c
ftype=getattr(c,ftype_c)
gtype=getattr(c,gtype_c)

ftype_p=c.POINTER(ftype)
gtype_p=c.POINTER(gtype)
c_ubyte_p=c.POINTER(c.c_ubyte)
