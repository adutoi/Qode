#!/bin/csh
#    (C) Copyright 2018 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#

set list = ( commute_____Fo_Ex.py commute_____Fv_Ex.py commute_____Dx_Ex.py commute___ExFo_Ex.py commute___ExFv_Ex.py commute___ExDx_Ex.py commute___FoFo_Ex.py commute___FvFv_Ex.py commute___FoFv_Ex.py commute___FoDx_Ex.py commute___FvDx_Ex.py commute___DxDx_Ex.py commute_____Fo_ExEx.py commute_____Fv_ExEx.py commute_____Dx_ExEx.py commute___ExDx_ExEx.py commute___FoFo_ExEx.py commute___FvFv_ExEx.py commute___FoFv_ExEx.py commute___FoDx_ExEx.py commute___FvDx_ExEx.py commute___DxDx_ExEx.py commute_ExExDx_Ex.py commute_ExFoFo_Ex.py commute_ExFvFv_Ex.py commute_ExFoDx_Ex.py commute_ExFvDx_Ex.py commute_ExFoDx_ExEx.py commute_ExFvDx_ExEx.py )

foreach i ( $list )
 python3.3 $i >& $i.std
end
