#    (C) Copyright 2019 Anthony D. Dutoi
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



# List of all of the element symbols such that the index of any symbol is its atomic number.  Gh is for "ghost"
symbols = [\
"Gh", \
"H",  \
"He", \
"Li", \
"Be", \
"B",  \
"C",  \
"N",  \
"O",  \
"F",  \
"Ne", \
"Na", \
"Mg", \
"Al", \
"Si", \
"P",  \
"S",  \
"Cl", \
"Ar", \
"K",  \
"Ca", \
"Sc", \
"Ti", \
"V",  \
"Cr", \
"Mn", \
"Fe", \
"Co", \
"Ni", \
"Cu", \
"Zn", \
"Ga", \
"Ge", \
"As", \
"Se", \
"Br", \
"Kr", \
"Rb", \
"Sr", \
"Y",  \
"Zr", \
"Nb", \
"Mo", \
"Tc", \
"Ru", \
"Rh", \
"Pd", \
"Ag", \
"Cd", \
"In", \
"Sn", \
"Sb", \
"Te", \
"I",  \
"Xe", \
"Cs", \
"Ba", \
"La", \
"Ce", \
"Pr", \
"Nd", \
"Pm", \
"Sm", \
"Eu", \
"Gd", \
"Tb", \
"Dy", \
"Ho", \
"Er", \
"Tm", \
"Yb", \
"Lu", \
"Hf", \
"Ta", \
"W",  \
"Re", \
"Os", \
"Ir", \
"Pt", \
"Au", \
"Hg", \
"Tl", \
"Pb", \
"Bi", \
"Po", \
"At", \
"Rn", \
"Fr", \
"Ra", \
"Ac", \
"Th", \
"Pa", \
"U",  \
"Np", \
"Pu", \
"Am", \
"Cm", \
"Bk", \
"Cf", \
"Es", \
"Fm", \
"Md", \
"No", \
"Lr", \
"Rf", \
"Db", \
"Sg", \
"Bh", \
"Hs", \
"Mt", \
"Ds", \
"Rg", \
"Cn", \
"Nh", \
"Fl", \
"Mc", \
"Lv", \
"Ts", \
"Og"  ]



# Make this more general so that it can also accept a string that resolves to an int (or a float?) ... and/or accepts and int/float as an argument?
def atomic_number(symbol):
    """ returns the atomic number associated with a given element symbol """
    return symbols.index(symbol)
