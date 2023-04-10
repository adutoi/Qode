#    (C) Copyright 2016 Yuhong Liu
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
# This code reads Q-Chem Input file and retrieves all the information for any calculations.
import sys

'''
_list_ is a datatype defined as:
[ [x, y, z, exponential, contraction coefficient],[...],[...],[...],... ]

'''

def read_lines(name):
    'read all lines of input file'
    f = open(name,'r')
    lines = f.readlines()
    f.close()
    return lines






def read_char_mul(name):
    'Read charge, multiplicity'
    # Read by the lines of q_input, the line (i+1) of $molecule line
    # contains the CHARGE and the MULTIPLICITY

    q_input = read_lines(name)
    for i in range(len(q_input)):
        tmp = q_input[i].split()
        if len(tmp) == 1 and tmp[0] =='$molecule':
            tmp_next_line = q_input[i+1].split()
            charge = int(tmp_next_line[0])
            multi = int(tmp_next_line[1])
    
    
    return charge,multi
 




def jobtype(name):
    ' Return exchange, basis, and correlation'
    # By default correlation is None.
    correlation='None'
    
    q_input = read_lines(name)
    # Read exchange, correlation and basis
    for item in q_input:
        tmp = item.split()
        if len(tmp) == 2:
            if tmp[0]=='exchange':
                exch = tmp[1]
            elif tmp[0] == 'basis':
                basis = tmp[1]
            elif tmp[0] == 'correlation':
                correlation = tmp[1]
    
    return exch,basis,correlation




def get_gaussians(input_tmp,basis):   
    'Open corresponding basis file to retrieve primitive gaussians.'
    element = input_tmp[0]
    x = input_tmp[1]
    y = input_tmp[2]
    z = input_tmp[3]


    
    try:
        f = open('/home/yhliu/Qode-reorganized/qode/SCF_diag_full/'+basis+'.bas','r')
        # To find the atom in .bas file
        read = True
        while read:
            line =f.readline()
            basis_tmp = line.split()
            
            if len(line) == 0: # SO NO BLANK LINE ALLOWED IN BASIS FILES!!
                # Atom not Found, Pop out this message.
                print('Atom Not Supported in this Basis!')
                sys.exit(1)
                
            elif len(basis_tmp)==2:
                # Atom Found, Stop reading, start new section.
                if basis_tmp[1]=='0':
                    if basis_tmp[0]==element:
                        read = False

        # To read the Gaussian Functions in this section
        # If Found ****, exit.
        # If Found Primer: Element, # of Gaussians, 1.00. Keep reading Gaussians.
        
        return_list = []
        no_gaussians=[]
        orbitype=[]
        read = True
        while read:
            line =f.readline()
            read_basis_tmp = line.split()
            if read_basis_tmp[0]=='****':
                read = False
                # If **** is Found, it's the END of the Gaussians.
            elif len(read_basis_tmp)==3 and read_basis_tmp[0] =='S':
                if read_basis_tmp[2]=='1.00':
                    # If THREE components, also has 1.00 it's the Primer.
                    # Now Start Reading int(read_basis_tmp[1]) lines.
                    no_gaussians += [int(read_basis_tmp[1])]
                    orbitype += [read_basis_tmp[0]]
                    for i in range(int(read_basis_tmp[1])):
                        newline= f.readline()
                        read_basis_tmp=newline.split()
                        return_list += [[float(x),float(y),float(z),float(read_basis_tmp[0]),float(read_basis_tmp[1])]]

            elif len(read_basis_tmp)==3 and read_basis_tmp[0] =='P':
                if read_basis_tmp[2]=='1.00':
                    # If THREE components, also has 1.00 it's the Primer.
                    # Now Start Reading int(read_basis_tmp[1]) lines.
                    p_list = []
                    for i in range(3):
                        no_gaussians += [int(read_basis_tmp[1])]
                        orbitype += ['P']
                    for i in range(int(read_basis_tmp[1])):
                        newline= f.readline()
                        read_basis_tmp=newline.split()
                        p_list += [[float(x),float(y),float(z),float(read_basis_tmp[0]),float(read_basis_tmp[1])]]
                    return_list += p_list + p_list + p_list

            elif len(read_basis_tmp)==3 and read_basis_tmp[0] =='D':
                if read_basis_tmp[2]=='1.00':
                    # If THREE components, also has 1.00 it's the Primer.
                    # Now Start Reading int(read_basis_tmp[1]) lines.
                    d_list = []
                    for i in range(5):
                        no_gaussians += [int(read_basis_tmp[1])]
                        orbitype += ['D']
                    for i in range(int(read_basis_tmp[1])):
                        newline= f.readline()
                        read_basis_tmp=newline.split()
                        d_list += [[float(x),float(y),float(z),float(read_basis_tmp[0]),float(read_basis_tmp[1])]]
                    for i in range(5):
                        return_list += d_list

            elif len(read_basis_tmp)==3 and read_basis_tmp[0] =='SP':
                ################  CHANGE HERE  CHANGE HERE   ##########################################################
                #print("Not Supported")
                if read_basis_tmp[2]=='1.00':
                    # If THREE components, also has 1.00 it's the Primer.
                    # Now Start Reading int(read_basis_tmp[1]) lines.
                    no_gaussians += [int(read_basis_tmp[1])]  # One S plus one P
                    orbitype += ['S']
                    for i in range(3):
                        no_gaussians += [int(read_basis_tmp[1])]
                        orbitype += ['P']
                    s_list = []
                    p_list = []
                    for i in range(int(read_basis_tmp[1])):
                        newline= f.readline()
                        read_basis_tmp=newline.split()
                        s_list += [[float(x),float(y),float(z),float(read_basis_tmp[0]),float(read_basis_tmp[1])]]
                        p_list += [[float(x),float(y),float(z),float(read_basis_tmp[0]),float(read_basis_tmp[2])]]
                    return_list += s_list + p_list + p_list + p_list

                ################  CHANGE HERE  CHANGE HERE   ##########################################################
                
                    
                    
        f.close()
        return return_list,no_gaussians,orbitype




    # Basis not supported: go to IOError
    except IOError:
        print(basis,'basis is currently unavailable, sorry!')
        sys.exit(1)





   
def build_pri_basis_list(name,basis):   
    'read input file, generate lists, positions, and info'

    q_input = read_lines(name)
    input_list = []
    no_gauss = []
    orbitals = []
    for item in q_input:
        tmp = item.split()
        if len(tmp) == 4:   # Only lines with len==4 have the coordinates
            list_element, no_gaussians,orbitype = get_gaussians(tmp,basis)
            input_list += list_element
            no_gauss += no_gaussians
            orbitals += orbitype
    
    #for i in _list_:
    #    print(i)
    #print(no_gauss)
    #print(orbitals)
    return input_list,no_gauss,orbitals



    

def make_atom_list(name):
    "Read Q-Chem input, return all atom symbols"
    "name: input file name"
    "atom: all elemental symbols"
    qinput = read_lines(name)
    atom = []
    for line in qinput:
        tmp = line.split()
        if len(tmp)==4:
            atom += [tmp[0]]
    return atom




def chk_atomic_num(atom,atomtable):
    "Loop over all symbols to seek the atomic number for 'atom'"
    "atom: One Symbol"
    "atomtable: H\t1\nHe\t2\nLi\t3\n..."
    for i in atomtable:
        tmp = i.split()
        if tmp[0] == atom:
            anum = tmp[1]
            break
    return anum





def make_atomic_number_list(inputatoms):
    "Read all symbols, call chk_atomic_num, return number list"
    "inputatom:['H','H']"
    "num_electron.in: H\t1\nHe\t2\nLi\t3\n..."
    "anum:'1','2','3',..."
    "atomtable: num_electron.in -> list object"
    atomtable = read_lines('/home/yhliu/Qode-clone/qode/SCF_diag_full/num_electron.in')
    num_e = []
    for item in inputatoms:
        anum = chk_atomic_num(item,atomtable)
        num_e += [int(anum)]
    return num_e

def cal_num_electron(atomiclist,charge):
    "Input total charge, number of electrons in each atom"
    "Output total electron number as summation"
    summation = 0
    for i in atomiclist:
        summation += i
    summation -= charge
    return summation



#================== MAIN FUNCTION ========================================   

def mol_main(name):

    charge,mul = read_char_mul(name)
    exch,basis,correlation = jobtype(name)
    # input_list,no_gauss,orbitype = build_pri_basis_list(name,basis)
    input_list,no_gauss,orbitype = [],[],[]
    atoms = make_atom_list(name)
    A_num_list = make_atomic_number_list(atoms)
    num_electron = cal_num_electron(A_num_list,charge)

    
    print('Charge:',charge,'\t','Multiplicity:',mul)
    print('Exchange:',exch,'\nBasis:',basis,'\nCorrelation:',correlation)

    for i in input_list:
        print(i)

    print('# of gaussians in each orbital:')
    print(no_gauss)
    print('Orbital Types:')
    print(orbitype)
    print('Input Symbols:')
    print(atoms)
    print('Number of electrons:')
    print(num_electron)
    
    return input_list,no_gauss,orbitype,atoms,num_electron,charge,mul,A_num_list,basis,exch,correlation


