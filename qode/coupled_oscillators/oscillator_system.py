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
import numpy as np
import copy
from qode.coupled_oscillators import build_coupling_mat

# The OscillatorSystem class supports:
#    1) Initializing ONE Harmonic Oscillator A TIME with a list of [ mass, k, and h bar ].
#    2) Combine SEVERAL oscillator SYSTEMS. ( The SYSTEM can be ONE oscillator or a list of HOs.)
#

class primitive_oscillator(object):
    def __init__( self,
                  m,
                  k
                ):
        self.m = m
        self.k = k


class group_of_oscillators(object):
    def __init__( self,
                  recursive_list,
  	              coupling
                 ):
        self.recursive_list = recursive_list
        self.coupling = coupling


class OscillatorSystem(object):
    def __init__( self,
                  argument=None 
                 ):
        """ parameters for single oscillator:
                argument = [ mass, k, h_bar ]
            parameters for HO groups:
                argument = [recursive_list, coupling ]
                             recursive_list = a list of OscillatorSystem objects ( single HO )
                                              OR SEVERAL complex oscillator SYSTEMS
                             coupling = a coupling matrix
        """
        # Lists of HOs case:
        if (  isinstance(argument, group_of_oscillators)  ):
            self.set_up_list( argument )
        # Single HO case:
        elif (  isinstance(argument, primitive_oscillator)  ):
            self.set_up_prim( argument )
        # NO OTHER CASES:
        else:
            raise ValueError('INVALID INPUT!')

    def set_up_list( self, argument ): 
        self.is_primitive = False
        self.group = copy.deepcopy( argument )

    def set_up_prim( self, argument ):
        self.is_primitive = True
        self.primitive = copy.deepcopy( argument )

    def list_of_force_const(self):
        force_constants = []
        if self.is_primitive:
            force_constants = [self.primitive.k]
        else:
            for i in self.group.recursive_list:
                force_constants.extend( i.list_of_force_const() )
        return force_constants

    def list_of_masses(self):
        masses = []
        if self.is_primitive:
            masses = [self.primitive.m]
        else:
            for i in self.group.recursive_list:
                masses.extend( i.list_of_masses() )
        return masses

    def top_level_coupling_mat(self):
        if self.is_primitive:
            coupling_mat = self.primitive.k
        else:
            coupling_mat = self.group.coupling 
        return coupling_mat

    def top_level_groups(self):
        return self.group.recursive_list
 

    def get_k_mat(self):
        if self.is_primitive:
            k_mat = [[ self.primitive.k ]]    # has to be made as a 2D matrix to fit in the correlate function 
        else:
            k_mat = build_coupling_mat.correlate(
 
                                       [ i.get_k_mat() for i in self.group.recursive_list ], 
                                       self.group.coupling 
 
                                                  ) 
        return k_mat


if __name__ == "__main__":
    # Test N group with ONLY ONE OSCILLATOR in each group
    ho1 = OscillatorSystem( primitive_oscillator( m = 1.0, k = 1.0 ) )
    ho2 = OscillatorSystem( primitive_oscillator( m = 1.0, k = 0.316 ) )
    
    ho3 = OscillatorSystem( primitive_oscillator( m = 1.0, k = 1.0 ) )
    ho4 = OscillatorSystem( primitive_oscillator( m = 1.0, k = 0.316 ) )
    
    mol1 = OscillatorSystem( group_of_oscillators( [ho1], [] ))
    mol2 = OscillatorSystem( group_of_oscillators( [ho2], [] ))
    mol3 = OscillatorSystem( group_of_oscillators( [ho3], [] ))
    mol4 = OscillatorSystem( group_of_oscillators( [ho4], [] ))
    coupling = [[0.0, -0.1, -0.1, -0.1],
                [0.0,  0.0, -0.1, -0.1],
                [0.0,  0.0,  0.0, -0.1],
                [0.0,  0.0,  0.0,  0.0]]

    total = OscillatorSystem( group_of_oscillators( [mol1,mol2,mol3,mol4], coupling ))
    print(total.list_of_force_const())
    print(total.list_of_masses())
    #printfunc.printline(total.get_k_mat(), '\t')
    
    print("Newest Test:")
    print(total.top_level_coupling_mat() )


'''
if __name__ == '__main__':
    import printfunc
    ho1 = OscillatorSystem( primitive_oscillator( m = 10.0, k = 1. ) )
    ho2 = OscillatorSystem( primitive_oscillator( m = 20.0, k = 2. ) )
    ho3 = OscillatorSystem( primitive_oscillator( m = 30.0, k = 3. ) )
    ho4 = OscillatorSystem( primitive_oscillator( m = 40.0, k = 4. ) )
    ho5 = OscillatorSystem( primitive_oscillator( m = 50.0, k = 5. ) )

    frag1 = OscillatorSystem( group_of_oscillators( recursive_list=[ho1,ho2], coupling=[[0.0, 0.12],
                                                                                        [0.12, 0.0]] ) )

    frag2 = OscillatorSystem( group_of_oscillators( recursive_list=[ho3,ho4,ho5], coupling=[[0.0, 0.34, 0.35],
                                                                                            [0.34, 0.0, 0.45],
                                                                                            [0.35, 0.45, 0.0]] ) )

    mol1 = OscillatorSystem( group_of_oscillators( recursive_list=[frag1,frag2], coupling=[[0.0, 12.345],
                                                                                           [12.345, 0.0]] ) )

    ho6 = OscillatorSystem( primitive_oscillator( m = 60.0, k = 6. ) )
    ho7 = OscillatorSystem( primitive_oscillator( m = 70.0, k = 7. ) )
    ho8 = OscillatorSystem( primitive_oscillator( m = 80.0, k = 8.) )
    ho9 = OscillatorSystem( primitive_oscillator( m = 90.0, k = 9.) )

    frag3 = OscillatorSystem( group_of_oscillators( recursive_list=[ho6,ho7], coupling=[[0.0, 0.67],
                                                                                        [0.67, 0.0]] ) )


    frag4 = OscillatorSystem( group_of_oscillators( recursive_list=[ho8,ho9], coupling=[[0.0, 0.89],
                                                                                        [0.89, 0.0]] ) )



    mol2 = OscillatorSystem( group_of_oscillators( recursive_list=[frag3,frag4], coupling=[[0.0, 67.89],
                                                                                           [67.89, 0.0]] ) )
    total = OscillatorSystem( group_of_oscillators( recursive_list=[mol1, mol2], coupling=[[0.0, 12345.6789],
                                                                                           [12345.6789, 0.0]] ) )
    print(total.list_of_force_const())
    print(total.list_of_masses())
    printfunc.printline(total.get_k_mat(), '\t')
    
    print("Newest Test:")
    print(total.top_level_coupling_mat() )'''
