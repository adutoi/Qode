#    (C) Copyright 2019 Simeng Zhang
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
import numpy

def energy(fov, Voovv, t1, t2):
    energy = numpy.tensordot( fov, t1, axes = ([0,1],[1,0]))
    energy += (1/4.)* numpy.tensordot( Voovv, t2, axes = ([0,1,2,3],[2,3,0,1]))
    tempp = numpy.tensordot( Voovv, t1, axes = ([2,0],[0,1]))
    energy += (1/2.)*numpy.tensordot( tempp, t1, axes = ([0,1],[1,0]))
    return energy



def amp_eqs1(foo, fov, fvo, fvv,
             Voooo,Vooov,Voovo,Vovoo,Vvooo,Voovv,Vovov,Vvoov,Vovvo,Vvovo,Vvvoo,Vovvv,Vvovv,Vvvov,Vvvvo, Vvvvv,
             t1, t2):

    print("t1 = \n",t1)

    omega1 = numpy.copy(fvo)   #0

    omega1 += numpy.tensordot( fvv, t1, axes = ([1],[0]))  #1

    temp9 = (-1)*numpy.tensordot( foo, t1, axes = ([0],[1]))   #2
    temp9 = numpy.moveaxis(temp9, [0,1], [1,0])
    omega1 += temp9

    omega1 += numpy.tensordot( Vovvo, t1, axes = ([0,2],[1,0]))   #3

    omega1 += numpy.tensordot( fov, t2, axes = ([0,1],[3,1]))    #4

    omega1 += (1/2.)*numpy.tensordot( Vovvv, t2, axes = ([0,2,3],[2,0,1]))   #5

    temp10 = (-1/2.)*numpy.tensordot( Voovo, t2, axes = ([0,1,2],[2,3,0]))   #6
    temp10 = numpy.moveaxis(temp10, [0,1], [1,0])
    omega1 += temp10

    temp1 = numpy.tensordot( fov, t1, axes = ([1],[0]))
    temp11 = (-1)*numpy.tensordot( temp1, t1, axes = ([0],[1]))   #7
    temp11 = numpy.moveaxis(temp11, [0,1], [1,0])
    omega1 += temp11

    temp2 = numpy.tensordot( Voovo, t1, axes = ([0,2],[1,0]))
    temp12 = (-1)*numpy.tensordot( temp2, t1, axes = ([0],[1]))   #8
    temp12 = numpy.moveaxis(temp12, [0,1], [1,0])
    omega1 += temp12

    temp3a = numpy.tensordot( Vovvv, t1, axes = ([0,2],[1,0]))
    temp3b = numpy.tensordot( temp3a, t1, axes = ([1],[0]))    #9
    omega1 += temp3b      ##############   modified!20190804

    
    temp4 = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp5 = numpy.tensordot( temp4, t1, axes = ([1],[0]))
    temp13 = (-1)*numpy.tensordot( temp5, t1, axes = ([0],[1]))   #10
    temp13 = numpy.moveaxis(temp13, [0,1], [1,0])
    omega1 += temp13


    temp6a = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp6b = numpy.tensordot( temp6a, t2, axes = ([0,1],[2,0]))   #11
    omega1 += temp6b               #############     modified!20190804


    temp7 = numpy.tensordot( Voovv, t2, axes = ([0,2,3],[2,0,1]))
    temp14 = (-1/2.)*numpy.tensordot( temp7, t1, axes = ([0],[1]))   #12
    temp14 = numpy.moveaxis(temp14, [0,1], [1,0])
    omega1 += temp14

    temp8 = numpy.tensordot( Voovv, t2, axes = ([0,1,2],[2,3,0]))
    omega1 += (-1/2.)*numpy.tensordot( temp9, t1, axes = ([0],[0]))   #13

    #print("omega1 = \n", omega1)

    return omega1


def amp_eqs2(foo, fov, fvo, fvv,
             Voooo,Vooov,Voovo,Vovoo,Vvooo,Voovv,Vovov,Vvoov,Vovvo,Vvovo,Vvvoo,Vovvv,Vvovv,Vvvov,Vvvvo, Vvvvv,
             t1, t2):
    omega2 = numpy.zeros(t2.shape)
    omega2ij = numpy.zeros(t2.shape)
    omega2ab = numpy.zeros(t2.shape)
    omega2abij = numpy.zeros(t2.shape)

    print("test:t2 = \n", t2)


    test1 = numpy.copy(t2)
    test2 = numpy.copy(t2)
    test1 += numpy.moveaxis(test1, [0,1,2,3], [1,0,2,3])
    test2 += numpy.moveaxis(test2, [0,1,2,3], [0,1,3,2])

    print("test: tab = ", numpy.linalg.norm(test1))
    print("test: tij = ", numpy.linalg.norm(test2))


    omega2 = numpy.copy(Vvvoo)  #1
 
    temp36 = numpy.tensordot( fvv, t2, axes = ([1],[1])) #2.1
    omega2 += numpy.moveaxis(temp36, [0,1,2,3], [1,0,2,3])
    omega2 += (-1)*numpy.tensordot( fvv, t2, axes = ([1],[1])) #2.2


    temp37 = (-1)*numpy.tensordot( foo, t2, axes = ([0],[3])) #3.1
    omega2 += numpy.moveaxis(temp37, [0,1,2,3], [3,0,1,2])
    temp38 = numpy.tensordot( foo, t2, axes = ([0],[3])) #3.2
    omega2 += numpy.moveaxis(temp38, [0,1,2,3], [2,0,1,3])

    temp39 = (1/2.)*numpy.tensordot( Voooo, t2, axes = ([0,1],[2,3])) #4
    omega2 += numpy.moveaxis(temp39, [0,1,2,3], [2,3,0,1])

    omega2 += (1/2.)*numpy.tensordot( Vvvvv, t2, axes = ([2,3],[0,1])) #5

    temp2a = numpy.tensordot( Voovv, t2, axes = ([2,3],[0,1])) 
    temp2b = (1/4.)*numpy.tensordot( temp2a, t2, axes = ([0,1],[2,3])) #10
    omega2 += numpy.moveaxis(temp2b, [0,1,2,3], [2,3,0,1])

#only ij swap
    temp33a = numpy.tensordot( Vvvvo, t1, axes = ([2],[0])) #7
    temp33a = numpy.moveaxis(temp33a, [0,1,2,3], [0,1,3,2])
    omega2ij = numpy.copy(temp33a)     ######  modified! 20190804

    temp4a = numpy.tensordot( Voovv, t2, axes = ([0],[3]))
    temp4b = (-1/2.)*numpy.tensordot( temp4a, t2, axes = ([0,1,2],[3,0,1])) #12
    omega2ij += temp4b

    temp6a = numpy.tensordot( Vvvvv, t1, axes = ([2],[0]))
    temp6b = (1/2.)*numpy.tensordot( temp6a, t1, axes = ([2],[0])) #14
    omega2ij += temp6b

    temp9a = numpy.tensordot( fov, t1, axes = ([1],[0]))
    temp9b = numpy.tensordot( temp9a, t2, axes = ([0],[3])) #17
    temp9b = numpy.moveaxis(temp9b, [0,1,2,3], [2,0,1,3])
    omega2ij += temp9b

    temp10a = numpy.tensordot( Voovo, t1, axes = ([0,2],[1,0]))
    temp10b = (-1)*numpy.tensordot( temp10a, t2, axes = ([0],[2])) #18
    temp10b = numpy.moveaxis(temp10b, [0,1,2,3], [2,0,1,3])
    omega2ij += temp10b

    temp14a = numpy.tensordot( Voovo, t1, axes = ([2],[0]))
    temp14b = (1/2.)*numpy.tensordot( temp14a, t2, axes = ([0,1],[2,3])) #22
    temp14b = numpy.moveaxis(temp14b, [0,1,2,3], [3,2,0,1])
    omega2ij += temp14b

    temp20 = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp21a = numpy.tensordot( temp20, t1, axes = ([1],[0]))
    temp21b = (-1)*numpy.tensordot( temp21a, t2, axes = ([0],[2])) #26
    temp21b = numpy.moveaxis(temp21b, [0,1,2,3], [2,0,1,3])
    omega2ij += temp21b

    temp24 = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp25a = numpy.tensordot( temp24, t1, axes = ([2],[0]))
    temp25b = (1/4.)*numpy.tensordot( temp25a, t2, axes = ([0,1],[2,3])) #28
    temp25b = numpy.moveaxis(temp25b, [0,1,2,3], [2,3,0,1])
    omega2ij += temp25b


    omega2ij1 = numpy.copy(omega2ij)
    omega2ij1 -= numpy.moveaxis(omega2ij, [0,1,2,3], [0,1,3,2])


#only ab swap

    temp34a = (-1)*numpy.tensordot( Vovoo, t1, axes = ([0],[1])) #8
    temp34a = numpy.moveaxis(temp34a, [0,1,2,3], [1,2,3,0])
    omega2ab = numpy.copy(temp34a)              ###### modified!

    temp3a = numpy.tensordot( Voovv, t2, axes = ([2],[1]))
    temp3b = (-1/2.)*numpy.tensordot( temp3a, t2, axes = ([0,1,2],[2,3,1])) #11
    temp3b = numpy.moveaxis(temp3b, [0,1,2,3], [0,2,3,1])
    omega2ab += temp3b

    temp5a = numpy.tensordot( Voooo, t1, axes = ([0],[1]))
    temp5b = (1/2.)*numpy.tensordot( temp5a, t1, axes = ([0],[1])) #13
    temp5b = numpy.moveaxis(temp5b, [0,1,2,3], [2,3,0,1])
    omega2ab += temp5b

    temp8a = numpy.tensordot( fov, t1, axes = ([0],[1]))
    temp8b = numpy.tensordot( temp8a, t2, axes = ([0],[1])) #16
    omega2ab += temp8b

    temp11a = numpy.tensordot( Vovvv, t1, axes = ([0,2],[1,0]))
    temp11b = numpy.tensordot( temp11a, t2, axes = ([1],[0])) #19
    omega2ab += temp11b

    temp15a = numpy.tensordot( Vovvv, t1, axes = ([0],[1]))
    temp15b = (-1/2.)*numpy.tensordot( temp15a, t2, axes = ([1,2],[0,1])) #23
    temp15b = numpy.moveaxis(temp15b, [0,1,2,3], [1,0,2,3])
    omega2ab += temp15b

    temp22 = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp23a = numpy.tensordot( temp22, t1, axes = ([0],[1]))
    temp23b = (-1)*numpy.tensordot( temp23a, t2, axes = ([0],[0])) #27
    omega2ab += temp23b

    temp26 = numpy.tensordot( Voovv, t1, axes = ([0],[1]))
    temp27a = numpy.tensordot( temp26, t1, axes = ([0],[1]))
    temp27b = (1/4.)*numpy.tensordot( temp27a, t2, axes = ([0,1],[0,1])) #29
    omega2ab += temp27b

    omega2ab1 = numpy.copy(omega2ab)
    omega2ab1 -= numpy.moveaxis(omega2ab, [0,1,2,3], [1,0,2,3])

#both ab and ij

    temp35a = numpy.tensordot( Vovvo, t2, axes = ([0,2],[3,1])) #6
    temp35b = numpy.moveaxis(temp35a, [0,1,2,3], [1,3,0,2])
    omega2abij = numpy.copy(temp35b)     ######   modified!


    temp1a = numpy.tensordot( Voovv, t2, axes = ([0,2],[3,1]))
    temp1b = (1/2.)*numpy.tensordot( temp1a, t2, axes = ([0,1],[2,0])) #9
    temp1b = numpy.moveaxis(temp1b, [0,1,2,3], [0,2,1,3])
    omega2abij += temp1b

    temp7a = numpy.tensordot( Vovov, t1, axes = ([0],[1]))
    temp7b = (-1)*numpy.tensordot( temp7a, t1, axes = ([2],[0])) #15
    temp7b = numpy.moveaxis(temp7b, [0,1,2,3], [1,2,0,3])
    omega2abij += temp7b

    temp12a = numpy.tensordot( Vvovv, t1, axes = ([2],[0]))
    temp12b = numpy.tensordot( temp12a, t2, axes = ([1,2],[3,1])) #20
    temp12b = numpy.moveaxis(temp12b, [0,1,2,3], [0,2,1,3])
    omega2abij += temp12b

    temp13a = numpy.tensordot( Vooov, t1, axes = ([1],[1]))
    temp13b = numpy.tensordot( temp13a, t2, axes = ([0,2],[3,1])) #21
    temp13b = numpy.moveaxis(temp13b, [0,1,2,3], [2,0,1,3])
    omega2abij += temp13b

    temp16 = numpy.tensordot( Vovvv, t1, axes = ([2],[0]))
    temp17a = numpy.tensordot( temp16, t1, axes = ([0],[1]))
    temp17b = (-1/2.)*numpy.tensordot( temp17a, t1, axes = ([1],[0])) #24
    temp17b = numpy.moveaxis(temp17b, [0,1,2,3], [1,2,0,3])
    omega2abij += temp17b

    temp18 = numpy.tensordot( Voovo, t1, axes = ([2],[0]))
    temp19a = numpy.tensordot( temp18, t1, axes = ([0],[1]))
    temp19b = (1/2.)*numpy.tensordot( temp19a, t1, axes = ([0],[1])) #25
    temp19b = numpy.moveaxis(temp19b, [0,1,2,3], [3,2,0,1])
    omega2abij += temp19b

    temp28 = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp29a = numpy.tensordot( temp28, t1, axes = ([1],[1]))
    temp29b = numpy.tensordot( temp29a, t2, axes = ([0,1],[2,1])) #30
    temp29b = numpy.moveaxis(temp29b, [0,1,2,3], [2,1,0,3])
    omega2abij += temp29b

    temp30 = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp31 = numpy.tensordot( temp30, t1, axes = ([0],[1]))
    temp32a = numpy.tensordot( temp31, t1, axes = ([1],[0]))
    temp32b = (1/4.)*numpy.tensordot( temp32a, t1, axes = ([0],[1])) #31
    temp32b = numpy.moveaxis(temp32b, [0,1,2,3], [2,0,3,1])
    omega2abij += temp32b

    omega2abij1 = numpy.copy(omega2abij)
    omega2abij2 = numpy.copy(omega2abij)
    omega2abij3 = numpy.copy(omega2abij)
    omega2abij4 = numpy.copy(omega2abij)
    omega2abij1 -= numpy.moveaxis(omega2abij2, [0,1,2,3], [0,1,3,2])
    omega2abij1 -= numpy.moveaxis(omega2abij3, [0,1,2,3], [1,0,2,3])
    omega2abij1 += numpy.moveaxis(omega2abij4, [0,1,2,3], [1,0,3,2])

    omega2 += omega2ij1 + omega2ab1 + omega2abij1
    
    #print("omega2 = \n", omega2)

    test1 = numpy.copy(omega2)
    test2 = numpy.copy(omega2)
    test1 += numpy.moveaxis(test1, [0,1,2,3], [1,0,2,3])
    test2 += numpy.moveaxis(test2, [0,1,2,3], [0,1,3,2])


    print("test: Oab = ", numpy.linalg.norm(test1))
    print("test: Oij = ", numpy.linalg.norm(test2))


    return omega2
    
