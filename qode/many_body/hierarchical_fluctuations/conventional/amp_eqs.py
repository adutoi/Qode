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


     omega1 = numpy.copy(fvo)  #1

     omega1 += numpy.tensordot( fvv, t1, axes = ([1],[0]))  #2

     temp1a = (-1)*numpy.tensordot( foo, t1, axes = ([0],[1]))  #3
     temp1b = numpy.moveaxis(temp1a, [0,1], [1,0])
     omega1 += temp1b

     omega1 += numpy.tensordot( Vovvo, t1, axes = ([0,2],[1,0]))  #4

     omega1 += numpy.tensordot( fov, t2, axes = ([0,1],[3,1]))  #5

     omega1 += (1/2.)*numpy.tensordot( Vovvv, t2, axes = ([0,2,3],[2,0,1]))  #6

     temp2a = (-1/2.)*numpy.tensordot( Voovo, t2, axes = ([0,1,2],[2,3,0]))  #7
     temp2b = numpy.moveaxis(temp2a, [0,1], [1,0])
     omega1 += temp2b

     temp3a = numpy.tensordot( fov, t1, axes = ([1],[0]))  #8
     temp3b = (-1)*numpy.tensordot( temp3a, t1, axes = ([0],[1]))
     temp3c = numpy.moveaxis(temp3b, [0,1], [1,0])
     omega1 += temp3c

     temp4a = numpy.tensordot( Voovo, t1, axes = ([0,2],[1,0]))  #9
     temp4b = (-1)*numpy.tensordot( temp4a, t1, axes = ([0],[1]))
     temp4c = numpy.moveaxis(temp4b, [0,1], [1,0])
     omega1 += temp4c

     temp5a = numpy.tensordot( Vovvv, t1, axes = ([0,2],[1,0]))  #10
     temp5b = numpy.tensordot( temp5a, t1, axes = ([1],[0]))
     omega1 += temp5b

     temp6a = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))  #11
     temp6b = numpy.tensordot( temp6a, t1, axes = ([1],[0]))
     temp6c = (-1)*numpy.tensordot( temp6b, t1, axes = ([0],[1]))
     temp6d = numpy.moveaxis(temp6c, [0,1], [1,0])
     omega1 += temp6d

     temp7a = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))  #12
     temp7b = numpy.tensordot( temp7a, t2, axes = ([0,1],[2,0]))
     omega1 += temp7b

     temp8a = numpy.tensordot( Voovv, t2, axes = ([0,2,3],[2,0,1]))  #13
     temp8b = (-1/2.)*numpy.tensordot( temp8a, t1, axes = ([0],[1]))
     temp8c = numpy.moveaxis(temp8b, [0,1], [1,0])
     omega1 += temp8c

     temp9a = numpy.tensordot( Voovv, t2, axes = ([0,1,2],[2,3,0]))  #14
     temp9b = (-1/2.)*numpy.tensordot( temp9a, t1, axes = ([0],[0]))
     omega1 += temp9b

     return omega1



def amp_eqs2(foo, fov, fvo, fvv,
             Voooo,Vooov,Voovo,Vovoo,Vvooo,Voovv,Vovov,Vvoov,Vovvo,Vvovo,Vvvoo,Vovvv,Vvovv,Vvvov,Vvvvo, Vvvvv,
             t1, t2):


    omega2 = numpy.zeros(t2.shape)
    omega2ij = numpy.zeros(t2.shape)
    omega2ab = numpy.zeros(t2.shape)
    omega2abij = numpy.zeros(t2.shape)



    omega2 = numpy.copy(Vvvoo)  #1

    temp2a = numpy.tensordot( fvv, t2, axes = ([1],[1])) #2.1
    temp2b = numpy.moveaxis(temp2a, [0,1,2,3], [1,0,2,3])
    omega2 += temp2b
    omega2 += (-1)*numpy.tensordot( fvv, t2, axes = ([1],[1])) #2.2


    temp3a = (-1)*numpy.tensordot( foo, t2, axes = ([0],[3])) #3.1
    temp3b = numpy.moveaxis(temp3a, [0,1,2,3], [3,0,1,2])
    omega2 += temp3b
    temp3c = numpy.tensordot( foo, t2, axes = ([0],[3])) #3.2
    temp3d = numpy.moveaxis(temp3c, [0,1,2,3], [2,0,1,3])
    omega2 += temp3d

    temp4a = (1/2.)*numpy.tensordot( Voooo, t2, axes = ([0,1],[2,3])) #4
    temp4b = numpy.moveaxis(temp4a, [0,1,2,3], [2,3,0,1])
    omega2 += temp4b

    omega2 += (1/2.)*numpy.tensordot( Vvvvv, t2, axes = ([2,3],[0,1])) #5

    temp10a = numpy.tensordot( Voovv, t2, axes = ([2,3],[0,1]))
    temp10b = (1/4.)*numpy.tensordot( temp10a, t2, axes = ([0,1],[2,3])) #10
    temp10c = numpy.moveaxis(temp10b, [0,1,2,3], [2,3,0,1])
    omega2 += temp10c
     
#only ij swap
    temp7a = numpy.tensordot( Vvvvo, t1, axes = ([2],[0])) #7
    temp7b = numpy.moveaxis(temp7a, [0,1,2,3], [0,1,3,2])
    omega2ij = numpy.copy(temp7b)   ###  a error here!

    temp12a = numpy.tensordot( Voovv, t2, axes = ([0],[3]))
    temp12b = (-1/2.)*numpy.tensordot( temp12a, t2, axes = ([0,1,2],[3,0,1])) #12
    omega2ij += temp12b

    temp14a = numpy.tensordot( Vvvvv, t1, axes = ([2],[0]))
    temp14b = (1/2.)*numpy.tensordot( temp14a, t1, axes = ([2],[0])) #14
    omega2ij += temp14b

    temp17a = numpy.tensordot( fov, t1, axes = ([1],[0]))
    temp17b = numpy.tensordot( temp17a, t2, axes = ([0],[3])) #17
    temp17c = numpy.moveaxis(temp17b, [0,1,2,3], [2,0,1,3])
    omega2ij += temp17c

    temp18a = numpy.tensordot( Voovo, t1, axes = ([0,2],[1,0]))
    temp18b = (-1)*numpy.tensordot( temp18a, t2, axes = ([0],[2])) #18
    temp18c = numpy.moveaxis(temp18b, [0,1,2,3], [2,0,1,3])
    omega2ij += temp18c

    temp22a = numpy.tensordot( Voovo, t1, axes = ([2],[0]))
    temp22b = (1/2.)*numpy.tensordot( temp22a, t2, axes = ([0,1],[2,3])) #22
    temp22c = numpy.moveaxis(temp22b, [0,1,2,3], [3,2,0,1])
    omega2ij += temp22c

    temp26a = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp26b = numpy.tensordot( temp26a, t1, axes = ([1],[0]))
    temp26c = (-1)*numpy.tensordot( temp26b, t2, axes = ([0],[2])) #26
    temp26d = numpy.moveaxis(temp26c, [0,1,2,3], [2,0,1,3])
    omega2ij += temp26d

    temp28a = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp28b = numpy.tensordot( temp28a, t1, axes = ([2],[0]))
    temp28c = (1/4.)*numpy.tensordot( temp28b, t2, axes = ([0,1],[2,3])) #28
    temp28d = numpy.moveaxis(temp28c, [0,1,2,3], [2,3,0,1])
    omega2ij += temp28d


    omega2ij1 = numpy.copy(omega2ij)
    omega2ij1 -= numpy.moveaxis(omega2ij, [0,1,2,3], [0,1,3,2])



#only ab swap

    temp8a = (-1)*numpy.tensordot( Vovoo, t1, axes = ([0],[1])) #8
    temp8b = numpy.moveaxis(temp8a, [0,1,2,3], [1,2,3,0])
    omega2ab = numpy.copy(temp8b)

    temp11a = numpy.tensordot( Voovv, t2, axes = ([2],[1]))
    temp11b = (-1/2.)*numpy.tensordot( temp11a, t2, axes = ([0,1,2],[2,3,1])) #11
    temp11c = numpy.moveaxis(temp11b, [0,1,2,3], [0,2,3,1])
    omega2ab += temp11c

    temp13a = numpy.tensordot( Voooo, t1, axes = ([0],[1]))
    temp13b = (1/2.)*numpy.tensordot( temp13a, t1, axes = ([0],[1])) #13
    temp13c = numpy.moveaxis(temp13b, [0,1,2,3], [2,3,0,1])
    omega2ab += temp13c

    temp16a = numpy.tensordot( fov, t1, axes = ([0],[1]))
    temp16b = numpy.tensordot( temp16a, t2, axes = ([0],[1])) #16
    omega2ab += temp16b

    temp19a = numpy.tensordot( Vovvv, t1, axes = ([0,2],[1,0]))
    temp19b = numpy.tensordot( temp19a, t2, axes = ([1],[0])) #19
    omega2ab += temp19b

    temp23a = numpy.tensordot( Vovvv, t1, axes = ([0],[1]))
    temp23b = (-1/2.)*numpy.tensordot( temp23a, t2, axes = ([1,2],[0,1])) #23
    temp23c = numpy.moveaxis(temp23b, [0,1,2,3], [1,0,2,3])
    omega2ab += temp23c

    temp27a = numpy.tensordot( Voovv, t1, axes = ([0,2],[1,0]))
    temp27b = numpy.tensordot( temp27a, t1, axes = ([0],[1]))
    temp27c = (-1)*numpy.tensordot( temp27b, t2, axes = ([0],[0])) #27
    omega2ab += temp27c

    temp29a = numpy.tensordot( Voovv, t1, axes = ([0],[1]))
    temp29b = numpy.tensordot( temp29a, t1, axes = ([0],[1]))
    temp29c = (1/4.)*numpy.tensordot( temp29b, t2, axes = ([0,1],[0,1])) #29
    omega2ab += temp29c

    omega2ab1 = numpy.copy(omega2ab)
    omega2ab1 -= numpy.moveaxis(omega2ab, [0,1,2,3], [1,0,2,3])

#both ab and ij

    temp6a = numpy.tensordot( Vovvo, t2, axes = ([0,2],[3,1])) #6
    temp6b = numpy.moveaxis(temp6a, [0,1,2,3], [1,3,0,2])
    omega2abij =numpy.copy(temp6b)


    temp9a = numpy.tensordot( Voovv, t2, axes = ([0,2],[3,1]))
    temp9b = (1/2.)*numpy.tensordot( temp9a, t2, axes = ([0,1],[2,0])) #9
    temp9c = numpy.moveaxis(temp9b, [0,1,2,3], [0,2,1,3])
    omega2abij += temp9c

    temp15a = numpy.tensordot( Vovov, t1, axes = ([0],[1]))
    temp15b = (-1)*numpy.tensordot( temp15a, t1, axes = ([2],[0])) #15
    temp15c = numpy.moveaxis(temp15b, [0,1,2,3], [1,2,0,3])
    omega2abij += temp15c

    temp20a = numpy.tensordot( Vvovv, t1, axes = ([2],[0]))
    temp20b = numpy.tensordot( temp20a, t2, axes = ([1,2],[3,1])) #20
    temp20c = numpy.moveaxis(temp20b, [0,1,2,3], [0,2,1,3])
    omega2abij += temp20c

    temp21a = numpy.tensordot( Vooov, t1, axes = ([1],[1]))
    temp21b = numpy.tensordot( temp21a, t2, axes = ([0,2],[3,1])) #21
    temp21c = numpy.moveaxis(temp21b, [0,1,2,3], [2,0,1,3])
    omega2abij += temp21c

    temp24a = numpy.tensordot( Vovvv, t1, axes = ([2],[0]))
    temp24b = numpy.tensordot( temp24a, t1, axes = ([0],[1]))
    temp24c = (-1/2.)*numpy.tensordot( temp24b, t1, axes = ([1],[0])) #24
    temp24d = numpy.moveaxis(temp24c, [0,1,2,3], [1,2,0,3])
    omega2abij += temp24d

    temp25a = numpy.tensordot( Voovo, t1, axes = ([2],[0]))
    temp25b = numpy.tensordot( temp25a, t1, axes = ([0],[1]))
    temp25c = (1/2.)*numpy.tensordot( temp25b, t1, axes = ([0],[1])) #25
    temp25d = numpy.moveaxis(temp25c, [0,1,2,3], [3,2,0,1])
    omega2abij += temp25d

    temp30a = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp30b = numpy.tensordot( temp30a, t1, axes = ([1],[1]))
    temp30c = numpy.tensordot( temp30b, t2, axes = ([0,1],[2,1])) #30
    temp30d = numpy.moveaxis(temp30c, [0,1,2,3], [2,1,0,3])
    omega2abij += temp30d

    temp31a = numpy.tensordot( Voovv, t1, axes = ([2],[0]))
    temp31b = numpy.tensordot( temp31a, t1, axes = ([0],[1]))
    temp31c = numpy.tensordot( temp31b, t1, axes = ([1],[0]))
    temp31d = (1/4.)*numpy.tensordot( temp31c, t1, axes = ([0],[1])) #31
    temp31e = numpy.moveaxis(temp31d, [0,1,2,3], [2,0,3,1])
    omega2abij += temp31e

    omega2abij1 = numpy.copy(omega2abij)
    omega2abij2 = numpy.copy(omega2abij)
    omega2abij3 = numpy.copy(omega2abij)
    omega2abij4 = numpy.copy(omega2abij)
    omega2abij1 -= numpy.moveaxis(omega2abij2, [0,1,2,3], [0,1,3,2])
    omega2abij1 -= numpy.moveaxis(omega2abij3, [0,1,2,3], [1,0,2,3])
    omega2abij1 += numpy.moveaxis(omega2abij4, [0,1,2,3], [1,0,3,2])


    omega2 += omega2ij1 + omega2ab1 + omega2abij1


    return omega2
