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
import  math
import cmath
from . import field_traits

# Biggest hole at present, if two functionable operators are added, should they not be functionable if in the same basis ... how to let the user implement
# that?  Would need a layer that does not immediately create a composite_op, but falls back to it if something more detailed not implemented.

# Also what about making vectors that are tensor products of vectors in different spaces ... how does that interact with dyadic notation

# How to make v+=A|w is all one operation? ... of course, the functions already defined below retain their behavior.
# Had thought about a thin wrapper vector class that stores a multiplier separately, for things like v+=c*w, but the problem is if something then changes the 
# underlying data in w ... about the only think I can think of is to make yet a deeper wrapper to the actual data that
# stores a list of vectors which have a lock on that data.  If that data should need to change in a way that is more fundamental than a scaling,
# then it forces a copy, then it forces that guy to copy the data ... still how to know when c*v falls out of scope (esp. if temporary).

def is_int_zero(vec):
    # This replaces the expression 'vec is 0' in if statements, which, as of python 3.8 began resulting in the following warning
    # SyntaxWarning: "is" with a literal. Did you mean "=="?
    # This is explained in the following article
    # https://adamj.eu/tech/2020/01/21/why-does-python-3-8-syntaxwarning-for-is-literal/
    # Summarized, some implementations of python do not guarantee that multiple instances of the integer 0 are all the same object.
    # This seems like an oversight in the language spec to me, but maybe its just really hard to enforce at the implementation level?
    # In any case, this function has the intented effect, where the explicit integer '0' is used as a signal/convenience in my code.
    return (isinstance(vec,int) and vec==0)

class _member(object):
    """\
    + An object of this class is returned by the function linear_inner_product_space.member.  It is thought of as a vector, but,
        strictly speaking could be some representation of a function in a function space.
    + The utility of this class is to define all the slick shorthands like "v = a*(a|w) + b*(b|w)" that automatically call
        user-implemented functions to resolve inner products, vector addition, etc.
    + The focus is on implementational simplicity, requiring the minimum amount of work from the user.
      Objects of this class are informed of the user implemented functions by storing a reference to the space (whose identity
        also serves as a check for consistent usage), which contains a traits object that stores references to the functions.
      This class ultimately calls the user-defined functions:
        - dot(v,w):   . . . . . . .  Returns the scalar product of v with w
	- add_to(v,w,c):  . . . . .  In-place addition of c*w to v (c is optional! default is 1), leaving w unchanged, No return value.
	- scale(n,v):   . . . . . .  In-place scaling of v by n.  No return value.
	- copy(v):    . . . . . . .  Instantiates and returns an independent (deep) copy of vector v.
        - back_act_on_vec(v,A):   .  Returns a *new vector* of the same type as v, which is the result of [v^+ A]^+ (same as A^+ v),
                                     where "^+" is interpreted as Hermitian conjugate and the "product" is interpreted as operator action.
                                     In this way, there is no distinction between left- and right-hand vectors, but consistency is
                                     maintained to the extent that "u = w|A ; n = (u|v)" gives the same value for n as "n = (w|A|v)"
                                     because conjugation of the left-hand argument of an inner product is part of the definition.
        - act_on_vec is not called by this class, since it is called by the operator classes directly, due to precedence.
      In case some of the assumptions below break the efficiency of the wrapped classes (seems unlikely), this interface can be
       replicated and other classes could be used in code that expects members of a space, as a last resort, I guess. (Once upon
       a time, "copy" took a scale argument which would scale the copy, rather than doing this in two steps here.  Is it possible that
       matters?)
    + For the most part, this class does not need to know anything about the insides of the user-implemented function, except that
        much higher level math does depend on the scalar field over which a space is defined.  The scalar field is the type returned
        by a general inner product of two arbitrary vectors.  Objects of this type are informed of the field type by the space member.
    + A non-obvious special feature of this class is that numbers may be interpreted as operators or scalars when their product
        is taken with a vector.  This means that (letting n be a number, and v a vector) "n*v" is the same as "v*n" and "n|v", but
        "v|n" might be different because n will be conjugated, according to the field type given, to maintain consistency (see above).
      This "scalar" promotion is useful because some operators might be defined as "1+A", where A is of operator class, and letting "1"
        be quietly promoted at resolution time allow us to skip checking for special cases elsewhere, which is error prone.
      There is a special case though; the integer "0" (always checked for by identity not value) may be interpreted as a vector,
        in addition to a scalar or operator.  This is safe, since the null vector must exist in all spaces (but the promotion of other
        scalars would be fundamentally ambiguous).  Any action that uses an identical zero and results in the null vector is "demoted"
        to the scalar "0" in its representation; this is not only efficient, but it is also valuable to know when a vector is identically zero.
        It also is necessary to maintain consistency since "0|v" is ambiguous in its interpretation;  if "0" is an operator (why would one do
        this anyway?), the result is a vector, but this expression might also be an inner product, but since the null vector is demoted,
        the result is the same, and the user will likely not even notice the ambiguity.
    """
    def __init__(self,v,space):
        self.v     = v
        self.space = space
        self.field = self.space.field
    def __or__(self,other):
        """ Handles (v|w) and v|A, in the case where w or A is either another vector an an operator, respectively,
            or a number interpreted as an operator, unless it is 0, when it could be interpreted as either a vector or operator. """
        if   isinstance(other,(float,int,complex)):  return self*field_traits.conjugate(other)		# calls __mul__ below (handles "0" case)
        elif isinstance(other,_operator_base):       return self.space.traits.back_act_on_vec(self,other)
        else:                                        return self.space.traits.dot(self,other)		# checks that both are _member class
    def __ror__(self,other):
        """ Handles n|v, where n is a number interpreted as an operator, unless it is 0, when it could be interpreted as either a vector or operator.
            If the left-hand operand would be of vector or operator class, then __or__ for that class would take precedence. """
        if isinstance(other,(float,int,complex)):  return self*other					# calls __mul__ below (handles "0" case)
        else:                                      raise Exception('Illegal inner product with unknown object')
    def __matmul__(self,other):
        return self|other
    def __rmatmul__(self,other):
        return other|self
    def __iadd__(self,other):
        """ Handles v+=w, where w is another vector (perhaps 0). """
        if not is_int_zero(other):  self.space.traits.add_to(self,other)		# checks that both are _member class
        return self
    def __add__(self,other):
        """ Handles v+w, where w is another vector (perhaps 0). """
        value = self.space.traits.copy(self)
        value += other		# calls __iadd__ above (handles "0" case)
        return value
    def __radd__(self,other):
        """ Handles 0+v.  If the left-hand vector would be of vector class, then __add__ above would take precedence. """
        if is_int_zero(other):  return self.space.traits.copy(self)
        else:                   raise Exception('Illegal linear combination with unknown object')
    def __isub__(self,other):
        """ Handles v-=w, where w is another vector (perhaps 0). """
        if not is_int_zero(other):  self.space.traits.add_to(self,other,-1)		# checks that both are _member class
        return self
    def __sub__(self,other):
        """ Handles v-w, where w is another vector (perhaps 0). """
        value = self.space.traits.copy(self)
        value -= other		# calls __isub__ above (handles "0" case)
        return value
    def __rsub__(self,other):
        """ Handles 0-v.  If the left-hand vector would be of vector class, then __sub__ above would take precedence. """
        if is_int_zero(other):  return -self		# calls __neg__ below (handles "0" case)
        else:                   raise Exception('Illegal linear combination with unknown object')
    def __imul__(self,n):
        """ Handles v*=n, where v is a vector and n is a scalar. """
        self.space.traits.scale(n,self)
        return self
    def __mul__(self,other):
        """ Handles v*n, where v is a vector and n is a scalar and v*w, where both v and w are vectors. """
        if is_int_zero(other):
            return 0
        elif isinstance(other,(float,int,complex)):
            value = self.space.traits.copy(self)
            self.space.traits.scale(other,value)
            return value
        else:
            return _dyadic_op(self,other)			# checks that both are _member class
    def __rmul__(self,n):
        """ Handles n*v, where v is a vector and n is a scalar. """
        if isinstance(n,(float,int,complex)):  return self*n		# calls __mul__ above
        else:   raise Exception('Illegal scalar product with unsupported object')
    def __truediv__(self,n):
        """ Handles v/n, where v is a vector and n is a scalar. """
        return self * (1./n)		# calls __mul__ above (implicit check that n is not zero)
    def __itruediv__(self,n):
        """ Handles v/=n, where v is a vector and n is a scalar. """
        self *= (1./n)			# calls __imul__ above (implicit check that n is not zero)
        return self
    def __neg__(self):
        """ Handles -v. """
        value = self.space.traits.copy(self)
        self.space.traits.scale(-1,value)
        return value



class _operator_base(object):
    """\
    + This is the base class for all objects returned by the function linear_inner_product_space.lin_op.  It is thought of as 
        representing a matrix, but it could be any linear mapping (also on non-discrete spaces).
    + The utility of this class is to define all the slick shorthands like "v = A|w" or "A = T + c*V" that automatically call
        user-implemented functions to resolve operator action or build internal representations of composite operators.
      Specifically this base class implements the two things we insist that all operators (composite, primitive or specialized) can do:
        - be added or contracted with other operators (i.e., form composites), which might be numbers (see comments for _member above).
        - act on a vector (backwards action is handled above in the _member class, due to precedence)
    + The focus is on implementational simplicity, requiring the minimum amount of work from the user.
      Objects of this class are informed of the user implemented functions by storing a reference to the space (whose identity
        also serves as a check for consistent usage), which contains a traits object that stores references to the functions.
      This class ultimately calls the user-defined functions:
        - act_on_vec(A,v)   . . . . Returns a *new vector* of the same type as v, which is the result of Av, where the "product" is
                                    interpreted as operator action.
        - The other functions a user will implement (dot, add_to, scale, copy, back_act_on_vec) are not called explicitly by objects of
            this class.
    + For the most part, this class does not need to know anything about the insides of the user-implemented function, except that
        much higher level math does depend on the scalar field over which a space is defined.  The scalar field is the type returned
        by a general inner product of two arbitrary vectors.  Objects of this type are informed of the field type by the space member.

    !!!
    + There is a big difference between the math of vectors and operators!  Operator math just makes wrappers.  If you want true
      operator math, you should make a space of operators and define how to do it there.  You are then at liberty to define how members
      of that space furthermore operate on vectors of another space, but you will have to realize that you will always have two versions
      of your operators, those that perform real math (as though they were vectors), and those that perform encapsulated abstracted math
      ("promoted" to operators).
    !!!

    """
    # Right now, multiplication * is interpreted identically to contraction |, which is ok for scalars, which is the intention, but
    # there will be errors if used between non-commuting operators, so should check that this does not happen.
    def __init__(self,space):
        self.space = space
        self.field = self.space.field
    def __or__(self,other):
        if   is_int_zero(other):                                    return 0
        elif isinstance(other,(_operator_base,float,int,complex)):  return _composite_op(self,other,'|')
        else:                                                       return self.space.traits.act_on_vec(self,other)	# checks that other is _member class
    def __ror__(self,number):       return   self | number
    def __matmul__(self,other):     return   self | other
    def __rmatmul__(self,other):    return   other| self
    def __mul__(self,number):       return   self | number
    def __rmul__(self,number):      return   self | number
    def __truediv__(self,number):   return   self | (self.field(1)/number)
    def __neg__(self):              return   self | (-1)
    def __add__(self,other):        return _composite_op(self,other,'+')
    def __radd__(self,number):      return   self  +  number
    def __sub__(self,other):        return   self  + (-other)
    def __rsub__(self,number):      return (-self) +  number	# Always expect first component of composite to be an operator
    # The following all rely on the fact that A+=B is exactly equivalent to A=A+B, etc.  That is, the return value, which is a composite
    # is simply assigned to the original label.  This is slightly different than the rule for vectors, in which the vector is modified
    # in place and returned.  You would never notice the difference unless you had another label that is attached to the object that is
    # either modified or not ...
    def __iadd__(self,other):       return   self + other
    def __isub__(self,other):       return   self - other
    def __imul__(self,number):      return   self * number
    def __itruediv__(self,number):  return   self / number
    def act_on_vec_block(self,v_block):       return self.space.traits.act_on_vec_block(self,v_block)
    def back_act_on_vec_block(self,v_block):  return self.space.traits.back_act_on_vec_block(v_block,self)

# Some checks advertised above are skipped below!

class _composite_op(_operator_base):
    """\
    + Class for holding added or contracted (multiplied) operators X (e.g., X=A|B or X=A+B), so that they have X|v and v|X on a vector v defined.
    """
    def __init__(self,A,B,algebra):
        _operator_base.__init__(self,A.space)
        if not isinstance(B,(float,int,complex)):
            if A.space is not B.space:  raise Exception('Illegal sequence of operators on different spaces, or unknown objects acting on operators')
        self.A       = A
        self.B       = B
        self.algebra = algebra

class _dyadic_op(_operator_base):
    """\
    + Class for holding and operator that is the outer product of two vectors v*w
    """
    def __init__(self,v,w):
        _operator_base.__init__(self,v.space)
        if v.space is not w.space:  raise Exception('Illegal outer product of vectors in different spaces')
        self.v = v
        self.w = w

class _primitive_op(_operator_base):
    """\
    + Holds the user-defined operator representation and reference (through) traits to how the dirty work is actually done
    """
    def __init__(self,op,space):
        _operator_base.__init__(self,space)
        self.op = op
    def diagonal(self):  return self.space.traits.diagonal(self)

class _fn_prim_op(_primitive_op):
    """\
    + Allows one to write things like 2*A and exp(A), assuming the implementor can implement member functions
        - copy
        - function_on_diags
    + These sort of assume the user is storing the operator in some diagonalized form, but it is not strictly necessary
    """
    def __init__(self,op,space):
        _primitive_op.__init__(self,op,space)
    def generic_fn(self,func):
        return self.space.traits.function_on_diags(func,self)
    def __rmul__(self,number):
        func = lambda x: number*x
        return self.generic_fn(func)
    def __mul__(self,number):
        func = lambda x: number*x
        return self.generic_fn(func)
    def __truediv__(self,number):
        func = lambda x: x/number
        return self.generic_fn(func)
    def __neg__(self):
        func = lambda x: -x
        return self.generic_fn(func)
    def __pow__(self,number):
        func = lambda x: x**number
        return self.generic_fn(func)



def abs(obj):
    """ Overrides the abs() function to seamlessly handle all numbers and allowed (functionable) operators. """
    if isinstance(obj,_fn_prim_op):  return obj.generic_fn(field_traits.abs)
    else:                            return field_traits.abs(obj)

def conjugate(obj):
    """ Overrides the conjugate() function to seamlessly handle all numbers and allowed (functionable) operators. """
    if isinstance(obj,_fn_prim_op):  return obj.generic_fn(field_traits.conjugate)
    else:                            return field_traits.conjugate(obj)

def sqrt(obj):
    """ Overrides the sqrt() function to seamlessly handle all numbers and allowed (functionable) operators. """
    if isinstance(obj,_fn_prim_op):  return obj.generic_fn(field_traits.sqrt)
    else:                            return field_traits.sqrt(obj)

def exp(obj):
    """ Overrides the exp() function to seamlessly handle all numbers and allowed (functionable) operators. """
    if isinstance(obj,_fn_prim_op):  return obj.generic_fn(field_traits.exp)
    else:                            return field_traits.exp(obj)





class _space_traits(object):
    """\
    among other things . . .
    + Maps X|v and v|X of composite X in terms of + and | for the stored primitive operators (see operator class for details).
    + the traits class taken as a constructor argument is meant to have some things like the following defined, but these are out of date ...
        - ops_are_compatible(A,B):  . . .  Returns a boolean; True if A and B act on the same vector space (user creativity needed
                                             for implementation, but, if type saftey is not desired, you can just return True and be
                                             careful).  The implementor may assume that both A and be are non-composite operators;
                                             the recursion needed to deal with testing composites is taken care of in this file
        - op_and_vec_are_compatible(A,v):  Returns a boolean; True if A acts on the vector space to which v belongs (user creativity
                                             needed for implementation, but, if type saftey is not desired, you can just return True
                                             and be careful).
        - act_on_vec(A,v):  . . . . . . .  Returns a *new vector* of the same type as v, which is the action of A onto v from the left
                                             (compatibility is checked).
        - back_act_on_vec(v,A):   . . . .  Returns a *new vector* of the same type as v, which is the action of A onto v from the right
                                             (compatibility is checked).  Note that this function is not called from the operator
                                      class but rather it is called from the vector class when something like v|A is parsed.
    """
    def __init__(self,traits):
        self.traits = traits
        self.field  = self.traits.field
    def dot(self,a,b):
        if a.space is not b.space:  raise Exception('Illegal inner product of members of different spaces')
        value = self.traits.dot(a.v,b.v)
        if a is b:
            if abs(value)!=0 and abs(value.imag/value.real)>1e-14:  raise Exception('Something seriously wrong.  Norm^2 is not real.')      # Danger, hardcoded threshold!
            return value.real
        else:
            return value
    # Hmmm, a lot of these functions use raw constructors which avoid check functions ... move check function calls into constructors (presently inside generator functions of space class)
    def add_to(self,a,b,n=1):
        if a.space is not b.space:  raise Exception('Illegal linear combination of members from different spaces')
        self.traits.add_to(a.v,b.v,n)
    def scale(self,n,a):
        self.traits.scale(n,a.v)
    def copy(self,a):
        return _member(self.traits.copy(a.v),a.space)
    def act_on_vec(self,X,a):
        if X.space is not a.space:  raise Exception('Illegal operator on member of different space')
        if isinstance(X,_composite_op):
            A = X.A
            B = X.B
            if   X.algebra=='+':  return (A|a) + (B|a)
            elif X.algebra=='|':  return (A|(B|a))
        elif isinstance(X,_dyadic_op):
            return (X.v)*(X.w|a)
        else:  return _member(self.traits.act_on_vec(X.op,a.v),a.space)
    def back_act_on_vec(self,a,X):
        if X.space is not a.space:  raise Exception('Illegal operator on member of different space')
        if isinstance(X,_composite_op):
            A = X.A
            B = X.B
            if   X.algebra=='+':  return (a|A) + (a|B)
            elif X.algebra=='|':  return ((a|A)|B)
        elif isinstance(X,_dyadic_op):
            return (X.w)*(X.v|a)
        else:  return _member(self.traits.back_act_on_vec(a.v,X.op),a.space)
    def function_on_diags(self,func,X):
        return _fn_prim_op(self.traits.function_on_diags(func,X.op),X.space)
    def diagonal(self,X):
        return _fn_prim_op(self.traits.diagonal(X.op),X.space)
    def check_member(self,v):      self.traits.check_member(v)
    def check_lin_op(self,op):     return self.traits.check_lin_op(op)
    #
    # Part of an incomplete progression towards also allowing blocks of vectors.
    # This assumes that vector blocks come in and out as lists, which has the advantage of being flexible (for examply in use with Lanczos that makes new vectors).
    # However, if the efficiency of acting on a block is tied to having them be contiguous in memory, then the user would be stuck on the other side of this doing
    # one of two things:  (1) copying vectors into contiguous block, acting and then copying back to individually indexable vectors or (2) implementing some rather
    # sophisticated back-end memory management and then passing around indexable "vectors" that really just point to the physical information.  More about my philosophy
    # about how to handle "blocks" of vectors can be found in vector_set.py
    #
    # The biggest gaping hole right now is that there is no slick syntax for calling this.  The user has to call this function by name directly on a list of vectors (or a vector_set)
    #
    def act_on_vec_block(self,X,a_block):
        if isinstance(X,(float,int,complex)):  return [ X|a for a in a_block ]		# This only needs to be here since act_on_vec_block called directly since not finished yet
        # Note return above . . . the remainder does not appy to scalar operators
        the_space = X.space
        for a in a_block:
            if a.space is not the_space:  raise Exception('Illegal operator on member of different space')
        if isinstance(X,_composite_op):
            A = X.A
            B = X.B
            if   X.algebra=='+':
                Aa_block = self.act_on_vec_block(A,a_block)
                Ba_block = self.act_on_vec_block(B,a_block)
                return [ Aa+Ba for Aa,Ba in zip(Aa_block,Ba_block) ]
            elif X.algebra=='|':
                return self.act_on_vec_block(A, self.act_on_vec_block(B,a_block))
        elif isinstance(X,_dyadic_op):  return [ X|a for a in a_block ]
        else:                           return [ _member(v,the_space) for v in self.traits.act_on_vec_block(X.op,[a.v for a in a_block]) ]
    def back_act_on_vec_block(self,a_block,X):
        if isinstance(X,(float,int,complex)):  return [ a|X for a in a_block ]		# This only needs to be here since act_on_vec_block called directly since not finished yet
        # Note return above . . . the remainder does not appy to scalar operators
        the_space = X.space
        for a in a_block:
            if a.space is not the_space:  raise Exception('Illegal operator on member of different space')
        if isinstance(X,_composite_op):
            A = X.A
            B = X.B
            if   X.algebra=='+':
                aA_block = self.back_act_on_vec_block(a_block,A)
                aB_block = self.back_act_on_vec_block(a_block,B)
                return [ aA+aB for aA,aB in zip(aA_block,aB_block) ]
            elif X.algebra=='|':
                return self.back_act_on_vec_block(self.back_act_on_vec_block(a_block,A),B)
        elif isinstance(X,_dyadic_op):  return [ a|X for a in a_block ]
        else:                           return [ _member(v,the_space) for v in self.traits.back_act_on_vec_block([a.v for a in a_block],X.op) ]
    # How woudl you call this if you really just did have to lists.  Global function?  space.dot()?
    def dot_vec_blocks(self,a_block,b_block):
        the_space = a_block[0].space
        for a in a_block:
            if a.space is not the_space:  raise Exception('Illegal inner product of members of different spaces')
        for b in b_block:
            if b.space is not the_space:  raise Exception('Illegal inner product of members of different spaces')
        values = self.traits.dot_vec_blocks([a.v for a in a_block],[b.v for b in b_block])
        for a,val_row in zip(a_block,values):
            for b,val in zip(b_block,val_row):
                if a is b:
                    if abs(val)!=0 and abs(val.imag/val.real)>1e-14:  raise Exception('Something seriously wrong.  Norm^2 is not real.')      # Danger, hardcoded threshold!
        return values



# I think the point of spinning off the traits above was just modularity?  Or is it a historical artifact?
# Leaning towards historical artifact, since it would be nice to have the traits functions above available from the outside (see hack for dot_vec_blocks)

class linear_inner_product_space(object):
    def __init__(self,traits):
        self.traits = _space_traits(traits)
        self.field  = self.traits.field
        self.aux    = traits.aux
    def member(self,v):
        self.traits.check_member(v)
        return _member(v,self)
    def lin_op(self,op):
        functionable = self.traits.check_lin_op(op)
        if functionable:  return _fn_prim_op(op,self)
        else:             return _primitive_op(op,self)
    def dot_vec_blocks(self,a_block,b_block):  return self.traits.dot_vec_blocks(a_block,b_block)




# The idea is that one can build a space traits class and instance from scratch or use this class to build an instance by just assigning functions
# to the members in __init__.  If there are some functions that are not needed, they can be safely left out and the code will crash if you
# try to use them.  If any of these (particularly the "check" routines) need access to any stored information, these functions can themselves
# be callable classes or generated by external functions that store that data.

def _default_dot(v,w):                    raise NotImplementedError
def _default_add_to(v,w,n=1):             raise NotImplementedError
def _default_scale(n,v):                  raise NotImplementedError
def _default_copy(v):                     raise NotImplementedError
def _default_act_on_vec(op,v):            raise NotImplementedError
def _default_back_act_on_vec(v,op):       raise NotImplementedError
def _default_function_on_diags(func,op):  raise NotImplementedError
def _default_diagonal(op):                raise NotImplementedError
def _default_check_member(v):             raise NotImplementedError
def _default_check_lin_op(op):            raise NotImplementedError

class traits(object):
    def __init__(self, field, dot=_default_dot, add_to=_default_add_to, scale=_default_scale, copy=_default_copy, act_on_vec=_default_act_on_vec, back_act_on_vec=_default_back_act_on_vec, function_on_diags=_default_function_on_diags, diagonal=_default_diagonal, check_member=_default_check_member, check_lin_op=_default_check_lin_op, aux=None):
        self.field             = field
        self.dot               = dot
        self.add_to            = add_to
        self.scale             = scale
        self.copy              = copy
        self.act_on_vec        = act_on_vec
        self.back_act_on_vec   = back_act_on_vec
        self.function_on_diags = function_on_diags
        self.diagonal          = diagonal
        self.check_member      = check_member
        self.check_lin_op      = check_lin_op
        self.aux               = aux
