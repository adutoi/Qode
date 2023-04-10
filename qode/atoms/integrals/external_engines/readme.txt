someday these should all have a uniform interface

right now, we force the return of a simple numpy tensor of AO integrals in bra-bra-ket-ket for the two-electron ints.
eventually, these will be of my special tensor class, where the internal storage is not specified.
the interface should remain bra-bra-ket-ket, however.
