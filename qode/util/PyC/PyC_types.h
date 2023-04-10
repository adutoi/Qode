/*   (C) Copyright 2019 Anthony D. Dutoi
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _PYC_TYPES_H_
#define _PYC_TYPES_H_

#include <unistd.h>  // Defines ssize_t



/*
 * Defines some C data types that are guaranteed to be compatible with properly communicated numpy types.
 * See comments in the sister file:  PyC_types.py
 */
typedef int     Int;
typedef ssize_t BigInt;
typedef double  Double;

/*
 * The types above to which native python types are automatically converted to and from.
 */
typedef BigInt PyInt;
typedef Double PyFloat;



#endif	// _PYC_TYPES_H_
