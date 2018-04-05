/* NDArray.hpp:  A STL-style wrapper for Array.h
   Copyright (C) 2018 Patrik Schoenfeldt, Karlsruhe Institute of Technology

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include "Array.h"

namespace Array
{

template <typename T>
using NDarray1 = Array::array1<T>;

} // namespace Array

#endif // NDARRAY_HPP
