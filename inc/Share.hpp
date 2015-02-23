/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2007-2009: Peter Schregle (Fixed Point Math Library)         */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#ifndef SHARE_HPP
#define SHARE_HPP

#include <algorithm>
#include <array>
#include <climits>
#include <list>
#include <vector>

namespace vfps
{

/**
 * @brief The Share class is a fixed number implementation
 * designed to represent numbers between 0 and 1 (including 0 and 1).
 *
 * The Share class has been designed to represent parts of a whole.
 * In Contrast to floating point numbers [0,1]
 * Share does not suffer from cancellation or absorption.
 * (So when "Share s = Share::ONE - Share(foo)",
 * "Share(foo) + s == Share::ONE" is guaranteed to be true.)
 */
class Share
{
public: // constructors and destructors
	/**
	 * @brief Share
	 */
	Share();

	/**
	 * @brief Share
	 * @param other
	 */
	Share(const Share& other);

	/**
	 * @brief Share
	 * @param share
	 */
	Share(int32_t share);

	/**
	 * @brief Share
	 * @param share
	 *
	 * Easiest way to use is Share(Share::ONE*numerator/denominator).
	 */
	Share(uint32_t share);

	/**
	 * @brief Share
	 * @param share number in the range [0,1]
	 */
	Share(double share);

	~Share();

public: // important definitions
	/**
	 * @brief ONE
	 */
	static constexpr uint32_t ONE = UINT_MAX/2 + 1;

public: // operators
	/**
	 * @brief operator=
	 * @param other
	 * @return
	 */
	Share& operator=(Share other);

	/**
	 * @brief operator+=
	 * @param other
	 * @return
	 */
	Share& operator+=(const Share& other);

	/**
	 * @brief operator -=
	 * @param other
	 * @return
	 */
	Share& operator-=(const Share& other);

	/**
	 * @brief operator*=
	 * @param other
	 * @return
	 */
	Share& operator*=(const Share& other);

	/**
	 * @brief operator /=
	 * @param rhs
	 * @return
	 */
	Share& operator/=(const uint32_t rhs);


	inline explicit operator double() const
		{ return double(__s)/ONE; }

	/**
	 * @brief operator float
	 *
	 * The conversion to float is designed to give an estimate
	 * of the share in an easy format.
	 * Doing calculations with the returned float instead of the Share itself,
	 * you will loose precision.
	 */
	inline explicit operator float() const
		{ return float(__s)/ONE; }

	/**
	 * @brief operator int
	 */
	inline explicit operator int32_t() const
		{ return static_cast<uint32_t>(*this); }

	/**
	 * @brief operator unsigned int
	 */
	inline explicit operator uint32_t() const
		{ return __s; }

	/**
	 * @brief operator long int
	 */
	inline explicit operator int64_t() const
		{ return static_cast<uint32_t>(*this); }

	/**
	 * @brief operator long unsigned int
	 */
	inline explicit operator uint64_t() const
		{ return static_cast<unsigned int>(*this); }

public: // other functions
	/**
	 * @brief renormalize rescales the arguments so that their sum is Share::ONE
	 * @param n number of arguments
	 * @param args pointer to first element in an array of Share
	 *
	 * Finds offset from Share::ONE and distributes it accoding to the
	 * largest remainder method.
	 */
	friend void renormalize(size_t n, Share* args);

	/**
	 * @brief swap
	 * @param first
	 * @param second
	 */
	friend inline void swap(Share& first, Share& second)
		{ std::swap(first.__s, second.__s); }

private:
	/**
	 * @brief __s
	 */
	uint32_t __s;
};

/**
 * @brief operator+
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator+(Share lhs, const Share rhs)
	{ return lhs += rhs; }

/**
 * @brief operator+
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator+(Share lhs, const double rhs)
	{ return lhs+Share(rhs); }

/**
 * @brief operator-
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const int32_t rhs)
	{ return lhs-=Share(Share::ONE*rhs); }

/**
 * @brief operator-
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const uint32_t rhs)
	{ return lhs-=Share(Share::ONE*rhs); }

inline Share operator-(const double lhs, const Share rhs)
	{ return Share(lhs)-=rhs; }

/**
 * @brief operator-
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const Share rhs)
	{ return lhs -= rhs; }

/**
 * @brief operator-
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const double rhs)
	{ return lhs-Share(rhs); }


/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline double operator*(const double lhs, const Share& rhs)
	{ return lhs*static_cast<uint32_t>(rhs)/Share::ONE; }

/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline double operator*(const Share& lhs, const double rhs)
	{ return rhs*lhs; }

/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline float operator*(const float lhs, const Share& rhs)
	{ return lhs*static_cast<uint32_t>(rhs)/Share::ONE; }

/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline float operator*(const Share& lhs, const float rhs)
	{ return rhs*lhs; }

/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline uint32_t operator*(const uint32_t lhs, const Share& rhs)
	{ return (lhs)*static_cast<uint64_t>(rhs)/Share::ONE; }

/**
 * @brief operator*
 * @param num
 * @return
 */
inline uint32_t operator*(const Share& lhs, const uint32_t rhs)
	{ return rhs*lhs; }

/**
 * @brief operator*
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator*(Share lhs, const Share& rhs)
	{ lhs *= rhs; return lhs; }

/**
 * @brief operator*=
 * @param lhs
 * @param rhs
 * @return
 */
inline unsigned int operator*=(uint32_t lhs, const Share& rhs)
	{ lhs = lhs*rhs; return lhs; }

inline Share operator/(Share lhs, const uint32_t rhs)
	{ lhs /= rhs; return lhs; }

/**
 * @brief operator==
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator==(const Share& lhs, const Share& rhs)
	{ return static_cast<uint32_t>(lhs) == static_cast<uint32_t>(rhs); }

/**
 * @brief operator!=
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator!=(const Share& lhs, const Share& rhs)
	{ return !operator==(lhs,rhs); }

/**
 * @brief operator<
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator< (const Share& lhs, const Share& rhs)
	{ return static_cast<unsigned int>(lhs) < static_cast<unsigned int>(rhs); }

/**
 * @brief operator>
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator> (const Share& lhs, const Share& rhs)
	{ return  operator< (rhs,lhs); }

/**
 * @brief operator<=
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator<=(const Share& lhs, const Share& rhs)
	{ return !operator> (lhs,rhs); }

/**
 * @brief operator>=
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator>=(const Share& lhs, const Share& rhs)
	{ return !operator< (lhs,rhs); }

/**
 * @brief renormalize
 * @param n
 * @param args
 * @todo optimize runtime when no renormalization is needed
 */
void renormalize(size_t n, Share* args);

}

#endif // SHARE_HPP
