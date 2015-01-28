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
	Share(int share);

	/**
	 * @brief Share
	 * @param share
	 *
	 * Easiest way to use is Share(Share::ONE*numerator/denominator).
	 */
	Share(unsigned int share);

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
	static constexpr unsigned int ONE = UINT_MAX/2 + 1;

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
	Share& operator/=(const unsigned int rhs);


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
	inline explicit operator int() const
		{ return static_cast<unsigned int>(*this); }

	/**
	 * @brief operator unsigned int
	 */
	inline explicit operator unsigned int() const
		{ return __s; }

	/**
	 * @brief operator long int
	 */
	inline explicit operator long int() const
		{ return static_cast<unsigned int>(*this); }

	/**
	 * @brief operator long unsigned int
	 */
	inline explicit operator long unsigned int() const
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
	friend void renormalize(unsigned int n, Share* args);

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
	unsigned int __s;
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
 * @brief operator -
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const int rhs)
	{ return lhs-=rhs; }

/**
 * @brief operator -
 * @param lhs
 * @param rhs
 * @return
 */
inline Share operator-(Share lhs, const unsigned int rhs)
	{ return lhs-=rhs; }

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
inline unsigned int operator*(const unsigned int lhs, const Share& rhs)
	{ return (lhs)*static_cast<unsigned int>(rhs)/Share::ONE; }

/**
 * @brief operator*
 * @param num
 * @return
 */
inline unsigned int operator*(const Share& lhs, const unsigned int rhs)
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
inline unsigned int operator*=(unsigned int lhs, const Share& rhs)
	{ lhs = lhs*rhs; return lhs; }

inline Share operator/(Share lhs, const unsigned int rhs)
	{ lhs /= rhs; return lhs; }

/**
 * @brief operator==
 * @param lhs
 * @param rhs
 * @return
 */
inline bool operator==(const Share& lhs, const Share& rhs)
	{ return static_cast<unsigned int>(lhs) == static_cast<unsigned int>(rhs); }

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

void renormalize(unsigned int n, Share* args);

}

#endif // SHARE_HPP
