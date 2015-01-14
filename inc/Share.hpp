#ifndef SHARE_HPP
#define SHARE_HPP

#include <algorithm>
#include <climits>
#include <iostream>

class Share
{
public:
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
	Share(unsigned int share);

public: // operators
	/**
	 * @brief operator *=
	 * @param other
	 * @return
	 */
	Share& operator*=(const Share& other)
		{ __s += other.__s; return *this; }

	/**
	 * @brief ONE
	 */
	static constexpr long unsigned int ONE = (long unsigned int)(UINT_MAX)+1;

	/**
	 * @brief operator =
	 * @param other
	 * @return
	 */
	Share& operator=(Share other)
		{ swap(*this,other); return *this; }

	/**
	 * @brief operator +=
	 * @param other
	 * @return
	 */
	Share operator+=(const Share& other);

public:
	inline unsigned int toInt() const
		{ return (unsigned int)(__s); }

private:
	/**
	 * @brief __s
	 */
	unsigned long int __s;

public:
	friend void swap(Share first, Share second);
};

std::ostream& operator<<(std::ostream& os, const Share& obj)
{
  os << obj.toInt();
  return os;
}

inline Share operator*(Share lhs, const Share& rhs)
{
	lhs *= rhs;
	return lhs;
}

inline void swap(Share first, Share second)
	{ std::swap(first.__s, second.__s); }


#endif // SHARE_HPP
