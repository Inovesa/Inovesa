#include "Share.hpp"

Share::Share() :
	__s(ONE)
{
}

Share::Share(const Share& other) :
	Share()
{
	operator=(other);
}

Share::Share(unsigned int share) :
	__s(share)
{
}

Share::~Share()
{
}

Share& Share::operator+=(const Share& other)
{
	__s += other.__s;
	return *this;
}

Share& Share::operator-=(const Share& other)
{
	__s -= other.__s;
	return *this;
}

Share& Share::operator*=(const Share& other)
{
	__s *= other.__s;
	__s /= ONE;
	return *this;
}

Share& Share::operator=(Share other)
{
	swap(*this,other);
	return *this;
}

std::istream&operator>>(std::istream& is, Share& obj)
{
	unsigned int val;
	is >> val;
	obj = val;
	if (false)
	{
		is.setstate(std::ios::failbit);
	}
	return is;
}

std::ostream&operator<<(std::ostream& os, const Share& obj)
{
	os << static_cast<float>(obj);
	return os;
}
