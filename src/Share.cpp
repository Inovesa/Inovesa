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
{}

Share Share::operator+=(const Share& other)
{
	__s += other.__s;
	return *this;
}
