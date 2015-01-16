#include "Share.hpp"

Share::Share() :
	__s(ONE)
{
}

Share::Share(const Share& other) :
	__s(other.__s)
{
}

Share::Share(unsigned int share) :
	__s(share)
{
}

Share::Share(double share) :
	__s(share*ONE)
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

void renormalize(unsigned int n, Share* args)
{
	std::vector<unsigned long int> tmpshare;
	tmpshare.reserve(n);
	unsigned long int sum=0;
	for( unsigned int i=0; i<n; i++)
	{
		sum += static_cast<unsigned int>(args[i]);
		tmpshare.push_back(static_cast<unsigned int>(args[i]));
	}
	std::list<std::array<unsigned int,2>> remainders;
	unsigned long int rensum=0;
	for( unsigned int i=0; i<n; i++)
	{
		long unsigned int val = tmpshare[i]*Share::ONE;
		val /= sum;
		unsigned int remainder = val%sum;
		rensum += val;
		args[i] = (unsigned int)(val);
		remainders.push_back({{remainder,i}});
	}

	remainders.sort();
	unsigned int rest = Share::ONE - rensum;
	while (rest > 0)
	{
		args[remainders.front()[1]] += 1U;
		remainders.pop_front();
		rest--;
	}
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
