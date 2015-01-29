#include "Share.hpp"

namespace vfps
{

Share::Share() :
	__s(ONE)
{
}

Share::Share(const Share& other) :
	__s(other.__s)
{
}

Share::Share(int32_t share) :
	Share(static_cast<uint32_t>(share))
{
}

Share::Share(uint32_t share) :
	__s(share)
{
}

Share::Share(double share)
{
	if (share >= 0) {
		__s = share*ONE;
	} else {
		__s = (1+share)*ONE;
	}

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
	uint64_t s = __s;
	s *= other.__s;
	s /= ONE;
	__s = s;
	return *this;
}

Share&Share::operator/=(const uint32_t rhs)
{
	__s /= rhs;
	return *this;
}

void renormalize(size_t n, Share* args)
{
	std::vector<unsigned long int> tmpshare;
	tmpshare.reserve(n);
	unsigned long int sum=0;
	for( unsigned int i=0; i<n; i++)
	{
		sum += static_cast<unsigned int>(args[i]);
		tmpshare.push_back(static_cast<unsigned int>(args[i]));
	}
	if (sum == 0)
	{
		for( unsigned int i=0; i<n; i++)
		{
			tmpshare[i] = 1UL;
		}
		sum = n;
	}
	std::list<std::array<unsigned int,2>> remainders;
	long int rensum=0;
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
	long int rest = static_cast<long int>(Share::ONE) - rensum;
	if (rest >= 0)
		while (rest > 0)
		{
			args[remainders.front()[1]] += 1U;
			remainders.pop_front();
			rest--;
		}
	else
	{
		while (rest < 0)
		{
			args[remainders.front()[1]] -= 1U;
			remainders.pop_front();
			rest++;
		}
	}
}

Share& Share::operator=(Share other)
{
	swap(*this,other);
	return *this;
}

}
