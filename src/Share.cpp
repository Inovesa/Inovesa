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
	std::vector<uint32_t> tmpshare;
	tmpshare.reserve(n);
	uint64_t sum=0;
	for( size_t i=0; i<n; i++) {
		sum += args[i].__s;
		tmpshare.push_back(args[i].__s);
	}

	if (sum == 0) {
		for( size_t i=0; i<n; i++) {
			tmpshare[i] = 1U;
		}
		sum = n;
	}
	std::list<std::pair<uint32_t,size_t>> remainders;
	uint32_t rensum=0;
	for( size_t i=0; i<n; i++) {
		uint64_t val = static_cast<uint64_t>(tmpshare[i])*Share::ONE;
		uint32_t remainder = val%sum;
		val /= sum;
		rensum += val;
		args[i].__s = val;
		remainders.push_back(std::pair<uint32_t,size_t>(remainder,i));
	}

	size_t rest = Share::ONE - rensum;
	if (rest > 0) {
		// sort remainders rescending
		remainders.sort([](	const std::pair<uint32_t,size_t> &lhs,
							const std::pair<uint32_t,size_t> &rhs)
						-> bool
						{ return lhs.first > rhs.first; }
		);
		do {
			args[remainders.front().second].__s += 1U;
			remainders.pop_front();
			rest--;
		} while (rest > 0);
	}
}

Share& Share::operator=(Share other)
{
	swap(*this,other);
	return *this;
}

}
