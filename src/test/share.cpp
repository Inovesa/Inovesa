#include <iostream>
#include <stdlib.h>
#include <vector>

#include "Share.cpp"

int main(int argc, char** argv)
{	
	size_t n;
	if (argc == 1) {
		n = 4;
	} else {
		n = argc -1;
	}
	std::vector<vfps::Share> a(n);
	
	if (argc == 1) {
		a[0] = 0.25;
		a[1] = 0.125;
		a[2] = 0.125;
		a[3] = 0;
	} else {
		for (int i=1; i<argc; i++) {
			a[i-1] = std::atof(argv[i]);
		}
	}
	vfps::Share sum(0);
	for (size_t i=0; i<n; i++) {
		std::cout << static_cast<float>(a[i]) << std::endl;
		sum += a[i];
	}
	std::cout << "_____" << std::endl
			  << static_cast<float>(sum) << std::endl << std::endl;


	std::cout << "Normalized to:" << std::endl;
	vfps::renormalize(n,a.data());

	sum = 0;
	for (size_t i=0; i<n; i++) {
		std::cout << static_cast<float>(a[i]) << std::endl;
		sum += a[i];
	}
	std::cout << "_____" << std::endl
			  << static_cast<float>(sum) << std::endl;

	std::cout << "Arithmetics:" << std::endl;
	vfps::Share diff = 1-a[0];
	std::cout << "1 - " << static_cast<float>(a[0]) << " = "
			  << static_cast<float>(diff) << std::endl;

	
	return 0;
}

