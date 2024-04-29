#include <iomanip>
#include <iostream>

using namespace std;

int main() {
	float f, d;
	f = 0.1;
	d = 0.2;
	std::cout << f << " " << d<< " " << setprecision(10) << " " << f+d << std::endl;
return 0;
}


