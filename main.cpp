#include <iostream>
#include <string>

#include "parallel-range.hpp"

using namespace std;

int main(int argc, char* args[]) {
	if (argc != 3) {
		cout << "Expected two arguments (input and output file)."
			<< endl;	
		return -1;
	}
	return 0;

}
