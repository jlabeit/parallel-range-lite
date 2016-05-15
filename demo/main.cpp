#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <parallel-range.h>

using namespace std;

typedef uint32_t uint;
int main(int argc, char* args[]) {
	if (argc != 2) {
		cout << "Expected one argument (input file)."
			<< endl;	
		return -1;
	}

	string text;
	{ // Read input file.
		ifstream input_file(args[1]);
		input_file.seekg(0, ios::end);   
		text.reserve(input_file.tellg());
		input_file.seekg(0, ios::beg);
		text.assign((istreambuf_iterator<char>(input_file)),
				istreambuf_iterator<char>());
	}
	uint size = text.size();
	uint *SA,*ISA,*textInt;
	{ // Allocate integer buffers.
		SA = new uint[size];
		ISA = new uint[size];
		textInt = new uint[size];
		for (int i = 0; i < size; ++i) {
			SA[i] = i;
			textInt[i] = ISA[i] = (uint)(text[i]+128); // Cannot be negative.
		}
	}
	auto start = chrono::steady_clock::now();
	parallelrangelite(ISA, SA, size);
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << "Parallel Range Light time: " <<
		chrono::duration <double, milli> (diff).count()<< " ms" << endl;
	if (sufcheck(textInt, SA, size, false)) {
		cout << "Sufcheck failed!" << endl;
		return -1;
	}
	return 0;
}
