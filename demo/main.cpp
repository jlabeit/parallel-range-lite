#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <parallel-range.h>

using namespace std;

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
	int size = text.size();
	int *SA,*ISA,*textInt;
	{ // Allocate integer buffers.
		SA = new int[size];		
		ISA = new int[size];
		textInt = new int[size];
		for (int i = 0; i < size; ++i) {
			SA[i] = i;
			textInt[i] = ISA[i] = ((int)text[i]) + 256; // Cannot be negative.
			if (ISA[i] < 0) cout << i << endl;
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
