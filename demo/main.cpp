#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

//#include <divsufsort.h>

#include <parallel-range.h>

using namespace std;

void printPos(int pos, string& text) {
	for (int i = 0; i < 10 && i + pos < (int)text.size(); ++i) {
		cout << text[i+pos];
	} cout << endl;
}

int main(int argc, char* args[]) {
	if (argc != 3) {
		cout << "Expected two arguments (input and output file)."
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
	text += '\0';
	int size = text.size();
	int *SA,*ISA;
	{ // Allocate integer buffers.
		SA = new int[size];		
		ISA = new int[size];
		for (int i = 0; i < size; ++i) {
			SA[i] = i;
			ISA[i] = text[i];
		}
	}
	auto start = chrono::steady_clock::now();
	parallelrangelight(ISA, SA, size);
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << "Parallel Range Light time: " <<
		chrono::duration <double, milli> (diff).count()<< " ms" << endl;
	{ // Write output.
		fstream output_file(args[2], std::ios::out);
		for (int i = 0; i < size; ++i)
			output_file << SA[i] << endl;
		cout << "Output written to " << args[2] << endl;
		output_file.close();
	}
	/*
	// Time divsufsort.
	start = chrono::steady_clock::now();
	divsufsort((unsigned char*)text.data(), ISA, size);
	end = chrono::steady_clock::now();
	diff = end - start;
	cout << "DivSufSort time: " <<
		chrono::duration <double, milli> (diff).count()<< " ms" << endl;
	for (int i = 0; i < size; i++) {
		if (ISA[i] != SA[i]) {
			cout << "SA[" << i << "] = " << ISA[i] <<
				" expected but was " << SA[i] << endl;
			printPos(ISA[i], text);
			printPos(SA[i], text);
			break;
		}	
	}
	*/
	return 0;

}
