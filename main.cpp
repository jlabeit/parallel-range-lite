#include <fstream>
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
	int *SA,*ISA;
	{ // Allocate integer buffers.
		SA = new int[size];		
		ISA = new int[size];
		for (int i = 0; i < size; ++i) {
			SA[i] = i;
			ISA[i] = text[i];
		}
	}
	cout << "Starting tr sort" << endl;
	paralleltrsort(ISA, SA, size);
	cout << "Done sorting..." << endl;
	{
		fstream output_file(args[2], ios::out);
		for (int i = 0; i < size; ++i) {
			output_file << SA[i] << endl;
		}
	}
	return 0;

}
