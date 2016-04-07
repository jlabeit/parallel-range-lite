// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "parallel.hpp"
#include "sequence.hpp"

#include <iostream>
#include <algorithm>
#include <functional>

#include "blockRadixSort.hpp"
#include "quicksort.hpp"

using namespace std;


//#define printInfo

#ifdef printInfo
#define nextTimeM(_str) nextTime(_str)
#else
#define nextTimeM(_str) 
#endif



template <class saidx_t>
struct cmp_offset {
	saidx_t *ISA;
	saidx_t n,offset;
	cmp_offset(saidx_t* ISA_, saidx_t n_, saidx_t offset_ ) : ISA(ISA_), n(n_), offset(offset_) {}
	// Return rank of the suffix offset characters after the suffix at position pos
	saidx_t operator() (const saidx_t& pos)const {
		return (pos + offset >= n) ? (n-pos-offset) : ISA[pos+offset]; 
	}
	bool operator()(const saidx_t& a, const saidx_t& b)const {
		return (*this)(a) < (*this)(b);
	}
};

const int32_t BLOCK_SIZE = 128*1024;

template <class saidx_t>
struct segment_info {
	// Text size.
	saidx_t n;
	saidx_t* SA;
	saidx_t* ISA;
	saidx_t num_blocks;
	vector<bool> bitvector;
	// Is set to fals if all suffixes are compared different.
	volatile bool done;

	segment_info(saidx_t n_, saidx_t* SA_, saidx_t* ISA_) {
		n = n_;
		SA = SA_;
		ISA = ISA_;
		num_blocks = n / BLOCK_SIZE + 1;
		bitvector.resize(n, false);
		bitvector[0] = true;
		bitvector[n-1] = true;
		done = false;
	}

	// Find next set bit after pos in block or end in a following block.
	// Return false if there is no following 1 bit or the the following 1
	// bit is a start of a segment starting in a new block.
	inline bool next_one(saidx_t& pos, saidx_t block) const {
		// TODO implement.
	}

	// Find first segement start in a block or return false if block is empty.
	inline bool find_first_open(saidx_t& pos, saidx_t block) const {
		// TODO implement.
	}

	// TODO skip empty blocks.
	void iterate_segments(std::function<void(saidx_t, saidx_t)> predicate) const {
		parallel_for(saidx_t b = 0; b < num_blocks; b++) {
			saidx_t start_segment, end_segment;
			if (find_first_open(start_segment, b)) {
				end_segment = start_segment;
				next_one(end_segment, b);
				predicate(start_segment, end_segment);
				start_segment = end_segment;
				while (next_one(start_segment, b)) {
					end_segment = start_segment;
					next_one(end_segment, b);
					predicate(start_segment, end_segment);
				}
			}
		}
	}

	// Lambda arguments: start, end, global_start, global_end.
	void iterate_blocked_segments(std::function<void(saidx_t, saidx_t, saidx_t, saidx_t)>
			predicate) const {
		iterate_segments([&predicate](saidx_t start, saidx_t end) {
				saidx_t bs = start / BLOCK_SIZE;
				saidx_t be = end / BLOCK_SIZE;
				parallel_for (saidx_t b = bs; b <= be; ++b) {
					saidx_t s = std::max(b * BLOCK_SIZE, start);	
					saidx_t e = std::min((b+1)* BLOCK_SIZE, end);
					predicate(s, e, start, end);
				}
				});
	}

	// Update additional data structure used to navigate segments.
	// TODO parallelize (use iterate_blocked_segments).
	void update_segments(saidx_t offset) {
		cmp_offset<saidx_t> F(ISA, n, offset);
		iterate_segments([&F, this](saidx_t start, saidx_t end) {
			saidx_t old_f, cur_f, new_f; 
			cur_f = F(SA[start]);
			new_f = F(SA[start+1]);
			// Beginning of segment.
			if (start == 0) {
				bitvector[start] = cur_f == new_f;			
			} else {
				old_f = F(SA[start-1]);
				bitvector[start] = (old_f != cur_f) ^ (cur_f != new_f); 
			}
			old_f = cur_f; cur_f = new_f;
			// Inner part.
			for (saidx_t i = start+1; i < end; i++) {
				new_f = F(SA[i + 1]);
				bitvector[start] = (old_f != cur_f) ^ (cur_f != new_f); 
				old_f = cur_f; cur_f = new_f;
			}
			// End of segment.
			bitvector[end] = old_f == cur_f;
			});
	}

	// Update 'names' (rank values in ISA). 
	void update_names() {
	}
};




//	if (l >= 256) 
//		intSort::iSort(SAi, l, n , cmp_offset<saidx_t>(SA, ISA, n, offset));
//	else
//		quickSort(SAi, l, cmp_offset<saidx_t>(SA, ISA, n, offset));



template <class saidx_t>
void paralleltrsort(saidx_t* ISA, saidx_t* SA, saidx_t n) {
	
	// segments = [0,n]
	segment_info<saidx_t> segs(n, SA, ISA);
	// make all comparisons
	segs.update_segments(0);
	saidx_t offset = 1;
	while (!segs.done) {
	 	parallel_for(saidx_t b = 0; b < segs.num_blocks; ++b) {	
	 		sort_segs_in_block(b, segs);	
		}
		segs.update_segments();
		segs.update_names();
	 	offset *= 2;
	}
}
