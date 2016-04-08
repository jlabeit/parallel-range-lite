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

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "blockRadixSort.hpp"
#include "parallel.hpp"
#include "sequence.hpp"
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


template <class saidx_t, int32_t BLOCK_SIZE = 128*1024>
struct segment_info {
	// Text size.
	saidx_t n;
	saidx_t* SA;
	saidx_t* ISA;
	saidx_t num_blocks;
	vector<bool> bitvector;
	vector<bool> write_bv;
	vector<bool> odd_prefix_sum;

	segment_info(saidx_t n_, saidx_t* SA_, saidx_t* ISA_) {
		n = n_;
		SA = SA_;
		ISA = ISA_;
		num_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
		bitvector.resize(n, false);
		write_bv.resize(n);
		bitvector[0] = true;
		bitvector[n-1] = true;
		odd_prefix_sum.resize(num_blocks, true);
		odd_prefix_sum[0] = false;
	}


	// Update precomputed structures to answere queries faster.
	void update_structure() {
		parallel_for(saidx_t b = 0; b < num_blocks; ++b) {
			saidx_t pos = b * BLOCK_SIZE - 1; // Do not skip first.
			saidx_t count = 0;
			while (next_one_in_block(pos, b)) {
				count++;
			}
			odd_prefix_sum[b] = count & 1;
		}
		bool sum = 0;
		bool tmp = 0;
		for (saidx_t b = 0; b < num_blocks; ++b) {
			tmp = sum;
			sum ^= odd_prefix_sum[b];			
			odd_prefix_sum[b] = tmp; // exlusive.
		}
	}

	inline bool next_one_in_block(saidx_t& pos, saidx_t block) const {
		saidx_t end = std::min((block + 1) * BLOCK_SIZE, n);
		++pos;
		while (pos < end) {
			if (bitvector[pos]) {
				return true;
			}
			++pos;
		}
		return false;
	}

	// Find next set bit after pos in block or end in a following block.
	// Return false if there is no following 1 bit or the the following 1
	// bit is a start of a segment starting in a new block.
	inline bool next_one(saidx_t& pos) const {
		saidx_t end = n; // TODO fast block skip implementation.
		++pos;
		while (pos < end) {
			if (bitvector[pos]) {
				return true;
			}
			++pos;
		}
		return false;
	}

	// Find previous set bit before pos in block or end in a following block.
	inline void previous_one(saidx_t& pos) const {
		// TODO fast block skip implementation.
		while (!bitvector[--pos]) 
			continue;
	}

	// Find first segement start in a block or return false if block is empty.
	// Note: Pos is only used for output.
	inline bool find_first_open_in_block(saidx_t& pos, saidx_t block) const {
		// TODO, maybe do faster precompute.
		pos = block * BLOCK_SIZE - 1; // -1 to not skip first bit.
		if (next_one_in_block(pos, block)) {
			if (odd_prefix_sum[block]) {
				return next_one_in_block(pos, block); // Skip end of segment.
			} else {
				return true;
			}
		}
		return false;
	}

	inline bool not_done() {
		saidx_t tmp = -1;
		return next_one(tmp);
	}

	// Important: Segments are processed in parallel, even in same block.
	// TODO skip empty blocks.
	void iterate_segments(function<void(saidx_t, saidx_t)> predicate) const {
		parallel_for(saidx_t b = 0; b < num_blocks; b++) {
			saidx_t start_segment, end_segment;
			if (find_first_open_in_block(start_segment, b)) {
				end_segment = start_segment;
				next_one(end_segment);
				predicate(start_segment, end_segment);
				start_segment = end_segment;
				while (next_one_in_block(start_segment, b)) {
					end_segment = start_segment;
					next_one(end_segment);
					predicate(start_segment, end_segment);
					start_segment = end_segment;
				}
			}
		}
	}

	// Lambda arguments: start, end, global_start, global_end.
	// Important: Parallel calls only between blocks not between segments
	// of the same block.
	void iterate_segments_blocked(std::function<void(saidx_t, saidx_t, saidx_t, saidx_t)>
			predicate) const {
		parallel_for(saidx_t b = 0; b < num_blocks; b++) {
			saidx_t start_segment, end_segment;
			saidx_t start_block = b * BLOCK_SIZE;
			saidx_t end_block = std::min((b+1)*BLOCK_SIZE, n) - 1;
			if (odd_prefix_sum[b]) { // Close segment.
				start_segment = start_block;
				previous_one(start_segment);
				end_segment = start_block-1;
				next_one(end_segment);
				predicate(start_block, std::min(end_block, end_segment),
					       start_segment, end_segment);
			} else {
				end_segment = start_block-1;
			}
			while (end_segment <= end_block) {
				start_segment = end_segment;	
				if (!next_one_in_block(start_segment, b)) 
					break;
				end_segment = start_segment;				
				next_one(end_segment);
				predicate(start_segment, std::min(end_block, end_segment),
					       start_segment, end_segment);

			}	
		}
	}

	// Update additional data structure used to navigate segments.
	// TODO parallelize (use iterate_blocked_segments).
	void update_segments(saidx_t offset) {
		parallel_for(saidx_t i = 0; i < n; i++) // TODO do this more efficient.
			write_bv[i] = false;
		cmp_offset<saidx_t> F(ISA, n, offset);
		iterate_segments_blocked([&F, this](saidx_t start, saidx_t end,
					saidx_t start_segment, saidx_t end_segment) {
			saidx_t old_f, cur_f, new_f; 
			// old_f and new_f can only compare as equal if they
			// are in the segment bounds.
			cur_f = start_segment < start ? F(SA[start-1]) : n;
			new_f = F(SA[start]);
			for (saidx_t i = start; i <= end; i++) {
				old_f = cur_f; cur_f = new_f;
				new_f = i < end_segment ? F(SA[i+1]) : n;
				write_bv[i] = (old_f == cur_f) ^ (cur_f == new_f); 
			}
			});
		update_names_1();
		swap(write_bv, bitvector);
		update_structure();
		update_names_2();		
	}

	// Assign to all suffixes in the current segments their position as ISA value.
	void update_names_1() {
		// TODO: Test whats faster blocked or not blocked update.
		iterate_segments([this](saidx_t start, saidx_t end) {
				parallel_for (saidx_t i = start; i <= end; ++i) {
					ISA[SA[i]] = i;
				}
			});
	}
	// Assign all suffixes in the current segments the same ISA value.
	void update_names_2() {
		// TODO: Test whats faster blocked or not blocked update.
		iterate_segments([this](saidx_t start, saidx_t end) {
				parallel_for (saidx_t i = start; i <= end; ++i) {
					ISA[SA[i]] = start;
				}
			});
	}

	void prefix_sort(saidx_t offset) {
		cmp_offset<saidx_t> F(ISA, n, offset); 	
		iterate_segments([F,offset, this](saidx_t start, saidx_t end) {
				saidx_t l = end-start+1;
				if (l >= 256)
					intSort::iSort(SA + start, l, n , F);
				else
					quickSort(SA + start, l, F);

				});
	}
};


template <class saidx_t>
void paralleltrsort(saidx_t* ISA, saidx_t* SA, saidx_t n) {
	// segments = [0,n]
	segment_info<saidx_t> segs(n, SA, ISA);
	// make all comparisons
	cout << "starting prefix sort " << n << endl;
	segs.prefix_sort(0); // Not necessary if already sorted by first character.
	cout << "updating segments" << endl;
	segs.update_segments(0);
	saidx_t offset = 1;
	while (segs.not_done()) {
		cout << "offset: " << offset << endl;
		segs.prefix_sort(offset);
		segs.update_segments(offset);
	 	offset *= 2;
	}
}
