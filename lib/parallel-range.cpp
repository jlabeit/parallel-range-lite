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
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

#include "blockRadixSort.hpp"
#include "parallel.hpp"
#include "sequence.hpp"
#include "quicksort.hpp"

#include "timer.hpp"


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

class BV {
	public:
	int64_t n;
	BV(int64_t n_): n(n_) {
		num_words = (n + 63) / 64;
		data = new uint64_t[num_words];
	}	
	~BV() {
		delete [] data;
	}

	void set(int64_t pos) {
		data[pos / 64] |= 1ULL << (pos % 64);  
	}
	void set(int64_t pos, bool value) {
		if (value)
			data[pos / 64] |= 1ULL << (pos % 64);  
		else 
			data[pos / 64] &= ~(1ULL << (pos % 64));  
	}

	void print() {
		for (int i = 0; i < n; i++)
			cout << is_set(i);
		cout << endl;
	}
	
	void fill(bool value) {
		if (value)
			parallel_for(int64_t i = 0; i < num_words; ++i)
				data[i] = ~0ULL;
		else
			parallel_for(int64_t i = 0; i < num_words; ++i)
				data[i] = 0ULL;
	}
	bool is_set(int64_t pos) const {
		return (data[pos/64] & (1ULL << (pos % 64))) != 0;
	}

	uint64_t first_set(int64_t start, int64_t end) const {
		uint64_t tmp;
		if ((tmp = __builtin_ffsll(data[start/64] >> (start % 64)))) {
			// Answer in first block.
			start += tmp -1;
			if (start > end)
				return end;
			return start;
		} 
		start = (start + 63) / 64 * 64; // Round to next word.
		while (start < end) {
			if ((tmp = __builtin_ffsll(data[start/64]))) {
				// Answer in this block.
				start += tmp -1;
				if (start > end)
					return end;
				return start; 
			}
			start += 64;
		}
		return end;
	}
	// Returns last set bit in interval (end, start].
	uint64_t reverse_first_set(int64_t start, int64_t end) const {
		uint64_t word;
		if ((word = data[start/64] << (63 - (start % 64)))) {
			start -= __builtin_clzll(word);
			if (start < end)
				return end;
			return start;
		}
		start = start / 64 * 64 - 1; // Round to next word end.
		while (start > end) {
			if ((word = data[start / 64])) {
				start -= __builtin_clzll(word);
				if (start < end)
					return end;
				return start;
			}
			start -= 64;
		}
		return end;
	}
	bool operator==(const BV& obj) const {
		if (obj.n != n)
			return false;
		for (int64_t i = 0; i < n; i++)
			if (obj.is_set(i) != is_set(i))
				return false;
		return true;
	}
	BV(const vector<bool>& obj) : BV(obj.size()) {
		init(obj);		
	}

	void init(const vector<bool>& obj) {
		assert(obj.size() == (uint64_t)n);
		for (int64_t i = 0; i < n; ++i)
			set(i, obj[i]);
	}

	void swap(BV& bv) {
		std::swap(n, bv.n);
		std::swap(num_words, bv.num_words);
		std::swap(data, bv.data);
	}

	private:
	int64_t num_words;
	uint64_t* data;
};


template <class saidx_t, int32_t BLOCK_SIZE = 128*1024>
struct segment_info {
	// Text size.
	saidx_t n;
	saidx_t* SA;
	saidx_t* ISA;
	saidx_t num_blocks;
	BV bitvector;
	BV write_bv;
	vector<uint64_t> popcount_sum;
	// Precompute arrays pointing to the next/previous set bit.
	// Values are stored for all block boundaries.
	vector<saidx_t> next_one_arr;
	vector<saidx_t> previous_one_arr;

	segment_info(saidx_t n_, saidx_t* SA_, saidx_t* ISA_) : n(n_), bitvector(n), write_bv(n) {
		SA = SA_;
		ISA = ISA_;
		num_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
		bitvector.fill(0);
		bitvector.set(0);
		bitvector.set(n-1);
		next_one_arr.resize(num_blocks, n-1); next_one_arr[0] = 0;
		previous_one_arr.resize(num_blocks, 0);
		popcount_sum.resize(num_blocks);
		popcount_sum[0] = 0;
		parallel_for(saidx_t i = 1; i < num_blocks; ++i) 
			popcount_sum[i] = 1;
	}


	// Update precomputed structures to answere queries faster.
	void update_structure() {
		parallel_for(saidx_t b = 0; b < num_blocks; ++b) {
			previous_one_arr[b] = -1;
			next_one_arr[b] = n;
			saidx_t pos = b * BLOCK_SIZE - 1; // Do not skip first.
			saidx_t count = 0;
			while (next_one_in_block(pos, b)) {
				next_one_arr[b] = min(next_one_arr[b], pos);
				previous_one_arr[b] = pos;
				count++;
			}
			popcount_sum[b] = count;
		}
		saidx_t sum = 0;
		saidx_t tmp = 0;
		for (saidx_t b = 0; b < num_blocks; ++b) {
			tmp = sum;
			sum += popcount_sum[b];			
			popcount_sum[b] = tmp; // exlusive.
			if (b > 0)
				previous_one_arr[b] = max(previous_one_arr[b], previous_one_arr[b-1]);
		}
		for (saidx_t b = num_blocks-2; b != -1; --b) {
			next_one_arr[b] = min(next_one_arr[b], next_one_arr[b+1]);
		}
	}

	inline bool next_one_in_block(saidx_t& pos, saidx_t block) const {
		saidx_t end = std::min((block + 1) * BLOCK_SIZE, n);
		++pos;
		pos = bitvector.first_set(pos, end);
		// pos = std::find(bitvector.begin() + pos, bitvector.begin() + end, true) - bitvector.begin();
		return pos < end;
	}

	// Find next set bit after pos in block or end in a following block.
	// Return false if there is no following 1 bit or the the following 1
	// bit is a start of a segment starting in a new block.
	inline bool next_one(saidx_t& pos) const {
		++pos;
		saidx_t b = pos / BLOCK_SIZE;
		if (pos % BLOCK_SIZE == 0) {
			pos = next_one_arr[b];
			return pos < n;
		}
		--pos;
		if (next_one_in_block(pos, b))
			return true;
		if (b+1 < num_blocks)
			pos = next_one_arr[b+1];
			return pos < n;
		return false;
	}

	// Find previous set bit before pos in block or end in a following block.
	inline void previous_one(saidx_t& pos) const {
		// Assuming: -	There is always a 1 set before pos.
		// 	     -	Pos > 0
		saidx_t start = pos-1;
		saidx_t end = pos / BLOCK_SIZE * BLOCK_SIZE-1;
		pos = bitvector.reverse_first_set(start, end);
		if (pos == end) {
			pos = previous_one_arr[end / BLOCK_SIZE];
		}
	}

	// Find first segement start in a block or return false if block is empty.
	// Note: Pos is only used for output.
	inline bool find_first_open_in_block(saidx_t& pos, saidx_t block) const {
		pos = next_one_arr[block];
		if (pos / BLOCK_SIZE == block && pos < n) {
			if (popcount_sum[block] % 2) {
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
	void iterate_segments(function<void(saidx_t, saidx_t)> predicate) const {
		parallel_for(saidx_t b = 0; b < num_blocks; b++) {
			saidx_t start_segment, end_segment;
			if (find_first_open_in_block(start_segment, b)) {
				end_segment = start_segment;
				assert(next_one(end_segment));
				predicate(start_segment, end_segment);
				start_segment = end_segment;
				while (next_one_in_block(start_segment, b)) {
					end_segment = start_segment;
					assert(next_one(end_segment));
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
			if (popcount_sum[b] % 2) { // Close segment.
				start_segment = start_block;
				previous_one(start_segment);
				end_segment = start_block-1;
				assert(next_one(end_segment));
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
				assert(next_one(end_segment));
				predicate(start_segment, std::min(end_block, end_segment),
					       start_segment, end_segment);

			}	
		}
	}

	// Update additional data structure used to navigate segments.
	void update_segments(saidx_t offset) {
		write_bv.fill(false);
		cmp_offset<saidx_t> F(ISA, n, offset);
		//Timer::start("write");
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
				if ((old_f == cur_f) ^ (cur_f == new_f))
					write_bv.set(i);
			}
			});
		//Timer::stop("write");
		//Timer::start("name1");
		update_names_1();
		//Timer::stop("name1");
		bitvector.swap(write_bv);
		//Timer::start("update_structure");
		update_structure();
		//Timer::stop("update_structure");
		//Timer::start("name2");
		update_names_2();		
		//Timer::stop("name2");
	}

	// Assign to all suffixes in the current segments their position as ISA value.
	void update_names_1() {
		iterate_segments_blocked([this](saidx_t start, saidx_t end,
					saidx_t s, saidx_t e) {
				for (saidx_t i = start; i <= end; ++i) {
					ISA[SA[i]] = i;
				}
			});
	}
	// Assign all suffixes in the current segments the same ISA value.
	void update_names_2() {
		iterate_segments_blocked([this](saidx_t start, saidx_t end,
					saidx_t s, saidx_t e) {
				for (saidx_t i = start; i <= end; ++i) {
					ISA[SA[i]] = s;
				}
			});
	}

	void prefix_sort(saidx_t offset) {
		//Timer::start("sort");
		cmp_offset<saidx_t> F(ISA, n, offset); 	
		iterate_segments([F,offset, this](saidx_t start, saidx_t end) {
				saidx_t l = end-start+1;
				if (l >= 256)
					intSort::iSort(SA + start, l, n , F);
				else
					quickSort(SA + start, l, F);

				});
		//Timer::stop("sort");
	}
};


template <class saidx_t>
void paralleltrsort(saidx_t* ISA, saidx_t* SA, saidx_t n) {
	// segments = [0,n]
	segment_info<saidx_t> segs(n, SA, ISA);
	// make all comparisons
	segs.prefix_sort(0); // Not necessary if already sorted by first character.
	segs.update_segments(0);
	saidx_t offset = 1;
	while (segs.not_done()) {
		segs.prefix_sort(offset);
		segs.update_segments(offset);
	 	offset *= 2;
	}
	Timer::print();
}

void parallelrangelight(int32_t* ISA, int32_t* SA, int32_t n) {
	paralleltrsort(ISA, SA, n);
}
void parallelrangelight(int64_t* ISA, int64_t* SA, int64_t n) {
	paralleltrsort(ISA, SA, n);
}
