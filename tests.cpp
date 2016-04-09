#include <gtest/gtest.h>

#include "parallel-range.hpp"

TEST(SegmentInfo, constructor) {
	segment_info<int, 5> segs(7, NULL, NULL);
	EXPECT_EQ(segs.n, 7);	
	EXPECT_EQ(segs.num_blocks, 2);	
	EXPECT_EQ(segs.bitvector, vector<bool>({1,0,0,0,0,0,1}));	
	EXPECT_EQ(segs.odd_prefix_sum, vector<bool>({0,1}));
}

TEST(SegmentInfo, next_one) {
	segment_info<int, 5> segs(7, NULL, NULL);
	segs.bitvector = vector<bool>({1,0,1,0,1,0,1});	
	segs.update_structure();
	int pos = 0;
	EXPECT_EQ(true, segs.next_one(pos));
	EXPECT_EQ(2, pos);
	EXPECT_EQ(true, segs.next_one(pos));
	EXPECT_EQ(4, pos);
	EXPECT_EQ(true, segs.next_one(pos));
	EXPECT_EQ(6, pos);
	EXPECT_EQ(false, segs.next_one(pos));
}

TEST(SegmentInfo, update_structure) {
	segment_info<int, 5> segs(9, NULL, NULL);
	segs.bitvector = vector<bool>({0,1,1,1,0,0,0,1,0});	
	segs.update_structure();
	EXPECT_EQ(vector<bool>({0,1}), segs.odd_prefix_sum);
	segs.bitvector = vector<bool>({0,1,1,0,0,0,0,0,0});	
	segs.update_structure();
	EXPECT_EQ(segs.odd_prefix_sum, vector<bool>({0,0}));
}

TEST(SegmentInfo, find_first_open) {
	segment_info<int, 5> segs(9, NULL, NULL);
	segs.bitvector = vector<bool>({0,1,1,1,0,0,0,1,0});	
	segs.update_structure();
	int pos = -42; // Should not matter.
	EXPECT_EQ(true, segs.find_first_open_in_block(pos, 0));
	EXPECT_EQ(1, pos);
	// Last block only contains end not start!.
	EXPECT_EQ(false, segs.find_first_open_in_block(pos, 1));
}


TEST(SegmentInfo, iterate_segments) {
	segment_info<int, 5> segs(9, NULL, NULL);
	segs.bitvector = vector<bool>({1,1,0,1,0,0,0,1,0});	
	segs.update_structure();
	vector<vector<int>> result;
	segs.iterate_segments([&result](int s, int e) { result.push_back({s,e});});
	EXPECT_EQ( vector<vector<int>>({{0,1},{3,7}}), result);
}

TEST(SegmentInfo, iterate_segments_blocked) {
	segment_info<int, 5> segs(9, NULL, NULL);
	segs.bitvector = vector<bool>({1,1,0,1,0,0,0,1,0});	
	segs.update_structure();
	vector<vector<int>> result;
	segs.iterate_segments_blocked([&result](int s, int e, int gs, int ge) {
			result.push_back({s,e, gs, ge});
			});
	EXPECT_EQ( vector<vector<int>>({{0,1,0,1},{3,4,3,7}, {5,7,3,7}}), result);
}

TEST(SegmentInfo, iterate_segments_blocked_without_update) {
	segment_info<int, 5> segs(10, NULL, NULL);
	segs.bitvector = vector<bool>({1,0,0,0,0,0,0,0,0,1});	
	vector<vector<int>> result;
	segs.iterate_segments_blocked([&result](int s, int e, int gs, int ge) {
			result.push_back({s,e, gs, ge});
			});
	EXPECT_EQ( vector<vector<int>>({{0,4,0,9},{5,9,0,9}}), result);
}

TEST(SegmentInfo, update_segments) {
	vector<int> ISA = {0,0,0,0,4,5,5,5};
	vector<int> SA = {0,1,2,3,4,5,6,7};
	segment_info<int, 5> segs(8, SA.data(), ISA.data());
	segs.update_segments(0);
	EXPECT_EQ(vector<bool>({1,0,0,1,0,1,0,1}), segs.bitvector);
	EXPECT_EQ(4, segs.ISA[4]);
}

TEST(SegmentInfo, update_names_1) {
	vector<int> ISA = {0,0,0,0,0,0,6,7};
	vector<int> SA = {0,1,2,3,4,5,6,7};
	segment_info<int, 5> segs(8, SA.data(), ISA.data());
	segs.bitvector = vector<bool>({0,1,0,1,1,1,0,0});
	segs.update_structure();
	segs.update_names_1();
	vector<int> result({0,1,2,3,4,5,6,7});
	for (size_t i = 0; i < result.size(); ++i) {
		EXPECT_EQ(result[i], segs.ISA[i]);	
	}
}

TEST(SegmentInfo, update_names_2) {
	vector<int> ISA = {0,0,0,0,0,0,6,7};
	vector<int> SA = {0,1,2,3,4,5,6,7};
	segment_info<int, 5> segs(8, SA.data(), ISA.data());
	segs.bitvector = vector<bool>({0,1,0,1,1,1,0,0});
	segs.update_structure();
	segs.update_names_2();
	vector<int> result({0,1,1,1,4,4,6,7});
	for (size_t i = 0; i < result.size(); ++i) {
		EXPECT_EQ(result[i], segs.ISA[i]);	
	}
}

TEST(SegmentInfo, suffix_sort) {
	vector<int> ISA = {3,2,1,0,0,0,7,6};
	vector<int> SA = {0,1,2,3,4,5,6,7};
	segment_info<int, 5> segs(8, SA.data(), ISA.data());
	segs.bitvector = vector<bool>({1,0,0,1,0,0,1,1});
	segs.update_structure();
	segs.prefix_sort(0);
	EXPECT_EQ(vector<int>({3,2,1,0,4,5,7,6}), SA);
}

//TEST(ParallelTrSort, simple_integration) {
//	vector<int> ISA = {0,0,0,0,0,1,0,1};
//	vector<int> SA = {0,1,2,3,4,5,6,7};
//	paralleltrsort(ISA.data(), SA.data(), (int)SA.size());
//	vector<int> expectSA = {0,1,2,3,4,5,6,7};
//	vector<int> expectISA = {0,1,2,3,4,5,6,7};
//	EXPECT_EQ(expectSA, SA);
//	EXPECT_EQ(expectISA, ISA);
//}

TEST(ParallelTrSort, error) {
	// 	   SA: 0123456789
	string text = "banananaaa";
	vector<int> ISA,SA;
	for (char c : text) {
		ISA.push_back(c);
		SA.push_back(SA.size());
	}
	paralleltrsort(ISA.data(), SA.data(), (int)SA.size());
	vector<int> expectSA = {9,8,7,5,3,1,0,6,4,2};
	EXPECT_EQ(expectSA, SA);
}
int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
