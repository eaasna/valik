add_api_test (consolidate_matches_test.cpp)

target_use_datasources (consolidate_matches_test FILES multi_seq_query.fasta)
target_use_datasources (consolidate_matches_test FILES 8bins50overlap_dream_all.gff)
target_use_datasources (consolidate_matches_test FILES 16bins50overlap_dream_all.gff)
target_use_datasources (consolidate_matches_test FILES 8bins50overlap_reference_metadata.bin)
target_use_datasources (consolidate_matches_test FILES 16bins50overlap_reference_metadata.bin)
target_use_datasources (consolidate_matches_test FILES stellar_truth.gff)

target_use_datasources (consolidate_matches_test FILES stellar_truth_num12_dis3.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num12_dis8.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num12_dis13.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num9_dis3.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num9_dis8.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num9_dis13.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num3_dis3.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num3_dis8.gff)
target_use_datasources (consolidate_matches_test FILES stellar_truth_num3_dis13.gff)

add_api_test (stellar_match_test.cpp)
target_use_datasources (stellar_match_test FILES 8bins50overlap_reference_metadata.bin)
