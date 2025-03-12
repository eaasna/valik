add_api_test (compute_bin_size_test.cpp)

target_use_datasources (compute_bin_size_test FILES ref_bin_0.fasta)
target_use_datasources (compute_bin_size_test FILES ref_bin_20.fasta)
target_use_datasources (compute_bin_size_test FILES ref_bin_40.fasta)
target_use_datasources (compute_bin_size_test FILES ref_bin_60.fasta)

target_use_datasources (compute_bin_size_test FILES s1001_ref.0.header)
target_use_datasources (compute_bin_size_test FILES s1001_ref.20.header)
target_use_datasources (compute_bin_size_test FILES s1001_ref.40.header)
target_use_datasources (compute_bin_size_test FILES s1001_ref.60.header)

target_use_datasources (compute_bin_size_test FILES s11111111_ref.0.header)
target_use_datasources (compute_bin_size_test FILES s11111111_ref.20.header)
target_use_datasources (compute_bin_size_test FILES s11111111_ref.40.header)
target_use_datasources (compute_bin_size_test FILES s11111111_ref.60.header)
