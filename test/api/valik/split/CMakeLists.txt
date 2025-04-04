add_api_test (write_seg_sequences_test.cpp)

target_use_datasources (write_seg_sequences_test FILES database.fasta)
target_use_datasources (write_seg_sequences_test FILES write_out_20_16_reference_metadata.bin)
target_use_datasources (write_seg_sequences_test FILES write_out_20_4_reference_metadata.bin)

add_api_test (trim_fasta_id_test.cpp)

add_api_test (metadata_test.cpp)
target_use_datasources (metadata_test FILES query.fasta)
target_use_datasources (metadata_test FILES ref.fasta)
target_use_datasources (metadata_test FILES bin_0.fasta)
target_use_datasources (metadata_test FILES bin_1.fasta)
target_use_datasources (metadata_test FILES bin_2.fasta)
target_use_datasources (metadata_test FILES bin_3.fasta)
target_use_datasources (metadata_test FILES bin_4.fasta)
target_use_datasources (metadata_test FILES bin_5.fasta)
target_use_datasources (metadata_test FILES bin_6.fasta)
target_use_datasources (metadata_test FILES bin_7.fasta)
target_use_datasources (metadata_test FILES various_chromosome_lengths.fasta)
target_use_datasources (metadata_test FILES 20overlap4bins.bin)
target_use_datasources (metadata_test FILES 20overlap16bins.bin)

add_api_test (adjust_bin_count_test.cpp)
