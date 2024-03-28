add_api_test (basics_test.cpp)
add_api_test (kmer_loss_test.cpp)
add_api_test (param_set_test.cpp)

add_api_test (find_test.cpp)
target_use_datasources (find_test FILES 150overlap4bins.bin)
target_use_datasources (find_test FILES 150overlap4bins.arg)
target_use_datasources (find_test FILES 150overlap16bins.bin)
target_use_datasources (find_test FILES 150overlap16bins.arg)