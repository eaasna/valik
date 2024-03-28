cmake_minimum_required (VERSION 3.16)

include (cmake/app_datasources.cmake)

declare_datasource (FILE database.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/split/database.fasta
                URL_HASH SHA256=7c7a8fcdd52a932cda76219f24024c1624292377103d9fd5a55abd288c6072be)
declare_datasource (FILE single_query.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single_query.fasta
                URL_HASH SHA256=e0b4924e4b9b47df8ecaed90c508e2786f27ec8d54ad80c059307c9f7ccbbb12)
declare_datasource (FILE single_reference.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single_reference.fasta
                URL_HASH SHA256=53d633474b01a68927d3ab1fd970b200e96403bb1fdcc53feb0367a2093be273)
declare_datasource (FILE various_chromosome_lengths.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/split/various_chromosome_lengths.fasta
                URL_HASH SHA256=7c7a8fcdd52a932cda76219f24024c1624292377103d9fd5a55abd288c6072be)


declare_datasource (FILE 150overlap16bins13window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins13window1errors.gff.out
                URL_HASH SHA256=118f6f235b77ffe9fcf91067cbf718a2ca7d3ce0d3406b9f83a8238be603e155)
declare_datasource (FILE 150overlap16bins13window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins13window.ibf
                URL_HASH SHA256=f774150b74f9c66f108b17bcccf8f7e7782c9940c6ad0faa8cc7910c4a397725)
declare_datasource (FILE 150overlap16bins15window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins15window1errors.gff.out
                URL_HASH SHA256=118f6f235b77ffe9fcf91067cbf718a2ca7d3ce0d3406b9f83a8238be603e155)
declare_datasource (FILE 150overlap16bins15window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins15window.ibf
                URL_HASH SHA256=fde51c6b696e5b0e1904fd20c9396385f0e1770beb24e9a6c63de57cdc9ae9e8)
declare_datasource (FILE 150overlap16bins.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins.arg
                URL_HASH SHA256=5abba452e93768fc2654f3ef3b72b74dfffda1061e25c2450e87b39bd5be6a0b)
declare_datasource (FILE 150overlap16bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins.bin
                URL_HASH SHA256=e12cdf696d2384083acb8cbfec40de6c350001f367ed231b978d2adfaece27e3)
declare_datasource (FILE 150overlap4bins13window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins13window1errors.gff.out
                URL_HASH SHA256=d45982f54310c2037e3b75da6cfff62179eff043d41e4987a575f0727415c4ea)
declare_datasource (FILE 150overlap4bins13window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins13window.ibf
                URL_HASH SHA256=690c52011ba874eb76aee594a66cf682726a5332eaea76e4272d5f40dfe12865)
declare_datasource (FILE 150overlap4bins15window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins15window1errors.gff.out
                URL_HASH SHA256=a12cad7ceae234ab1dd051d08e9a8463ab32793980ee42d383cf63fcf890405b)
declare_datasource (FILE 150overlap4bins15window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins15window.ibf
                URL_HASH SHA256=f38d6f8c9bbeb9eaf55690e979a32517b7963673a8b23433308db8b92333bf4a)
declare_datasource (FILE 150overlap4bins.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins.arg
                URL_HASH SHA256=5abba452e93768fc2654f3ef3b72b74dfffda1061e25c2450e87b39bd5be6a0b)
declare_datasource (FILE 150overlap4bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins.bin
                URL_HASH SHA256=751fb56b2f26371e1b672de8bd304051f937ff21b9919bc77c3a2df53f4e88ee)


declare_datasource (FILE 20overlap16bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/multi/20overlap16bins.bin
                URL_HASH SHA256=89f8a14fe18fa7a64fffb46516eadbde46d4eb0992305587a22be3602f454b77)
declare_datasource (FILE 20overlap4bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/multi/20overlap4bins.bin
                URL_HASH SHA256=87d789f02607c0bda9b881d3a69beee8cf72e4d1a0836401d0956caa60b3c82c)


declare_datasource (FILE write_out_20_16_reference_metadata.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_16/reference_metadata.arg
                URL_HASH SHA256=111ebf21c6bfa60f37392d13f67e81580f6d40ae332e6ffd547b624cca797645)
declare_datasource (FILE write_out_20_16_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_16/reference_metadata.bin
                URL_HASH SHA256=41a02debb90f108a779e42079a0637ac0c7aebea40766f98095b060ebe37a599)
declare_datasource (FILE write_out_20_4_reference_metadata.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_4/reference_metadata.arg
                URL_HASH SHA256=111ebf21c6bfa60f37392d13f67e81580f6d40ae332e6ffd547b624cca797645)
declare_datasource (FILE write_out_20_4_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_4/reference_metadata.bin
                URL_HASH SHA256=41a02debb90f108a779e42079a0637ac0c7aebea40766f98095b060ebe37a599)
declare_datasource (FILE 8bins19window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/build/8bins19window.ibf
                URL_HASH SHA256=3a13c890650bf857770816244ed9420295ad8bbe681dac335f687863fc79a603)
declare_datasource (FILE 8bins23window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/build/8bins23window.ibf
                URL_HASH SHA256=250578b9e0c47df929ed628931441ada4569ab7df193a5ecb5c069e6339bd56a)
declare_datasource (FILE bin_0.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_0.fasta
                URL_HASH SHA256=f9836f233fe459f8e387f8723dc030a10e44f3490cc1c89bed36222742bd6c35)
declare_datasource (FILE bin_1.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_1.fasta
                URL_HASH SHA256=782cb2eb50722e4a4fb2b5ca82b39817bf78305146a0ae086bab3273997e9237)
declare_datasource (FILE bin_2.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_2.fasta
                URL_HASH SHA256=5ff43c19f3b2d0d7cd113efd5941ff204a5e31820c6edc49c7c47f2d16388e77)
declare_datasource (FILE bin_3.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_3.fasta
                URL_HASH SHA256=0fc4021828129d1752dc2d9a1c163cfb642547ca8e889e969b119e6a2942a39f)
declare_datasource (FILE bin_4.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_4.fasta
                URL_HASH SHA256=d62dbba326f03da7fd4f0d61b460efcc73aee52582fc0903f5358978f33d1a9f)
declare_datasource (FILE bin_5.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_5.fasta
                URL_HASH SHA256=6f16ee171f262e6d1c621b94563adad42a5aab27afee48340d9c2273e67e396f)
declare_datasource (FILE bin_6.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_6.fasta
                URL_HASH SHA256=c4efc778c575e103041b300e0f110677d740c169e55421b9d16650fe92a8c872)
declare_datasource (FILE bin_7.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_7.fasta
                URL_HASH SHA256=45063104427a48892ba4ccb45cc295dc94760f08c1db90e6f3a73744591ada68)
declare_datasource (FILE bin_meta.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_meta.arg
                URL_HASH SHA256=d3440ae40cbff43dee6941cee09ad3c4b9bc860c516c76afd3d2e4ae2319aac3)
declare_datasource (FILE bin_meta.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_meta.bin
                URL_HASH SHA256=303da9e67ae926a711e1ff86950f63aac5e5a57db7d4cb103a2163240f26de08)
declare_datasource (FILE bin_paths.txt
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_paths.txt
                URL_HASH SHA256=614e23263b689c7b4cc0ae41e99aeb5b43b351f865b4604f892320f2cc4377c7)


declare_datasource (FILE 8bins19window0error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window0error100pattern1query_every.gff.out
                URL_HASH SHA256=af90795ca7b27c762c2b9c8d2fa45b130fb528d038020279a1dea690f31bf8cd)
declare_datasource (FILE 8bins19window0error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window0error100pattern50query_every.gff.out
                URL_HASH SHA256=af90795ca7b27c762c2b9c8d2fa45b130fb528d038020279a1dea690f31bf8cd)
declare_datasource (FILE 8bins19window0error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window0error50pattern1query_every.gff.out
                URL_HASH SHA256=cc70823e384e508f353c15e651dea6982f93210849fe795d7f0f4e64978a7f71)
declare_datasource (FILE 8bins19window0error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window0error50pattern50query_every.gff.out
                URL_HASH SHA256=cc70823e384e508f353c15e651dea6982f93210849fe795d7f0f4e64978a7f71)
declare_datasource (FILE 8bins19window1error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error100pattern1query_every.gff.out
                URL_HASH SHA256=5e6059b76815b454926439a477cbb1e48a3e611a75c50bd5421aae258a9888b0)
declare_datasource (FILE 8bins19window1error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error100pattern50query_every.gff.out
                URL_HASH SHA256=7d2ebcdf575567501c959e462c4b2c881b0f501eeb401d673285c3fd29269464)
declare_datasource (FILE 8bins19window1error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error50pattern1query_every.gff.out
                URL_HASH SHA256=e87d9177c94fe52fb156dd09ff7d074dc88672d36e25ef625a08800fa94525a2)
declare_datasource (FILE 8bins19window1error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error50pattern50query_every.gff.out
                URL_HASH SHA256=3b8ab3cf3ddd8bde12610265053b6e388008548a684d3ce9ffe86cc0aded9dd0)
declare_datasource (FILE 8bins23window0error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error100pattern1query_every.gff.out
                URL_HASH SHA256=af90795ca7b27c762c2b9c8d2fa45b130fb528d038020279a1dea690f31bf8cd)
declare_datasource (FILE 8bins23window0error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error100pattern50query_every.gff.out
                URL_HASH SHA256=af90795ca7b27c762c2b9c8d2fa45b130fb528d038020279a1dea690f31bf8cd)
declare_datasource (FILE 8bins23window0error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error50pattern1query_every.gff.out
                URL_HASH SHA256=cc70823e384e508f353c15e651dea6982f93210849fe795d7f0f4e64978a7f71)
declare_datasource (FILE 8bins23window0error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error50pattern50query_every.gff.out
                URL_HASH SHA256=4247fe19cc73d47140462af68ef447892e0b937236f5f8c86797f4aa26438c50)
declare_datasource (FILE 8bins23window1error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error100pattern1query_every.gff.out
                URL_HASH SHA256=921225eae9d8e54963d2423e39ee349158a2e39fc52c8926707749bcda1163a4)
declare_datasource (FILE 8bins23window1error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error100pattern50query_every.gff.out
                URL_HASH SHA256=583d848d9c2e5c082404283f0f4fa04ff55e99a2ee50418c4ced6f7673fbaa9e)
declare_datasource (FILE 8bins23window1error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error50pattern1query_every.gff.out
                URL_HASH SHA256=3b8ab3cf3ddd8bde12610265053b6e388008548a684d3ce9ffe86cc0aded9dd0)
declare_datasource (FILE 8bins23window1error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error50pattern50query_every.gff.out
                URL_HASH SHA256=3b8ab3cf3ddd8bde12610265053b6e388008548a684d3ce9ffe86cc0aded9dd0)
declare_datasource (FILE query.fq
                URL ${CMAKE_SOURCE_DIR}/test/data/search/query.fq
                URL_HASH SHA256=65fbd58c14ca2b4c2274f44fff14dbfce54dc04a89cf6759e1a69cecba933130)


declare_datasource (FILE 16bins50overlap_dream_all.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/16bins50overlap_dream_all.gff
                URL_HASH SHA256=3ecd0e55e704cfd71442bd16805fe587de720e646cbe2d13675ed6249b30045f)
declare_datasource (FILE 16bins50overlap_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/16bins50overlap_reference_metadata.bin
                URL_HASH SHA256=d135d1e85cf27a7b1f256e6fff90d8b112c214abf0969bff3a7d7e05d8b40a40)
declare_datasource (FILE 8bins50overlap_dream_all.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/8bins50overlap_dream_all.gff
                URL_HASH SHA256=f03c93ddb758a7e9ef89e3243dca9fb49b97a1ba3239408ffc30f8ff486982c8)
declare_datasource (FILE 8bins50overlap_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/8bins50overlap_reference_metadata.bin
                URL_HASH SHA256=07a8cb9b040f7dd8aec1e897907ba8cc3263db07452c41dd7e85f36353d9cd64)
declare_datasource (FILE multi_seq_query.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/multi_seq_query.fasta
                URL_HASH SHA256=6d1260dc701802924487453c19715044cbeaa0f42021be8f2996bbaac1f24c58)
declare_datasource (FILE multi_seq_ref.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/multi_seq_ref.fasta
                URL_HASH SHA256=a37f4be29ec99b66efb6ac235224b3145aadf8d5d6ff2c5c2f5324890170ce92)
declare_datasource (FILE stellar_truth.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth.gff
                URL_HASH SHA256=95bf8d1fcbcfde4dfea3cd2055bcf35be161e6b21519e2c936b3595ef37c017b)
declare_datasource (FILE stellar_truth_num12_dis13.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num12_dis13.gff
                URL_HASH SHA256=95bf8d1fcbcfde4dfea3cd2055bcf35be161e6b21519e2c936b3595ef37c017b)
declare_datasource (FILE stellar_truth_num12_dis3.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num12_dis3.gff
                URL_HASH SHA256=3ba650aa865f8190f21704275c9006e31a75c23b8adc27b1572d0f8a370285e6)
declare_datasource (FILE stellar_truth_num12_dis8.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num12_dis8.gff
                URL_HASH SHA256=18203e3ca735c4416fcc2444013640f13326babff0eea4456de604f0e1882e3e)
declare_datasource (FILE stellar_truth_num3_dis13.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num3_dis13.gff
                URL_HASH SHA256=315f57cedc484dcf074c4d8f0c199d73480fe0d9d5a58e5072df28c8e7d1694f)
declare_datasource (FILE stellar_truth_num3_dis3.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num3_dis3.gff
                URL_HASH SHA256=3ba650aa865f8190f21704275c9006e31a75c23b8adc27b1572d0f8a370285e6)
declare_datasource (FILE stellar_truth_num3_dis8.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num3_dis8.gff
                URL_HASH SHA256=f11cb143f32804c9fb51ee5c865fcb08fa1862b2fb67557ba95d17b38a4fb0e5)
declare_datasource (FILE stellar_truth_num9_dis13.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num9_dis13.gff
                URL_HASH SHA256=95bf8d1fcbcfde4dfea3cd2055bcf35be161e6b21519e2c936b3595ef37c017b)
declare_datasource (FILE stellar_truth_num9_dis3.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num9_dis3.gff
                URL_HASH SHA256=3ba650aa865f8190f21704275c9006e31a75c23b8adc27b1572d0f8a370285e6)
declare_datasource (FILE stellar_truth_num9_dis8.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/stellar_truth_num9_dis8.gff
                URL_HASH SHA256=18203e3ca735c4416fcc2444013640f13326babff0eea4456de604f0e1882e3e)


declare_datasource (FILE 1error.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/1error.gff
                URL_HASH SHA256=cc85d4bad448839816a359b43f25b9cedd840d1a9970026d7d2ea876310fa7d7)
declare_datasource (FILE 1error.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/1error.gff.out
                URL_HASH SHA256=f44f0c78c791566e52635a0f26b3356abd4c1a7a8a72243dbed3d364fd7c6f0f)
declare_datasource (FILE 1error_stellar.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/1error_stellar.gff
                URL_HASH SHA256=01993f28b0973e612a7cc3e84abdbe551c9b47a6cc7507106ff4b48071c21613)
declare_datasource (FILE 2error.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error.gff
                URL_HASH SHA256=f2b7196410d5d9dfc2fac9ba5fe0e92c48eac7197913f6ea6f4219e303f43a78)
declare_datasource (FILE 2error.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error.gff.out
                URL_HASH SHA256=9ad3fff98e94f85a2c07b0d673327bb2dce8a3ea82f5e1e642de7e7a70ceebf3)
declare_datasource (FILE 2error_stellar.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error_stellar.gff
                URL_HASH SHA256=85cef348cd514cd12c442c71c7b44a6e44da4d3ea99e56752b3cc88d6f8d5b3f)
declare_datasource (FILE dummy_reads.fastq
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/dummy_reads.fastq
                URL_HASH SHA256=f1aa9ca0fb0b87393923848f0389cc3fb5cfd4841566afaf72e6c55829b64d73)
declare_datasource (FILE meta_param_deduction.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/meta_param_deduction.arg
                URL_HASH SHA256=99630eace6da93c95750a3f060a2f2de847848a54180fcac75b5c11d4dea2946)
declare_datasource (FILE meta_param_deduction.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/meta_param_deduction.bin
                URL_HASH SHA256=9f103c5a09fad5ce02b8066c8106758fa3f27f6322d73f84a559e7f2cb7fc992)
declare_datasource (FILE query.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/query.fasta
                URL_HASH SHA256=40246c2de99c2d41508f3eec7fbfb1006001e6132c4b2fbdbe9ba500ce8be887)
declare_datasource (FILE ref.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/ref.fasta
                URL_HASH SHA256=47f808d207c4c90afebbe1c8ab28990ec0e3e777c75ec787099279005428f3da)
