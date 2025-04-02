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
                URL_HASH SHA256=3c0588d0f9f39bf4fee504b13999696aad3681699f88966ab031d58664828b5a)
declare_datasource (FILE 150overlap16bins15window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins15window1errors.gff.out
                URL_HASH SHA256=e076d1239a1a906265f1d81b6e2566fd8ccfae47c715a41a034288f42f13e7c6)
declare_datasource (FILE 150overlap16bins15window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins15window.ibf
                URL_HASH SHA256=228a243796aa2e2f98cb8e0783c3806d8e75f32a0026c4219a4ff79665436394)
declare_datasource (FILE 150overlap16bins.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins.arg
                URL_HASH SHA256=5abba452e93768fc2654f3ef3b72b74dfffda1061e25c2450e87b39bd5be6a0b)
declare_datasource (FILE 150overlap16bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap16bins.bin
                URL_HASH SHA256=065abd432e34fed34fe074b7fc078700823be510b76a37620ec59a635ea01798)
declare_datasource (FILE 150overlap4bins13window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins13window1errors.gff.out
                URL_HASH SHA256=d45982f54310c2037e3b75da6cfff62179eff043d41e4987a575f0727415c4ea)
declare_datasource (FILE 150overlap4bins13window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins13window.ibf
                URL_HASH SHA256=6122ec0ca7e0c7ee9c5671d85eb9693bd1b2eb2bce09f82d22a8f5865b592920)
declare_datasource (FILE 150overlap4bins15window1errors.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins15window1errors.gff.out
                URL_HASH SHA256=588fd82cfcdd32215f4b1d96d1b8ccc31ed96547e8ca75ea306913b99eef37ee)
declare_datasource (FILE 150overlap4bins15window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins15window.ibf
                URL_HASH SHA256=518131b888ba0e2192aecce2eed79c4cc679dcc2413e792237e6bd564455bb43)
declare_datasource (FILE 150overlap4bins.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins.arg
                URL_HASH SHA256=5abba452e93768fc2654f3ef3b72b74dfffda1061e25c2450e87b39bd5be6a0b)
declare_datasource (FILE 150overlap4bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/single/150overlap4bins.bin
                URL_HASH SHA256=17ab4ce8126c191cf827faf87bb477058226fd00ca3798989ae926a9e0e351a6)


declare_datasource (FILE 20overlap16bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/multi/20overlap16bins.bin
                URL_HASH SHA256=07c233ad32a88fb4475fc797351721546e5ac7ac246579a64c333b5f0d84487b)
declare_datasource (FILE 20overlap4bins.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/multi/20overlap4bins.bin
                URL_HASH SHA256=4a68c6f7c58394d507d86ba3dc8d35f6c4f3947607a118b653108185f18842de)


declare_datasource (FILE write_out_20_16_reference_metadata.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_16/reference_metadata.arg
                URL_HASH SHA256=d05af1f21f5da0837ad321f814dfaa39d6ad2ac4b02099144f0d99ff8875f3c9)
declare_datasource (FILE write_out_20_16_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_16/reference_metadata.bin
                URL_HASH SHA256=f16bd951afbde25ebc7ccceca828d0e98db8efb62e85f3f79c9cbbb7cd26641d)
declare_datasource (FILE write_out_20_4_reference_metadata.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_4/reference_metadata.arg
                URL_HASH SHA256=d05af1f21f5da0837ad321f814dfaa39d6ad2ac4b02099144f0d99ff8875f3c9)
declare_datasource (FILE write_out_20_4_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/split/write_out_20_4/reference_metadata.bin
                URL_HASH SHA256=f16bd951afbde25ebc7ccceca828d0e98db8efb62e85f3f79c9cbbb7cd26641d)
declare_datasource (FILE 8bins19window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/build/8bins19window.ibf
                URL_HASH SHA256=b1acbb08d6199e1f0ae1593fe8847437f2e1a336bbb814409b82d6047a6876fe)
declare_datasource (FILE 8bins23window.ibf
                URL ${CMAKE_SOURCE_DIR}/test/data/build/8bins23window.ibf
                URL_HASH SHA256=8f04312ce36498e6e3768206ba8e7c83c0e3a77242773730f5fb5eddb42fcf1c)
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
                URL_HASH SHA256=2ed9939a1815cd5f85ad1068296506ecad861a69f450b67e90ce3bf136d1448c)
declare_datasource (FILE bin_meta.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/build/bin_meta.bin
                URL_HASH SHA256=f4a3dec17b486d55319a385da189709e6678c1cc03b4224d98dd9cb9c3b247f0)
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
                URL_HASH SHA256=5e6059b76815b454926439a477cbb1e48a3e611a75c50bd5421aae258a9888b0)
declare_datasource (FILE 8bins19window1error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error50pattern1query_every.gff.out
                URL_HASH SHA256=cc29ddfb27b9327c4ccc7249f6b802720ac84af4c7695bee028ed4a404f9cba5)
declare_datasource (FILE 8bins19window1error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins19window1error50pattern50query_every.gff.out
                URL_HASH SHA256=3b8ab3cf3ddd8bde12610265053b6e388008548a684d3ce9ffe86cc0aded9dd0)
declare_datasource (FILE 8bins23window0error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error100pattern1query_every.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 8bins23window0error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error100pattern50query_every.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 8bins23window0error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error50pattern1query_every.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 8bins23window0error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window0error50pattern50query_every.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 8bins23window1error100pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error100pattern1query_every.gff.out
                URL_HASH SHA256=0bcbfdc35e39f3a7649aae50741220e438633c8b0d2f9aeb9e0c397f8d706f0f)
declare_datasource (FILE 8bins23window1error100pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error100pattern50query_every.gff.out
                URL_HASH SHA256=0bcbfdc35e39f3a7649aae50741220e438633c8b0d2f9aeb9e0c397f8d706f0f)
declare_datasource (FILE 8bins23window1error50pattern1query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error50pattern1query_every.gff.out
                URL_HASH SHA256=5824c908e91f53f9c9575191dba9a5491e5a0d5eb46fb3ffdcb13338c4016bc8)
declare_datasource (FILE 8bins23window1error50pattern50query_every.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/search/8bins23window1error50pattern50query_every.gff.out
                URL_HASH SHA256=1a8aa9aaff1ca312473b83c43c940978d0cd8f1765a61199c96d1f6d7e7b52ed)
declare_datasource (FILE query.fq
                URL ${CMAKE_SOURCE_DIR}/test/data/search/query.fq
                URL_HASH SHA256=65fbd58c14ca2b4c2274f44fff14dbfce54dc04a89cf6759e1a69cecba933130)


declare_datasource (FILE 16bins50overlap_dream_all.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/16bins50overlap_dream_all.gff
                URL_HASH SHA256=3ecd0e55e704cfd71442bd16805fe587de720e646cbe2d13675ed6249b30045f)
declare_datasource (FILE 16bins50overlap_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/16bins50overlap_reference_metadata.bin
                URL_HASH SHA256=416170a9f77e366f4d0c2c74f8e8417a0945b3b82cf609e9ceb7c30fc8577a2f)
declare_datasource (FILE 8bins50overlap_dream_all.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/8bins50overlap_dream_all.gff
                URL_HASH SHA256=f03c93ddb758a7e9ef89e3243dca9fb49b97a1ba3239408ffc30f8ff486982c8)
declare_datasource (FILE 8bins50overlap_reference_metadata.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/consolidate/8bins50overlap_reference_metadata.bin
                URL_HASH SHA256=6d6855c8bb964b619dcb0065e3a8c16ec5bde2c4ec99c345d37af8e2d9e79660)
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
                URL_HASH SHA256=320bd247df63f0af37b069914a690e168460666f0c1b6040d3138bf387e2cd4a)
declare_datasource (FILE 1error.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/1error.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 1error_stellar.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/1error_stellar.gff
                URL_HASH SHA256=01993f28b0973e612a7cc3e84abdbe551c9b47a6cc7507106ff4b48071c21613)
declare_datasource (FILE 2error.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error.gff
                URL_HASH SHA256=ae95fcdba6f8082424bb9a554439fd3d3335de6a24f89fb75ca039a9abfe8255)
declare_datasource (FILE 2error.gff.out
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error.gff.out
                URL_HASH SHA256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855)
declare_datasource (FILE 2error_stellar.gff
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/2error_stellar.gff
                URL_HASH SHA256=85cef348cd514cd12c442c71c7b44a6e44da4d3ea99e56752b3cc88d6f8d5b3f)
declare_datasource (FILE dummy_reads.fastq
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/dummy_reads.fastq
                URL_HASH SHA256=f1aa9ca0fb0b87393923848f0389cc3fb5cfd4841566afaf72e6c55829b64d73)
declare_datasource (FILE meta_param_deduction.arg
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/meta_param_deduction.arg
                URL_HASH SHA256=c4fb0c22033073d875aee2aed25a020949aee52ce662e6e595ee629f8808a997)
declare_datasource (FILE meta_param_deduction.bin
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/meta_param_deduction.bin
                URL_HASH SHA256=19d9f4fd273cf03754e3f2dda4b3e767d6c96f8a227f5a5c5231262b18296f13)
declare_datasource (FILE query.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/query.fasta
                URL_HASH SHA256=40246c2de99c2d41508f3eec7fbfb1006001e6132c4b2fbdbe9ba500ce8be887)
declare_datasource (FILE ref.fasta
                URL ${CMAKE_SOURCE_DIR}/test/data/dream/ref.fasta
                URL_HASH SHA256=47f808d207c4c90afebbe1c8ab28990ec0e3e777c75ec787099279005428f3da)
