#!/bin/bash
cd "$(dirname "$0")"
set -Eeuo pipefail

echo -e "cmake_minimum_required (VERSION 3.16)\n" > datasources.cmake

echo -e "include (cmake/app_datasources.cmake)\n" >> datasources.cmake

cd split

for file in *
do
    [[ $file == *.sh ]] && continue # don't add script as datasource
    [[ -d $file ]] && continue      # skip folders
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/split/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

echo -e "\n" >> ../datasources.cmake

cd single

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/split/single/${file}
                URL_HASH SHA256=" >> ../../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../../datasources.cmake
    echo ")" >> ../../datasources.cmake
done

echo -e "\n" >> ../../datasources.cmake

cd ../multi

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/split/multi/${file}
                URL_HASH SHA256=" >> ../../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../../datasources.cmake
    echo ")" >> ../../datasources.cmake
done

echo -e "\n" >> ../../datasources.cmake

cd ../

for relative_path in write_out_*_*/*
do
    [[ -d $file ]] && continue      # skip folders
    filename=${relative_path/"/"/"_"}
    echo -n "declare_datasource (FILE ${filename}
                URL \${CMAKE_SOURCE_DIR}/test/data/split/${relative_path}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $relative_path))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

cd ../build

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    [[ $file == *.sh ]] && continue
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/build/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

echo -e "\n" >> ../datasources.cmake

cd ../search

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    [[ $file == *.sh ]] && continue
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/search/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

echo -e "\n" >> ../datasources.cmake

cd ../consolidate

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    [[ $file == *.sh ]] && continue
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/consolidate/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

echo -e "\n" >> ../datasources.cmake

cd ../dream

for file in *
do
    [[ -d $file ]] && continue      # skip folders
    [[ $file == *.sh ]] && continue
    echo -n "declare_datasource (FILE ${file}
                URL \${CMAKE_SOURCE_DIR}/test/data/dream/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

cd ../prepare

for shape in "11111111" "1001"; do
    file="s$shape.bin"
    echo -n "declare_datasource (FILE s${shape}.bin
                URL \${CMAKE_SOURCE_DIR}/test/data/prepare/${file}
                URL_HASH SHA256=" >> ../datasources.cmake

    sha=($(shasum -a 256 $file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake
done

for b in $(seq 0 7); do
    
    file="ref_$b.fasta"
    echo -n "declare_datasource (FILE ref_bin_${b}.fasta
                URL \${CMAKE_SOURCE_DIR}/test/data/prepare/bins/${file}
                URL_HASH SHA256=" >> ../datasources.cmake
    
    sha=($(shasum -a 256 bins/$file))
    echo -n $sha >> ../datasources.cmake
    echo ")" >> ../datasources.cmake

    for shape in "11111111" "1001"; do
        file="ref.$b.header"
        echo -n "declare_datasource (FILE s${shape}_${file}
                    URL \${CMAKE_SOURCE_DIR}/test/data/prepare/s$shape/${file}
                    URL_HASH SHA256=" >> ../datasources.cmake

        sha=($(shasum -a 256 s$shape/$file))
        echo -n $sha >> ../datasources.cmake
        echo ")" >> ../datasources.cmake

        file="ref.$b.minimiser"
        echo -n "declare_datasource (FILE s${shape}_${file}
                    URL \${CMAKE_SOURCE_DIR}/test/data/prepare/s$shape/${file}
                    URL_HASH SHA256=" >> ../datasources.cmake

        sha=($(shasum -a 256 s$shape/$file))
        echo -n $sha >> ../datasources.cmake
        echo ")" >> ../datasources.cmake
    done
done
