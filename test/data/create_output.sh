#!/bin/bash
cd "$(dirname "$0")"
set -Eeuo pipefail

if ! which valik 2> /dev/null; then
    echo "valik not found. please compile it and add its path to the PATH variable"
    exit 255
fi

export VALIK_TMP=tmp/valik/my_dir
export VALIK_STELLAR=echo
export VALIK_MERGE=echo

execs=(valik stellar)
for exec in "${execs[@]}"; do
    if ! which ${exec} &>/dev/null; then
        echo "${exec} is not available"
        echo ""
        echo "make sure \"${execs[@]}\" are reachable via the \${PATH} variable"
        echo ""

        # trying to do some guessing here:
        paths=(../../build/bin)
        paths+=(../../../stellar3/build/bin)

        p=""
        for pp in ${paths[@]}; do
            p=${p}$(realpath -m $pp):
        done
        echo "you could try "
        echo "export PATH=${p}\${PATH}"

        exit 127
    fi
done

echo "### Running valik split ###"
./split/cli_test_output.sh

echo "### Running valik build ###"
./build/cli_test_output.sh

echo "### Running valik search ###"
./search/cli_test_output.sh

export VALIK_STELLAR=stellar
export VALIK_MERGE=cat

echo "### Running distributed DREAM-Stellar ###"
./dream/cli_test_output.sh

echo "### Running sequential Stellar ###"
./consolidate/api_test_output.sh

echo "### Finished ###"
