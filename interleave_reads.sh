#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "No read 1 file specified"
    exit 1
fi
if [ -z "$2" ]; then
    echo "No read 2 file specified"
    exit 1
fi
if [ ! -f "$1" ]; then
    echo "Read 1 file does not exist"
    exit 1
fi
if [ ! -f "$2" ]; then
    echo "Read 2 file does not exist"
    exit 1
fi

{
    # Paired reads
    paste <(zcat -f "$1" | paste - - - -) <(zcat -f "$2" | paste - - - -) \
    | tr '\t' '\n'
    shift
    shift

    # Unpaired reads (if present)
    while test ${#} -gt 0; do
        zcat -f "$1"
        shift
    done

} | gzip -c
