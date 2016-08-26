#!/bin/bash

TMPDIR=~/.tmp

if [ $# -lt 1 ]; then
    echo -e "\tusage: $0 <FASTA FILE 1> .. <FASTA FILE N>"
    exit 1
fi

INPUT_FILES=$@

OUT_FILE_NAME=""
for f in $INPUT_FILES; do
    if [ ! -f $f ]; then
        echo File $f does not exist. Exiting.
        exit 1
    fi
    fb=$(basename $f)
    if [ -z $OUT_FILE_NAME ]; then
        OUT_FILE_NAME=${fb%.*}
    else
        OUT_FILE_NAME=${OUT_FILE_NAME}_${fb%.*}
    fi
done


PARAMS="--output=$OUT_FILE_NAME.xmfa --output-guide-tree=$OUT_FILE_NAME.tree --backbone-output=$OUT_FILE_NAME.backbone --scratch-path-1=$TMPDIR $INPUT_FILES"

echo progressiveMauve $PARAMS > mauve_params.log

progressiveMauve $PARAMS 2>&1 | tee mauve.log
