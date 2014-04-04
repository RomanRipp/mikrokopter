#! /bin/bash

suf="_download.zip"
pfd=$(pwd)
dnm=${pfd##*/}

cd ..

for D in $(find . -maxdepth 1 -type d)
do
    name=${D##*/}
 
    if [ -n "$name"  -a  "$name" != "$dnm"  -a  "$name" != "." ]; then
 
       echo "$D --> $pfd/$name$suf"

       zip -r  $pfd/$name$suf $D -x \*_download.zip

    fi

done


cd $pfd
