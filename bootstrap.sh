#!/bin/bash
#autoheader
aclocal -I m4 --install

unamestr=`uname`

if [[ "$unamestr" == 'Linux' ]];
    then
        libtoolize --automake

elif [[ "$unamestr" == 'Darwin' ]];
    then
        glibtoolize --automake
fi

automake --add-missing
autoconf
