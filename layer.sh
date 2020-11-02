#!/bin/bash

x="`echo "$1+1.7" | bc`d0"

sed -i.bak "31s/.*/$1d0/g" BOXtemp 
sed -i.bak "32s/.*/$x/g" BOXtemp 
sed -i.bak "33s/.*/$2d0/g" BOXtemp 
rm BOXtemp.bak
echo "BOXtemp modified"
