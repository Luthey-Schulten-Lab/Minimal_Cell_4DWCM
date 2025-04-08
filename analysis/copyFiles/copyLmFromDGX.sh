#!/bin/bash

# Loop over indices from 1 to 10
index=1
while [ $index -lt 12 ]; do
  dir="Mar21_$index"
  echo $dir
  mkdir $dir
  sshpass -p 'yourPasswordHere' scp "username@zandgx.scs.illinois.edu:/raid/andrew/4DWCM/Data/Mar21/$dir/*.lm" "$dir"
  echo $((index+=1))
done
