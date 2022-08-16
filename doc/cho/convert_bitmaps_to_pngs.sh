#!/usr/bin/env bash

for bmp in *.bmp ; do

  echo $bmp
  basename=${bmp%.bmp}
  ## Mac OS X command: sips
  sips -s format png $bmp --out $basename.png
  rm -f $bmp

done
