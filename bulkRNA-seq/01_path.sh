#!/bin/bash

tar=/data/pingjiale/05_Result/01_qz_pbmc/01_bulkRNA

cd $tar
mkdir -p 02_trim 03_map 04_count
for path in `ls $tar`
do
mkdir -p  $path/scripts
mkdir -p $path/logs
done
