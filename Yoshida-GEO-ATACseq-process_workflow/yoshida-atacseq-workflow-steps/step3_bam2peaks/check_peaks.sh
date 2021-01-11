#!/bin/bash

base_dir="s3://yoshida-atacseq/peaks/"
all_cell_dirs=( $(aws s3 ls $base_dir | awk '{print $2}') )
echo ${all_cell_dirs[@]}

if [ -f check_output_file.txt ]
then
  rm check_output_file.txt
fi
touch check_output_file.txt

for cd in ${all_cell_dirs[@]}:
do
  echo $cd
  echo $cd >> check_output_file.txt
  aws s3 ls "$base_dir$cd" >> check_output_file.txt
done
