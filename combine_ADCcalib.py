#!/usr/bin/python3
import os

base_dir = input("combine file path(til/): ")
output_file = base_dir+"combined_ADCcalib.dat"

combined_data = []

for root, dirs, files in os.walk(base_dir):
    for dir_name in dirs:
        file_path = os.path.join(root, dir_name, 'ADCcalib.dat')
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                combined_data.extend(file.readlines())

# 統合されたデータを新しいファイルに書き込む
with open(output_file, 'w') as output:
    output.writelines(combined_data)

print(f"Completed!: {output_file}")