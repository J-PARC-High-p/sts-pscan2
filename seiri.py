#!/usr/bin/python3
import json
import os
import shutil
import re  # 正規表現モジュールを追加

# EfuseMap.json のパス
efuse_map_path = './EfuseMap.json'
# ファイルを分類したいディレクトリのパス
<<<<<<< HEAD
directory_a_path = '../data/2023_pfar/pscan/20231121_pscan_Scan_after_Run113'
=======
#directory_a_path ='/home/ryamada/pscan_data/20231120_pscan_104PN_107P_108PN'
>>>>>>> edcde05db87044102d0957169e96222ecec989b3

# EfuseMap.json を読み込む
with open(efuse_map_path, 'r') as file:
    efuse_map = json.load(file)

# Module_ID と識別子のマッピングを作成
id_to_module = {}
for key, value in efuse_map.items():
    id_to_module[key] = value.get('Module_ID', 'unknown')

# ディレクトリAの下のファイルを処理
for filename in os.listdir(directory_a_path):
    # ファイル名から識別子を抽出（ファイル名フォーマットに基づく）
    match = re.search(r'XA-\d{3}-\d{2}-\d{3}-\d{3}-\d{3}-\d{3}-\d{2}', filename)
    if match:
        identifier = match.group(0)
        # 対応するModule_IDを取得
        module_id = id_to_module.get(identifier, 'unknown')
    else:
        module_id = 'unknown'  # 識別子が見つからない場合

    # 対象のModule_IDディレクトリへのパス
    target_dir = os.path.join(directory_a_path, module_id)
    # ディレクトリが存在しなければ作成
    os.makedirs(target_dir, exist_ok=True)
    # ファイルを移動
    shutil.move(os.path.join(directory_a_path, filename), os.path.join(target_dir, filename))

