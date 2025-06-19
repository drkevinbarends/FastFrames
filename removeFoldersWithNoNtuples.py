import os
import shutil

input_dir = '/eos/user/k/kebarend/tWZ/FastFrames/input/'

for subfolder in os.listdir(input_dir):
    subfolder_path = os.path.join(input_dir, subfolder)
    if os.path.isdir(subfolder_path):
        # Check if the folder is empty (no files or subfolders)
        if not os.listdir(subfolder_path):
            print(f"Removing empty folder: {subfolder_path}")
            shutil.rmtree(subfolder_path)
