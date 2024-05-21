import os
import shutil
from zipfile import ZipFile


# ########################################################################################
def compress_files(localdir="./", dirbase="./", zipname="gaussian_inputs.zip"):

    # Empty path file paths list
    file_paths = []
    prev_dir = os.getcwd()+"/"

    if dirbase[-1] != "/":
        dirbase += "/"

    # Walking through directory and subidrectories
    for root, directories, files in os.walk(prev_dir):
        for filename in files:
            # join the two strings in order to form the full filepath.
            if filename.find("zip") != -1 or filename.find("full_send") != -1:
                continue
            filepath = os.path.join(root, filename)
            file_paths.append(localdir+filepath.replace(prev_dir, ''))

    # writing files to a zipfile
    with ZipFile(zipname, 'w') as zipp:
        # writing each file
        for file in file_paths:
            zipp.write(file)


# ########################################################################################
def delete_folder(folder):

    shutil.rmtree(folder, ignore_errors=False)
