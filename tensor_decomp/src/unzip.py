"""## unzipping and untar-ing data files"""

data_path = '/home/jjzhong/tensor_decomp/data'

import glob
import os
import re
import tarfile
import zipfile

with zipfile.ZipFile(os.path.join(data_path, 'harris_2023.zip'), 'r') as zip_ref:
    zip_ref.extractall(path=os.path.join(data_path,'harris_2023'))

for file in glob.glob(os.path.join(data_path,'harris_2023/*.tar')):
  match = re.search(r'^/(.+/)*(.+)\.(.+)$', file)
  new_dir = match.group(2)

  with tarfile.open(file) as tar:
    tar.extractall(path=os.path.join(data_path, 'raw_data', new_dir), filter='fully_trusted')

print('files unzipped!')