


inputDir="/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/Top_300_From_Each_Group"
rootdir = inputDir

import os

from pathlib import Path

for path in Path(rootdir).iterdir():
    if path.is_dir():
        print(path)
