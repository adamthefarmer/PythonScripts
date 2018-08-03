################################################################################
#
# Title: AppendFileName
# Purpose: Search 'targetDir' and all subdirectories for filenames which contain
#          'rmvStr' and remove that string from the filename.
#
# Author: Adam Farmer
# Date: 6/28/2018
#
################################################################################

# Import necessary libraries.
import os

# Set the target directory here.
targetDir = "C:\\Users\\afarm\\Desktop\\3-D Printer"
# Set the string in the filename that you want to remove.
rmvStr = " (2018_01_11 05_17_57 UTC)"

for root, dirs, filenames in os.walk(targetDir):
    if not filenames:
        print("no files in this directory!")
        continue
    for file in filenames:
        print(file)
        if rmvStr in file:
            print("found a file")
            os.rename(os.path.join(root, file), os.path.join(root,file.replace(rmvStr,"")))
        else:
            continue
