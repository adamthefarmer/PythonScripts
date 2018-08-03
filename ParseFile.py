################################################################################
#
# Title: ParseFile
# Purpose: Open a .txt file. Parse through the file line by line searching for a
#          a target string 'targetStr'.
#
# Author: Adam Farmer
# Date: 6/13/2018
#
################################################################################

# Import necessary libraries. Not needed currently.


# Set the target string, target directory and target file.
targetStr = "ref_mrp"
targetDir = "C:\\Users\\afarm\\Desktop\\Hard_New_Code_ADCSverif\\Generated ECEF Pos Vel and Corresponding MRPs from STK\\3hr Astrogator vs CDH\HW ITL\\Different Orbit\\"
targetFileBase = "ACS_sunsoak.txt"

# Set the output file name and path.
outputDir = "C:\\Users\\afarm\\Desktop\\Hard_New_Code_ADCSverif\\Generated ECEF Pos Vel and Corresponding MRPs from STK\\3hr Astrogator vs CDH\HW ITL\\Different Orbit\\"
outputFileBase = "mrp_sun.txt"

# Open the output file.
OFID = open(outputDir + outputFileBase,'w')

# Open the file.
counter = 1
tempLine = [" "," "," "]
with open(targetDir + targetFileBase, 'r') as IFID:
    for line in IFID:
        if targetStr in line:
            if counter == 1:
                tempLine[0] = line[44:60]
                counter = counter + 1
            elif counter == 2:
                tempLine[1] = line[44:60]
                counter = counter + 1
            else:
                tempLine[2] = line[44:60]
                lineToFile = tempLine[0] + "," + tempLine[1] + "," + tempLine[2]
                print(lineToFile)
                OFID.write(lineToFile + '\n')
                counter = 1
        
    ### DON'T FORGET TO CLOSE THE FILES ###

# Note that the input file is closed automatically when opened with 'with open'
OFID.close()



