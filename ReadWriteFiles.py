################################################################################
#
# Title: Read and Write Text File
#
# Purpose: Basic example of how to read a text file and write to various file
#          types.
#
#
# Author: Adam Farmer
#
#
################################################################################



    ### OPENING AND READING FROM A TEXT FILE ###

# Get the input file's directory and base name.
inputFileDir = ""
inputBaseFileName = "" # Include file extension.

# Open the file with read permission.
# Don't do this:  inFileID = open(inputFileDir + inputBaseFileName, "r")
with open(inputFileDir + inputBaseFileName, "r") as inFileID:
    for line in inFileID:
        doSomething()


    ### CREATING AN OUTPUT TEXT FILE ###

# Get the output file's directory and base name.
outputFileDir = ""
outputBaseFileName = "" # Include file extension.

# Open the output file with write permission.
with open(outputFileDir + outputBaseFileName, "w") as outFileID:



    ### DON'T FORGET TO CLOSE THE FILES ###

# Note that the input file is closed automatically when opened with 'with open'
outFileID.close()
