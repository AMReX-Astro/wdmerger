import sys
import string

inputFileName  = sys.argv[1]
inputFile      = open(inputFileName, 'r')
outputFileName = string.join([inputFileName,".flash"],'')
outputFile     = open(outputFileName, 'w')

# Read in the first line and determine the number of points in the model

firstLine = inputFile.readline()

nPts = firstLine[10:14]

# Write out the column headers and the number of points

outputFile.write("# radius       dens         temp         c12           ne22")
outputFile.write("\n")
outputFile.write(nPts)
outputFile.write("\n")

# Skip the next few lines

for i in range(0,8):
  inputFile.readline()

# Read in the radius, density, temp and carbon content and then write them to file

line = inputFile.readline()

while line:
  outputFile.write("  ")
  outputFile.write('%e' % float(line[6:23]))    # radius
  outputFile.write(" ")
  outputFile.write('%e' % float(line[31:54]))   # density
  outputFile.write(" ")
  outputFile.write('%e' % float(line[59:77]))   # temperature
  outputFile.write(" ")
  outputFile.write('%e' % float(line[139:158])) # C12
  outputFile.write(" ")
  outputFile.write(" 0.00000000e+00")   # Ne22
  outputFile.write("\n")
  line = inputFile.readline()


