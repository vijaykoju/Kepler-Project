#!/usr/bin/env python
# File:            unfits_0.9.py
# Author:          Chad Williamson
# Created:         28 March 2011
# Last modified:   6 April 2011
#
#      unfits is a program for converting between file formats frequently
# used in astronomy research, in particular FITS files. It is intended to
# be run from the command line in a UNIX-like environment such as Linux,
# Mac OS X, Solaris, or FreeBSD. Run this file with the -h option or see
# the display_help() function below for further details.
#      You may use, adapt, and redistribute this software subject to the
# terms of the Creative Commons Attribution-ShareAlike 2.5 Generic License.
#   --> See http://creativecommons.org/licenses/by-sa/2.5/ for details.

import sys, os, pyfits, png

# Return the substring of file_name resulting from removing the last
# period and everything that follows it. For example, trunk('foo.bar')
# returns 'foo'.
def trunk(file_name):
  for i in range(1, len(file_name)):
    if '.' in file_name[len(file_name) - i:len(file_name)]:
      return file_name[0:len(file_name) - i]
  return file_name

# Return the substring of file_name following the last period. For
# example, extension('foo.bar') returns 'bar' (forces lowercase).
def extension(file_name):
  for i in range(1, len(file_name)):
    if '.' in file_name[len(file_name) - i:len(file_name)]:
      return file_name[len(file_name) - i + 1:len(file_name)].lower()
  return ''

# Return the file name (everything after the last slash) for a path
def name_only(path):
  for i in range(len(path)):
    if '/' in path[len(path) - i:len(path)]:
      return path[len(path) - i + 1:len(path)]
  return path

# Check a file path. If OK it returns True, otherwise it prints an
# appropriate error message and returns False.
def good_path(path):
  if not os.access(path, os.F_OK):
    print "Error: \'" + path + "\' is not a valid filesystem path."
    return False
  if os.path.isdir(path):
    print "Error: \'" + path + "\' is a directory."
    return False
  return True

# Determine whether a given string is one of the valid output types
def valid_type(output_type):
  if output_type in ['csv', 'txt', 'dat', 'png']:
    return True
  return False

# Return the type of file at path
def classify(path):
  if not good_path(path):
    return 'bad'
  elif extension(path) in ['fits', 'fts', 'fit']:
    return 'fits'
  elif extension(path) == 'csv':
    return 'csv'
  elif extension(path) == 'txt':
    return 'txt'
  elif extension(path) == 'dat':
    return 'dat'
  elif extension(path) == 'png':
    return 'png'
  elif extension(path) == '':
    print "Error: Unrecognized or unsupported file type."
  else:
    print "Error: Unrecognized or unsupported file with extension: \'",
    print "\b." + extension(path) + "\'"
  return 'unrec'

# Return an array containing data from the given FITS file HDU
def array_from_fits(path, hdu):
  hdulist = pyfits.open(path)
  if hdu == 'default':
    if len(hdulist) > 1:
      hdu = 1
    else:
      hdu = 0
  array = []
  for row in hdulist[hdu].data:
    line = []
    for entry in row:
      line.append(entry)
    array.append(line)
  hdulist.close()
  return array

# Return an array containing data from the given CSV or TXT file
def array_from_textfile(path):
  textfile = open(path, 'r')
  linelist = textfile.readlines()
  textfile.close()
  array = []
  for line in linelist:
    row = []
    entry = ''
    if line[0] == ' ':
      gap = True
    else:
      gap = False
    for char in line:
      if char in '0123456789.Ee-infNa':
        entry += char
        gap = False
      elif not gap:
        row.append(entry)
        entry = ''
        gap = True
      else:
        continue
    array.append(row)
  return array

# Write an array of values to a plaintext table delimited by the given delimiter
def array_to_txt(array, path, delimiter):
  output = []
  for row in array:
    line = ''
    for entry in row:
      line += str(entry) + delimiter
    output.append(line[0:len(line)-len(delimiter)] + '\n')
  txt_file = open(path, 'w')
  txt_file.writelines(output)
  txt_file.close()

# Write an array of values to a PNG image at the specified path
def array_to_png(array, path, floor, ceiling):
  using_min, using_max = False, False
  if floor == 'min':
    using_min = True
    floor = float(array[0][0])
    for row in array:
      for entry in row:
        if float(entry) < floor:
          floor = float(entry)
  if ceiling == 'max':
    using_max = True
    ceiling = float(array[0][0])
    for row in array:
      for entry in row:
        if float(entry) > ceiling:
          ceiling = float(entry)
  if not (using_min and using_max) and ceiling <= floor:
    print "failed."
    if using_max and not using_min:
      print "  Error: The floor you specified is at or above the ceiling of", ceiling, '\b.'
      return 1
    elif using_min and not using_max:
      print "  Error: The ceiling you specified is at or below the floor of", floor, '\b.'
      return 2
    else:
      print "  Error: The ceiling must be greater than the floor."
      print "Aborting further operations..."
      exit(1) # in this case all files will fail so we shouldn't keep going
  image = []
  for row in array:
    line = []
    for entry in row:
      value = float(entry)
      if value > ceiling:
        value = ceiling
      if value < floor:
        value = floor
      value -= floor
      value *= 255./(ceiling - floor)
      line.append(int(value))
    image.append(line)
  image_file = open(path, 'wb')
  image_writer = png.Writer(len(image[0]), len(image), greyscale=True)
  image_writer.write(image_file, image)
  image_file.close()
  return 0

# Process a file
def process(path, output_type, hdu, floor, ceiling, delimiter):
  if output_type == 'dat':
    if delimiter == 'unspecified':
      print "Error: To use the DAT output type you must specify a delimiter with the"
      print "       \'-d\' option (and remember to put your delimiter in quotes)."
      return 1
    if delimiter == '':
      print "WARNING: Empty delimiter string. Output columns will not be separated"
      print "         by anything!"
  
  input_type = classify(path)
  print "Loading data from \'" + name_only(path) + "\' into RAM...",
  if input_type in ['bad', 'unrec', 'png']:
    print "failed."
    return 1
  elif input_type == 'fits':
    array = array_from_fits(path, hdu)
  elif input_type in ['csv','txt','dat']:
    array = array_from_textfile(path)
  print "done."
  
  path = trunk(path) + '.' + output_type
  print "Writing data to \'" + name_only(path) + "\'...",
  if output_type == 'csv':
    array_to_txt(array, path, ', ')
  elif output_type == 'txt':
    array_to_txt(array, path, ' ')
  elif output_type == 'dat':
    array_to_txt(array, path, delimiter)
  elif output_type == 'png':
    if array_to_png(array, path, floor, ceiling) > 0:
      return 1
  print "done."
  return 0

# Print details about a file
def info_gen(path):
  file_type = classify(path)
  if file_type in ['bad', 'unrec']:
    return 1
  print "Information for \'" + path + "\':"
  if file_type == 'fits':
    print "  Type:                FITS (NASA's [F]lexible [I]mage [T]ransport [S]ystem)"
    hdulist = pyfits.open(path)
    print "  Header Data Units:   " + str(len(hdulist))
    for i in range(len(hdulist)):
      print "  HDU", i, "\b:"
      if type(hdulist[i]) == pyfits.core.PrimaryHDU:
        print "    Type:              Primary HDU"
      elif type(hdulist[i]) == pyfits.core.BinTableHDU:
        print "    Type:              Binary Table"
        print "    Rows:              " + str(len(hdulist[i].data))
        print "    Columns:           " + str(len(hdulist[i].data[0]))
      elif type(hdulist[i]) == pyfits.core.ImageHDU:
        print "    Type:              Image"
        print "    Width:             " + str(len(hdulist[i].data[0]))
        print "    Height:            " + str(len(hdulist[i].data))
      else:
        print "    Type:              Unknown/Unsupported"
    hdulist.close()
  elif file_type == 'png':
    print "  Type:        PNG ([P]ortable [N]etwork [G]raphics image)"
  elif file_type == 'csv':
    print "  Type:        CSV ([C]omma [S]eparated [V]alue plaintext spreadsheet)"
    test_array = array_from_textfile(path)
    print "  Rows:        " + str(len(test_array))
    print "  Columns:     " + str(len(test_array[0]))
  elif file_type in ['txt','dat']:
    print "  Type:        plaintext (assuming a tab/space delimited table)"
    test_array = array_from_textfile(path, True)
    print "  Rows:        " + str(len(test_array))
    print "  Columns:     " + str(len(test_array[0]))
  return 0
  
# Display detailed help
def display_help():
  print "USAGE: unfits [OPTIONS] [FILES]"
  print "   Ex: unfits -i vega_123.fits"
  print "       > Print information about the structure of the file \'vega_123.fits\'."
  print "   Ex: unfits -o csv -u 2 /my/data/*.fits"
  print "       > Take any files ending in \'.fits\' in the directory \'/my/data/\' and"
  print "       > write HDU 2 of that fits file to a CSV spreadsheet."
  print "   Ex: unfits -o dat -d \"; \" foo.txt"
  print "       > Write the data in foo.txt to a dat file where columns are delimited"
  print "       > by semicolons and some whitespace."
  print "Available options:"
  print "   -h              Print this page."
  print "   -i              Print information about files that follow."
  print "   -o [type]       Place output in file of type [type]. Available types are"
  print "                   csv, txt, dat, and png (case insensitive)."
  print "   -d [delimiter]  Specify a custom column delimiter string (in quotes) for"
  print "                   use with the output type dat (which requires this)."
  print "   -u [number]     Specify which HDU to use (for FITS files only). Numbering"
  print "                   starts at 0. Default is 1 unless only 0 is available."
  print "   -f [value]      Specify an input floor for pixel values (for PNG output"
  print "                   only). Any data value at or below this value will be mapped"
  print "                   to a pixel value of 0 in the resulting PNG image. Default"
  print "                   is the minimum of all pixel values."
  print "   -c [value]      Specify an input ceiling for pixel values (for PNG output"
  print "                   only). Any data value at or above this value will be mapped"
  print "                   to a pixel value of 255 in the resulting PNG image. Default"
  print "                   is the maximum of all pixel values."
  print "                   > Values between the floor and ceiling are scaled linearly"
  print "                   > to fit those values at the endpoints."
  print "Supported input filetypes:"
  print "   FITS            Flexible Image Transport System. Must end in \'.fts\',"
  print "                   \'.fit\', or \'.fits\'."
  print "   CSV             Comma Separated Value spreadsheet. Must end in \'.csv\'"
  print "   TXT             Space or tab delimited textual table. Must end in \'.txt\'"
  print "   DAT             Generic textual table delimited by spaces, tabs, semicolons,"
  print "                   colons, absolute value bars, commas, or anything else that is"
  print "                   not a number, letter, period or hyphen. Must end in \'.dat\'."

# Display usage information
def display_usage():
  print "USAGE: unfits [OPTIONS] [FILES]"
  print "Type \'unfits -h\' for detailed help."

# Initialize flags for parsing arguments...
title_flag = True
info_flag = False
output_flag = False
dat_flag = False
header_flag = False
floor_flag = False
ceiling_flag = False
empty = True

# Initialize our other variables...
output_type = 'none'
delimiter = 'unspecified'
header = 'default'
floor = 'min'
ceiling = 'max'

# Parse arguments...
for argument in sys.argv:
  if title_flag:
    title_flag = False
    continue
  if argument == '-h':
    empty = False
    display_help()
    break
  if argument == '-i':
    info_flag = True
    continue
  if info_flag:
    empty = False
    info_gen(argument)
    continue
  if argument == '-o':
    output_flag = True
    continue
  if output_flag:
    output_type = argument.lower()
    if valid_type(output_type):
      output_flag = False
      continue
    else:
      print "Error: \'" + output_type + "\' is not a valid output type."
      print "Supported output types are: txt, csv, and png"
      print "No action taken due to errors."
      empty = False
      break
  if argument == '-d':
    dat_flag = True
    continue
  if dat_flag:
    delimiter = argument
    dat_flag = False
    continue
  if argument == '-u':
    header_flag = True
    continue
  if header_flag:
    header = int(argument)
    header_flag = False
    continue
  if argument == '-f':
    floor_flag = True
    continue
  if floor_flag:
    floor = float(argument)
    floor_flag = False
    continue
  if argument == '-c':
    ceiling_flag = True
    continue
  if ceiling_flag:
    ceiling = float(argument)
    ceiling_flag = False
    continue
  if output_type == 'none':
    break
  empty = False
  process(argument, output_type, header, floor, ceiling, delimiter)

if empty:
  display_usage()

