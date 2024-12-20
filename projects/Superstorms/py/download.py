#!/usr/bin/env python3
# download.py: Download geomagnetic data from Edinburgh GIN
#
# This script is designed to download geomagnetic data
# from the Edinburgh INTERMAGNET Geomagnetic Information Node.
# To use the script, follow these steps:
#
# 1.) Create a new directory.
# 2.) Move to the new directory and copy the script there.
# 3.) Execute the script by typing 'python3 download.py'.
#
# The script keeps track of its progress through the list of files
# to download. If it fails at any time, you can restart it and it will
# resume at the point where it failed.
#
#
import os
import sys
import shutil
import re
from urllib import request
from urllib.error import URLError
from contextlib import ExitStack


# Configurable parameters - these variables control use of proxy
# servers at the users site - there's also an option to use
# authentication - change them as required (though the
# defaults should work OK)
gin_username = ''
gin_password = ''
proxy_address = ''
n_retries = 4
#
#
#----------------------------------------------------------------------
# upd_op_co:         Increment and serialize the operation count
# Input parameters:  op_count the current operation number
# Returns:           op_count is incremented and returned (and written to disk)
def upd_op_co(op_count):
  op_count = op_count + 1
  with open('counter.dat', 'w') as f:
    f.write('%d' % op_count)
  return op_count

#----------------------------------------------------------------------
# safemd:            Safely create a folder (no error if it already exists)
# Input parameters:  folder the directory to create
#                    op_number the operation number for this call
#                    op_count the current operation number
# Returns:           op_count is incremented and returned (and written to disk)
def safemd (folder, op_number, op_count):
  if op_number >= op_count:
    if op_number == 0:
      print ('Creating directories...')
    try:
      os.makedirs (folder, exist_ok=True)
    except OSError:
      print ('Error: unable to create directory: ' + str(folder))
      sys.exit (1)
    op_count = upd_op_co (op_count)
  return op_count

#----------------------------------------------------------------------
# getfile:           Download a file from a web server
# Input parameters:  url URL to download from
#                    local_file local file to download to
#                    n_retries number of retries to do
#                    op_number the operation number for this call
#                    gin_username the username of the GIN (or empty string)
#                    gin_password the username of the GIN (or empty string)
#                    proxy_address address of proxy server (or empty string)
#                    n_folders the number of folders to create
#                    n_downloads the number of files to download
#                    op_count the current operation number
# Returns:           op_count is incremented and returned (and written to disk)
def getfile (url, local_file, n_retries, op_number, gin_username,
             gin_password, proxy_address, n_folders,
             n_downloads, op_count, rename=False):
  if op_number >= op_count:
    # tell the user what's going on
    percent = ((op_number - n_folders) * 100) / n_downloads
    print ('%d%% - downloading file: %s' % (percent, local_file))
    
    # remove any existing file
    try:
      os.remove (local_file)
    except FileNotFoundError:
      pass
    except OSError:
      print ('Error: unable to remove file: ' + str(local_file))
      sys.exit (1)
    
    # handle authentication and proxy server
    proxy = auth = None
    if len (proxy_address) > 0:
      proxy = request.ProxyHandler({'http': proxy_address, 'https': proxy_address})
    if len (gin_username) > 0:
      pwd_mgr = request.HTTPPasswordMgrWithPriorAuth()
      pwd_mgr.add_password (None,
                            'https://imag-data.bgs.ac.uk/GIN_V1',
                            gin_username,
                            gin_password,
                            is_authenticated=True)
      auth = request.HTTPBasicAuthHandler(pwd_mgr)
    if url.startswith ('https'):
      default_handler = request.HTTPSHandler
    else:
      default_handler = request.HTTPHandler
    if auth and proxy:
      opener = request.build_opener(proxy, auth, default_handler)
    elif auth:
      opener = request.build_opener(auth, default_handler)
    elif proxy:
      opener = request.build_opener(proxy, default_handler)
    else:
      opener = request.build_opener(default_handler)
    
    # download the file
    success = False
    while (not success) and (n_retries > 0):
      try:
        with opener.open (url) as f_in:
          with open (local_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out, 4096)
        success = True
      except (URLError, IOError, OSError):
        n_retries -= 1
    if not success:
      print ('Error: cannot download ' + local_file)
      sys.exit (1)
    
    # rename IAGA-2002 files
    if rename:
      dt = None
      try:
        with open(local_file, 'r') as f:
          for line in f.readlines():
            if re.search ('^ Data Type', line):
              dt = line[24:25].lower()
      except (IOError, OSError):
        pass
      if dt:
        if not dt.isalpha():
          dt = None
      if dt:
        new_local_file = local_file[:len(local_file) - 7] + dt + local_file[len(local_file) - 7:]
        try:
          os.remove (new_local_file)
        except (FileNotFoundError, OSError):
          pass
        try:
          os.rename (local_file, new_local_file)
        except (IOError, OSError):
          print ('Warning: unable to rename ' + local_file + ' to ' + new_local_file)
      else:
        print ('Warning: unable to determine data type for renaming of ' + local_file)
      
      op_count = upd_op_co (op_count)
  return op_count

if __name__ == '__main__':
  # are we restarting a failed operation or starting from new
  try:
    with open ('counter.dat') as f:
      op_count = int(f.read())
      print ('Information: resuming download after previous failure')
  except (IOError, ValueError, OSError):
    op_count = 0
  n_folders = 1
  n_downloads = 2

  folder = os.path.join('dataset', 'May2024')
  op_count = safemd (folder, 0, op_count)
  stn = "frd"
  local_file = os.path.join(folder, f'{stn}20240510min.min')
  op_count = getfile (
    f'https://imag-data.bgs.ac.uk/GIN_V1/GINServices?Request=GetData&format=IAGA2002&testObsys=0&observatoryIagaCode={stn.upper()}&samplesPerDay=1440&orientation=Native&publicationState=adj-or-rep&recordTermination=UNIX&dataStartDate=2024-05-10&dataDuration=2', 
    local_file, n_retries, 1, gin_username, gin_password, proxy_address, n_folders, n_downloads, op_count
  )

  # tidy up
  print ('100% - data download complete')
  os.remove ('counter.dat')

