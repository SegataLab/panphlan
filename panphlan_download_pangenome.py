#!/usr/bin/env python

"""
Downloading PanPhlAn pangenome files 
"""

import os, subprocess, sys, time, bz2
import re
import argparse as ap
from urllib.request import urlretrieve, urlcleanup

from misc import info

author__ = 'Leonard Dubois and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.0'
__date__ = '24 April 2020'

DOWNLOAD_URL = "https://www.dropbox.com/s/1gxpwk8ba0rmopp/panphlan_pangenomes_links.tsv?dl=1"


# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-i','--input_name', type=str, default = None,
                    help='')
    p.add_argument('-o', '--output', type = str, default = ".",
                    help='')
    p.add_argument('-v', '--verbose', action='store_true',
                    help='Show progress information')
    return p.parse_args()


# ------------------------------------------------------------------------------
#   DOWNLOAD FILES
# ------------------------------------------------------------------------------
def byte_to_megabyte(byte):
    """Convert byte value to megabyte
    """
    return (byte / 1048576)
    
class ReportHook():

    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """Print download progress message
        """
        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                info("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            info(status)


def download(url, download_file, overwrite=False, verbose=False):
    """Download a file from a url
    """
    urlcleanup()
    if (not os.path.isfile(download_file)) or overwrite:
        try:
            if verbose:
                info('Downloading "{}" to "{}"\n'.format(url, download_file))
            urlretrieve(url, download_file, reporthook=ReportHook().report)
            info('\n')
        except EnvironmentError:
            info('unable to download "{}"'.format(url), exit=True)
    else:
        info('File "{}" already present\n'.format(download_file))


# ------------------------------------------------------------------------------
#   DATABASE MAPPING FILE 
# ------------------------------------------------------------------------------


def find_url(query, verbose):
    if verbose:
        print("Retrieving mapping file...")
    mapping_file = os.path.basename(DOWNLOAD_URL).replace('?dl=1', '')
    download(DOWNLOAD_URL, mapping_file, overwrite=True, verbose=False)    
    IN = open(mapping_file, mode='r')
    
    url = ""
    for line in IN:
        line = line.strip()
        if re.match(query, line):
            filename, url = line.split('\t')
            if verbose:
                print("Pangenome found : {}".format(filename))
            break
    IN.close()
    if url == "":
        print("Pangenome not found.\nEither the species is not available or the name provided is incorrect")
        sys.exit()

    url = url.replace('?dl=0', '?dl=1')          
    return(url, filename)


def extract_pangenome(archive_name, output_path):
    print("Extracting archive")
    process = subprocess.Popen(['tar', 'jxvf', os.path.join(output_path, archive_name), '-C', output_path],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()


# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()

    url, filename = find_url(args.input_name, args.verbose)
    download(url, os.path.join(args.output, filename) )
    extract_pangenome(filename, args.output)
    

if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
