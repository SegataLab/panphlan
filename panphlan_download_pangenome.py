#!/usr/bin/env python

"""
Downloading PanPhlAn pangenome files
"""

import os, subprocess, sys, time, bz2
import hashlib
import re
import argparse as ap
from urllib.request import urlretrieve, urlcleanup

from misc import info

author__ = 'Leonard Dubois and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '3.1'
__date__ = '12 Jan 2021'

DOWNLOAD_URL = "https://www.dropbox.com/s/c6fkhz4g42w4pf2/panphlan_pangenomes_links_md5.tsv?dl=1"


# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser()
    required = p.add_argument_group('required arguments')
    required.add_argument('-i','--input_name', type = str, default = None,
                        help='Name of species to download', required=True)
    required.add_argument('-o', '--output', type = str, default = ".",
                        help='output location', required=True)
    p.add_argument('-v', '--verbose', action='store_true',
                    help='Show progress information')
    p.add_argument('--retry', type = int, default = 5,
                    help='Number of retry in pangenome download. Default is 5')
    p.add_argument('--wait', type = int, default = 30,
                    help='Number of second spend waiting between download retries. Default 60')
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
            info('unable to download "{}"'.format(url), exit = True, exit_value = 1)
    else:
        info('File "{}" already present\n'.format(download_file))


# ------------------------------------------------------------------------------
#   DATABASE MAPPING FILE
# ------------------------------------------------------------------------------


def find_url(query, verbose, output_path):
    if verbose:
        print("Retrieving mapping file...")
    mapping_file = os.path.basename(DOWNLOAD_URL).replace('?dl=1', '')
    download(DOWNLOAD_URL, os.path.join(output_path, mapping_file), overwrite=False, verbose=False)
    IN = open(os.path.join(output_path, mapping_file), mode='r')

    url = ""
    for line in IN:
        line = line.strip()
        if re.match(query, line):
            filename, url, true_md5 = line.split('\t')
            if verbose:
                print("Pangenome found : {}".format(filename))
            break
    IN.close()
    if url == "":
        print("Pangenome not found.\nEither the species is not available or the name provided is incorrect")
        sys.exit(2)

    url = url.replace('?dl=0', '?dl=1')
    return(url, filename, true_md5)


def extract_pangenome(archive_name, output_path, verbose):
    print("Extracting archive")                     
    try:
        proc = subprocess.check_output(['tar', 'jxvf', os.path.join(output_path, archive_name), '-C', output_path], 
                        stderr=subprocess.STDOUT)
        if verbose:
            sys.stdout.write('[I] Archive extracted !\n')
        os.remove(os.path.join(output_path, archive_name))    
    except subprocess.CalledProcessError:
        print("Error in extracting pangenome archive \n")
        sys.exit(3)



# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()

    if not os.path.exists(args.output) :
        os.mkdir(args.output)

    url, filename, true_md5 = find_url(args.input_name, args.verbose, args.output)
    download(url, os.path.join(args.output, filename) )
    current_hash = hashlib.md5(open( os.path.join(args.output, filename),'rb').read()).hexdigest()
    retries_done = 0
    while current_hash != true_md5 and retries_done < args.retry :
        sys.stdout.write('[W] Incorrect MD5 sum. PanPhlAn will try to re-dowload the file...\n')
        retries_done += 1
        time.sleep(args.wait)
        download(url, os.path.join(args.output, filename), overwrite=True )
        current_hash = hashlib.md5(open( os.path.join(args.output, filename),'rb').read()).hexdigest()
        
    if current_hash != true_md5:
        info("Incorrect file integrity after retries", exit=True, exit_value = 3)
    
    sys.stdout.write('[I] File downloaded ! MD5 checked\n')
    extract_pangenome(filename, args.output, args.verbose)
    


if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')
