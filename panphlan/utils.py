#!/usr/bin/env python3

import os
import sys
import time
import fnmatch

def end_program(total_time):
    print('\n[TERMINATING...] ' + __file__ + ', ' + str(round(total_time / 60.0, 2)) + ' minutes.\n')

def show_interruption_message():
    sys.stderr.flush()
    sys.stderr.write('\r')
    sys.stderr.write(INTERRUPTION_MESSAGE)

def show_error_message(error):
    sys.stderr.write('[E] Execution has encountered an error!\n')
    sys.stderr.write('    ' + str(error) + '\n')

def time_message(start_time, message):
    current_time = time.time()
    print(' [I] ' + message + ' Execution time: ' + str(round(current_time - start_time, 2)) + ' seconds.')
    return current_time

def find(pattern, path):
    '''
    Find all the files in the path whose name matches with the specified pattern
    '''
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                target = os.path.join(root, name)
                result.append(target)
    return result
