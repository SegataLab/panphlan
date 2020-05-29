#!/usr/bin/env python


import os, subprocess, sys, time, bz2
from random import randint


"""Check if bowtie2 is installed. Stops programm if not"""
def check_bowtie2():

    platform = sys.platform.lower()[0:3]
    try:
        if platform == 'win':
            bowtie2 = subprocess.Popen(['where', 'bowtie2'], stdout=subprocess.PIPE)
        else: # Linux, Mac, ...
            bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE)
        bowtie2 = bowtie2.communicate()[0].decode('utf-8')
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE)
        bowtie2_version = bowtie2_version.communicate()[0].decode('utf-8').split()[2]
        print('[I] Bowtie2 is installed')
        print('    version: ' + str(bowtie2_version) + ', path: ' + str(bowtie2).strip())

    except OSError as err:
        sys.stderr.write('\n[E] Execution has encountered an error!\n')
        sys.stderr.write('    ' + str(err) + '\n')
        print('\n[E] Please, install Bowtie2.\n')
        print('    Bowtie2 is used to generate the .bt2 index files used for mapping\n')
        sys.exit()


"""Get a random color for plotting coverage curves"""
def random_color(used):
    reset = False
    total = ['#ff0000', '#800000', '#ffff00', '#808000', '#00ff00',
            '#008000', '#00ffff', '#008080', '#008080', '#0000ff',
            '#000080', '#ff00ff', '#800080', '#fa8072', '#ffa07a',
            '#dc143c', '#b22222', '#8b0000', '#ff69b4', '#ff1493',
            '#c71585', '#ff7f50', '#ff4500', '#ffa500', '#ffd700',
            '#bdb76b', '#9400d3', '#4b0082', '#483d8b', '#6a5acd',
            '#7fff00', '#32cd32', '#00fa9a', '#2e8b57', '#006400',
            '#20b2aa', '#4682b4', '#4169e1', '#ffdead', '#f4a460',
            '#d2691e', '#a52a2a', '#a0522d', '#b8860b', '#000000']
    available = [c for c in total if c not in used]
    # If we have no other available colors, than repeat the picking
    if len(available) == 0:
        available = total
        reset = True
    return (available[randint(0, len(available) - 1)], reset)


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')
    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()
    if exit:
        sys.exit(exit_value)