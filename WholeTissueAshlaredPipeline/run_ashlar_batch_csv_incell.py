from __future__ import print_function
import csv
from subprocess import call
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib
import argparse
import os
import datetime

def text_to_bool(text):
    return False \
        if not str(text) \
        else str(text).lower() in '1,yes,y,true,t'
def path_to_date(path):
    return os.path.getmtime(str(path))

parser = argparse.ArgumentParser(
    description='Read a csv config file and run ashlar'
)
parser.add_argument(
    'csv_filepath', metavar='CSVFILE',
    help='a csv file with header row: Directory, Correction, Pyramid'
)
parser.add_argument(
    '-f', '--from-dir', dest='from_dir', default=0, type=int, metavar='FROM',
    help=('starting directory; numbering starts at 0')
)
parser.add_argument(
    '-t', '--to-dir', dest='to_dir', default=None, type=int, metavar='TO',
    help=('ending directory; numbering starts at 0')
)
args = parser.parse_args()

csv_path = pathlib.Path(args.csv_filepath)
if not csv_path.is_file() or csv_path.suffix != '.csv':
    print('A csv file is required to proceed.')
    parser.print_usage()
    parser.exit()

with open(str(csv_path)) as exp_config:
    exp_config_reader = csv.DictReader(exp_config)
    exps = [dict(row) for row in exp_config_reader]


for exp in exps[args.from_dir : args.to_dir]:
    path_exp = pathlib.Path(exp['Directory'])
    files_exp = sorted(path_exp.glob('*/*xdce'))
    files_exp.sort(key=path_to_date)

    if len(files_exp) == 0:
        print('No xdce files found in', str(path_exp))
        continue

   # print('Processing files in', str(path_exp))
    print('Processing files') 
    print(datetime.datetime.now())
    print()
    if text_to_bool(exp['Correction']):
        lambda_flat = '0.1'
        lambda_dark = '0.01'
        ffp_list = []
        dfp_list = []
        for j in files_exp:
            print('Run BaSiC')
            # print('\r    ' + 'Generating ffp and dfp for ' + j.name)
            ffp_file_name = j.name.replace('.xdce', '-ffp-basic.tif')
            dfp_file_name = j.name.replace('.xdce', '-dfp-basic.tif')
            if (j.parent / ffp_file_name).exists() and (j.parent / dfp_file_name).exists():
                # print('\r        ' + ffp_file_name + ' already exists')
                # print('\r        ' + dfp_file_name + ' already exists')
                print('already exists')
            else:
                call("C:\Users\Public\Downloads\Fiji.app\ImageJ-win64.exe --ij2 --headless --run C:\Users\Public\Downloads\Fiji.app\plugins\imagej_basic_ashlar.py \"filename='%s', output_dir='%s', experiment_name='%s', lambda_flat=%s, lambda_dark=%s\"" %(str(j), str(j.parent), j.name.replace('.xdce', ''), lambda_flat, lambda_dark))
                # print('\r        ' + ffp_file_name + ' generated')
                # print('\r        ' + dfp_file_name + ' generated')
                pass
            ffp_list.append(str(j.parent / ffp_file_name))
            dfp_list.append(str(j.parent / dfp_file_name))

    print('Run ashlar')
    xdce_files = ' '.join([str(f) for f in files_exp])
    command = 'ashlar ' + xdce_files + ' -o ' + str(path_exp)

    if text_to_bool(exp['Pyramid']):
        command += ' --pyramid -f ' + path_exp.name + '.ome.tif'
    
    if text_to_bool(exp['Correction']):
        ffps = ' '.join(ffp_list)
        dfps = ' '.join(dfp_list)
        command += ' --ffp ' + ffps + ' --dfp ' + dfps
    
    print(command)
    call(command)