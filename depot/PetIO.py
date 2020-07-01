import sys
import os
import shutil
import gzip
import subprocess

# Check if program is in path

def check_programs(*arg):
    error_list = []
    for program in arg:
        if which(program) is False:
            error_list.append("\t[!] {} not found! Please install and add to it PATH variable".format(program))
    if error_list:
        print("\n".join(error_list))
        sys.exit()

def which(program):
    if shutil.which(program):
        return program
    else:
        return False

# Directory checking

def get_outdir(out_directory, add_dir=""):
    """generates output directory in case it does not exist."""
    if type(out_directory) != str:
        print("\t[!] {} is NOT a directory! Please specify an output directory".format(out_directory))
        sys.exit()
    elif os.path.isfile(out_directory):
        print("\t[!] {} is a File! Please specify an output directory".format(out_directory))
        sys.exit()
    elif not os.path.exists(os.path.join(out_directory, add_dir)):
        os.mkdir(os.path.join(out_directory, add_dir))
        return os.path.abspath(os.path.join(out_directory, add_dir))
    else:
        return os.path.abspath(os.path.join(out_directory, add_dir))

def check_indir(input_dir):
    if not os.path.exists(input_dir):
        print("\t[!] FATAL ERROR: '{}' directory not found".format(input_dir))
        sys.exit()
    elif not os.path.isdir(input_dir):
        print("\t[!] FATAL ERROR: '{}' is not a directory".format(input_dir))
        sys.exit()
    else:
        return input_dir

# Can open gzip files

def open_file(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname,'rt')
    else:
        return open(fname,'r')

# Run FastTree
def run_fasttree(t_list):
    fasttree_params = t_list[0]
    fasttree_threads = t_list[1]
    os.environ['OMP_NUM_THREADS'] = str(fasttree_threads)
    FNULL = open(os.devnull, 'w')
    subprocess.call("FastTreeMP " + fasttree_params, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

# Run iqtree
def run_iqtree(iqtree_params):
    FNULL = open(os.devnull, 'w')
    subprocess.call("iqtree " + iqtree_params,shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

# Fix IQtree result
def fix_iqtree(fname):
    file = open(fname,'r')
    for line in file:
        line = line.rstrip('\n')
        line_splitted = line.split(":")
    good_tree_list = []
    for item in line_splitted:
        item=item[::-1].replace("_","@",1)[::-1]
        good_tree_list.append(item)
    file.close()
    file = open(fname,'w')
    file.write(":".join(good_tree_list)+"\n")
    file.close()

# Progress bar woohoo!

def progress(iteration, steps, max_value, no_limit=False):
    if int(iteration) == max_value:
        if no_limit == True:
            sys.stdout.write('\r')
            print ("[x] \t%d%%" % (100), end='\r')
        else:
            sys.stdout.write('\r')
            print ("[x] \t%d%%" % (100))
    elif int(iteration) % steps == 0:
        sys.stdout.write('\r')
        print ("[x] \t%d%%" % (float(int(iteration) / int(max_value)) * 100), end='\r')
        sys.stdout.flush()
    else:
        pass
