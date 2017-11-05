import os
from subprocess import call
import argparse
import sys

'''Run this script to download all samples written in SRR_Acc_List.txt file '''
ap = argparse.ArgumentParser()
ap.add_argument("-c", "--cmd", required = False, default="fastq-dump", help="Path to fastq-dump. Default assumes fastq-dump is in user's bin file")
args = vars(ap.parse_args())

if args["cmd"] != "fastq-dump":
    cmd_path = os.path.abspath(args["cmd"])
else:
    cmd_path = "fastq-dump"
    
sample = list()
sampleNumbersTxt = open("SRR_Acc_List.txt", "r")

for samp in sampleNumbersTxt:
    sample.append(samp[:-1])    #Without \n character

if not os.path.exists("data"):
    os.mkdir("data")
    
os.chdir("data")

for samp in sample:
    if not os.path.exists(samp):
        os.mkdir(samp)
    os.chdir(samp)
    print("~~~~~ Downloading sample {} ~~~~~~\n".format(samp))
    try:
        call([cmd_path, "--split-3", samp])
    except OSError:
        print("Error in downloading. Please check the path to fastq-dump is correct.")
        sys.exit(-1)
#    os.system("fastq-dump --split-3 {}".format(samp))
    print("~~~~~ Downloaded sample {} ~~~~~\n".format(samp))
    os.chdir("..")  
    
print("**** All samples downloaded ****")
os.chdir("..")
    
    