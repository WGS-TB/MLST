import os

'''Run this script to download all samples written in SRR_Acc_List.txt file '''

sample = list()
sampleNumbersTxt = open("SRR_Acc_List.txt", "r")

for samp in sampleNumbersTxt:
    sample.append(samp[:-1])    #Without \n character
    
os.chdir("data")

for samp in sample:
    if not os.path.exists(samp):
        os.mkdir(samp)
    os.chdir(samp)
    print("~~~~~ Downloading sample {} ~~~~~~\n".format(samp))
    os.system("fastq-dump --split-3 {}".format(samp))
    print("~~~~~ Downloaded sample {} ~~~~~\n".format(samp))
    os.chdir("..")
    
print("**** All samples downloaded ****")
os.chdir("..")
    
    