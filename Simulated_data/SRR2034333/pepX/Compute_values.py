#!/usr/bin/python
import subprocess
import sh
import math
import random
import os
random.seed(a=1994) #set random seed

for x in range(1,41): 
        k =random.randint(2,7) #generate a random integer k between 2 and 7
	#generate k random fractions that sum up to 1
        r = [random.random() for j in range(k)] 
        s = sum(r)
        r = [ i/s for i in r ]
	#run a bash sub command to sort the variants text file and return the firt k variants
        variants = (sh.head(sh.sort("variants.txt", "-R"), "-n", k))
        variants = list(variants) #convert output from runningCommand to a list
        print variants #for testing purposes
        ratios = [] #list to store the proprotions of reads to generate
        for i in range(len(r)):
                string = str(variants[i]) 
                string = string.rstrip() #remove the "\n" character that is returned by bash
                string1= string.split(">")
                sim_name = string1[1]+"_"+str(x)+"_"
                file_name = sim_name+"reference.fa"
                val = math.ceil(r[i]*197) #compute the proportions 
                ratios.append(int(val)) 
                temp=sh.grep(string,"linear.txt","-w","-A1") #use bash to extract the variant sequence
                temp = str(temp) 
                cmd = "art_illumina -ss HS25 -sam -i "+file_name+" -p -l 76 -c "+str(ratios[i])+" -m 200 -s 30 -o "+sim_name #the ART command to generate the simulated data for a variant
                #write the variant sequence to a text file
		with open(file_name, "w") as text_file:
                        text_file.write(temp)
                print sim_name #testing purposes
                print file_name #testing purposes
                print cmd #testing purposes
                os.system(cmd) #run the ART command 
        new_cmd = "cat *_1.fq > PepX_all_"+str(x)+"_1.fa" #append all the first of the pairs together
        new_cmd2 ="cat *_2.fq > PepX_all_"+str(x)+"_2.fa" #append all the second of the pairs together
        os.system(new_cmd) #run the command
        os.system(new_cmd2) #run the command
        os.system("rm pepX_*") #remove unneccessary files for the next iteration.
        print ratios #for testing purposes
