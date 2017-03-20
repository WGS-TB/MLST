#!/usr/bin/python
import subprocess
import sh
import math
import random
import os
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-g", "--gene", required = True,  help="name of gene")
args = vars(ap.parse_args())
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
	gene = args["gene"]
	w = gene[0].capitalize()
	gene = gene.replace(gene[0],w)
	temp2 = ""
	file_name2 = gene+"_"+str(x)+".fas"
        for i in range(len(r)):
                string = str(variants[i]) 
                string = string.rstrip() #remove the "\n" character that is returned by bash
                string1= string.split(">")
		print string1
                sim_name = string1[1]+"_"+str(x)+"_"
                file_name = sim_name+"reference.fa"
                val = math.ceil(r[i]*197) #compute the proportions 
                ratios.append(int(val)) 
                temp=sh.grep(string,"linear.txt","-w","-A1") #use bash to extract the variant sequence
		temp = temp.rstrip()
                temp = str(temp)
		temp2 = temp2 + temp
		print temp2 
                cmd = "art_illumina -ss HS25 -sam -i "+file_name+" -p -l 76 -c "+str(ratios[i])+" -m 200 -s 30 -o "+sim_name #the ART command to generate the simulated data for a variant
                #write the variant sequence to a text file
		with open(file_name, "w") as text_file:
                        text_file.write(temp)
                print sim_name #testing purposes
                print file_name #testing purposes
                print cmd #testing purposes
                os.system(cmd) #run the ART command
	with open(file_name2, "w") as text_file2:
		text_file2.write(temp2)
        new_cmd = "cat *_1.fq > "+gene+"_"+str(x)+"_1.fa" #append all the first of the pairs together
        new_cmd2 ="cat *_2.fq > "+gene+"_"+str(x)+"_2.fa" #append all the second of the pairs together
	new_cmd3 = "bash /home/elijah/Desktop/SRA_bowtie/Alignments/scripts/temp.sh "+ gene+"_"+str(x) + " " + gene + "_" + str(x)
	print new_cmd3 
        os.system(new_cmd) #run the command
        os.system(new_cmd2) #run the command
	os.system(new_cmd3)
        os.system("rm "+args["gene"]+"*") #remove unneccessary files for the next iteration.
        print ratios #for testing purposes
