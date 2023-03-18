#!/usr/bin/env python3
import os
import sys
import random
from collections import Counter
from math import exp
import argparse
from collections import defaultdict
from operator import itemgetter
import matplotlib.pyplot as plt
from skbio.diversity import alpha_diversity
from statistics import mean
from statistics import stdev
##################################################
PARSER = argparse.ArgumentParser(
    prog="Staph_rarefaction_curve.py",
    description="Calculating shanon, simpson and Hill Numbers (q=1,2)",)
PARSER.add_argument("input", type=str, help="ST Frequency")
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
##################################################
QUERY = ARGS.input
QUERY = open(ARGS.input,'r')
data = []
ids = []
e = 2.718281
for line in QUERY:
    line_list = line.rstrip().split(',')
    data.append(line_list[1:])
    ids.append(line_list[0])
##############with replacement#############
#strain_names_list = list(range(100))
op_file = open('random_samples_with_replacement.csv','w')
op_file2 = open('diversity_metrics_with_replacement.csv','w')
op_file2.write('Period,Sample_size,Shanon Index Mean,Stdev,N,Simpson Index Mean,Stdev,N,Hill numbers (q=1) Mean,Stdev,N,Hill numbers (q=2) Mean,Stdev,N\n')
for (j, b) in zip(ids, data):
    for i in [10,20,30,40,50,60,70,80,90,100]:
        #try:
        #strains_random_list1 = list(random.sample(j, k=i) for x in range(100)) #without replacement
        #print(strains_random_list1)
        strains_choice_list1 = list(random.choices(b, k=i) for x in range(100)) #with replacement
        div_sh = []
        div_si = []
        div_h1 = []
        div_h2 = []
        for choices_lst in strains_choice_list1:
            op_file.write('{},{}\n'.format(j,','.join(choices_lst)))
            CCs_counter = Counter(choices_lst)
            CCs_count = [CC[1] for CC in CCs_counter.most_common()]
            shannon_CC = alpha_diversity('shannon',CCs_count,base=e)
            simpson_CC = alpha_diversity('simpson',CCs_count)
            Hill_1 = exp(shannon_CC[0])
            Hill_2 = 1/(1-simpson_CC[0])
            div_sh.append(shannon_CC[0])
            div_si.append(simpson_CC[0])
            div_h1.append(Hill_1)
            div_h2.append(Hill_2)
        op_lst = [mean(div_sh),stdev(div_sh),100,mean(div_si),stdev(div_si),100,
                        mean(div_h1),stdev(div_h1),100,mean(div_h2),stdev(div_h2),100]
        op_lst2 = [str(i) for i in op_lst]
            #print(div_si)
            #print(div_h1)
            #print(div_h2)
        #print(len(div_sh))
        op_file2.write('{},{},{}\n'.format(j,i,','.join(op_lst2)))
        #print(len(strains_random_list1))
        #print(len(strains_choice_list1))
        #except:
        #    continue

#print(strains_random_list1[0])
#print(strains_random_list1[1])
#print(strains_choice_list1[0])
#print(list(strains_random_list1))
##############without replacement#############
op_file = open('random_samples_without_replacement.csv','w')
op_file2 = open('diversity_metrics_without_replacement.csv','w')
op_file2.write('Period,Sample_size,Shanon Index Mean,Stdev,N,Simpson Index Mean,Stdev,N,Hill numbers (q=1) Mean,Stdev,N,Hill numbers (q=2) Mean,Stdev,N\n')
for (j, b) in zip(ids, data):
    for i in [10,20,30,40,50,60,70,80,90,100]:
        try:
            strains_random_list1 = list(random.sample(b, k=i) for x in range(100)) #without replacement
            #print(strains_random_list1)
            #strains_choice_list1 = list(random.choices(b, k=i) for x in range(100)) #with replacement
            div_sh = []
            div_si = []
            div_h1 = []
            div_h2 = []
            for random_lst in strains_random_list1:
                op_file.write('{},{}\n'.format(j,','.join(random_lst)))
                CCs_counter = Counter(random_lst)
                CCs_count = [CC[1] for CC in CCs_counter.most_common()]
                shannon_CC = alpha_diversity('shannon',CCs_count,base=e)
                simpson_CC = alpha_diversity('simpson',CCs_count)
                Hill_1 = exp(shannon_CC[0])
                Hill_2 = 1/(1-simpson_CC[0])
                div_sh.append(shannon_CC[0])
                div_si.append(simpson_CC[0])
                div_h1.append(Hill_1)
                div_h2.append(Hill_2)
            op_lst = [mean(div_sh),stdev(div_sh),100,mean(div_si),stdev(div_si),100,
                            mean(div_h1),stdev(div_h1),100,mean(div_h2),stdev(div_h2),100]
            op_lst2 = [str(i) for i in op_lst]
                #print(div_si)
                #print(div_h1)
                #print(div_h2)
            #print(len(div_sh))
            op_file2.write('{},{},{}\n'.format(j,i,','.join(op_lst2)))
        except:
            op_lst2 = ['NA']*12
            op_file2.write('{},{},{}\n'.format(j,i,','.join(op_lst2)))
####################
