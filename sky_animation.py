### hard coding section
ant1 = 136
ant2 = 140
tGoal = 0.7854 # almost exactly 3 hours

###

import os

os.chdir(r"./skyflux/simulations")
#print(os.getcwd())

try:
    exec(open("generate_wedge.py").read())
except:
    print("Failure to import g_w")

os.chdir(r"../..")

try:
    exec(open("index_test.py").read())
except:
    print("Failure to import i_t")

print("Importing skyflux.")

import skyflux as sf

print("Skyflux imported.")

srcs = []

for src in sf.catalog.srcs:
    if src.alpha == src.alpha:
        srcs.append(src)

print("Full-alpha catalog assembled. Length is " + str(len(srcs)))
        
import numpy as np

srcs = np.array(srcs)

index = len(srcs) - 1

while index >= 0:
    ### todo
    # just show one plot at each plotting step

    label = "computer_" + str(len(srcs) - 1 - index)
    
    print("Starting a new source.")
    auto_wedge(srcs[np.array([index])], label=label, ptitle=label)
    
    print("Printing said source.")
    load_wedge_sim(label + ".pickle", ant1, ant2, tGoal)
    
    if index == len(srcs) - 1:
        continue
        
    previous_label = "computer" + str(len(srcs) - index)
    merge_files(previous_label, label, label, label)
    print("Printing the entire story so far.")
    load_wedge_sim(label + ".pickle", ant1, ant2, tGoal)
    
    index -= 1

