# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:00:50 2016

@author: j.cosfeld
"""

linelength=10
filelength=30
Sea_entry=5.E-1
Sea_file=zeros((10,linelength))
Sea_file[:,:]=Sea_entry

text_file = open("Output.txt", "w")
for i in range(1,filelength+1):
    text_file.write("5.E-1 5.E-1 5.E-1 5.E-1 5.E-1 5.E-1 5.E-1 5.E-1 5.E-1 5.E-1\n")
    
text_file.close()