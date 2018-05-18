#!/usr/bin/env python2

import os,time,sys,argparse
from shutil import copyfile

#python code to create a user defined case

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--case_root', dest="case_root", metavar='root', type=str, nargs=1, default=["./"],
  help='root at which to search for models')

parser.add_argument('--case_name', dest="case_name", metavar='name', type=str, nargs=1, default=["sbetr_run"],
  help='nase of the case ')


args = parser.parse_args()
mdir=os.getcwd()
directory=cwd = mdir + '/'+args.case_root[0]+args.case_name[0]

#create case directory
if not os.path.exists(directory):
    os.makedirs(directory)
#copy run files
copyfile(mdir+'/templates/reaction.jar.sbetr.nl',directory+'/reaction.jar.sbetr.nl')
copyfile(mdir+'/templates/reaction.1d.sbetr.nl',directory+'/reaction.1d.sbetr.nl')
