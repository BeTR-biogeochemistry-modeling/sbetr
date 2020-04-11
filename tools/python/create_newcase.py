#!/usr/bin/env python

import os,time,sys,argparse
from shutil import copyfile
import re
#python code to create a user defined case

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--case_root', dest="case_root", metavar='root', type=str, nargs=1, default=["./"],
  help='root at which to search for models')

parser.add_argument('--case_name', dest="case_name", metavar='name', type=str, nargs=1, default=["sbetr_run"],
  help='nase of the case ')


args = parser.parse_args()
#obtain the root directory of betr

dir_path = os.path.dirname(os.path.realpath(__file__))
print (dir_path)

d=dir_path.split('/sbetr')
mdir=''
nn=len(d)-1
for j in range(nn):
    if j==0:
        mdir=d[j]
    else:
        mdir=mdir+'/sbetr'+d[j]
mdir=mdir+'/'
if not re.search('/sbetr/',mdir):
    mdir=mdir+'/sbetr/'
if args.case_root[0][0]=='/':
    directory=args.case_root[0]+args.case_name[0]
else:
    directory= mdir + '/'+args.case_root[0]+args.case_name[0]
print ('The case is built in directory: '+directory)

#create case directory
if not os.path.exists(directory):
    os.makedirs(directory)
#copy run files
copyfile(mdir+'/templates/reaction.jar.sbetr.nl',directory+'/reaction.jar.sbetr.nl')
copyfile(mdir+'/templates/reaction.1d.sbetr.nl',directory+'/reaction.1d.sbetr.nl')

#create build script
build_file=directory+'/case_make.py'

file = open(build_file,'w')
file.write("#!/usr/bin/env python\n")
file.write("#build script\n")
file.write("import os,sys,argparse\n")
file.write("parser = argparse.ArgumentParser(description=__doc__)\n")
file.write("parser.add_argument('--task', dest='task', metavar='task', type=str, nargs=1, default=['none'],"+
    "help='task for case_make.py: config, install or distclean')\n")

file.write("parser.add_argument('--CC', dest='CC', metavar='CC', type=str, nargs=1, default=['gcc8'],"+
    "help='C compiler')\n")
file.write("parser.add_argument('--CXX', dest='CXX', metavar='CXX', type=str, nargs=1, default=['g++8'],"+
    "help='C++ compiler')\n")
file.write("parser.add_argument('--FC', dest='FC', metavar='FC', type=str, nargs=1, default=['gf8'],"+
    "help='Fortran compiler')\n")
file.write("parser.add_argument('--extra', dest='extra', metavar='extra', type=str, nargs=1, default=[''],"+
    "help='other options')\n")
file.write("args = parser.parse_args()\n")
file.write("task =args.task[0]\n")
file.write("path='"+mdir+"'\n")
file.write("os.chdir(path)\n")
file.write("if task=='clean':\n")
file.write("    os.system('make distclean')\n")
file.write("elif task=='config':\n")
file.write("    command='make config CC='+args.CC[0]+' CXX='+args.CXX[0]+' FC='+args.FC[0]+' '+args.extra[0]\n")
file.write("    os.system(command)\n")
file.write("elif task=='install':\n")
file.write("    command='make install CC='+args.CC[0]+' CXX='+args.CXX[0]+' FC='+args.FC[0]+' '+args.extra[0]\n")
file.write("    os.system(command)\n")
file.write("    command='cp local/bin/sbetr "+directory+"/'\n")
file.write("    os.system(command)\n")
file.write("    command='cp local/bin/jarmodel "+directory+"/'\n")
file.write("    os.system(command)\n")
file.close()
#os.system(command)


#os.system(command)
