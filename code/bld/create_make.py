#
#create_makef.py
#This script generates a Makefile based on the contents in the src code directory
#however some modification is necessary to have the correct cross-dependency figured out in the source files

print "list all contents in src"

import os
def list_files(mypath):
    file_paths =[]
    #Walk through the directories
    for root, directories, files in os.walk(mypath):
        for filename in files:
            filepath=os.path.join(root,filename)
            file_paths.append(filepath)
    return file_paths

curdir=os.path.dirname(os.path.realpath(__file__))
print "Current directory is:"
print curdir
print ' '
#list all source codes
print "source codes:"
mypath="../src"
full_file_paths=list_files(mypath)
objs=[]
file_lowercase=[]
for file in full_file_paths:
    stuff=file.split('/')
    file_new=file.replace("..","$SRCDIR")
    file_lowercase.append(file.lower())
    obj=stuff[-1].replace(".F90",".o")
    objs.append(obj)
    print "%s : %s"  % (obj,file_new)
    print "\t$(F90C) %s " % file_new

#define the grep function

def grep_use(file):
    bashcommand="grep use "+file
    import subprocess
    process=subprocess.Popen(bashcommand.split(),stdout=subprocess.PIPE)
    output=process.communicate()[0]
    print output

#now create the cross-dependency using the backward tracing method
#First: set the main file
#file1="sbetr.o"

#def get_nextnodes(file1, srcfiles_lower, srcfiles):
    #get the actual location of the file
    #get stem of the obj file
#    stuff=file1.split('.')
#    for file in srcfiles_lower:
#        loc=file.find(stuff[1])
#        if loc > 0:
#            print "do something"
    #grep all use commands
#    print "do sth else"
    #remove duplicating files

    #return

#grep all use commands in the source file and create a non-duplicating file list

#walk through the new list of files, if there's cross-dependency, put it as next, repeat untill the end of the tree

#find all files that have not been included in the tree

#rerverse and output

#now creates the actual makefile,
