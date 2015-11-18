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
    outputs=output.split('\n')
    mods=[]
    for line in outputs:
        tline=line.strip()
        loc=tline.find("use")
        if loc == 0:
            stuff=tline.split(' ')
            stuff_low=stuff[1].lower()
            #determine if it is a netcdf lib
            loc1=stuff_low.find("netcdf")
            loc2=stuff_low.find("intrinsic")
            if loc1 == -1 and loc2==-1 :
                mods.append(stuff_low.replace(",","")+".o")
    return mods


#get the stem of the obj file
def get_objstem(obj):
    stuff=obj.split(".")
    stem=stuff[0]
    return stem

#now create the cross-dependency using the backward tracing method
#First: set the main file
obj_main="sbetr.o"
def find_fname(obj, srcfiles):
    stem=get_objstem(obj)
    srcfl=""
    for srcf in srcfiles:
        srcf_low=srcf.lower()
        loc=srcf_low.find(stem)
        if loc > 0:
            srcfl=srcf
            break
    if len(srcfl) == 0:
        import sys
        sys.exit("source file not found for " + obj)
    return srcfl
#remove duplicates from a list of strings
def remove_listdup(inlist):
    list1=set(inlist)
    oulist=[]
    seen=set()
    for item in list1:
        if item not in seen:
            seen.add(item)
            oulist.append(item)
    return oulist



#start from the main file
objs=[]
objs.append(obj_main)
var=1
while var == 1:
    lens=len(objs)
    if lens==0:
        break
    #get all use mods
    all_mods=[]
    for obj in objs:
        obj_file=find_fname(obj,full_file_paths)
        print "grep mod for %s" % obj
        use_mods=grep_use(obj_file)
        all_mods.extend(use_mods)
    #remove duplicate mods
    all_mods=remove_listdup(all_mods)

    #determine if an obj is in the next layer
    result=[]
    all_mods=set(all_mods)
    for obj in objs:
        if obj not in all_mods:
            result.append(obj)

    objs=list(all_mods)
    print "RESULT"
    print result
    print "NEXT"
    print objs
    import time
    time.sleep(1)
    #go to the next loop


obj_file=find_fname(obj_main,full_file_paths)
use_mods=grep_use(obj_file)
print use_mods
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
