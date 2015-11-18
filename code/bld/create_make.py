#
#create_makef.py
#This script generates a Makefile based on the contents in the src code directory
#however some modification is necessary to have the correct cross-dependency figured out in the source files

#open the Makefile to write
ofile=open("Makefile","w")
#empty the file
ofile.truncate()

#write headers
ofile.write("#This is an automatically generated Makefile from the script create_make.py\n")
ofile.write("shell=/bin/sh\n")
ofile.write("include Macros\n")

import os
def list_files(mypath):
    file_paths =[]
    #Walk through the directories
    for root, directories, files in os.walk(mypath):
        for filename in files:
            filepath=os.path.join(root,filename)
            file_paths.append(filepath)
    return file_paths

curdir=os.path.dirname(os.path.realpath("../src/"))
print 'source code directory is'
ofile.write("#define source code directory\n")
ofile.write("SRCDIR="+curdir+"\n")

#define F90CC
ofile.write("F90CC=$(F90C) -I $(SRCDIR)/"+"src/esmf_wrf_timemgr\n")
#list all source codes
print "source codes:"
mypath="../src"

#list all directories

full_file_paths=list_files(mypath)
objs=[]
file_lowercase=[]

for file in full_file_paths:
    stuff=file.split('/')
    file_new=file.replace("..","$(SRCDIR)")
    file_lowercase.append(file.lower())
    obj=stuff[-1].replace(".F90",".o")
    objs.append(obj)




#define the grep function

def grep_use(file):
    bashcommand="grep use "+file
    import subprocess
    process=subprocess.Popen(bashcommand.split(),stdout=subprocess.PIPE)
    output=process.communicate()[0]
    outputs=output.split('\n')
    #try upper case

    bashcommand="grep USE "+file
    process=subprocess.Popen(bashcommand.split(),stdout=subprocess.PIPE)
    output=process.communicate()[0]
    outputs1=output.split('\n')
    outputs.extend(outputs1)

    mods=[]
    for line in outputs:
        tline=line.strip()
        tline=tline.lower()
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
tobjs=[]
tobjs.append(obj_main)
var=1
objsf=[]
while var == 1:
    lens=len(tobjs)
    if lens==0:
        break
    #get all use mods
    all_mods=[]
    for obj in tobjs:
        obj_file=find_fname(obj,full_file_paths)
#        print "grep mod for %s" % obj
        use_mods=grep_use(obj_file)
        all_mods.extend(use_mods)
    #remove duplicate mods
    all_mods=remove_listdup(all_mods)

    #determine if an obj is in the next layer
    result=[]
    all_mods=set(all_mods)
    for obj in tobjs:
        if obj not in all_mods:
            result.append(obj)
            objsf.append(obj)

    tobjs=list(all_mods)

#now add files that are not in objsf
objs1=set(objsf)
objs2=objsf
objs3=[]
for obj in objs:
    if obj.lower() not in objs1:
        objs2.append(obj)
        objs3.append(obj)
ofile.write("\n\n\n")
ofile.write("# define OBJS\n")
ofile.write("OBJS=")
cnt=0
for obj in reversed(objs2):
    loc=obj.find(".o")
    if loc < 0:
        continue
    #find its corresponding entry in objs
    for objj in objs:
        if obj == objj.lower():
            ofile.write(" "+objj)
            break
    cnt=cnt+1
    if cnt%7 == 0:
        ofile.write(" \\")
        ofile.write("\n")
        ofile.write("\t")

ofile.write("\n\n\n")
ofile.write("# define dependence\n")
for file in full_file_paths:
    stuff=file.split('/')
    file_new=file.replace("..","$(SRCDIR)")
    obj=stuff[-1].replace(".F90",".o")
    ofile.write(obj+": "+file_new+"\n")
    ofile.write("\t$(F90CC) "+file_new+"\n")

ofile.write("#define executable\n")
ofile.write(".PHONY : clean\n")
ofile.write("all: sbetr.exe\n")
ofile.write("sbetr.exe: $(OBJS)\n")
ofile.write("\t$(F90L) sbetr.exe $(OBJS)\n")
ofile.write("clean:\n")
ofile.write("\t@rm -rf build/*.o build/*.mod build/*.exe")
ofile.close()



#now creates the actual makefile,
