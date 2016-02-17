#
#create_makef.py
#This script generates a Makefile based on the contents in the src code directory
#however some modification is necessary to have the correct cross-dependency figured out in the source files
#The algorithm works as following
#first a search is done to find main files, which have the key word program, the main files are excluded
#from attaching to the lowest layer after the dependence tree is constructred.
#then for each program, a dependence tree is constructred as following
#search for use/USE key words in the mother layer, find the dependent modules in the child layer,
#look for source code in the mother layer and determine if it belongs to the child layer, if so remove it
#from the mother layer and send it to the child layer. Next take the child layer as the new monther, then
#repeat. By the end of iteration, determine the residual as files that are not included in the tree.
#if the residual is not empty, the add the residual as the root, and redo the iteration. The residual from next iteration
#will be smaller. Repeat untill the residual is empty.
#Circular dependency is detected when in an iteraction all members in the mother layer are removed.



#do an early stop for debugging
from sys import exit

#define all internal functions

#list all files in the direcotry
import os
def list_files(mypath):
    file_paths =[]
    #Walk through the directories
    for root, directories, files in os.walk(mypath):
        for filename in files:
            loc1=root.find("betr_alm_coupler")
            loc2=root.find("betr_clm_coupler")
            if loc1 < 0 and loc2 < 0 :
                filepath=os.path.join(root,filename)
                loc3=filename.find("DS_Store")
                if loc3 < 0 :
                    file_paths.append(filepath)
    return file_paths

#define the grep function
def grep_func(file, words):
    bashcommand="grep "+words+" "+file
    import subprocess
    process=subprocess.Popen(bashcommand.split(),stdout=subprocess.PIPE)
    output=process.communicate()[0]
    outputs=output.split('\n')

    return outputs


#define the use grep function
def grep_mods(file):

    #try upper case
    outputs=grep_func(file,"use")
    outputs1=grep_func(file,"USE")
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
#for xxx.o the stem is xxx
def get_objstem(obj):
    stuff=obj.split(".")
    stem=stuff[0].lower()
    return stem

#find the corresponding source code for the named .o file
def find_fname(obj, srcfiles):
    stem=get_objstem(obj)
    srcfl=""
    for srcf in srcfiles:
        srcf_low=srcf.lower()
        loc=srcf_low.find(stem)
        if loc >= 0:
            srcfl=srcf
            break
    if len(srcfl) == 0:
        import sys
        sys.exit("source file not found for " + obj)
    return srcfl

#remove duplicates from a list of strings
def remove_listdup(inlist):
    list1=inlist
    oulist=[]
    seen=set()
    for item in list1:
        if item not in seen:
            seen.add(item)
            oulist.append(item)
    return oulist

#determine if the given file is a main file
def find_main(file):

    outputs=grep_func(file,"program")
    outputs1=grep_func(file,"PROGRAM")

    outputs.extend(outputs1)
    yesno=0
    for outstr in outputs:
        tstr=outstr.strip()
        tstr=tstr.lower()
        stuff=outstr.split(' ')
        if stuff[0] == "program":
             yesno=1

    return yesno

#define a class for value returning
class treeValue(object):
    def _init_(self,tree,residual):
        self.tree=tree
        self.residual=residual

#create dependence tree, return both the tree and residual
def create_deptree(mainobjs, full_file_paths, objs):
    #start from the main file
    tobjs=mainobjs
    var=1
    objsf=[]  #stores the dependency tree
    while var == 1:
        lens=len(tobjs)
        if lens==0:
            break
        #get all use mods
        all_mods=[]
        for obj in tobjs:
            obj_file=find_fname(obj,full_file_paths)
            #print "grep mod for %s" % obj
            use_mods=grep_mods(obj_file)
            #print use_mods
            all_mods.extend(use_mods)
            #remove duplicate mods
            all_mods=remove_listdup(all_mods)

        result=[] #stores the mother layer
        #determine if an obj is in the next layer
        all_mods=set(all_mods)
        for obj in tobjs:
            if obj not in all_mods:
                result.append(obj)
                objsf.append(obj)
        if not result:
            print "cross-dependency detected, check following files"
            print tobjs
            exit(0)

        #use tobjs for next iteration
        tobjs=list(all_mods)

    #now determine the residual
    objs1=set(objsf)
    objs2=objsf
    residual=[]
    for obj in objs:
        if obj.lower() not in objs1:
            residual.append(obj.lower())
    tree_ret=treeValue()
    tree_ret.tree=objsf
    tree_ret.residual=residual
    return tree_ret



#
mypath="../src"

#list files in all directories
full_file_paths=list_files(mypath)

#create the obj list
obj_ls=[]
for file in full_file_paths:
    #make sure it is an .F90 file
    loc1=file.find(".F90")
    loc2=file.find(".f90")
    if loc1 > 0 or loc2 >0 :
        stuff=file.split('/')
        file_new=file.replace("..","$(SRCDIR)")
        obj=stuff[-1].replace(".F90",".o")
        obj_ls.append(obj)

mainobj_ls=[]
cnt=0
for file in full_file_paths:
    yesno=find_main(file)
    if yesno == 1:
        mainobj_ls.append(obj_ls[cnt])
    cnt=cnt + 1

#do iterations to contruct the dependence tree
var=1
while var==1:
    temp_tree=create_deptree(mainobj_ls, full_file_paths, obj_ls)

    if len(temp_tree.residual)<=1:
        temp_tree.tree.extend(temp_tree.residual)
        break
    mainobj_ls.extend(temp_tree.residual)

#remove duplicates from a list of strings
ttree=reversed(temp_tree.tree)

tree_dep=remove_listdup(ttree)

#now write the Makefile

#open the Makefile to write
ofile=open("Makefile","w")
#empty the file
ofile.truncate()

#write headers
ofile.write("#This is an automatically generated Makefile from the script create_make.py\n")
ofile.write("shell=/bin/sh\n")
ofile.write("include Macros\n")

curdir=os.path.dirname(os.path.realpath("../src/"))

ofile.write("#define source code directory\n")
ofile.write("SRCDIR="+curdir+"\n")

#define F90CC as a dirty hack for the inc problem
ofile.write("F90CC=$(F90C) -I$(SRCDIR)/src/esmf_wrf_timemgr -I$(SRCDIR)/src/shr/\n")

ofile.write("\n\n\n")
ofile.write("# define OBJS\n")
ofile.write("OBJS=")

cnt=0
for obj in tree_dep:
    loc=obj.find(".o")
    if loc < 0:
        continue
    #find its corresponding entry in objs
    for objj in obj_ls:
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
ofile.write("\t@rm -rf *.o *.*mod *.exe")
ofile.close()



#now creates the actual makefile,
