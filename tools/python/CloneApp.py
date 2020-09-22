#!/usr/bin/env python


import os,time,sys,argparse
import AppTools

def copy_file(fromf,tof, fromapp, toapp):
    """
    copy file fromf to tof and replace string fromapp to toapp
    """
    f2w=open(tof,"w")
    with open(fromf) as f:
        for line in f:
            newline=line.replace(fromapp,toapp)
            f2w.write(newline.replace(fromapp.upper(),toapp.upper()))
    f2w.close()

def copy_dir(dirf,dir2,fromapp,toapp):
    """
    copy files in directory dirf to dir2
    """
    import glob
    fls=glob.glob(dirf+"/*.F90")
    for f in fls:
        nf=f.replace(fromapp,toapp)
        copy_file(f,nf,fromapp,toapp)
    #copy the cmake file
    fromf=dirf+'/CMakeLists.txt'
    tof=dir2+'/CMakeLists.txt'
    copy_file(fromf,tof,fromapp,toapp)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--fromapp', dest="fromapp", metavar='fromapp', type=str, nargs=1, default=[""],
  help='the original app to be cloned')

parser.add_argument('--toapp', dest="toapp", metavar='toapp', type=str, nargs=1, default=[""],
  help='the new app to be created')


args = parser.parse_args()

fromapp=args.fromapp[0]

toapp=args.toapp[0]

dir_path = os.path.dirname(os.path.realpath(__file__))

mdir,chdir=dir_path.split('/sbetr')
sfarm_dir=mdir+'/sbetr/src/Applications/soil-farm'

#check existence of the app to cloned

if os.path.exists(sfarm_dir+'/'+fromapp):
    if toapp:
        #create directory
        os.system('mkdir -p '+sfarm_dir+'/'+toapp+"/"+toapp+'Para')
        os.system('mkdir -p '+sfarm_dir+'/'+toapp+"/"+toapp+'1layer')
        os.system('mkdir -p '+sfarm_dir+'/'+toapp+"/"+toapp+'Nlayer')
        #copy directory 1
        dir2=sfarm_dir+'/'+toapp+"/"+toapp+'Para'
        dirf=sfarm_dir+'/'+fromapp+'/'+fromapp+'Para'
        copy_dir(dirf,dir2,fromapp,toapp)
        #copy directory 2
        dir2=sfarm_dir+'/'+toapp+"/"+toapp+'1layer'
        dirf=sfarm_dir+'/'+fromapp+'/'+fromapp+'1layer'
        copy_dir(dirf,dir2,fromapp,toapp)
        #copy directory 3
        dir2=sfarm_dir+'/'+toapp+"/"+toapp+'Nlayer'
        dirf=sfarm_dir+'/'+fromapp+'/'+fromapp+'Nlayer'
        copy_dir(dirf,dir2,fromapp,toapp)
        #copy CMakeLists.txt
        fromf=sfarm_dir+'/'+fromapp+'/CMakeLists.txt'
        tof=sfarm_dir+'/'+toapp+"/CMakeLists.txt"
        copy_file(fromf, tof, fromapp, toapp)
        AppTools.add_case_file(mdir,toapp)
        print ("Clone finished!")
        print ("Check app in "+sfarm_dir+'/'+toapp)

    else:
        print ("please define new app name")
else:
    print ("app "+fromapp+" does not exist.")
    print ("Clone failed.")
