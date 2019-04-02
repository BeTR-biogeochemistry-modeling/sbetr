#!/usr/bin/env python2


import os,time,sys,argparse
import AppTools


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--newapp', dest="newapp", metavar='newapp', type=str, nargs=1, default=[""],
  help='the app to be added')


args = parser.parse_args()

newapp=args.newapp[0]

dir_path = os.path.dirname(os.path.realpath(__file__))

mdir,chdir=dir_path.split('/sbetr')
sfarm_dir=mdir+'/sbetr/src/Applications/soil-farm'

#check existence of the app to cloned


if newapp:
    AppTools.add_case_file(mdir,newapp)
    print "Clone finished!"
    print "Check app in "+sfarm_dir+'/'+newapp

else:
    print "please define new app name"
