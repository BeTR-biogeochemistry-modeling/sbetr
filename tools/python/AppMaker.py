#!/usr/bin/env python


import os,time,sys,argparse
import AppMakerPara
import AppMaker1layer
import AppMakerNlayer

import AppMakerParadef
import AppMaker1layerdef
import AppMakerNlayerdef

import shutil


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--appname', dest="appname", metavar='newapp', type=str, nargs=1, default=["newApp"],
  help='the name of the new app')

parser.add_argument('--action', dest="action", metavar='action', type=str, nargs=1, default=["add"],
  help='action needs to take for the code generator, add or remove')

parser.add_argument('--type', dest="type", metavar='type', type=str, nargs=1, default=["empty"],
  help='type of template to be used, empty or soilbgc')


args = parser.parse_args()

app_name=args.appname[0]


dir_path = os.path.dirname(os.path.realpath(__file__))

print dir_path

mdir,chdir=dir_path.split('/sbetr')
sfarm_dir=mdir+'/sbetr/src/Applications/soil-farm'

#mkdir directories
print ("Create directories")

if args.action[0] == 'add':
    print "add a new case"
    os.system('mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'Para')
    os.system('mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'1layer')
    os.system('mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'Nlayer')
    print ""
elif args.action[0] == 'remove':
    yourvar = input('really to remove the app '+app_name+': yes/no')
    print('you entered: ' + yourvar)
    if yourvar =='yes':
        try:
            print 'app: '+app_name +' is  removed'
            os.system('rm -rf '+sfarm_dir+'/'+app_name)
        except:
            pass
    else:
        print 'AppMaker did nothing'
    quit()
else:
    print 'AppMaker did nothing'
    quit()
print ("Create files")


if args.type[0] == 'soilbgc':
    #in app_namePara create two files
    #CMakeLists.txt
    AppMakerPara.MakePara(sfarm_dir, app_name)

    #in app_name1layer write three files
    #CMakeLists.txt
    #app_nameBGCIndexType.F90
    #app_nameBGCType.F90
    AppMaker1layer.Make1layer(sfarm_dir, app_name)


    AppMakerNlayer.MakeNlayer(sfarm_dir, app_name)
    #in app_nameNlayer write three files
    #CMakeLists.txt
    #app_nameBGCReactionsType.F90
    #app_namePlantSoilBGCType.F90

else:
    AppMakerParadef.MakePara(sfarm_dir, app_name)

    #in app_name1layer write three files
    #CMakeLists.txt
    #app_nameBGCIndexType.F90
    #app_nameBGCType.F90
    AppMaker1layerdef.Make1layer(sfarm_dir, app_name)

    AppMakerNlayerdef.MakeNlayer(sfarm_dir, app_name)

fcmake=open(sfarm_dir+'/'+app_name+"/CMakeLists.txt","w")
fcmake.write("set(BETR_LIBRARIES "+app_name+'Para '+app_name+'1layer '+ app_name+'Nlayer;${BETR_LIBRARIES} PARENT_SCOPE)\n')
fcmake.write('add_subdirectory('+app_name+'Para)\n')
fcmake.write('add_subdirectory('+app_name+'1layer)\n')
fcmake.write('add_subdirectory('+app_name+'Nlayer)\n')
fcmake.close()
