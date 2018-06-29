#!/usr/bin/env python2

import os
import AppMakerPara
import AppMaker1layer
import AppMakerNlayer
dir_path = os.path.dirname(os.path.realpath(__file__))
print dir_path

mdir,chdir=dir_path.split('sbetr')
sfarm_dir=mdir+'src/Applications/soil-farm'

#mkdir directories
print ("Create directories")
app_name='newApp'
print 'mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'Para'
print 'mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'1layer'
print 'mkdir -p '+sfarm_dir+'/'+app_name+"/"+app_name+'Nlayer'
print ""

print ("Create files")

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
