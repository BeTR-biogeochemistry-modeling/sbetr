
import shutil

def add_case_cmake1(mdir, model):
    """
    add new model to camke files
    """
    fromf=mdir+'/sbetr/src/Applications/app_util/CMakeLists.txt'
    tof=mdir+'/sbetr/src/Applications/app_util/CMakeLists.txt1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    with open(fromf) as f:
        for line in f:
            if model in line:
                model_found=True
            if "end_appadd" in line:
                if model_found:
                    print (model+' already exist')
                else:
                    line1='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)\n'
                    line2='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnp1layer)\n'
                    line3='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpNlayer)\n'
                    newline1=line1.replace('ecacnp',model)
                    newline2=line2.replace('ecacnp',model)
                    newline3=line3.replace('ecacnp',model)
                    f2w.write(newline1)
                    f2w.write(newline2)
                    f2w.write(newline3)
                sec_start=False

            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True

    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)

def add_case_cmake2(mdir, model):
    """
    add new model to camke files
    """
    fromf=mdir+'/sbetr/src/Applications/soil-farm/CMakeLists.txt'
    tof=mdir+'/sbetr/src/Applications/soil-farm/CMakeLists.txt1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    with open(fromf) as f:
        for line in f:
            if model in line:
                model_found=True
            if "end_appadd" in line:
                if model_found:
                    print (model+' already exist')
                else:
                    line1='add_subdirectory(ecacnp)\n'
                    newline1=line1.replace('ecacnp',model)
                    f2w.write(newline1)
                sec_start=False

            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True
    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)

def add_case_cmake2(mdir, model):
    """
    add new model to camke files
    """
    fromf=mdir+'/sbetr/src/Applications/soil-farm/CMakeLists.txt'
    tof=mdir+'/sbetr/src/Applications/soil-farm/CMakeLists.txt1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    with open(fromf) as f:
        for line in f:
            if model in line:
                model_found=True
            if "end_appadd" in line:
                if model_found:
                    print (model+' already exist')
                else:
                    line1='add_subdirectory(ecacnp)\n'
                    newline1=line1.replace('ecacnp',model)
                    f2w.write(newline1)
                sec_start=False

            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True
    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)

def add_case_cmake3(mdir, model):
    """
    add new model to camke files
    """
    fromf=mdir+'/sbetr/src/jarmodel/driver/CMakeLists.txt'
    tof=mdir+'/sbetr/src/jarmodel/driver/CMakeLists.txt1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    with open(fromf) as f:
        for line in f:
            if model in line:
                model_found=True
            if "end_appadd" in line:
                if model_found:
                    print (model+' already exist')
                else:
                    line1='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)\n'
                    line2='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnp1layer)\n'
                    newline1=line1.replace('ecacnp',model)
                    newline2=line2.replace('ecacnp',model)
                    f2w.write(newline1)
                    f2w.write(newline2)
                sec_start=False

            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True
    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)

def add_case_cmake4(mdir, model):
    """
    add new model to camke files
    """
    fromf=mdir+'/sbetr/src/jarmodel/forcing/CMakeLists.txt'
    tof=mdir+'/sbetr/src/jarmodel/forcing/CMakeLists.txt1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    with open(fromf) as f:
        for line in f:
            if model in line:
                model_found=True
            if "end_appadd" in line:
                if model_found:
                    print (model+' already exist')
                else:
                    line1='include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)\n'
                    newline1=line1.replace('ecacnp',model)
                    f2w.write(newline1)
                sec_start=False

            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True
    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)


def add_case_fortranf(mdir, model,fromf):
    """
    add new model to camke files
    """
    tof=fromf+'1'
    f2w=open(tof,"w")
    sec_start=False
    model_found=False
    use_sec=False
    kk=0
    with open(fromf) as f:
        for line in f:
            if sec_start:
                #determine if it is a use section
                if 'use' in line:
                    use_sec=True
                    #read ecacnp line
                    if 'ecacnp' in line:
                        line0=line.replace('ecacnp',model)
                else:
                    if 'ecacnp' in line:
                        if kk == 0:
                            line0=line.replace('ecacnp',model)
                            kk=kk+1
                        elif kk == 1:
                            line1=line.replace('ecacnp',model)
                            kk=0
            if model in line:
                model_found=True
            if 'end_appadd' in line:
                if model_found:
                    print (model+' already exist')
                else:
                    if use_sec:
                        f2w.write(line0)
                        use_sec=False
                    else:
                        f2w.write(line0)
                        f2w.write(line1)
                sec_start=False
            f2w.write(line)
            if "begin_appadd" in line:
                sec_start=True
    f2w.close()
    shutil.move(tof,fromf)
    if not model_found:
        print ('Added new model '+model)

def further_instructions(mdir):
    """
    """
    print ("To turn on the new model, remember to accordingly modify the following files:")
    print (mdir+'/sbetr/src/dirver/main/HistBGCMod.F90')
    print (mdir+'/sbetr/src/dirver/main/sbetrDriverMod.F90')
    print (mdir+'/sbetr/src/dirver/standalone/BeTRSimulationStandalone.F90')
    print (mdir+'/sbetr/src/dirver/alm/BeTRSimulationALM.F90')
    print (mdir+'/sbetr/src/jarmodel/forcing/SetJarForcMod.F90')

def add_case_file(mdir,model):
    add_case_cmake1(mdir,model)
    add_case_cmake2(mdir,model)
    add_case_cmake3(mdir,model)
    add_case_cmake4(mdir,model)
    fromf=mdir+'/sbetr/src/jarmodel/driver/JarModelFactory.F90'
    add_case_fortranf(mdir,model,fromf)
    fromf=mdir+'/sbetr/src/Applications/app_util/ApplicationsFactory.F90'
    add_case_fortranf(mdir,model,fromf)
    fromf=mdir+'/sbetr/src/jarmodel/forcing/SetJarForcMod.F90'
    add_case_fortranf(mdir,model,fromf)

    further_instructions(mdir)
