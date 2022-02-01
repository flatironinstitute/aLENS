import os
import multiprocessing

##################################
#######   USER SWITCHES ##########
##################################

install = True  # compile only or compile + install
use_openblas = False   # use OpenBLAS by default, set to False if you use MKL
enable_packages = ['boost', 'eigen', 'trilinos', 'vtk']

# your installation destination, use ABSOLUTE PATH
os.environ["SFTPATH"] = os.environ['HOME'] + \
    '/envs/alens_intel/'

os.environ["CXXFLAGS"] = '-march=native -O3 -DNDEBUG'
os.environ["OPENMP_CXX_FLAGS"] = '-fopenmp'  # -qopenmp for intel compiler


##################################
####### END USER SWITCHES ########
##################################

os.environ["CFLAGS"] = os.environ["CXXFLAGS"]
os.environ["OPENMP_C_FLAGS"] = os.environ["OPENMP_CXX_FLAGS"]

print("install destination:", os.environ["SFTPATH"])
print("enable packages:", enable_packages)
print("flags:")
print(os.environ["CXXFLAGS"])
print(os.environ["CFLAGS"])
print(os.environ["OPENMP_CXX_FLAGS"])
print(os.environ["OPENMP_C_FLAGS"])

if use_openblas:
    pass
else:
    # use mkl
    print("Enable MKL")
    if 'MKLROOT' in os.environ:
        os.environ['MKL_INCLUDE_DIRS'] = os.environ['MKLROOT']+'/include'
        os.environ['MKL_LIB_DIRS'] = os.environ['MKLROOT']+'/lib/intel64'
    elif 'MKL_INCLUDE_DIRS' in os.environ and 'MKL_LIB_DIRS' in os.environ:
        pass
    else:
        msg = "must set either MKL_INCLUDE_DIRS/MKL_LIB_DIRS or MKLROOT in config.yaml\n"
        print(msg)
        exit()
    os.system('env | grep MKL')


make_jobs = multiprocessing.cpu_count()//2
if make_jobs <= 0:
    make_jobs = 4
print("make_jobs: ", make_jobs)


k = input("Press Y to continue, else to quit...  ")
if k != 'y' and k != 'Y':
    exit()


def run(cmd):
    print(cmd)
    os.system(cmd)
    return


depwd = os.getcwd()
log = depwd+'/compile.log'
err = depwd+'/compile.err'
os.system('date >'+log)
os.system('date >'+err)

if 'boost' in enable_packages:
    boost = 'boost_1_78_0'
    os.chdir(depwd)
    run('tar xf {}.tar.bz2'.format(boost))
    os.chdir('{}'.format(boost))
    run('./bootstrap.sh --prefix=' +
        os.environ["SFTPATH"]+' >> '+log+'  2>>'+err)
    if install:
        run('./b2 install'+' >> '+log+'  2>>'+err)

if 'eigen' in enable_packages:
    eigen = 'eigen-3.4.0'
    os.chdir(depwd)
    run('rm -rf ./build_eigen && mkdir ./build_eigen')
    run('tar xf {}.tar.bz2'.format(eigen))
    os.chdir('build_eigen')
    run('bash ../cmake-Eigen.sh && make -j' +
        str(make_jobs)+'  >> '+log+'  2>>'+err)
    if install:
        run('make install'+' >> '+log+'  2>>'+err)

if 'vtk' in enable_packages:
    vtk = 'VTK-9.1.0'
    os.chdir(depwd)
    run('rm -rf ./build_vtk && mkdir ./build_vtk')
    run('tar xf {}.tar.gz'.format(vtk))
    os.chdir('build_vtk')
    run('bash ../cmake-vtk.sh && make -j' +
        str(make_jobs)+'  >> '+log+'  2>>'+err)
    if install:
        run('make install'+' >> '+log+'  2>>'+err)


# Trilinos
if 'trilinos' in enable_packages:
    trilinos = 'trilinos-release-12-18-1'
    os.chdir(depwd)
    run('tar xf {}.tar.gz'.format(trilinos))
    if use_openblas:
        folder = 'build_trilinos_openblas'
        run('rm -rf ./{} && mkdir ./{}'.format(folder, folder))
        os.chdir(folder)
        run('bash ../cmake-Trilinos-OpenBLAS.sh && make -j' +
            str(make_jobs)+'  >> '+log+'  2>>'+err)
        if install:
            run('make install'+' >> '+log+'  2>>'+err)
    else:
        folder = 'build_trilinos_MKL'
        run('rm -rf ./{} && mkdir ./{}'.format(folder, folder))
        os.chdir(folder)
        run('bash ../cmake-Trilinos-MKL.sh && make -j' +
            str(make_jobs)+'  >> '+log+'  2>>'+err)
        if install:
            run('make install'+' >> '+log+'  2>>'+err)
