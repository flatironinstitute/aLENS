import glob
import re
import os


def getFrameNumber_lambda(filename): return int(
    re.search('_([^_.]+)(?:\.[^_]*)?$', filename).group(1))


def result2pvd(name, ext):
    FileList = glob.glob('./result*/'+name+'_*.'+ext)
    # sort as numerical order
    FileList.sort(key=getFrameNumber_lambda)
    file = open(name+ext+'.pvd', "w")
    file.write(
        '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64"> \n')
    file.write('<Collection>\n')
    for i in range(len(FileList)):
        file.write('<DataSet timestep="'+str(i)+'"  file="'+FileList[i]+'"/>\n')
    file.write('</Collection>\n')
    file.write('</VTKFile> \n')


result2pvd('Sylinder', 'pvtp')
result2pvd('Protein', 'pvtp')
result2pvd('ConBlock', 'pvtp')
