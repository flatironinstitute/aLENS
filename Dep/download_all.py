import urllib.request
import os


def download(url):
    urllib.request.urlretrieve(url, os.path.basename(url))
    return


download('https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2')
download(
    'https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.bz2')
download('https://www.vtk.org/files/release/9.1/VTK-9.1.0.tar.gz')
download('https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-12-18-1.tar.gz')
