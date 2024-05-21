import subprocess
import sys
import os
import site
import glob
import logging
import urllib.request
import tarfile
import shutil
import tempfile
from distutils.ccompiler import new_compiler
from datetime import datetime
from distutils.sysconfig import customize_compiler
from setuptools import setup, Extension


# Formatter for the logger
class CustomFormatter(logging.Formatter):

    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# Install packages from pip ==============================================================
def install_with_pip(pack, vers=None, log=None, namepkg=None):

    # Update pip
    p = subprocess.Popen([sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()

    # sys.executable gives the path of the python interpreter
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if vers is None:
        mm = "{}: ** {}: Installing {}".format(now, namepkg, pack)
        print(mm) if log is None else log.info(mm)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}".format(pack)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if len(err) != 0:
            mm = "\t ERROR: installing {} from pip.\n".format(pack)
            mm += "\t {}".format(err)
            print(mm) if log is None else log.error(mm)
            exit()
        mm = "\t{}".format(out)
        print(mm) if log is None else log.info(mm)
    else:
        mm = "{}: ** {}: Installing {}=={}".format(now, namepkg, pack, vers)
        print(mm) if log is None else log.info(mm)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers), " &>install.log"])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if len(err) != 0:
            mm = "\t ERROR: installing {} from pip.\n".format(pack)
            mm += "\t {}".format(err)
            print(mm) if log is None else log.error(mm)
            exit()
        mm = "\t{}".format(out)
        print(mm) if log is None else log.info(mm)


# Install indigox-bond software ======================================================================================
def install_indigox_bond(log=None, namepkg=None):
    """
    Installing the indigo-bond library if is not present in the python enviorement.
    """

    import git

    giturl = 'https://github.com/allison-group/indigo-bondorder.git'
    install_dir = 'thirdparty/indigo-bondorder'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    mm = "\n\t\t ============== COMPILING & INSTALLING INDIGO-BONDORDER PACKAGE ==============\n\n"

    try:
        import indigox as ix
        mm += "{}: ** {}: indigo-bondorder is already installed in your system. {}".format(now, namepkg, giturl)
        print(mm) if log is None else log.info(mm)
    except ModuleNotFoundError:
        mm += "{}: ** {}: indigo-bondorder is not installed in your system\n".format(now, namepkg)
        mm += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        print(mm) if log is None else log.info(mm)

        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            mm = "================= ERROR INSTALL ================\n"
            mm += "** {}: Cannot find CMake executable needed \n".format(namepkg)
            mm += "** {}: for indigo-bondorder compilation.\n".format(namepkg)
            mm += "** {}: Instal CMake in your Linux\n".format(namepkg)
            mm += "** {}: The installation is aborted\n".format(namepkg)
            mm += "================= ERROR INSTALL ================"
            print(mm) if log is None else log.info(mm)
            exit()

        # Look at thirdparty directory
        if os.path.isdir("thirdparty"):
            pass
        else:
            os.makedirs("thirdparty")

        # Share data in indigox
        envdir = None
        for ipath in site.getsitepackages():
            g = glob.glob(os.path.join(ipath))
            if g:
                envdir = g[0]
                break

        fullpathlib_cmake = os.path.abspath(install_dir)
        fullpathdata_cmake = os.path.abspath(envdir+"/indigox/share")
        # Check if exists a distribution of indigox in the thirdparty directory
        # git clone https://github.com/allison-group/indigo-bondorder.git
        if os.path.isdir("thirdparty/indigo-bondorder"):
            pass
        else:
            try:
                git.Repo.clone_from(giturl, install_dir)
            except git.GitCommandError:
                if not os.path.isdir(install_dir):
                    mm = "================= ERROR INSTALL ================"
                    mm += "** {}: The github repository for indigo-bondorder is not valid " \
                          "or not exists.!!!\n".format(namepkg)
                    mm += "** {}: giturl     : {}\n".format(namepkg, giturl)
                    mm += "** {}: install_dir: {}\n".format(namepkg, install_dir)
                    mm += "** {}: Indigo-bondorder cannot be installed\n".format(namepkg)
                    mm += "** {}: The installation is aborted\n".format(namepkg)
                    mm += "================= ERROR INSTALL ================"
                    print(mm) if log is None else log.info(mm)
                    exit()
                else:
                    pass

            subprocess.call(["rm", "-rf", "thirdparty/indigo-bondorder/build"])
            subprocess.call(["mkdir", "thirdparty/indigo-bondorder/build"])
            os.chdir("thirdparty/indigo-bondorder/build")
            cmake_arguments = ["-DCMAKE_INSTALL_PREFIX={}".format(fullpathlib_cmake),
                               "-DCMAKE_INSTALL_DATAROOTDIR={}".format(fullpathdata_cmake)]
            subprocess.check_call(["cmake", "{}".format(fullpathlib_cmake)] + cmake_arguments)
            subprocess.call("make")
            subprocess.call(["make", "install"])
            os.chdir("../../")
            subprocess.call(["tar", "cvfz", "indigo-bondorder.tar.gz", "indigo-bondorder"])
            subprocess.call(["rm", "-rf", "./indigo-bondorder"])
            os.chdir("..")

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "** {}: {}\n".format(namepkg, now)
        mm += "** {}: envdir={}\n".format(namepkg, envdir)
        mm += "** {}: The *.so library has been installed in: {}/indigox/" \
              "pyindigox.cpython-36m-x86_64-linux-gnu.so\n".format(namepkg, envdir)
        mm += "                                        {envdir}/indigox/__init__.py\n"
        mm += "** {}: Share library for indigo in {}\n".format(namepkg, envdir+"/indigox/share")
        print(mm) if log is None else log.info(mm)

    try:
        import indigox as ix
        mm = "\n{}: ** {}: indigo-bondorder is correctly imported. {}".format(now, namepkg, giturl)
        print(mm) if log is None else log.info(mm)
    except (ModuleNotFoundError, ImportError):
        mm = "================= ERROR INSTALL ================\n"
        mm += "{}: ** {}: indigo-bondorder libray cannot be imported as:\n".format(now, namepkg)
        mm += "{}: ** {}: \timport indigox as ix\n".format(now, namepkg)
        mm += "{}: ** {}: Something wrong during compilation.\n".format(now, namepkg)
        mm += "================= ERROR INSTALL ================"
        print(mm) if log is None else log.info(mm)
        exit()


# Install eigen library software ======================================================================================
def install_eigen(log=None, namepkg=None):

    """
    Installing the eigen library which is needed in openbabel.
    """

    # Version eigen 3.4.0 might not work. Error obtained compiler too old.
    # giturl = 'https://gitlab.com/libeigen/eigen.git'
    giturl = 'https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz'
    tar_download_file = os.path.basename(giturl)
    install_dir = 'thirdparty/eigen-3.3.9'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    mm = "\n\t\t ============== COMPILING & INSTALLING EIGEN PACKAGE ==============\n\n"
    print(mm) if log is None else log.info(mm)

    if not os.path.isdir(install_dir+"/eigen_library/include"):
        mm = "{}: ** {}: eigen is not installed in your system\n".format(now, namepkg)
        mm += "https://eigen.tuxfamily.org/index.php?title=Main_Page\n"
        mm += "** {}: Installing version 3.3.9 from git... {}".format(namepkg, giturl)
        print(mm) if log is None else log.info(mm)

    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        mm = "================= ERROR INSTALL ================\n"
        mm += "** {}: Cannot find CMake executable\n".format(namepkg)
        mm += "** {}: The installation is aborted\n".format(namepkg)
        mm += "================= ERROR INSTALL ================\n"
        print(mm) if log is None else log.info(mm)
        sys.exit()

    # Look at thirdparty directory
    if os.path.isdir("thirdparty"):
        pass
    else:
        os.makedirs("thirdparty")

    fullpath_cmake = os.path.abspath(install_dir)

    # Check if exists a distribution of indigox in the thirdparty directory
    if os.path.isdir(install_dir+"/eigen_library/include"):
        mm = "{}: ** {}: eigen_library seems to be already compiled in your system. " \
            "{}".format(now, namepkg, install_dir+"/eigen_library/include")
        print(mm) if log is None else log.info(mm)
    else:
        # git clone https://gitlab.com/libeigen/eigen.git
        try:
            # git.Repo.clone_from(giturl, install_dir)
            urllib.request.urlretrieve(giturl, "thirdparty/"+tar_download_file)
            tar = tarfile.open("thirdparty/"+tar_download_file)
            tar.extractall(path="./thirdparty/")
            tar.close()
        except (Exception,) as e:
            if not os.path.isdir(install_dir):
                mm = "================= ERROR INSTALL ================\n"
                mm += "** {}: The github repository for eigen is not valid or not exists.!!!\n".format(namepkg)
                mm += "** {}: giturl     : {}\n".format(namepkg, giturl)
                mm += "** {}: install_dir: {}\n".format(namepkg, install_dir)
                mm += "** {}: eigen cannot be installed\n".format(namepkg)
                mm += "** {}: The installation is aborted\n".format(namepkg)
                mm += "** {}\n".format(e)
                mm += "================= ERROR INSTALL ================"
                print(mm) if log is None else log.info(mm)
                exit()
            else:
                pass

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "{}: ** {}: Compiling eigen\n".format(now, namepkg)
        print(mm) if log is None else log.info(mm)
        subprocess.call(["rm", "-rf", install_dir+"/build"])
        subprocess.call(["mkdir",  install_dir+"/build"])
        subprocess.call(["mkdir", install_dir+"/eigen_library"])
        os.chdir(install_dir+"/build")
        cmake_arguments1 = ["-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake+"/eigen_library")]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake)]+cmake_arguments1)
        subprocess.call("make")
        subprocess.call(["make", "install"])
        os.chdir("../../")
        os.chdir("..")

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "{}: ** {}: The eigen library has been installed in: thirdparty/eigen/eigen_library\n".format(now, namepkg)
        print(mm) if log is None else log.info(mm)


# Install indigox-bond software ======================================================================================
def install_openbabel(log=None, namepkg=None):

    """
    Installing the openbabel library if is not present in the python enviorement.
    """

    import git

    giturl = 'https://github.com/openbabel/openbabel.git'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    mm = "\n\t\t ============== COMPILING & INSTALLING OPENBABEL PACKAGE ==============\n\n"
    print(mm) if log is None else log.info(mm)

    # Look at thirdparty directory
    if not os.path.isdir("thirdparty"):
        os.makedirs("thirdparty")

    install_dir = 'thirdparty'
    eigen_dir = 'thirdparty/eigen-3.3.9/eigen_library/include/eigen3'
    fullpath_cmake = os.path.abspath(install_dir)
    fullpath_eigen = os.path.abspath(eigen_dir)

    try:
        from openbabel import openbabel as ob
        mm = "{}: ** {}: openbabel is already installed in your system. " \
            "{}".format(now, namepkg, fullpath_cmake + "/openbabel")
        print(mm) if log is None else log.info(mm)
    except (ModuleNotFoundError, ImportError):
        if os.path.isdir(os.path.join(fullpath_cmake, "openbabel")):
            subprocess.call(["rm", "-rf", os.path.join(fullpath_cmake, "openbabel")])
        mm = "{}: Downloading: ... openbabel-3-1-1(Wait for one minute)\n".format(now)
        print(mm) if log is None else log.info(mm)
        git.Repo.clone_from(giturl, fullpath_cmake + "/openbabel")

    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        mm = "================= ERROR INSTALL ================\n"
        mm += "** {}: Cannot find CMake executable\n".format(namepkg)
        mm += "** {}: The installation is aborted\n".format(namepkg)
        mm += "================= ERROR INSTALL ================\n"
        print(mm) if log is None else log.info(mm)
        exit()

    # Check if swig is installed in the system. This is needed in order to build the python bindings
    error = subprocess.call(["which", "swig"])
    if error != 0:
        mm = "================= ERROR INSTALL ================\n"
        mm += "** {}: Cannot find Swig executable\n".format(namepkg)
        mm += "** {}: Try to install swig in your system (Ubuntu: apt-get install swig)\n".format(namepkg)
        mm += "** {}: The installation is aborted\n".format(namepkg)
        mm += "================= ERROR INSTALL ================\n"
        print(mm) if log is None else log.info(mm)
        exit()

    if not os.path.isdir(fullpath_cmake+"/openbabel/build"):
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "{}: ** {}: Compiling openbabel\n".format(now, namepkg)
        print(mm) if log is None else log.info(mm)
        subprocess.call(["rm", "-rf", "thirdparty/openbabel/build"])
        subprocess.call(["mkdir", "thirdparty/openbabel/build"])
        os.chdir("thirdparty/openbabel/build")
        cmake_arguments1 = ["-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake+"/openbabel")]
        cmake_arguments2 = ["-DPYTHON_BINDINGS=ON"]
        cmake_arguments3 = ["-DRUN_SWIG=ON"]
        cmake_arguments4 = ["-DEIGEN3_INCLUDE_DIR={}".format(fullpath_eigen)]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake+"/openbabel")] +
                              cmake_arguments1+cmake_arguments2+cmake_arguments3+cmake_arguments4)
        subprocess.call(["make", "-j4"])
        subprocess.call(["make", "install"])
        os.chdir("../../")
        os.chdir("..")

        # Copy the library to the root site-engines of the python distribution
        dir_env_python = site.getsitepackages()[0]
        label = str(sys.version_info.major)+"."+str(sys.version_info.minor)
        dir_openbabel_installed = fullpath_cmake+"/openbabel/lib/python{}/site-packages/openbabel".format(label)
        subprocess.call(["rm", "-rf", "{}".format(os.path.join(dir_env_python, "openbabel"))])
        subprocess.call(["mv", "-f", "{}".format(dir_openbabel_installed), "{}".format(dir_env_python)])
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "{}: ** {}: Open babel move from {} to {}\n".format(now, namepkg,
                                                                 os.path.join(dir_env_python, "openbabel"),
                                                                 dir_env_python)
        print(mm) if log is None else log.info(mm)

    try:
        from openbabel import openbabel as ob
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "\n{}: ** {}: openbabel is correctly imported. {}".format(now, namepkg, giturl)
        print(mm) if log is None else log.info(mm)
    except (ModuleNotFoundError, ImportError) as e:
        mm = "================= ERROR INSTALL ================\n"
        mm += "{}: ** {}: openbabel libray cannot be imported as:\n".format(now, namepkg)
        mm += "{}: ** {}: \tfrom openbabel import openbabel as ob\n".format(now, namepkg)
        mm += "{}: ** {}: Something wrong during compilation.\n".format(now, namepkg)
        mm += "Error: {}\n".format(e)
        mm += "================= ERROR INSTALL ================"
        print(mm) if log is None else log.info(mm)
        exit()


# Install Dock RMSD software =========================================================================================
def install_dockrmsd(log=None, namepkg=None):

    """
    Downloading, compiling and installing the DockRMSD software

    E.W. Bell, Y. Zhang. DockRMSD: an Open-Source Tool for Atom Mapping and
    RMSD calculation of Symmetric Molecules through Graph Isomorphism. Journal of Cheminformatics, 11:40 (2019)

    https://zhanglab.ccmb.med.umich.edu/DockRMSD/
    """

    import urllib.request

    install_dir = 'thirdparty/dock_rmsd'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    mm = "\n\t\t ============== COMPILING & INSTALLING DockRMSD PACKAGE ==============\n\n"

    if not os.path.isfile("thirdparty/dock_rmsd/dockrmsd.x"):
        mm += "{}: ** {}: DockRMSD is not installed in your system\n".format(now, namepkg)
        mm += "\t\thttps://zhanglab.ccmb.med.umich.edu/DockRMSD/\n"
        mm += "\t\t** {}: Installing from web... {}".format(namepkg, "https://zhanglab.ccmb.med.umich.edu/DockRMSD/\n")

        print(mm) if log is None else log.info(mm)
    else:
        mm += "{}: ** {}: DockRMSD is already installed in your system\n".format(now, namepkg)
        print(mm) if log is None else log.info(mm)
        return True

    # Look at thirdparty directory
    if os.path.isdir("thirdparty"):
        pass
    else:
        os.makedirs("thirdparty")

    fullpath_cmake = os.path.abspath(install_dir)+"/"

    url1 = "https://zhanglab.ccmb.med.umich.edu/DockRMSD/DockRMSD.h"
    url2 = "https://zhanglab.ccmb.med.umich.edu/DockRMSD/DockRMSD.c"

    if not os.path.isdir(fullpath_cmake):
        os.mkdir(fullpath_cmake)
    urllib.request.urlretrieve(url1, fullpath_cmake+"DockRMSD.h")
    urllib.request.urlretrieve(url2, fullpath_cmake+"DockRMSD.c")

    subprocess.call(["gcc", fullpath_cmake+"DockRMSD.c", "-o", fullpath_cmake+"dockrmsd.x", "-lm", "-O3"])


# Install CCLIB software =============================================================================================
def install_cclib(log=None, namepkg=None):
    """
    Installing cclib from github instead of pip. The more recent versions are in github
    """

    import git

    giturl = 'https://github.com/cclib/cclib.git'
    install_dir = 'thirdparty/cclib'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    mm = "\n\t\t ============== CLONING & INSTALLING CCLIB LIBRARY ==============\n\n"

    try:
        import cclib
    #     mm += "{}: ** {}: indigo-bondorder is already installed in your system. {}".format(now, namepkg, giturl)
    #     print(mm) if log is None else log.info(mm)
    except ModuleNotFoundError:
        mm += "{}: ** {}: cclib is not installed in your system\n".format(now, namepkg)
        mm += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        print(mm) if log is None else log.info(mm)

        # Look at thirdparty directory
        if os.path.isdir("thirdparty"):
            pass
        else:
            os.makedirs("thirdparty")

        # Share data in indigox
        envdir = None
        for ipath in site.getsitepackages():
            g = glob.glob(os.path.join(ipath))
            if g:
                envdir = g[0]
                break

        fullpathlib_cmake = os.path.abspath(install_dir)
        fullpathdata_cmake = os.path.abspath(envdir+"/indigox/share")
        # Check if exists a distribution of indigox in the thirdparty directory
        # git clone https://github.com/allison-group/indigo-bondorder.git
        if os.path.isdir("thirdparty/cclib"):
            pass
        else:
            try:
                git.Repo.clone_from(giturl, install_dir)
            except git.GitCommandError:
                if not os.path.isdir(install_dir):
                    mm = "================= ERROR INSTALL ================"
                    mm += "** {}: The github repository for cclib is not valid " \
                          "or not exists.!!!\n".format(namepkg)
                    mm += "** {}: giturl     : {}\n".format(namepkg, giturl)
                    mm += "** {}: install_dir: {}\n".format(namepkg, install_dir)
                    mm += "** {}: cclib cannot be installed\n".format(namepkg)
                    mm += "** {}: The installation is aborted\n".format(namepkg)
                    mm += "================= ERROR INSTALL ================"
                    print(mm) if log is None else log.info(mm)
                    exit()
                else:
                    pass

            os.chdir("thirdparty/cclib")
            cmake_arguments = ["setup.py", "install"]
            os.system("python setup.py install")
            os.chdir("../")
            subprocess.call(["tar", "cvfz", "cclib.tar.gz", "cclib"])
            subprocess.call(["rm", "-rf", "./cclib"])
            os.chdir("..")

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        mm = "** {}: {}\n".format(namepkg, now)
        mm += "** {}: envdir={}\n".format(namepkg, envdir)
        mm += "** {}: The cclib library has been installed\n"
        print(mm) if log is None else log.info(mm)

    # try:
    #     import cclib
    #     mm = "\n{}: ** {}: cclib is correctly imported. {}".format(now, namepkg, giturl)
    #     print(mm) if log is None else log.info(mm)
    # except (ModuleNotFoundError, ImportError):
    #     mm = "================= ERROR INSTALL ================\n"
    #     mm += "{}: ** {}: cclib libray cannot be imported as:\n".format(now, namepkg)
    #     mm += "{}: ** {}: \timport cclib\n".format(now, namepkg)
    #     mm += "{}: ** {}: Something wrong during compilation.\n".format(now, namepkg)
    #     mm += "================= ERROR INSTALL ================"
    #     print(mm) if log is None else log.info(mm)
    #     exit()


# Disabling-output-when-compiling-with-distutil =================================================
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as fout:
                if include is not None:
                    fout.write('#include {0!s}\n'.format(include))
                fout.write('int main(void) {\n')
                fout.write('    {0!s};\n'.format(funcname))
                fout.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except (Exception,):
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


# Does this compiler support OpenMP parallelization?""" ==============================================================
def detect_openmp():
    print("GECOS2: Attempting to autodetect OpenMP support... ", end="")
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("GECOS2: Compiler supports OpenMP")
    else:
        print("GECOS2: Did not detect OpenMP support.")

    return hasopenmp


# Setup external extensions ==============================================================
def setup_external_extensions(debug_cflags=False, use_openmp=True):
    has_openmp = detect_openmp()

    # parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    mathlib = ['m']
    define_macros = []
    extra_compile_args = ['-std=c99', '-ffast-math', '-O3', '-funroll-loops', '-Wno-cpp']
    if debug_cflags:
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    extensions_install = [
        Extension("ext_libc.c_distances_openmp",
                  ["gecos/ext_libc/c_distances_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),
    ]

    return extensions_install


# Main setup
if __name__ == '__main__':

    # Creating the logger to install.log file ===================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Starting installation!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)

    # Wheel must be needed for install
    try:
        import wheel
    except ModuleNotFoundError:
        nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "Install wheel in your system\n"
        m += "================= ERROR INSTALL ================\n"
        m += "{}: ** {}: wheel libray cannot be imported as:\n".format(nowm, "GECOS")
        m += "{}: ** {}: \timport wheel\n".format(nowm, "GECOS")
        m += "{}: ** {}: Install wheel in your system:\n".format(nowm, "GECOS")
        m += "{}: ** {}: \tpython -m pip install wheel\n".format(nowm, "GECOS")
        m += "================= ERROR INSTALL ================"
        print(m) if logger is None else logger.info(m)
        exit()

    # Print sys path ===================================
    m1 = "\t\t ============== SYS PATH ==============\n"
    for item in sys.path:
        m1 += item + "\n"
    m1 += "\n\t\t ============== INSTALLING PIP PACKAGES ==============\n"
    print(m1) if logger is None else logger.info(m1)

    # Install requirements ===================================
    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, vers=version, log=logger, namepkg="GECOS")
        except ValueError:
            pkg = ipack
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, log=logger, namepkg="GECOS")

    # Install third-party software ===========================
    import git
    install_indigox_bond(log=logger, namepkg="GECOS")
    install_eigen(log=logger, namepkg="GECOS")
    install_openbabel(log=logger, namepkg="GECOS")
    install_dockrmsd(log=logger, namepkg="GECOS")
    install_cclib(log=logger, namepkg="GECOS")     # Install cclib from github (more recent versions)

    # Setup Gecos ===========================================
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t ============== RUNNING SETUP FROM SETUPTOOLS {} ==============\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    # Setup TOPOLOGY ===========================================
    from Cython.Build import cythonize
    import Cython
    import numpy
    print(Cython.__version__)

    # Extensions
    extensions = setup_external_extensions()
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t ============== RUNNING SETUP FROM SETUPTOOLS {} ==============\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    setup(
        ext_modules=cythonize(extensions,
                              compiler_directives={'language_level': sys.version_info[0]}),
        include_dirs=[numpy.get_include()]
    )

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Installation Done!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
