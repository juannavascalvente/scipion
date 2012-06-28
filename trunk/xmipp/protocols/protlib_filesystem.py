'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

#---------------------------------------------------------------------------
# Filesystem utilities
#---------------------------------------------------------------------------
import os
from os.path import join, exists, dirname, basename
from protlib_utils import printLog 
from shutil import copyfile
#from xmipp import *

# The following are Wrappers to be used from Protocols
# providing filesystem utilities

def createDir(log, path):
    """ Create directory, does not add workingdir"""
    from distutils.dir_util import mkpath
    from distutils.errors import DistutilsFileError
    try:
        if not exists(path):
            mkpath(path, 0777, True)
            printLog("Created dir " + path, log)
    except DistutilsFileError, e:
        printLog("Couldn't create dir: '%(path)s': %(e)s" % locals(), 
                 log, err=True, isError=True)

def changeDir(log, path):
    """ Change to Directory """
    try:
        os.chdir(path)
        printLog("Changed to dir " + path,log)
    except os.error, (errno, errstr):
        printLog("Could not change to directory '%s': Error (%d): %s" % (path, errno, errstr),log,err=True,isError=True)

def deleteDir(log, path):
    from distutils.dir_util import remove_tree
    if exists(path):
        remove_tree(path, True)
        printLog("Deleted directory " + path, log)
           
def deleteFile(log, filename, verbose=True):
    if exists(filename):
        os.remove(filename)
        if verbose:
            printLog('Deleted file %s' % filename,log)
    else:
        if verbose:
            printLog('Do not need to delete %s; already gone' % filename,log)
            
def copyFile(log, source, dest):
    try:
        copyfile(source, dest)
        printLog("Copied '%s' to '%s'" % (source, dest), log)
    except Exception, e:
        printLog("Could not copy '%s' to '%s'. Error: %s" % (source, dest, str(e)), log, err=True, isError=True)
    
def renameFile(log, source, dest):
    try:
        os.rename(source, dest)
        printLog("Renamed '%s' to '%s'" % (source, dest))
    except Exception, e:
        printLog("Could not rename '%s' to '%s'. Error: %s" % (source, dest, str(e)), log, err=True, isError=True)
    
def copyDir(log, source, dest):
    try:
        from shutil import copytree
        if exists(dest):
            deleteDir(log,dest)
        copytree(source, dest, symlinks=True)
        printLog("Copied '%s' to '%s'" % (source, dest))
    except Exception, e:
        printLog("Could not copy '%s' to '%s'. Error: %s" % (source, dest, str(e)), log, err=True, isError=True)

def deleteFiles(log, filelist, verbose):
    for file in filelist:
        deleteFile(log, file, verbose)

def createLink(log, source, dest):
    try:
        if os.path.exists(dest):
            os.remove(dest)
        destDir=os.path.split(dest)[0]
        os.symlink(os.path.relpath(source,destDir),dest)
        printLog("Linked '%s' to '%s'" % (source, dest))
    except Exception, e:
        printLog("Could not link '%s' to '%s'. Error: %s" % (source, dest, str(e)), log, err=True, isError=True)

def createLink2(log, filename, dirSrc, dirDest):
    createLink(log,os.path.join(dirSrc,filename),os.path.join(dirDest,filename))

def uniqueFilename(file_name):
    ''' Create a unique filename (not file handler)
       this approach is insecure but good enough for most purposes'''
    counter = 1
    file_name_parts = os.path.splitext(file_name) # returns ('/path/file', '.ext')
    while os.path.isfile(file_name):
        file_name = file_name_parts[0] + '_' + str(counter) + file_name_parts[1]
        counter += 1
    return file_name 

def random_alphanumeric(limit):
    import random
    #ascii alphabet of all alphanumerals
    r = (range(48,58) + range(65,91) + range(97,123))
    random.shuffle(r)
    return reduce(lambda i,s: i + chr(s),r[:limit],"")

def uniqueRandomFilename(file_name,randomLength=5):
    ''' Create a unique filename (not file handler)
       this approach is insecure but good enough for most purposes'''
    counter = 1
    file_name_parts = os.path.splitext(file_name) # returns ('/path/file', '.ext')
    file_name = file_name_parts[0] + '_' + random_alphanumeric(randomLength) + file_name_parts[1]
    while os.path.isfile(file_name):
        file_name = file_name_parts[0] + '_' + str(counter) + file_name_parts[1]
        counter += 1
    return file_name 

def replaceBasenameExt(filename, new_ext):
    '''Take the basename of the filename 
    and replace the extension by a new one'''
    return replaceFilenameExt(basename(filename), new_ext)

def removeBasenameExt(filename):
    '''Take the basename of the filename and remove extension'''
    return replaceBasenameExt(filename, '')

def removeFilenameExt(filename):
    return replaceFilenameExt(filename, '')

def removeFilenamePrefix(filename):
    if '@' in filename:
        return filename.split('@')[1]
    return filename

def replaceFilenameExt(filename, new_ext):
    ''' Replace the current filename extension by a new one'''
    return os.path.splitext(filename)[0] + new_ext
    
def findFilePath(filename, *pathList):
    '''Search recursively in path to find filename path(excluding filename)
    None is returned if not found'''
    for path in pathList:
        for root, dirs, files in os.walk(path):
            if filename in files:
                return root
    return None

#--------------------------- Xmipp specific tools ---------------------------------
def getXmippPath(*subpath):
    '''Return the path the the Xmipp installation folder
    if a subfolder is provided, will be concatenated to the path'''
    if os.environ.has_key('XMIPP_HOME'):
        return join(os.environ['XMIPP_HOME'], *subpath)  
    else:
        raise Exception('XMIPP_HOME environment variable not set') 

def includeProtocolsDir():
    protDir = getXmippPath('protocols')
    import sys
    sys.path.append(protDir)

def getProtocolTemplate(prot):
    protDir = getXmippPath('protocols')
    srcProtName = '%s.py' % prot.key
    #srcProtDir = getXmippPath('protocols')
    srcProtAbsPath = join(protDir, srcProtName)
    return srcProtAbsPath

def findProjectInPathTree(filename):
    found=False
    filename = os.path.abspath(filename)
    while filename!="/" and not found:
        if exists(join(filename,".project.sqlite")):
            found=True
        else:
            filename = dirname(filename)
    if found:
        return filename
    else:
        return None

def fixPath(filename, *pathList):
    if os.path.isabs(filename):
        return filename
    for path in pathList:
        filepath = join(path, filename)
        if exists(filepath):
            return filepath
    return None

def findRealFile(path, recursive=True):
    '''This function behaves like readlink with -f in shell'''
    from os import readlink
    if not recursive:
        return readlink(path) 
    from os.path import islink
    while islink(path):
        path = readlink(path)
    return path

def xmippExists(path):
    from xmipp import FileName
    return FileName(path).exists()

acquisitionInfoFileName="acquisition_info.xmd"
def linkAcquisitionInfoIfPresent(log,InputFile,dirDest):
    dirSrc=os.path.dirname(InputFile)
    fnAcquisitionIn=os.path.join(dirSrc,acquisitionInfoFileName)
    if os.path.exists(fnAcquisitionIn):
        createLink2(log, acquisitionInfoFileName, dirSrc, dirDest)
        
def AcquisitionInfoExists(InputFile):
    dirSrc=os.path.dirname(removeFilenamePrefix(InputFile))
    fnAcquisitionIn=os.path.join(dirSrc,acquisitionInfoFileName)
    if os.path.exists(fnAcquisitionIn):
        return True
    else:
        return False
    
def AcquisitionInfoGetSamplingRate(InputFile):
    from xmipp import MetaData,MDL_SAMPLINGRATE
    dirSrc=os.path.dirname(removeFilenamePrefix(InputFile))
    fnAcquisitionIn=os.path.join(dirSrc,acquisitionInfoFileName)
    if os.path.exists(fnAcquisitionIn):
        mD =MetaData(fnAcquisitionIn)
        id = mD.firstObject()
        return(mD.getValue(MDL_SAMPLINGRATE,id))
    else:
        return -1.
