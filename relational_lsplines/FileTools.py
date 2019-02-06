# -*- coding: utf-8 -*-
"""
Created on Tue Jan 07 09:49:19 2014

@author: Luke.McCulloch
"""

import sys
import os


def fileTools():
    read = 'r'
    write = 'w'
    access_mode = [read,write]
    buffering = -1
    return access_mode, buffering


def GetFile(rootFldr, defaultFileName ):
    """Find and open a file"""
    if not defaultFileName:
        File = raw_input("Please enter Desired filename: %s"%defaultFileName + chr(8)*4)
        #if not File:
        #    File = defaultFileName 
    else:
        File = defaultFileName
    for root, dirs, files in os.walk(rootFldr): 
        for fCntr in files: 
            if fCntr == File: 
                rootFldr = root 
    return rootFldr, File
    

def GetLines(directory,filename):
    """ 
        inputs  :   string(directory)
                    string(filename)
                    
        outputs :   list[lines of the file]
    """
    access_mode, buffering = fileTools()    
    rootFldr = directory
    defaultFileName = filename
    FilePath, FileName = GetFile(rootFldr, defaultFileName )
    fileHandle = open(FilePath+'//'+FileName, access_mode[0], buffering)
    lines = fileHandle.readlines()
    fileHandle.close()
    return lines
    
def WriteLines(directory, filename, lines):
    access_mode, buffering = fileTools()    
    rootFldr = directory
    defaultFileName = filename
    FilePath, FileName = GetFile(rootFldr, defaultFileName )
    fileHandle = open(FilePath+'\\'+FileName, access_mode[1], buffering)
    for line in lines:
        fileHandle.write(line)
    fileHandle.close()
    return
    
def main():
    pre = '/home/luke/Documents/computational_naval_architecture/projects/relational_hull_design'
    directory = pre+'/relational_lsplines/relational_lsplines'
    filename =   'TemplateTable.tex'
    testlines = GetLines(directory,filename)
    for line in testlines:
        print line
    return testlines
    
if __name__=='__main__':
   lines =  main()
