import argparse
import numpy as np
import os
import re
from shutil import copyfile
from shutil import SameFileError

def checkForRestartFile(dir, file):
    section_ctr = 0
    line_ctr = 0
    restartFile = False
    restartFileName = ''
    with open(dir+'/'+file, 'r') as content:
        lines = content.readlines()
        while section_ctr<3 and line_ctr<len(lines):
            current_line = lines[line_ctr]
            if "=====" in current_line:
                section_ctr = section_ctr + 1
            elif "restartFile" in current_line:
                restartFile = True
                restartFileName = current_line.split('=')[-1].split("/")[-1][:-2]
            line_ctr = line_ctr + 1
    return restartFile, restartFileName

def sort_files(dir_in, dir_out, machine):
    vtkFiles = []
    logFiles = []
    logFilePattern = re.compile("(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})")
    for root, _, files in os.walk(dir_in):
        if not root.startswith('./'+dir_out):
            for file in files:
                if file.endswith('.vtk'):
                    vtkFiles.append(os.path.join(root,file))
                elif logFilePattern.match('_'.join(file.split('_')[:2])) and 'LInfError' not in file:
                    logFiles.append(os.path.join(root,file))
    for logFile in logFiles:
        if logFile.count('_') == 1 and os.stat(logFile + '_moments').st_size > 0:
            try:
                methodName, vtkBaseName = parse_logfile(logFile)
            except IndexError:
                continue
            dst = os.path.join(dir_out, methodName, machine, vtkBaseName)
            if not os.path.exists(dst):
                os.makedirs(dst)                
            vtkFilesToCopy = [file for file in vtkFiles if ''.join(file.split('/')[-1]).startswith(vtkBaseName)]            
            logFilesToCopy = [file for file in logFiles if file.startswith(logFile)]
            for file in vtkFilesToCopy + logFilesToCopy:
                try:
                    copyfile(file, os.path.join(dst, ''.join(file.split('/')[-1])))
                except SameFileError:
                    continue

def parse_logfile(file):
    section_ctr = 0
    line_ctr = 0
    with open(file, 'r') as content:
        lines = content.readlines()
        closure = ''
        nQuad = -1
        nMoments = -1
        maxIter = -1
        method = ''
        outputFile = ''        
        restartFile = False
        while section_ctr<3:
            current_line = lines[line_ctr]
            if "#" not in current_line:
                if "=====" in current_line:
                    section_ctr = section_ctr + 1
                elif "closure" in current_line:
                    closure = current_line.split("=")[1].strip().strip('"').strip('\n')
                elif "moments" in current_line:
                    try:
                        nMoments = np.asarray(current_line.split("=")[1][1:-2].split('],')[1][1:-1].split(','), dtype=int)
                    except IndexError:
                        nMoments = np.array([int(current_line.split("=")[1])])
                elif "quadPoints" in current_line:
                    nQuad = int(current_line.split("=")[1])
                elif "maxIterations" in current_line:
                    maxIter = int(current_line.split("=")[1])
                elif "restartFile" in current_line:
                    restartFile = True
                elif "outputFile" in current_line:
                    outputFile = current_line.split("=")[1].strip()[1:-1]
            line_ctr = line_ctr + 1
        if closure == 'StochasticGalerkin' and nQuad == nMoments and 'psc' in outputFile: method = 'PSC'          
        elif closure == 'StochasticGalerkin' and nQuad == nMoments: method = 'SC'        
        elif closure == 'StochasticGalerkin': method = 'SG'
        elif closure == 'Euler2D' and restartFile and maxIter == 1: method = 'CAOS-IPM'
        elif closure == 'Euler2D' and restartFile: method = 'CA-IPM'
        elif closure == 'Euler2D' and maxIter == 1: method = 'OS-IPM'
        elif closure == 'Euler2D': method = 'IPM'
        elif closure == 'L2Filter': method = 'L2Filter'
        elif closure == 'LassoFilter': method = 'LassoFilter'
        elif closure == 'BoundedBarrier': method = 'BoundedBarrier'
        else: method = 'unkown'
        return method, outputFile

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", "-m", type=str, required=True)
    parser.add_argument("--input" , "-i", type=str, required=False, default='.')
    parser.add_argument("--output", "-o", type=str, required=False, default='output')
    args = parser.parse_args()
    dir_in  = args.input
    dir_out = args.output
    machine = args.machine
    sort_files(dir_in, dir_out, machine)