import matplotlib.pyplot as plt
import argparse
import numpy as np
import datetime as dt
import os
import pandas as pd
import re
from io import StringIO

def time_to_float(time):
    t = dt.datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f")
    epoch = dt.datetime.utcfromtimestamp(0)
    return (t-epoch).total_seconds()

def checkForRestartFile(dir, file):
    section_ctr = 0
    line_ctr = 0    
    restartFile = False
    restartFileName = ''
    with open(dir+'/'+file, 'r') as content:
        lines = content.readlines()        
        while section_ctr<3 and line_ctr<len(lines):                     
            if "=====" in lines[line_ctr]:
                section_ctr = section_ctr + 1
            elif "restartFile" in lines[line_ctr] and not "#" in lines[line_ctr]:
                restartFile = True
                restartFileName = lines[line_ctr].split('=')[-1].split("/")[-1][:-2]
            line_ctr = line_ctr + 1 
    return restartFile, restartFileName

def create_label(dir, file):    
    section_ctr = 0
    line_ctr = 0    
    with open(dir+'/'+file, 'r') as content:
        lines = content.readlines()        
        closure = ''
        nQuad = -1
        nMoments = -1   
        maxIter = -1
        restartFile = False
        while section_ctr<3:            
            if "#" in lines[line_ctr]:
                line_ctr = line_ctr + 1
                continue
            if "=====" in lines[line_ctr]:
                section_ctr = section_ctr + 1
            if "closure" in lines[line_ctr]:
                closure = lines[line_ctr].split("=")[1].strip().strip('"').strip('\n')
            elif "moments" in lines[line_ctr]:
                nQuad = int(lines[line_ctr].split("=")[1])
            elif "quadPoints" in lines[line_ctr]:
                nMoments = int(lines[line_ctr].split("=")[1])
            elif "maxIterations" in lines[line_ctr]:
                maxIter = int(lines[line_ctr].split("=")[1])
            elif "restartFile" in lines[line_ctr] and not "#" in lines[line_ctr]:
                restartFile = True                
            line_ctr = line_ctr + 1        
        if closure == 'StochasticGalerkin' and nQuad == nMoments:
            return 'SC'
        elif closure == 'StochasticGalerkin':
            return 'SG'
        elif closure == 'Euler2D' and restartFile and maxIter == 1:
            return 'caos-IPM'            
        elif closure == 'Euler2D' and restartFile:
            return 'ca-IPM'
        elif closure == 'Euler2D' and maxIter == 1:
            return 'os-IPM'
        elif closure == 'Euler2D':
            return 'IPM'
        elif closure == 'L2Filter':
            return 'L2Filter'
        elif closure == 'LassoFilter':
            return 'LassoFilter'         
        elif closure == 'BoundedBarrier':
            return 'BoundedBarrier'         
        else:
            return 'unkown'

def parse_logfile(dir, file):
    with open(dir+'/'+file, 'r') as content:
        print("Parsing:\t" + dir + '/' + file)
        lines = content.readlines()        
        header = ['date', 'time', 'delim', 'rho', 'rhoU_x', 'rhoU_y', 'rhoE'] 
        df = pd.read_csv(StringIO("\n".join(lines)), header=None, names=header, delim_whitespace = True)
        df = df.drop(columns=['delim'])
        df = df.dropna()
        t0 = time_to_float(df['date'][0]+' '+df['time'][0])
        df['runtime'] = np.zeros(len(df['time']))
        for index, row in df.iterrows():
            df.loc[index,'runtime'] = time_to_float(df['date'][index]+' '+df['time'][index]) - t0
        df = df.drop(columns=['date','time'])
        return(df)

def create_plot(dir, configFiles, files, title, type):   
    nStates = 4
    stateLabel = ['rho', 'rhoU_x', 'rhoU_y', 'rhoE']
    data = []
    label = []
    for f in range(len(files)):
        label.append(create_label(dir, configFiles[f]))
        isRestartFile, prevFileName = checkForRestartFile(dir, configFiles[f])
        if isRestartFile and os.path.isfile(dir+'/'+prevFileName):
            ext = ''
            if title == 'L1': ext = 'L1Error'
            elif title == 'L2': ext = 'L2Error'
            elif title == 'L-infinity': ext = 'LInfError'
            if type == 'E': ext += 'Mean'
            elif type == 'Var': ext += 'Var'
            df_prev = parse_logfile(dir, prevFileName+'_'+ext)
            df = parse_logfile(dir, files[f])
            df['runtime'] += df_prev['runtime'].iloc[-1]
            data.append(pd.concat([df_prev, df]))
        else:
            data.append(parse_logfile(dir, files[f]))
    for s in range(nStates):
        plt.cla()
        plt.clf()
        plt.xlabel('Time [s]')
        plt.ylabel('Residual')
        plt.title(title + ' error ' + type + '[' + stateLabel[s] + ']')
        for d in data:
            plt.semilogy(d['runtime'], d[stateLabel[s]])
        plt.legend(label)
        plt.savefig(title + '_error_' + type + '[' + stateLabel[s] + '].pdf')
        print('Saved:\t\t' + title + '_error_' + type + '[' + stateLabel[s] + '].pdf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", "-d", type=str, required=True)
    args = parser.parse_args()
    data = []
    labels = []
    files = os.listdir(args.dir)
    ignoredFiles = []
    configFiles = []
    L1ErrMeanFiles = []
    L2ErrMeanFiles = []
    LInfErrMeanFiles = []
    L1ErrVarFiles = []
    L2ErrVarFiles = []
    LInfErrVarFiles = []
    for file in files:
        isRestartFile, prevFileName = checkForRestartFile(args.dir, file)
        if isRestartFile:
            if prevFileName not in ignoredFiles: ignoredFiles.append(prevFileName)
            if prevFileName not in files: print("WARNING: restart file found but previous file was not found")    
    for file in files:
        for ifile in ignoredFiles:
            if file.startswith(ifile):
                if file not in ignoredFiles: ignoredFiles.append(file)
    for ifile in ignoredFiles:
        files.remove(ifile)
    for file in files:      
        if re.search(r'\d+$', file) is not None:
            configFiles.append(file)
        if file.endswith('_L1ErrorMean'):
            L1ErrMeanFiles.append(file)
        elif file.endswith('_L2ErrorMean'):
            L2ErrMeanFiles.append(file)
        elif file.endswith('_LInfErrorMean'):
            LInfErrMeanFiles.append(file)
        elif file.endswith('_L1ErrorVar'):
            L1ErrVarFiles.append(file)
        elif file.endswith('_L2ErrorVar'):
            L2ErrVarFiles.append(file)
        elif file.endswith('_LInfErrorVar'):            
            LInfErrVarFiles.append(file)
    #create_plot(args.dir, configFiles, L1ErrMeanFiles, 'L1', 'E')
    create_plot(args.dir, configFiles, L2ErrMeanFiles, 'L2', 'E')
    #create_plot(args.dir, configFiles, LInfErrMeanFiles, 'L-infinity', 'E')
    #create_plot(args.dir, configFiles, L1ErrVarFiles, 'L1', 'Var')
    #create_plot(args.dir, configFiles, L2ErrVarFiles, 'L2', 'Var')
    #create_plot(args.dir, configFiles, LInfErrVarFiles, 'L-infinity', 'Var')