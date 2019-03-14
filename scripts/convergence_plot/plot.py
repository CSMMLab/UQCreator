import matplotlib.pyplot as plt
import argparse
import numpy as np
import datetime as dt
import os
import pandas as pd
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
        while section_ctr<3:
            if "=====" in lines[line_ctr]:
                section_ctr = section_ctr + 1
            elif "restartFile" in lines[line_ctr] and not "#" in lines[line_ctr]:
                restartFile = True
                restartFileName = lines[line_ctr].split('=')[-1].split("/")[-1][:-2]
            line_ctr = line_ctr + 1
    return restartFile, restartFileName

def getTEnd(dir,file):
    section_ctr = 0
    line_ctr = 0
    tEnd = 0.0
    with open(dir+'/'+file, 'r') as content:
        lines = content.readlines()
        while section_ctr<3:
            if "=====" in lines[line_ctr]:
                section_ctr = section_ctr + 1
            elif "tEnd" in lines[line_ctr] and not "#" in lines[line_ctr]:
                tEnd = float(lines[line_ctr].split('=')[-1][:-1])
            line_ctr = line_ctr + 1
    return tEnd

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
    section_ctr = 0
    line_ctr = 0
    with open(dir+'/'+file, 'r') as content:
        lines = content.readlines()
        while section_ctr<3:
            if "=====" in lines[line_ctr]:
                section_ctr = section_ctr + 1
            line_ctr += 1
        line_ctr += 1
        while "PE" in lines[line_ctr] or "residual" in lines[line_ctr]:
            line_ctr += 1
        lines = lines[line_ctr:]
        lines_to_del = []
        for i in range(len(lines)):
            if "Finished!" in lines[i]: 
                for j in range(i,len(lines)):
                    lines_to_del.append(j)    
                break
            if any(x in lines[i] for x in ['a','i','o','u']): #sketchy but works
                lines_to_del.append(i)
        for line in reversed(lines_to_del):
            del lines[line]
        header = ['date', 'time', 'delim', 't', 'residual_scaled', 'residual_abs']
        df = pd.read_csv(StringIO("\n".join(lines)), header=None, names=header, delim_whitespace = True)
        df = df.drop(columns=['delim'])
        df = df.dropna()
        t0 = time_to_float(df['date'][0]+' '+df['time'][0])
        df['runtime'] = np.zeros(len(df['time']))
        for index, row in df.iterrows():
            df.loc[index,'runtime'] = time_to_float(df['date'][index]+' '+df['time'][index]) - t0
        df = df.drop(columns=['date','time'])

        isRestartFile, prevFileName = checkForRestartFile(dir, file)
        if isRestartFile and os.path.isfile(dir+'/'+prevFileName) :
            df_prev = parse_logfile(dir, prevFileName)
            df['runtime'] += df_prev['runtime'].iloc[-1]
            tEnd_prev = getTEnd(dir, prevFileName)
            df['t'] += df_prev['t'].iloc[-1] - tEnd_prev
            return pd.concat([df_prev, df])
        else:
            return(df)

def create_runtime_residual_plots(data, labels, threshold=np.inf):
    plt.cla()
    plt.clf()
    for df in data:
        df = df[df['residual_scaled'] < threshold]
        plt.semilogy(df['runtime'],df['residual_scaled'])
    plt.xlabel('Runtime [s]')
    plt.ylabel('Residual')
    plt.legend(labels)
    plt.savefig("convergence_runtime_residual.pdf")
    print("Plot saved as: 'convergence_runtime_residual.pdf'")

def create_time_residual_plots(data, labels, tEnd=np.inf):
    plt.cla()
    plt.clf()
    for df in data:
        df = df[df['t'] < tEnd]
        plt.semilogy(df['t'],df['residual_scaled'])
    plt.xlabel('Time [s]')
    plt.ylabel('Residual')
    plt.legend(labels)
    plt.savefig("convergence_time_residual.pdf")
    print("Plot saved as: 'convergence_time_residual.pdf'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", "-d", type=str, required=True)
    parser.add_argument("--threshold", "-t", type=str, required=False)
    parser.add_argument("--tEnd", "-e", type=str, required=False)
    args = parser.parse_args()
    data = []
    labels = []
    files = os.listdir(args.dir)
    ignored_files = []
    for file in files:
        if file.endswith('_moments') or file.endswith('_duals'):
            ignored_files.append(file)
            continue
        isRestartFile, prevFileName = checkForRestartFile(args.dir, file)
        if isRestartFile:
            if prevFileName not in ignored_files: ignored_files.append(prevFileName)
            if prevFileName not in files: print("WARNING: restart file found but previous file was not found")              
    for file in ignored_files:
        if file in files: files.remove(file)
    for file in files:
        print("Processing '" + file + "'...")
        data.append(parse_logfile(args.dir, file))
        labels.append(create_label(args.dir, file))
    if(args.threshold):
        create_runtime_residual_plots(data, labels, args.threshold)
        create_time_residual_plots(data, labels, args.tEnd)
    else:
        create_runtime_residual_plots(data, labels)
        create_time_residual_plots(data, labels)
