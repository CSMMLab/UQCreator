import matplotlib.pyplot as plt
import argparse
import numpy as np
import datetime as dt
import os, time

parser = argparse.ArgumentParser()
parser.add_argument("--file", "-f", type=str, required=True)
parser.add_argument("--compareTo", "-c", type=str, required=False)
parser.add_argument("--live", "-l", action='store_true', required=False)
args = parser.parse_args()

def splitLine(line, splitchar):
    list = line.split(splitchar)
    list = [" ".join(x.split()) for x in list]
    return list

def time_to_float(time):
    t = dt.datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f")
    epoch = dt.datetime.utcfromtimestamp(0)
    return (t-epoch).total_seconds()

def plot(x,y,xlabel,ylabel):
    plt.cla()
    plt.clf()
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(xlabel+"_"+ylabel+".pdf")

def create_plots():
    t_arr = []
    residual_arr = []
    iter_arr = []
    timestamp_arr = []
    with open(args.file, 'r') as f:
        for line in f:
            [timestamp, data] = splitLine(line, '|')
            data_list = splitLine(data, ' ')
            try:
                for d in data_list:
                    d = float(d)
            except ValueError:
                continue
            data_list = [float(x) for x in data_list]
            t_arr.append(data_list[0])
            residual_arr.append(data_list[1])
            timestamp_arr.append(time_to_float(timestamp))
    for t in timestamp_arr:
        t = t-timestamp_arr[0]
    iter_arr = np.linspace(1, len(t_arr), len(t_arr))
    plot(iter_arr, residual_arr, 'iteration', 'residual')
    plot(iter_arr, timestamp_arr, 'iteration', 'walltime')
    plot(timestamp_arr, residual_arr, 'walltime', 'residual')
    plot(timestamp_arr, t_arr, 'walltime', 'time')


if(args.live):
    last_modified = os.stat(args.file).st_mtime
    while True:
        time.sleep(10)
        if(last_modified != os.stat(args.file).st_mtime):
            create_plots()
            last_modified = os.stat(args.file).st_mtime
else:
    create_plots()
