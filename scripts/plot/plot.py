import matplotlib.pyplot as plt
import argparse
import numpy as np
import datetime as dt
import os
import pandas as pd
import re
from io import StringIO
import vtk

def extract_intervals(sequence, num):
    length = float(len(sequence))
    for i in range(num):
        yield sequence[int(np.ceil(i * length / num))]

def time_to_float(time):
    t = dt.datetime.strptime(time, "%Y-%m-%d %H:%M:%S.%f")
    epoch = dt.datetime.utcfromtimestamp(0)
    return (t-epoch).total_seconds()

def sortByMethod(dir, configFiles, files):
    labels = []
    for file in configFiles:
        labels.append(create_label(dir,file))
    labels = [l.upper() for l in labels]
    order = np.argsort(labels)
    return [files[idx] for idx in order]

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
        label = ''
        restartFile = False
        while section_ctr<3:
            current_line = lines[line_ctr]
            if "#" not in current_line:
                if "=====" in current_line:
                    section_ctr = section_ctr + 1
                elif "closure" in current_line:
                    closure = current_line.split("=")[1].strip().strip('"').strip('\n')
                elif "moments" in current_line:
                    nMoments = int(current_line.split("=")[1])
                elif "quadPoints" in current_line:
                    nQuad = int(current_line.split("=")[1])
                elif "maxIterations" in current_line:
                    maxIter = int(current_line.split("=")[1])
                elif "restartFile" in current_line:
                    restartFile = True
            line_ctr = line_ctr + 1
        if closure == 'StochasticGalerkin' and nQuad == nMoments:
            label += 'SC'
            label += '$^{'+str(nQuad)+'}$'
        elif closure == 'StochasticGalerkin':
            label += 'SG'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'Euler2D' and restartFile and maxIter == 1:
            label += 'caos-IPM'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'Euler2D' and restartFile:
            label += 'ca-IPM'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'Euler2D' and maxIter == 1:
            label += 'os-IPM'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'Euler2D':
            label += 'IPM'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'L2Filter':
            label += 'L2Filter'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'LassoFilter':
            label += 'LassoFilter'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        elif closure == 'BoundedBarrier':
            label += 'BoundedBarrier'
            label += '$^{'+str(nQuad)+'}'
            label += '_{'+str(nMoments)+'}$'
        else:
            label += 'unkown'
        return label


def parse_error_logfile(dir, file):
    with open(dir+'/'+file, 'r') as content:
        print("Parsing:\t" + dir + '/' + file)
        lines = content.readlines()
        lines = extract_intervals(lines,1000)
        header = ['date', 'time', 'delim', 'rho', 'rhoU_x', 'rhoU_y', 'rhoE']
        df = pd.read_csv(StringIO("\n".join(lines)), header=None, names=header, delim_whitespace=True)
        df = df.drop(columns=['delim'])
        df = df.dropna()
        t0 = time_to_float(df['date'][0]+' '+df['time'][0])
        df['runtime'] = np.zeros(len(df['time']))
        for index, row in df.iterrows():
            df.loc[index,'runtime'] = time_to_float(df['date'][index]+' '+df['time'][index]) - t0
        df = df.drop(columns=['date','time'])
        return(df)

def parse_logfile(dir, file):
    with open(dir+'/'+file, 'r') as content:
        print("Parsing:\t" + dir + '/' + file)
        lines = content.readlines()
        pattern = re.compile("(\d{4}\-\d{2}\-\d{2}\s\d{2}\:\d{2}\:\d{2}.\d{6}\s\|){1}(\s+[+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+-]?\d+)?){3}")
        valid_lines = []
        for line in lines:
            if pattern.match(line):
                valid_lines.append(line)
        lines = extract_intervals(valid_lines,1000)
        header = ['date', 'time', 'delim', 't', 'residual_scaled', 'residual_abs']
        df = pd.read_csv(StringIO("\n".join(lines)), header=None, names=header, delim_whitespace = True)
        df = df.drop(columns=['delim'])
        df = df.dropna()
        t0 = time_to_float(df['date'][0]+' '+df['time'][0])
        df['runtime'] = np.zeros(len(df['time']))
        for index, _ in df.iterrows():
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

def create_error_plot(dir, configFiles, files, title, type):
    outputdir = 'plots/error'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
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
            df_prev = parse_error_logfile(dir, prevFileName+'_'+ext)
            df = parse_error_logfile(dir, files[f])
            df['runtime'] += df_prev['runtime'].iloc[-1]
            data.append(pd.concat([df_prev, df]))
        else:
            data.append(parse_error_logfile(dir, files[f]))
    for s in range(nStates):
        plt.cla()
        plt.clf()
        plt.xlabel('Time [s]')
        plt.ylabel('Error')
        plt.title(title + ' error ' + type + '[' + stateLabel[s] + ']')
        for d in data:
            plt.semilogy(d['runtime'], d[stateLabel[s]])
        plt.legend(label)
        plt.savefig(outputdir + '/' + title + '_error_' + type + '[' + stateLabel[s] + '].pdf')
        print('Plot created:\t' + outputdir + '/' + title + '_error_' + type + '[' + stateLabel[s] + '].pdf')

def create_convergence_plots(dir, configFiles):
    outputdir = 'plots/convergence'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    data = []
    labels = []
    for file in configFiles:
        data.append(parse_logfile(dir, file))
        labels.append(create_label(dir, file))
    plotname = 'convergence_runtime_residual.pdf'
    plt.cla()
    plt.clf()
    for df in data:
        plt.semilogy(df['runtime'],df['residual_scaled'])
    plt.xlabel('Runtime [s]')
    plt.ylabel('Residual')
    plt.legend(labels)
    plt.savefig(outputdir + '/' + plotname)
    print('Plot created:\t' + outputdir + '/' + plotname)

    plotname = 'convergence_time_residual.pdf'
    plt.cla()
    plt.clf()
    for df in data:
        plt.semilogy(df['t'],df['residual_scaled'])
    plt.xlabel('Time [s]')
    plt.ylabel('Residual')
    plt.legend(labels)
    plt.savefig(outputdir + '/' + plotname)
    print('Plot created:\t' + outputdir + '/' + plotname)

def create_vtk_plots(dir, vtkFiles, rescale):
    outputdir = 'plots/vtk'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    fields = ["E(ρ)", "E(ρU)", "E(ρE)", "Var(ρ)", "Var(ρU)", "Var(ρE)"]
    vtk_label = [''.join('${0}$'.format(f)).replace('ρ', '\\rho ') for f in fields]
    vtk_label_dict = {}
    for f in range(len(fields)):
        vtk_label_dict[fields[f]] = vtk_label[f]
    minVal = {}
    maxVal = {}
    minErr = {}
    maxErr = {}
    for field_name in fields:
        minVal[field_name] = -np.inf
        maxVal[field_name] = np.inf
        minErr[field_name] = -np.inf
        maxErr[field_name] = np.inf
    if rescale:
        for file in vtkFiles:
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(dir+'/'+file)
            reader.ReadAllScalarsOn()
            reader.ReadAllVectorsOn()
            reader.Update()
            output = reader.GetOutput()
            for field_name in fields:
                valRange = output.GetCellData().GetArray(field_name).GetRange()
                if file.endswith('_errors.vtk'):
                    minErr[field_name] = np.maximum(minErr[field_name], valRange[0])
                    maxErr[field_name] = np.minimum(maxErr[field_name], valRange[1])
                else:
                    minVal[field_name] = np.maximum(minVal[field_name], valRange[0])
                    maxVal[field_name] = np.minimum(maxVal[field_name], valRange[1])

    for file in vtkFiles:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(dir+'/'+file)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        output = reader.GetOutput()

        numvals = 1024
        ctf = vtk.vtkColorTransferFunction()
        ctf.SetColorSpaceToRGB()
        ctf.AddRGBPoint(-1, 1, 1, 0)
        ctf.AddRGBPoint(0, 0, 0, 0)
        ctf.AddRGBPoint(0.293069/1.72393, 0, 0, 1)
        ctf.AddRGBPoint(0.586138/1.72393, 0, 1, 1)
        ctf.AddRGBPoint(0.861967/1.72393, 0, 1, 0)
        ctf.AddRGBPoint(1.155040/1.72393, 1, 1, 0)
        ctf.AddRGBPoint(1.448100/1.72393, 1, 0, 0)
        ctf.AddRGBPoint(1.723930/1.72393, 0.87843, 0, 1)
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(numvals)
        lut.Build()
        for i in range(0,numvals):
            rgb = list(ctf.GetColor(float(i)/numvals))+[1]
            lut.SetTableValue(i,rgb)

        camera = vtk.vtkCamera()
        camera.SetPosition(0.5, 0.3, 3);
        camera.SetFocalPoint(0.5, 0.3, 0);

        for field_name in fields:
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(output)
            mapper.ScalarVisibilityOn()
            mapper.SetColorModeToMapScalars()
            mapper.SetLookupTable(lut)
            mapper.SetScalarModeToUsePointFieldData()
            mapper.SelectColorArray(field_name)
            if rescale:
                if file.endswith('_errors.vtk'):
                    mapper.SetScalarRange(minErr[field_name], maxErr[field_name])
                else:
                    mapper.SetScalarRange(minVal[field_name], maxVal[field_name])
            else:
                mapper.SetScalarRange(output.GetCellData().GetArray(field_name).GetRange())

            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(mapper.GetLookupTable())
            scalarBar.SetTitle(vtk_label_dict[field_name])
            scalarBar.SetOrientationToHorizontal()
            scalarBar.SetPosition(0.1,0)
            scalarBar.SetWidth(0.8)
            scalarBar.SetHeight(0.1)
            scalarBar.SetNumberOfLabels(4)
            scalarBar.SetMaximumNumberOfColors(numvals)
            scalarBar.SetTitleRatio(0.6)
            labelprop = scalarBar.GetLabelTextProperty()
            labelprop.ShadowOff()
            labelprop.BoldOff()
            if 'E(' in field_name and not file.endswith('_errors.vtk'):
                labelprop.SetColor(0,0,0)
                titleprop = scalarBar.GetTitleTextProperty()
                titleprop.SetColor(0,0,0)
                scalarBar.SetTitleTextProperty(titleprop)
            scalarBar.SetLabelTextProperty(labelprop)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            #actor.GetProperty().EdgeVisibilityOn()
            #actor.GetProperty().SetLineWidth(1.0)

            renderer = vtk.vtkRenderer()
            renderer.AddActor(actor)
            renderer.AddActor2D(scalarBar)
            renderer.UseFXAAOn()
            renderer.SetBackground(1, 1, 1)
            renderer.SetActiveCamera(camera)

            render_window = vtk.vtkRenderWindow()
            render_window.SetOffScreenRendering(True)
            render_window.AddRenderer(renderer)
            render_window.SetSize(800,800)
            render_window.Render()

            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(render_window)
            windowToImageFilter.ReadFrontBufferOff()
            windowToImageFilter.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(outputdir+'/'+os.path.splitext(file)[0]+"_"+field_name+".png")
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()
            print('Plot created:\t' + outputdir + '/' + os.path.splitext(file)[0]+"_"+field_name+".png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--logdir", "-l", type=str, required=False, default='logs')
    parser.add_argument("--vtkdir", "-v", type=str, required=False, default='.')
    parser.add_argument("--rescale", "-r", type=bool, required=False, default=True)
    args = parser.parse_args()
    logdir = args.logdir
    vtkdir = args.vtkdir
    rescale = args.rescale
    logdirfiles = os.listdir(logdir)
    vtkdirfiles = os.listdir(vtkdir)
    ignoredFiles = []
    configFiles = []
    vtkFiles = []
    L1ErrMeanFiles = []
    L2ErrMeanFiles = []
    LInfErrMeanFiles = []
    L1ErrVarFiles = []
    L2ErrVarFiles = []
    LInfErrVarFiles = []
    for file in vtkdirfiles:
        if file.endswith('.vtk'):
            vtkFiles.append(file)
    for file in logdirfiles:
        isRestartFile, prevFileName = checkForRestartFile(logdir, file)
        if isRestartFile:
            if prevFileName not in ignoredFiles: ignoredFiles.append(prevFileName)
            if prevFileName not in logdirfiles: print("WARNING: restart file found but previous file was not found")
    for file in logdirfiles:
        for ifile in ignoredFiles:
            if file.startswith(ifile):
                if file not in ignoredFiles: ignoredFiles.append(file)
    for ifile in ignoredFiles:
        logdirfiles.remove(ifile)
    for file in logdirfiles:
        if re.search(r'\d+$', file) is not None:
            configFiles.append(file)
        elif file.endswith('_L1ErrorMean'):
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
    configFiles.sort()
    L1ErrMeanFiles.sort()
    L2ErrMeanFiles.sort()
    LInfErrMeanFiles.sort()
    L1ErrVarFiles.sort()
    L2ErrVarFiles.sort()
    LInfErrVarFiles.sort()
    if any(L1ErrMeanFiles): L1ErrMeanFiles = sortByMethod(logdir, configFiles, L1ErrMeanFiles)
    if any(L2ErrMeanFiles): L2ErrMeanFiles = sortByMethod(logdir, configFiles, L2ErrMeanFiles)
    if any(LInfErrMeanFiles): LInfErrMeanFiles = sortByMethod(logdir, configFiles, LInfErrMeanFiles)
    if any(L1ErrVarFiles): L1ErrVarFiles = sortByMethod(logdir, configFiles, L1ErrVarFiles)
    if any(L2ErrVarFiles): L2ErrVarFiles = sortByMethod(logdir, configFiles, L2ErrVarFiles)
    if any(LInfErrVarFiles): LInfErrVarFiles = sortByMethod(logdir, configFiles, LInfErrVarFiles)
    if any(configFiles): configFiles = sortByMethod(logdir, configFiles, configFiles)

    if any(configFiles): create_convergence_plots(logdir, configFiles)
    if any(vtkFiles): create_vtk_plots(vtkdir, vtkFiles, rescale)
    if any(L1ErrMeanFiles): create_error_plot(logdir, configFiles, L1ErrMeanFiles, 'L1', 'E')
    if any(L2ErrMeanFiles): create_error_plot(logdir, configFiles, L2ErrMeanFiles, 'L2', 'E')
    if any(LInfErrMeanFiles): create_error_plot(logdir, configFiles, LInfErrMeanFiles, 'L-infinity', 'E')
    if any(L1ErrVarFiles): create_error_plot(logdir, configFiles, L1ErrVarFiles, 'L1', 'Var')
    if any(L2ErrVarFiles): create_error_plot(logdir, configFiles, L2ErrVarFiles, 'L2', 'Var')
    if any(LInfErrVarFiles): create_error_plot(logdir, configFiles, LInfErrVarFiles, 'L-infinity', 'Var')
