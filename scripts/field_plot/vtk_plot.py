import vtk

file_name = 'vtks/euler2D_nacaCoarse_caosipm_n10_nq15_s05_aoa.vtk'
field_name = "E(ρ)"
#field_name = "Var(ρ)"

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_name)
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

mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(output)
mapper.ScalarVisibilityOn()
mapper.SetColorModeToMapScalars()
mapper.SetLookupTable(lut)
mapper.SetScalarModeToUsePointFieldData()
mapper.SelectColorArray(field_name)
mapper.SetScalarRange(output.GetPointData().GetArray(field_name).GetRange())

scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(mapper.GetLookupTable())
#scalarBar.SetTitle(field_name)
scalarBar.SetOrientationToHorizontal()
scalarBar.SetPosition(0.1,0)
scalarBar.SetWidth(0.8)
scalarBar.SetHeight(0.05)
scalarBar.SetNumberOfLabels(4)

actor = vtk.vtkActor()
actor.SetMapper(mapper)
#actor.GetProperty().EdgeVisibilityOn()
#actor.GetProperty().SetLineWidth(1.0)

camera = vtk.vtkCamera()
camera.SetPosition(0.5, 0.3, 3);
camera.SetFocalPoint(0.5, 0.3, 0);

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.AddActor2D(scalarBar)
renderer.SetBackground(1, 1, 1) 
renderer.SetActiveCamera(camera)

render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window.SetSize(800,800)
render_window.Render()

windowToImageFilter = vtk.vtkWindowToImageFilter()
windowToImageFilter.SetInput(render_window)
#windowToImageFilter.SetInputBufferTypeToRGBA()
windowToImageFilter.ReadFrontBufferOff()
windowToImageFilter.Update()
  
writer = vtk.vtkPNGWriter()
writer.SetFileName("screenshot.png");
writer.SetInputConnection(windowToImageFilter.GetOutputPort());
writer.Write();