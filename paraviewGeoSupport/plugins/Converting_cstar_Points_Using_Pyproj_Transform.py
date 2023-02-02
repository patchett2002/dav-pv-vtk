import vtk
import pyproj
from vtkmodules.vtkGeovisCore import vtkGeoProjection, vtkGeoTransform

# Get the point coordinates from the cstar data set
cstar = "/Users/linneapalmstrom/Documents/LANL ISTI Projects/Combined_dataset_files/cstar_-180_to_180_lon.vts"

reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName(cstar)
reader.Update()
cstarPts = reader.GetOutput()
inputPts = cstarPts.GetPoints()

numPoints = inputPts.GetNumberOfPoints()
print(numPoints)

print("Transforming the latitude and longitude points using pyproj.transform (with the original z-coordinates set to zero):")

# Print out what the transformed coordinates are using pyproj.transform
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='lonlat', ellps='WGS84', datum='WGS84')

for i in range(0,numPoints):
    coord = inputPts.GetPoint(i)
    x,y,z = coord[:3]
    newX,newY,newZ = pyproj.transform(lla, ecef, x, y, 0, radians=False)
    print("Original point: ",x,",",y,sep="")
    print("Transformed point: ",newX,",",newY,",",newZ,sep="")
