import vtk
import pyproj
from vtkmodules.vtkGeovisCore import vtkGeoProjection, vtkGeoTransform

# Point coordinates to test pyproj.transform and vtkGeoTransform
#x0,y0 = 50,20
x0,y0 = 0.661,-2.137
x1,y1 = 50,50
x2,y2 = 50,70
x3,y3 = 80,20
x4,y4 = 80,50
x5,y5 = 80,70

print("Transforming the latitude and longitude points using pyproj.transform (with the original z-coordinates set to zero):")

# Print out what the transformed coordinates are using pyproj.transform
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='lonlat', ellps='WGS84', datum='WGS84')
newX0,newY0,newZ0 = pyproj.transform(lla, ecef, x0, y0, -20.531, radians=False)

print("Original first point: ",x0,",",y0,sep="")
print("Transformed first point: ",newX0,",",newY0,",",newZ0,sep="")

newX1,newY1,newZ1 = pyproj.transform(lla, ecef, x1, y1, 0, radians=False)

print("Original second point: ",x1,",",y1,sep="")
print("Transformed second point: ",newX1,",",newY1,",",newZ1,sep="")

newX2,newY2,newZ2 = pyproj.transform(lla, ecef, x2, y2, 0, radians=False)

print("Original third point: ",x2,",",y2,sep="")
print("Transformed third point: ",newX2,",",newY2,",",newZ2,sep="")

newX3,newY3,newZ3 = pyproj.transform(lla, ecef, x3, y3, 0, radians=False)

print("Original fourth point: ",x3,",",y3,sep="")
print("Transformed fourth point: ",newX3,",",newY3,",",newZ3,sep="")

newX4,newY4,newZ4 = pyproj.transform(lla, ecef, x4, y4, 0, radians=False)

print("Original fifth point: ",x4,",",y4,sep="")
print("Transformed fifth point: ",newX4,",",newY4,",",newZ4,sep="")

newX5,newY5,newZ5 = pyproj.transform(lla, ecef, x5, y5, 0, radians=False)

print("Original sixth point: ",x5,",",y5,sep="")
print("Transformed sixth point: ",newX5,",",newY5,",",newZ5,sep="")


print("\n\n")
print("Transforming the latitude and longitude points using vtkGeoTransform (with the original z-coordinates set to zero):")

# Print out what the transformed coordinates are using vtkGeoTransform
geo = vtkGeoTransform()
ps = vtkGeoProjection()
pd = vtkGeoProjection()
pd.SetPROJ4String("+proj=cart +ellps=WGS84")
ps.SetPROJ4String("+proj=lonlat +ellps=WGS84")
geo.SetSourceProjection(ps)
geo.SetDestinationProjection(pd)

# Create a vtkPoints array with the six original point coordinates to transform
points = vtk.vtkPoints()
pointsList = [(x0,y0),(x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5)]
for i in range(6):
    x,y = pointsList[i]
    points.InsertNextPoint(x,y,0)

newPoints = vtk.vtkPoints()

geo.TransformPoints(points, newPoints)

newX0,newY0,newZ0 = newPoints.GetPoint(0)
newX1,newY1,newZ1 = newPoints.GetPoint(1)
newX2,newY2,newZ2 = newPoints.GetPoint(2)
newX3,newY3,newZ3 = newPoints.GetPoint(3)
newX4,newY4,newZ4 = newPoints.GetPoint(4)
newX5,newY5,newZ5 = newPoints.GetPoint(5)

print("Original first point: ",x0,",",y0,sep="")
print("Transformed first point: ",newX0,",",newY0,",",newZ0,sep="")

print("Original second point: ",x1,",",y1,sep="")
print("Transformed second point: ",newX1,",",newY1,",",newZ1,sep="")

print("Original third point: ",x2,",",y2,sep="")
print("Transformed third point: ",newX2,",",newY2,",",newZ2,sep="")

print("Original fourth point: ",x3,",",y3,sep="")
print("Transformed fourth point: ",newX3,",",newY3,",",newZ3,sep="")

print("Original fifth point: ",x4,",",y4,sep="")
print("Transformed fifth point: ",newX4,",",newY4,",",newZ4,sep="")

print("Original sixth point: ",x5,",",y5,sep="")
print("Transformed sixth point: ",newX5,",",newY5,",",newZ5,sep="")
