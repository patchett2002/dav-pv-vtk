from vtkmodules.vtkCommonDataModel import vtkDataSet
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa

# new module for ParaView-specific decorators.
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

from paraview import vtk
import numpy as np # needed for interpolation and pi
import sys          # needed to get a command line argument
import math         # needed to calculate the sine and cosine of the longitude

@smproxy.filter(label="VTS to Stereographic VTS")
@smproperty.input(name="Input")
class VTStoStereographicVTS(VTKPythonAlgorithmBase):

    # Earth radius from https://github.com/nsidc/polarstereo-lonlat-convert-py
    #EARTH_RADIUS_KM = 6378.137

    EARTH_RADIUS_KM = 6378.137
    EARTH_ECCENTRICITY = 0.01671
    TRUE_SCALE_LATITUDE = 90.0
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)
        self.hemisphereToProject = ""
        # Create a list of the two hemispheres that can be shown in a stereographic projection
        self._availableHemispheres = ["Northern Hemisphere", "Southern Hemisphere"]

    def FillInputPortInformation(self, port, info):
        info.Set(vtk.vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
        return 1

    def FillOutputPortInformation(self, port, info):
        info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkStructuredGrid")
        return 1

    @smproperty.stringvector(name="AvailableHemispheres", information_only="1")
    def GetAvailableArrays(self):
        return(self._availableHemispheres)

    @smproperty.stringvector(name="Hemisphere", number_of_elements="1")
    @smproperty.xml("""
        <StringListDomain name="hemisphereChoice">
            <RequiredProperties>
                <Property name="AvailableHemispheres"
                    function="HemisphereSelection"/>
            </RequiredProperties>
        </StringListDomain>
        """)
    def SetHemisphere(self, val):
        print("Setting ", val)
        self.hemisphereToProject = val
        self.Modified()
    # Taken from https://github.com/nsidc/polarstereo-lonlat-convert-py/blob/main/polar_convert/polar_convert.py
    def polar_lonlat_to_xy(self, longitude, latitude, true_scale_lat, re, e, hemisphere):
        """Convert from geodetic longitude and latitude to Polar Stereographic
        (X, Y) coordinates in km.
        Args:
            longitude (float): longitude or longitude array in degrees
            latitude (float): latitude or latitude array in degrees (positive)
            true_scale_lat (float): true-scale latitude in degrees
            re (float): Earth radius in km
            e (float): Earth eccentricity
            hemisphere ('north' or 'south'): Northern or Southern hemisphere
        Returns:
            If longitude and latitude are scalars then the result is a
            two-element list containing [X, Y] in km.
            If longitude and latitude are numpy arrays then the result will be a
            two-element list where the first element is a numpy array containing
            the X coordinates and the second element is a numpy array containing
            the Y coordinates.
        """
        #hemisphere = validate_hemisphere(hemisphere)
        #hemi_direction = _hemi_direction(hemisphere)
        # We are projecting the northern hemisphere, so the hemi_direction is 1
        # (otherwise it is -1 for the southern hemisphere)
        #hemi_direction = 1
        if (hemisphere == "Northern Hemisphere"):
            hemi_direction = 1
        else:
            hemi_direction = -1
        lat = np.abs(latitude) * np.pi / 180
        lon = longitude * np.pi / 180
        slat = true_scale_lat * np.pi / 180
        e2 = e * e
        # Snyder (1987) p. 161 Eqn 15-9
        t = np.tan(np.pi / 4 - lat / 2) / \
            ((1 - e * np.sin(lat)) / (1 + e * np.sin(lat))) ** (e / 2)
        if np.abs(90 - true_scale_lat) < 1e-5:
            # Snyder (1987) p. 161 Eqn 21-33
            rho = 2 * re * t / np.sqrt((1 + e) ** (1 + e) * (1 - e) ** (1 - e))
        else:
            # Snyder (1987) p. 161 Eqn 21-34
            tc = np.tan(np.pi / 4 - slat / 2) / \
                ((1 - e * np.sin(slat)) / (1 + e * np.sin(slat))) ** (e / 2)
            mc = np.cos(slat) / np.sqrt(1 - e2 * (np.sin(slat) ** 2))
            rho = re * mc * t / tc
        x = rho * hemi_direction * np.sin(hemi_direction * lon)
        y = -rho * hemi_direction * np.cos(hemi_direction * lon)
        #print(x, y)
        return [x, y]
    def RequestData(self, request, inInfo, outInfo):
        print("Request data")
        EARTH_RADIUS_KM = 6378.137
        EARTH_ECCENTRICITY = 0.01671
        #TRUE_SCALE_LATITUDE = 90.0
        if (self.hemisphereToProject == "Northern Hemisphere"):
            TRUE_SCALE_LATITUDE = 90.0
        else:
            TRUE_SCALE_LATITUDE = -90.0
        #TRUE_SCALE_LATITUDE = self.trueScaleLatitude
        # get the input data set
        inputDataSet0 = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))
        newPoints = vtk.vtkPoints()
        numPoints = inputDataSet0.GetNumberOfPoints()
        num_arrays = inputDataSet0.GetPointData().GetNumberOfArrays()
        print("Number of arrays:", num_arrays)
        # Get the dimensions of the input dataset
        input_dimensions = inputDataSet0.GetDimensions()
        print(input_dimensions[0])   # should be 1025
        print(input_dimensions[1])   # should be 512
        print(input_dimensions[2])   # should be 1
        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]
        # populate newPoints with stereographic points
        for i in range(0, numPoints):
            coord = inputDataSet0.GetPoint(i)
            x0, y0, z0 = coord[:3]
            if (self.hemisphereToProject == "Northern Hemisphere"):
                if y0 >= 0:
                    x,y = self.polar_lonlat_to_xy(x0, y0, TRUE_SCALE_LATITUDE, EARTH_RADIUS_KM, EARTH_ECCENTRICITY, self.hemisphereToProject)
            else:
                if y0 <= 0:
                    x,y = self.polar_lonlat_to_xy(x0, y0, TRUE_SCALE_LATITUDE, EARTH_RADIUS_KM, EARTH_ECCENTRICITY, self.hemisphereToProject)
            newPoints.InsertPoint(i,x,y,0)
        # Get the dimensions of the newPoints vtk points data set that only
        # contains one hemisphere's points
        num_new_points = newPoints.GetNumberOfPoints()
        # Since the x-dimension is the same as the original input data set and the
        # z-dimension is flat with only one layer of points, then I can find the
        # y-dimension of newPoints by dividing the number of points in the newPoints
        # data set by the original x-dimension
        new_y_dimension = num_new_points / x_dimension
        # Create the output data set
        outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)
        # Loop through each of the scalar arrays in the dataset
        # and add the array to the output dataset
        for j in range(0, num_arrays):
            ivals = inputDataSet0.GetPointData().GetArray(j)
            ca = vtk.vtkFloatArray()
            ca.SetName(ivals.GetName())
            ca.SetNumberOfComponents(1)
            ca.SetNumberOfTuples(num_new_points)
            # Print out the number of points
            print("Number of points: ", newPoints.GetNumberOfPoints())
            # add to output
            #outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)
            print("Output data set: ", outputDataSet)
            #add the new array to the output
            outputDataSet.GetPointData().AddArray(ca)
            #copy the values over element by element
            indx = 0    # Index for the each scalar value in the ca vtk float array
            for i in range(0, numPoints):
                coord = inputDataSet0.GetPoint(i)
                x0, y0, z0 = coord[:3]
                if (self.hemisphereToProject == "Northern Hemisphere"):
                    if y0 >= 0:
                        ca.SetValue(indx, ivals.GetValue(i))
                        indx = indx + 1
                else:
                    if y0 <= 0:
                        ca.SetValue(indx, ivals.GetValue(i))
                        indx = indx + 1
            # Try printing out the value at the 1025th point
            print("Value at point 1024: ", ivals.GetValue(1024))
            print("Value at point 1024 in ca: ", ca.GetValue(1024))
        outputDataSet.SetDimensions(x_dimension,new_y_dimension,z_dimension)
        return 1
