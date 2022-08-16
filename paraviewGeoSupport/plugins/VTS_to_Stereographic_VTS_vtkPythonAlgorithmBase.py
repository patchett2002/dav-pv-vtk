from vtkmodules.vtkCommonDataModel import vtkDataSet
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
# new module for ParaView-specific decorators.
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain
from paraview import vtk
import numpy as np # needed for interpolation and pi
import sys          # needed to get a command line argument
import math         # needed to calculate the sine and cosine of the longitude
import copy
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
        # Set the default columnAtEnd value to 0
        self.columnAtEnd = 0
        # Create a list of the two hemispheres that can be shown in a stereographic projection
        self._availableHemispheres = ["Northern Hemisphere", "Southern Hemisphere"]
        # Initialize the input dimensions
        self.input_x_dimension = 0
        self.input_y_dimension = 0
        self.input_z_dimension = 0
        # Initialize the output dimensions
        self.output_x_dimension = 0
        self.output_y_dimension = 0
        self.output_z_dimension = 0
    def FillInputPortInformation(self, port, info):
        info.Set(vtk.vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
        return 1
    def FillOutputPortInformation(self, port, info):
        info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkStructuredGrid")
        return 1

    @smproperty.stringvector(name="AvailableHemispheres", information_only="1")
    def GetAvailableHemispheres(self):
        return(self._availableHemispheres)

    @smproperty.stringvector(name="Hemisphere", number_of_elements="1")
    @smdomain.xml("""
        <StringListDomain name="hemisphereChoice">
            <RequiredProperties>
                <Property name="AvailableHemispheres"
                    function="AvailableHemispheres"/>
            </RequiredProperties>
        </StringListDomain>
        """)
    def SetHemisphere(self, val):
        print("Setting ", val)
        self.hemisphereToProject = val
        self.Modified()
    @smproperty.xml("""
        <IntVectorProperty name="AddColumnToOneEnd"
            number_of_elements="1"
            default_values="0"
            command="SetColumnAtEnd">
            <BooleanDomain name="bool" />
            <Documentation>If on, copies the values in the first column of the
            input data set and adds that column to the end of the data set.
            </Documentation>
        </IntVectorProperty>""")
    def SetColumnAtEnd(self, x):
        self.columnAtEnd = x
        print("Set Column At End: ", self.columnAtEnd)
        self.Modified()

    def GetColumnAtEnd(self):
        print("Get Column At End: ", self.columnAtEnd)
        return self.columnAtEnd
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
            val_to_check_true_scale_lat = 90 - true_scale_lat
        else:
            hemi_direction = -1
            val_to_check_true_scale_lat = 90 + true_scale_lat
        lat = np.abs(latitude) * np.pi / 180
        lon = longitude * np.pi / 180
        slat = true_scale_lat * np.pi / 180
        e2 = e * e
        # Snyder (1987) p. 161 Eqn 15-9
        t = np.tan(np.pi / 4 - lat / 2) / \
            ((1 - e * np.sin(lat)) / (1 + e * np.sin(lat))) ** (e / 2)
        # Check that the true_scale_lat is about 90 degrees or -90 degrees depending
        # on whether you are projecting the northern or southern hemisphere
        #val_to_check_true_scale_lat = 50
        #if (hemisphere == "Northern Hemisphere")
        #if np.abs(90 - true_scale_lat) < 1e-5:
        if np.abs(val_to_check_true_scale_lat) < 1e-5:
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
    # This function will take in an array and find and return one higher
    # longitude value (as a whole degree) than the last element in the array
    def FindNextHighestLonValue(self, arr):
        # Get the highest value in the array arr and truncate the decimal places
        max_whole_lon_val = math.trunc(max(arr))

        # Return the truncated highest longitude value in the array plus one
        #return max_whole_lon_val + 1
        # for some reason, returned 181, so I made it return max_whole_lon_val
        # without the plus one for now
        return max_whole_lon_val


    def AddColumnToEnd(self, inputDataSet):
        # Create the new input data set
        newInputDataSet = vtk.vtkStructuredGrid()

        newPoints = vtk.vtkPoints()
        numPoints = inputDataSet.GetNumberOfPoints()
        print("input data Number of points", numPoints)

        num_arrays = inputDataSet.GetPointData().GetNumberOfArrays()
        print("Number of arrays:", num_arrays)

        # Get the dimensions of the input dataset
        input_dimensions = inputDataSet.GetDimensions()

        print("Dimensions:")
        print(input_dimensions[0])   # should be 1025
        print(input_dimensions[1])   # should be 512
        print(input_dimensions[2])   # should be 1
        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]

        original_lon_points = []
        for i in range(0, numPoints):
            original_coord = inputDataSet.GetPoint(i)
            x,y = original_coord[:2]
            original_lon_points.append(x)

        next_lon_value = self.FindNextHighestLonValue(original_lon_points)
        # Should print 360 for the next_lon_value
        print("Next highest longitude value: ", next_lon_value)

        latPoints = []
        lonPoints = []

        for i in range(0, numPoints):
            coord = inputDataSet.GetPoint(i)
            x0, y0 = coord[:2]
            latPoints.append(y0)
            lonPoints.append(x0)

        latPoints = copy.deepcopy(np.array(latPoints))
        lonPoints = copy.deepcopy(np.array(lonPoints))

        # Test out adding a column of y_dimension values of 360 degrees
        # to my longitude array and add a copy of the first
        # column of my latitude array to the end of the array


        # Add the first column of the latitude array to the end
        # of the latitude array
        latPoints=latPoints.reshape(y_dimension,x_dimension)
        ALat = copy.deepcopy(latPoints[0:y_dimension,:])
        BLat = copy.deepcopy(latPoints[0:y_dimension,0:1])
        NewLatArray = np.hstack((ALat, BLat))


        latPoints=latPoints.reshape(y_dimension * x_dimension,)
        NewLatArray=NewLatArray.reshape(y_dimension * (x_dimension + 1),)


        # Add a column of 360 degrees to the end of my longitude
        # array
        #lonPoints=lonPoints.reshape(512,1024)
        lonPoints=lonPoints.reshape(y_dimension,x_dimension)
        #ALon = lonPoints[0:512,:]
        #ALon = copy.deepcopy(lonPoints[0:512,:])
        ALon = copy.deepcopy(lonPoints[0:y_dimension,:])
        #BLon = np.full((512, 1), 360.0)
        #BLon = np.full((y_dimension, 1), 360.0)
        BLon = np.full((y_dimension, 1), next_lon_value)
        NewLonArray = np.hstack((ALon, BLon))


        #lonPoints=lonPoints.reshape(524288,)
        lonPoints=lonPoints.reshape(y_dimension * x_dimension,)
        #NewLonArray=NewLonArray.reshape(524800,)
        NewLonArray=NewLonArray.reshape(y_dimension * (x_dimension + 1),)


        # Check some of the points in the NewLonArray
        #print(NewLonArray[1023])
        #print(NewLonArray[1024])
        print(NewLonArray[x_dimension])
        #print(NewLonArray.shape)

        #for i in range(0,numPoints):
        #for i in range(0,524800):
        for i in range(0,y_dimension * (x_dimension + 1)):
            #mynewcoord = (lonPoints[i], latPoints[i], 0)
            mynewcoord = (NewLonArray[i], NewLatArray[i], 0)
            newPoints.InsertNextPoint(mynewcoord)


        # Loop through each of the scalar arrays in the dataset, add the first
        # column to the end and add the modified array to the new input dataset
        for j in range(0, num_arrays):
            dat = vtk.vtkFloatArray()
            rawdat = []
            ivals = inputDataSet.GetPointData().GetArray(j)
            # populates rawdat list with input points
            for i in range(0, numPoints):
                rawdat.append(ivals.GetValue(i))
            print("rawdat number of points", len(rawdat))
            # Input the values from the altered numpy array into a new vtk
            # float array
            dat.SetName(ivals.GetName())
            newNumPoints = 0
            # Test out the reshape function to see if that is messing up
            # my new input (it seems to work fine after running this code)
            #nprawdat = np.array(rawdat)
            nprawdat = copy.deepcopy(np.array(rawdat))
            nprawdat=nprawdat.reshape(y_dimension,x_dimension)
            nprawdat=nprawdat.reshape(y_dimension * x_dimension,)


            # Test out adding a new column to my rawdat array made up of
            # 512 ones to see if it messes up my code
            nprawdat=nprawdat.reshape(y_dimension,x_dimension)
            A = copy.deepcopy(nprawdat[0:y_dimension,:])
            B = nprawdat[0:y_dimension,0:1]
            newnprawdat = np.hstack((A, B))


            nprawdat=nprawdat.reshape(y_dimension * x_dimension,)
            newnprawdat=newnprawdat.reshape(y_dimension * (x_dimension + 1),)

            # populate newpoints vtkFloatArray
            for i in newnprawdat:
                dat.InsertNextValue(i)


            # Add the array to the new input dataset
            newInputDataSet.GetPointData().AddArray(dat)


        # Check some of the points in the vtk point array
        #print(newPoints.GetPoint(1024))
        print(newPoints.GetPoint(x_dimension))
        #print(newPoints.GetPoint(524799))


        #newInputDataSet.GetPointData().AddArray(dat)
        #newInputDataSet.SetDimensions(1024,512,1)
        #newInputDataSet.SetDimensions(1025,512,1)
        newInputDataSet.SetDimensions((x_dimension + 1),y_dimension,z_dimension)
        newInputDataSet.SetPoints(newPoints)

        return newInputDataSet
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
        # Check that the Add Column To One End checkbox is checked, and if it
        # is, then alter the input data set so that it has a copy of the first
        # column of values added to the end
        if (self.GetColumnAtEnd() == True):
            newDataSet = self.AddColumnToEnd(inputDataSet0)
        else:
            #newDataSet = copy.deepcopy(inputDataSet0)
            newDataSet = inputDataSet0
        newPoints = vtk.vtkPoints()
        numPoints = newDataSet.GetNumberOfPoints()
        num_arrays = newDataSet.GetPointData().GetNumberOfArrays()
        print("Number of arrays:", num_arrays)
        # Get the dimensions of the input dataset
        input_dimensions = newDataSet.GetDimensions()
        #print(input_dimensions[0])   # should be 1024
        #print(input_dimensions[1])   # should be 512
        #print(input_dimensions[2])   # should be 1
        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]
        print("RD Input Dimensions:", x_dimension, y_dimension, z_dimension)
        #y_vals = ""
        # populate newPoints with stereographic points
        indx = 0    # Index for each point in the newPoints vtk points data set
        for i in range(0, numPoints):
            coord = newDataSet.GetPoint(i)
            x0, y0, z0 = coord[:3]
            # Set initial values for x and y
            #x = x0
            #y = y0
            if (self.hemisphereToProject == "Northern Hemisphere"):
                if y0 >= 0:
                    #print("Populating newPoints with points in the Northern Hemisphere")
                    x,y = self.polar_lonlat_to_xy(x0, y0, TRUE_SCALE_LATITUDE, EARTH_RADIUS_KM, EARTH_ECCENTRICITY, self.hemisphereToProject)
                    newPoints.InsertPoint(indx,x,y,0)
                    indx = indx + 1
            else:
                if y0 <= 0:
                    x,y = self.polar_lonlat_to_xy(x0, y0, TRUE_SCALE_LATITUDE, EARTH_RADIUS_KM, EARTH_ECCENTRICITY, self.hemisphereToProject)
                    newPoints.InsertPoint(indx,x,y,0)
                    #y_vals = y_vals + " " + str(y0)
                    #print("Index of Southern Hemisphere point:", indx)
                    indx = indx + 1
        #print("Y-values in Southern hemisphere projection:", y_vals)
        # Get the dimensions of the newPoints vtk points data set that only
        # contains one hemisphere's points
        num_new_points = newPoints.GetNumberOfPoints()
        # Since the x-dimension is the same as the original input data set and the
        # z-dimension is flat with only one layer of points, then I can find the
        # y-dimension of newPoints by dividing the number of points in the newPoints
        # data set by the original x-dimension
        new_y_dimension = num_new_points // x_dimension
        print("New y dimension:", new_y_dimension)
        # Create the output data set
        outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)
        # Loop through each of the scalar arrays in the dataset
        # and add the array to the output dataset
        for j in range(0, num_arrays):
            ivals = newDataSet.GetPointData().GetArray(j)
            ca = vtk.vtkFloatArray()
            ca.SetName(ivals.GetName())
            ca.SetNumberOfComponents(1)
            ca.SetNumberOfTuples(num_new_points)
            # Print out the number of points
            print("Number of points: ", newPoints.GetNumberOfPoints())
            # add to output
            #outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)
            #print("Output data set: ", outputDataSet)
            #add the new array to the output
            outputDataSet.GetPointData().AddArray(ca)
            #copy the values over element by element
            indx = 0    # Index for each scalar value in the ca vtk float array
            for i in range(0, numPoints):
                coord = newDataSet.GetPoint(i)
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
        outputDataSet.SetPoints(newPoints)
        executive = self.GetExecutive()
        outInfo = executive.GetOutputInformation(0)

        self.output_x_dimension = x_dimension
        self.output_y_dimension = new_y_dimension
        self.output_z_dimension = z_dimension

        print("X-dimension:", x_dimension)

        outInfo.Set(executive.WHOLE_EXTENT(), 0, (self.output_x_dimension - 1), 0, (self.output_y_dimension - 1), 0, (self.output_z_dimension - 1))
        outputDataSet.SetDimensions(self.output_x_dimension,self.output_y_dimension,self.output_z_dimension)
        return 1
    def RequestInformation(self, request, inInfo, outInfo):
        print("I am running RequestInformation")
        #print("I am running RequestInformation")

        # get the first input.
        inputDataSet0 = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))

        # Get the list of array names from the input data set and assign that list to the
        # availableArrays variable
        num_arrays = inputDataSet0.GetPointData().GetNumberOfArrays()

        # Create a list of the array names inside of the input data set
        array_list = []
        for i in range(num_arrays):
            array_name = inputDataSet0.GetPointData().GetArray(i).GetName()
            array_list.append(array_name)
        self._availableArrays = array_list
        
        # Get the dimensions of the input dataset
        inputData = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))
        input_dimensions = inputData.GetDimensions()
        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]

        self.input_x_dimension = x_dimension
        self.input_y_dimension = y_dimension
        self.input_z_dimension = z_dimension

        # Set the output dimensions to equal the input dimensions so that the output
        # dimensions are not equal to 0, 0, 0 the first time the filter is applied
        # and nothing shows up in the output. These dimensions will later be changed
        # to the correct values in the RequestData section of the code.
        if (self.GetColumnAtEnd() == True):
            self.output_x_dimension = self.input_x_dimension + 1
        else:
            self.output_x_dimension = self.input_x_dimension
        if (self.output_y_dimension == 0):
            self.output_y_dimension = self.input_y_dimension
        if (self.output_z_dimension == 0):
            self.output_z_dimension = self.input_z_dimension
        print("Output dimensions:", self.output_x_dimension, self.output_y_dimension, self.output_z_dimension)

        executive = self.GetExecutive()
        outInfo = executive.GetOutputInformation(0)
        #outInfo.Set(executive.WHOLE_EXTENT(), 0, x_dimension, 0, (y_dimension - 1), 0, 0)
        outInfo.Set(executive.WHOLE_EXTENT(), 0, (self.output_x_dimension - 1), 0, (self.output_y_dimension - 1), 0, (self.output_z_dimension - 1))
        #outInfo.Set(executive.WHOLE_EXTENT(), 0, self.input_x_dimension, 0, (self.output_y_dimension - 1), 0, (self.output_z_dimension - 1))
        return 1
    def RequestUpdateExtent(self, request, inInfo, outInfo):
        print("I am running RequestUpdateExtent")

        
        # Get the dimensions of the input dataset
        inputData = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))
        input_dimensions = inputData.GetDimensions()
        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]


        executive = self.GetExecutive()
        inInfo = executive.GetInputInformation(0, 0)
        #inInfo.Set(executive.UPDATE_EXTENT(), 0, (x_dimension - 1), 0, (y_dimension - 1), 0, 0)
        inInfo.Set(executive.UPDATE_EXTENT(), 0, (self.input_x_dimension - 1), 0, (self.input_y_dimension - 1), 0, (self.input_z_dimension - 1))
        return 1
