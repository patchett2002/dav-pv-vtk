from vtkmodules.vtkCommonDataModel import vtkDataSet
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa

# new module for ParaView-specific decorators.
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

from paraview import vtk
import numpy as np #needed for interpolation and pi

#def createModifiedCallback(anobject):
    #import weakref
    #weakref_obj = weakref.ref(anobject)
    #anobject = None
    #def _markmodified(*args, **kwars):
        #o = weakref_obj()
        #if o is not None:
            #o.Modified()
    #return _markmodified


@smproxy.filter(label="VTS to Robinson VTS")
@smproperty.input(name="Input")
class VTStoRobinsonVTS(VTKPythonAlgorithmBase):

    # Arrays for interpolating Robinson Coordinates
    degrees=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]
    X=[1,0.9986,0.9954,0.99,0.9822,0.973,0.96,0.9427,0.9216,0.8962,\
   
       0.8679,0.835,0.7986,0.7597,0.7186,0.6732,0.6213,0.5722,0.5322]
    Y=[0,0.062,0.124,0.186,0.248,0.31,0.372,0.434,0.4958,0.5571,0.6176,\
   
       0.6769,0.7346,0.7903,0.8435,0.8936,0.9394,0.9761,1]

    def __init__(self):
        #VTKPythonAlgorithmBase.__init__(self)
        #VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType='vtkStructuredGrid')
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)

        #from vtkmodules.vtkCommonCore import vtkDataArraySelection
        #self._arrayselection = vtkDataArraySelection()
        #self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))
        #self._arrayselection.AddArray('One')
        #self._arrayselection.AddArray('Two')

        # Set the default realMeridian value to 0
        self.realMeridian = 0

    def FillInputPortInformation(self, port, info):
        info.Set(vtk.vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
        return 1

    def FillOutputPortInformation(self, port, info):
        info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkStructuredGrid")
        return 1

    #@smproperty.dataarrayselection(name="Arrays")
    #def GetDataArraySelection(self):
        #return self._arrayselection

    #@smproperty.intvector(name="PhiResolution", default_values=0)
    #@smdomain.intrange(min=0, max=1)
    @smproperty.xml("""
        <IntVectorProperty name="CentralMeridianAtZero"
            number_of_elements="1"
            default_values="0"
            command="SetCentralMeridian">
            <BooleanDomain name="bool" />
            <Documentation>If on, sets the central meridian of the Robinson
            projection to zero degrees longitude.</Documentation>
        </IntVectorProperty>""")
    def SetCentralMeridian(self, x):
        self.realMeridian = x
        self.Modified()

    def GetCentralMeridian(self):
        return self.realMeridian


    def GetRobinsonPoint(self,x,y,mylambda0=0):
        '''x,y is real longitude,latitude - Return Robinson point'''
        R = np.pi    # This number simply scales the output
        mylambda = x * 0.01745329252 # longitude of point to plot in Radians
        #mylambda0 = np.pi            # Central Meridian in Radians 180 degrees == pi
        #mylambda0 = 0            # Central Meridian in Radians 0 degrees == 0
        RobinsonX = np.interp([np.abs(y)],self.degrees,self.X) # interpolate the Robinson X
        RobinsonY = np.interp([np.abs(y)],self.degrees,self.Y) # interpolate the Robinson Y
        # calculate the robinson coordinates
        new_x = 0.8487 * R * RobinsonX * (mylambda - mylambda0)
        new_y = 1.3523 * R * RobinsonY
        # Check if this was in southern Hemisphere, interpolation didn't deal with negatives
        if y < 0:
            new_y = new_y * -1
        return (new_x[0], new_y[0])

    # This function will take in an array and the number of elements inside of it
    # and find and return the index of the element in the array that is equal to
    # or approximately equal to 180
    def FindIndxAt180(self, arr, arr_dimension):
        for i in range(arr_dimension):
            # Check that the element in the array is less than 0.35 away from 180
            if (abs((180 - arr[i])) < 0.35):
                return i
        
        #print(arr[512])
        return -1

    def CenterAtZero(self, inputDataSet):
        # Create the new input data set
        #newInputDataSet = vtk.vtkStructuredGrid.New()
        newInputDataSet = vtk.vtkStructuredGrid()

        newPoints = vtk.vtkPoints()
        numPoints = inputDataSet.GetNumberOfPoints()

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

        lon_indx_180 = self.FindIndxAt180(original_lon_points, x_dimension)
        # Should print 512 for the lon_indx_180
        print("Index of the x-dimension with a longitude of 180: ", lon_indx_180)

        latPoints = []
        lonPoints = []

        for i in range(0, numPoints):
            coord = inputDataSet.GetPoint(i)
            x0, y0 = coord[:2]
            if x0 > 180:
                x0 = x0 - 360
            latPoints.append(y0)
            lonPoints.append(x0)


        lonPoints = np.array(lonPoints)
        latPoints = np.array(latPoints)
        latPoints=latPoints.reshape(y_dimension, x_dimension)
        lonPoints=lonPoints.reshape(y_dimension, x_dimension)
        
        #A = lonPoints[0:512,0:513]    # [Rows, Columns]
        A = lonPoints[0:y_dimension,0:int(lon_indx_180 + 1)]    # [Rows, Columns]
        #B = lonPoints[0:512,512:]
        #B = lonPoints[0:512,513:]    # Right hand side of the dataset
        B = lonPoints[0:y_dimension,int(lon_indx_180 + 1):]    # Right hand side of the dataset
        NewLonArray = np.hstack((B,A))
        NewLonArray = NewLonArray.reshape(y_dimension * x_dimension,)


        #ALat = latPoints[0:512,0:512]
        ALat = latPoints[0:y_dimension,0:int(lon_indx_180)]
        #BLat = latPoints[0:512,512:]
        BLat = latPoints[0:y_dimension,int(lon_indx_180):]
        NewLatArray = np.hstack((BLat, ALat))
        NewLatArray = NewLatArray.reshape(y_dimension * x_dimension,)


        # Take the x, y and z coordinates from the input data's vtk
        # point array and convert them into a Robinson projection
        # before putting the new coordinate values into a list
        RobinsonPoints = []
        for x0,y0 in zip(NewLonArray, NewLatArray):
            RobinsonPoints.append((x0,y0,0))


        # Put the Robinson projection coordinates into the empty vtk
        # point array newPoints
        count = 0
        for i in RobinsonPoints:
            newPoints.InsertNextPoint(i)
            count += 1


        # Loop through each of the scalar arrays in the dataset
        # and add the array to the output dataset
        for j in range(0, num_arrays):
            # Convert the scalar values into a list
            rawdat = []
            ivals = inputDataSet.GetPointData().GetArray(j)
            for i in range(0, numPoints):
                rawdat.append(ivals.GetValue(i))


            # Create a numpy array from the rawdat list and move the right
            # half of the array (the scalar values found at longitudes
            # greater than 180 degrees) over to the left
            nprawdat = np.array(rawdat)
            nprawdat=nprawdat.reshape(y_dimension, x_dimension)


            #A = nprawdat[0:512,0:512]
            A = nprawdat[0:y_dimension,0:int(lon_indx_180)]
            #B = nprawdat[0:512,512:]
            B = nprawdat[0:y_dimension,int(lon_indx_180):]
            NewArray = np.hstack((B,A))
            NewArray = NewArray.reshape(y_dimension * x_dimension,)


            # Input the values from the altered numpy array into a new vtk
            # float array
            dat = vtk.vtkFloatArray()
            dat.SetName(ivals.GetName())
            for i in NewArray:
                dat.InsertNextValue(i)


            ca = vtk.vtkFloatArray()
            ca.SetName(ivals.GetName())
            ca.SetNumberOfComponents(1)
            ca.SetNumberOfTuples(numPoints)


            #add the new array to the newInputDataSet
            newInputDataSet.GetPointData().AddArray(ca)


            #copy the values over element by element
            for i in range(0, numPoints):
                ca.SetValue(i, dat.GetValue(i))



        newInputDataSet.SetDimensions(x_dimension,y_dimension,z_dimension)
        newInputDataSet.SetPoints(newPoints)
        
        return newInputDataSet

    def RequestData(self, request, inInfo, outInfo):
        # get the first input.
        inputDataSet0 = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0]))

        # compute a value.
        #data = inputDataSet0.PointData["Sea Level Change (m)"] / 2.0
        #pointX, pointY, pointZ = inputDataSet0.GetPoint(0)

        newPoints = vtk.vtkPoints()
        numPoints = inputDataSet0.GetNumberOfPoints()
        
        num_arrays = inputDataSet0.GetPointData().GetNumberOfArrays()
        print("Number of arrays:", num_arrays)
        
        # Get the dimensions of the input dataset
        input_dimensions = inputDataSet0.GetDimensions()


        print("Dimensions:")
        print(input_dimensions[0])   # should be 1025
        print(input_dimensions[1])   # should be 512
        print(input_dimensions[2])   # should be 1

        x_dimension = input_dimensions[0]
        y_dimension = input_dimensions[1]
        z_dimension = input_dimensions[2]

        print("Real Meridian: ", self.GetCentralMeridian())

        # Check that the Central Meridian At Zero checkbox is checked, and if it
        # is, then alter the input data set so that it is centered at zero
        # degrees longitude
        if (self.GetCentralMeridian() == True):
            inputDataSet0 = self.CenterAtZero(inputDataSet0)
            #inputDataSet1 = self.CenterAtZero(inputDataSet0)

        # Get the longitude value at the first and last columns and calculate
        # central meridian to use for the Robinson projection
        first_col = inputDataSet0.GetPoint(0)
        last_col = inputDataSet0.GetPoint(x_dimension - 1)

        first_lon = first_col[0]
        last_lon = last_col[0]

        mid_lon = first_lon + ((last_lon - first_lon)/2)

        # Convert the central meridian from degrees to radians
        mid_lon_rad = mid_lon * np.pi/180


        #print("Real Meridian: ", self.realMeridian)
        #print("Real Meridian: ", self.GetCentralMeridian())


        # populate newPoints with Robinson points
        for i in range(0, numPoints):
            coord = inputDataSet0.GetPoint(i)
            x0, y0, z0 = coord[:3]
            #x,y = self.GetRobinsonPoint(x0,y0)
            x,y = self.GetRobinsonPoint(x0,y0,mid_lon_rad)
            newPoints.InsertPoint(i,x,y,0)

        # Create the output data set
        outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)


        # Loop through each of the scalar arrays in the dataset
        # and add the array to the output dataset
        for j in range(0, num_arrays):
            ivals = inputDataSet0.GetPointData().GetArray(j)
            ca = vtk.vtkFloatArray()
            ca.SetName(ivals.GetName())
            ca.SetNumberOfComponents(1)
            ca.SetNumberOfTuples(numPoints)

            # Print out the number of points
            print("Number of points: ", newPoints.GetNumberOfPoints())

            # add to output
            #outputDataSet = vtk.vtkStructuredGrid.GetData(outInfo)
        

            print("Output data set: ", outputDataSet)



            #add the new array to the output
            outputDataSet.GetPointData().AddArray(ca)


            #copy the values over element by element
            for i in range(0, numPoints):
                ca.SetValue(i, ivals.GetValue(i))

            # Try printing out the value at the 1025th point
            print("Value at point 1024: ", ivals.GetValue(1024))
            print("Value at point 1024 in ca: ", ca.GetValue(1024))

        outputDataSet.SetDimensions(x_dimension,y_dimension,z_dimension)
        outputDataSet.SetPoints(newPoints)

        print(outputDataSet)


        return 1
