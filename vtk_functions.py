import vtk
print vtk.VTK_VERSION

class image_wrapper(object):
    """
    """

    def __init__(self, filename, array_name):
        """
        """
        self.filename = filename
        self.array_name = array_name
        self.read_image()

    def read_image(self):
        """
        """

        reader = vtk.vtkNIFTIImageReader()
        reader.SetFileName(self.filename)
        reader.Update()

        #print reader.GetNIFTIHeader()
        #print reader.GetQFormMatrix()
        #print reader.GetOutput().GetDimensions()
        #print reader.GetOutput().GetCenter()
        #print reader.GetOutput().GetSpacing()

        ox = reader.GetQFormMatrix().GetElement(0, 3)
        oy = reader.GetQFormMatrix().GetElement(1, 3)
        oz = reader.GetQFormMatrix().GetElement(2, 3)

        self.image = reader.GetOutput()
        self.image.SetOrigin((ox, oy, oz))
        # print self.image.GetPointData().GetArray(0).GetNumberOfComponents()
        # print self.image.GetPointData().GetArray(0).GetDataTypeAsString()
        # print self.image.GetPointData().GetArray(0).GetDataType()
        # print self.image.GetPointData().GetArray(0).GetTuple(0)

    def probe_points_old(self, point_list):
        """
        """

        vtk_points_to_probe = points_to_vtk_points(point_list)

        self._probe = vtk.vtkProbeFilter()
        self._probe.SetSourceData(self.image)
        self._probe.SetInputData(vtk_points_to_probe)
        self._probe.Update()

        self.probed_points = self._probe.GetOutput()
        self.probed_points.GetPointData().GetArray(0).SetName(self.array_name)

    def probe_points(self, point_list):
        """
        """
#       print self.image.GetPointData().GetArray(0).GetNumberOfComponents()
#       print self.image.GetPointData().GetArray(0).GetDataTypeAsString()
#       print self.image.GetPointData().GetArray(0).GetDataType()
#       print self.image.GetPointData().GetArray(0).GetTuple(0)
        arrays = {"unsigned char": "vtkUnsignedIntArray",
                  "unsigned short": "vtkUnsignedShortArray",
                  "unsigned int": "vtkUnsignedIntArray",
                  "float": "vtkFloatArray"}

        a_type = self.image.GetPointData().GetArray(0).GetDataTypeAsString()

        self._probed_array = getattr(vtk, arrays[a_type])()
        self._probed_array.SetName(self.array_name)
        self._probed_array.SetNumberOfComponents(
            self.image.GetPointData().GetArray(0).GetNumberOfComponents())
        self._probed_array.SetNumberOfTuples(len(point_list))

        for (idx, point) in enumerate(point_list):

            x, y, z = point
            img_point_id = self.image.FindPoint(x, y, z)
            value = self.image.GetPointData().GetArray(0).GetTuple(img_point_id)
            self._probed_array.SetTuple(idx, value)

    def get_array(self):
        """
        """

        cached_array = self._probed_array
        temp_array = cached_array.NewInstance()
        temp_array.DeepCopy(cached_array)
        temp_array.SetName(cached_array.GetName())

        return temp_array

    def get_array_old(self):
        """
        """
        cached_array = self.probed_points_old.GetPointData().GetArray(self.array_name)
        temp_array = cached_array.NewInstance()
        temp_array.DeepCopy(cached_array)
        temp_array.SetName(cached_array.GetName())

        return temp_array

def points_to_vtk_points(points_list):
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    for i in range(len(points_list)):
        id_ = points.InsertNextPoint(points_list[i])
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(id_)
    point = vtk.vtkPolyData()
    point.SetPoints(points)
    point.SetVerts(vertices)
    return point
