# -*- coding: utf-8 -*-
#*************************************************************************
#***File Name: binary_io.py
#***Author: Zhonghai Zhao
#***Mail: zhaozhonghi@126.com 
#***Created Time: 2018年03月25日 星期日 14时39分05秒
#*************************************************************************
class binary_io(object):
    '''
    This class contains some functions to read data from sdf file and write data to anothet file.
    '''
    # initialization
    def __init__(self):
        pass
    # write particle information
    def particle_write(self, filenumber=0, subset='tracer_p', particle='tracer_pro', prefix='3'):
        '''
        This function is used to write particle data to file.
        Parameters:
            filenumber    - sdf file number, an integer or an integer list.
            subset        - particle subset to read an write.
            particle      - particle name.
            prefix        - sdf file prefix.
        Returns:
            None.
        Raises:
            KeyError.
        '''
        #import numpy as np
        import struct
        from sdf_class import sdf_class as mysdf
        sc = mysdf()
        namelist = sc.get_list(filenumber)
        field = ['Px', 'Py', 'Pz']
        position = ['X', 'Y', 'Z']
        for i in range(len(namelist)):
            data_dict = sc.get_data(namelist[i], prefix=prefix)
            # read momentum
            p = []
            for j in range(len(field)):
                data = data_dict['Particles_' + field[j] + '_subset_' + subset + '_' + particle].data
                p = p + [data]
            # read grid
            grid = data_dict['Grid_Particles_subset_' + subset + '_' + particle].data
            weight = data_dict['Particles_Weight_subset_' + subset + '_' + particle].data
            # write data to file
            n = len(p[0])
            files = open(str(namelist[i]) + '.dat', 'wb')
            for j in range(n):
                for k in range(len(field)):
                    files.write(struct.pack('d', grid[k][j]))
                    files.write(struct.pack('d', p[k][j]))
                files.write(struct.pack('d', weight[j]))
            files.close()
    # write field information
    def field_write(self, filenumber=0, field='bx', prefix='1'):
        '''
        This function is used to write field data to file.
        Parameters:
            filenumber    - sdf file number, an integer or an integer list.
            field         - field information.
            prefix        - sdf file prefix.
        Returns:
            None.
        Raises:
            KeyError.
        '''
        import struct
        import numpy as np
        from sdf_class import sdf_class as mysdf
        sc = mysdf()
        namelist = sc.get_list(filenumber)
        data_dict = sc.get_data(namelist[0], prefix=prefix)
        keys = data_dict.keys()
        fields = sc.get_field(field=field, keys=keys)
        narray = []
        for i in range(len(namelist)):
            data_dict = sc.get_data(namelist[i], prefix=prefix)
            data = data_dict[fields].data
            narray.append(np.transpose(data))
        shape = narray[0].shape
        dims = len(shape)
        for each in range(len(namelist)):
            files = open(field + '_' + str(namelist[each]).zfill(4) + '.dat', 'wb')
            if(dims == 2):
                for j in range(shape[1]):
                    for i in range(shape[0]):
                        files.write(struct.pack('d', narray[each][i, j]))
            elif(dims == 3):
                for k in range (shape[2]):
                    for j in range(shape[1]):
                        for i in range(shape[0]):
                            files.write(struct.pack('d', narray[each][i, j, k]))
            else:
                pass
            files.close()
    # read binary file
    def binary_read(self, filename=r'density.dat', shape=[300, 200, 200], axis='x', index=100):
        '''
        This function is used to read binary file.
        Parameters:
            filename      - binary file name.
            shape         - array shape.
            axis          - if 3d, axis is slice.
            index         - if 3d, slice at index
        Returns:
            array
        Raises:
            KeyError.
        '''
        import numpy as np
        import struct
        length = len(shape)
        files = open(filename, 'rb')
        if (length == 1):
            files.seek(0, 0)
            array = np.zeros(shape, np.float)
            for i in range(shape[0]):
                element = struct.unpack('d', files.read(8))
                array[i] = element[0]
        elif(length == 2):
            files.seek(0, 0)
            array = np.zeros(shape, np.float)
            for j in range(shape[1]):
                for i in range(shape[0]):
                    element = struct.unpack('d', files.read(8))
                    array[i, j] = element[0]
        elif(length == 3):
            if (axis == 'x'):
                array = np.zeros([shape[1], shape[2]], np.float)
                files.seek((index-1)*8, 0)
                for k in range(shape[2]):
                    for j in range(shape[1]):
                        element = struct.unpack('d', files.read(8))
                        array[j, k] = element[0]
                        files.seek((shape[0]-1)*8, 1)
            elif (axis == 'y'):
                array = np.zeros([shape[0], shape[2]], np.float)
                files.seek(shape[0]*(index-1)*8, 0)
                for k in range(shape[2]):
                    for i in range(shape[0]):
                        element = struct.unpack('d', files.read(8))
                        array[i, k] = element[0]
                    files.seek(shape[0]*(shape[1]-1)*8, 1)
            else:
                array = np.zeros([shape[0], shape[1]], np.float)
                files.seek(shape[0]*shape[1]*(index-1)*8, 0)
                for j in range(shape[1]):
                    for i in range(shape[0]):
                        element = struct.unpack('d', files.read(8))
                        array[i, j] = element[0]
            #for k in range(shape[2]):
            #    for j in range(shape[1]):
            #        for i in range(shape[0]):
            #            element = struct.unpack('d', files.read(8))
            #            array[i, j, k] = element[0]
        else:
            pass
        return array
    def get_extent(self, axis='z', shape=[300, 200, 200], constant=1.0):
        '''
        This function is used to get extent.
        Parameters:
            axis          - slice axis.
            shape         - array shape.
            constant      - dimensionless constants.
        Returns:
            None.
        Raise:
            KeyError.
        '''
        dimen = len(shape)
        x = binary_io.binary_read(self, filename='x.dat', shape=[shape[0]])
        y = binary_io.binary_read(self, filename='y.dat', shape=[shape[1]])
        if (dimen == 3):
            z = binary_io.binary_read(self, filename='z.dat', shape=[shape[2]])
        extent = []
        if(axis == 'x'):
            extent.append(y[0]/constant)
            extent.append(y[shape[1]-1]/constant)
            extent.append(z[0]/constant)
            extent.append(z[shape[2]-1]/constant)
        elif (axis == 'y'):
            extent.append(x[0]/constant)
            extent.append(x[shape[0]-1]/constant)
            extent.append(z[0]/constant)
            extent.append(z[shape[2]-1]/constant)
        else:
            extent.append(x[0]/constant)
            extent.append(x[shape[0]-1]/constant)
            extent.append(y[0]/constant)
            extent.append(y[shape[1]-1]/constant)
        return extent
    def polar2line(self, filename='energy', axis='x', index=30, shape=[50, 2200, 2200]):
        '''
        This function is used to average a polar plane to 1d line.
        Parameters:
            filename      - binary file name.
            axis          - slice axis.
            index         - slice index.
            shape         - array shape.
        Returns:
            None.
        Raise:
            KeyError.
        '''
        import numpy as np
        data = binary_io.binary_read(self, filename=filename, shape=shape, axis=axis, index=index)
        n = shape[2]/2
        center = [n, n]
        array = np.zeros(n, np.float)
        number = np.zeros(n, np.int)
        for i in range(n):
            for j in range(n):
                if (data[i, j] != 0.0):
                    length = int(round(np.sqrt((i-n)**2 + (j-n)**2)))
                    if (length < n):
                        array[length] = array[length] + data[i,j]
                        number[length] = number[length] + 1
        for i in range(n):
            if(number[i] != 0):
                array[i] = array[i]/np.float(number[i])
        return array
