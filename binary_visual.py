# -*- coding: utf-8 -*-
#*************************************************************************
#***File Name: binary_visual.py
#***Author: Zhonghai Zhao
#***Mail: zhaozhonghi@126.com 
#***Created Time: 2018年03月25日 星期日 14时39分05秒
#*************************************************************************
class binary_visual(object):
    '''
    This class contains some function to visualize binary file.
    '''
    # initialization
    def __init__(self):
        pass
    def plot_array(self, filenumber=0, filename='density', case='radiography', factor=1.0, iflist=False, axis='z', index=100, vlim=1.0, arti_v=False, shape=[300, 200, 200], display=True, dpix=10, dpiy=5, ifaverage=False, nslices=1):
        '''
        This function is used to visualize 2d field data, from binary file.
        Parameters:
            filenumber    - binary file number list.
            filename      - binary file name.
            case          - case, to chose dimensionless constants.
            factor        - vmax multifactor.
            iflist        - if a list of binary files.
            axis          - slice axis.
            index         - slice index.
            vlim          - limit range.
            arti_v        - limit range artificially.
            shape         - array shape.
            display       - if display figue, or save it.
            dpix          - figure size x.
            dpiy          - figure size y.
            ifaverage     - if average on axis.
            nslice        - n slices over axis.         
        Returns:
            None.
        Raises:
            KeyError.
        '''
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from binary_io import binary_io as bio
        from sdf_class import sdf_class as mysdf
        bio = bio()
        sc = mysdf()
        dimen = len(shape)
        # get namelist
        files = []
        if (iflist == True):
            n = len(filenumber)
            for i in range(n):
                string = filename + '_' + str(filenumber[i]).zfill(4) + '.dat'
                files.append(string)
        else:
            string = filename + '.dat'
            files.append(string)
        # read array
        narray = []
        for i in range(len(files)):
            if (ifaverage == False):
                array = bio.binary_read(filename=files[i], shape=shape, axis=axis, index=index)
            else:
                array = 0
                for j in range(index-nslices, index+nslices+1):
                    array0 = bio.binary_read(filename=files[i], shape=shape, axis=axis, index=j)
                    array = array + array0
                array = array/(2*nslices+1.0)
            #if (dimen == 3):
            #    if (axis == 'x'):
            #        array = data[index, :, :]
            #    elif (axis == 'y'):
            #        array = data[:, index, :]
            #    else:
            #        array = data[:, :, index]
            #else:
            #    array = data
            narray = narray + [np.transpose(array)]
        # differrent case
        if (case == 'radiography'):
            #from constants import proton_radiography as const
            from constants import proton_benchmark as const
            dict_const = {"density":const.n0, "energy":const.ek, "length":const.di, \
                    "bx":const.B0, "by":const.B0, "bz":const.B0}
            constants = dict_const[filename]
            dx = dict_const["length"]
            labels = sc.select_label(case='radiography', axis=axis)
            title = []
            path = []
            if (iflist == True):
                for i in range(n):
                    title = title + [filename + ' ' + 'time = ' + str(filenumber[i]) + ' ' + '$ T_0 $']
                    path = path + [filename + '_' + str(filenumber[i]) + '_' + axis + '_' + str(index) + '.png']
            else:
                title = title + [filename]
                path = path + [filename + '_' + axis + '_' + str(index) + '.png']
        else:
            pass
        # find max
        if (arti_v == False):
            total = np.abs(np.array(narray))
            vmax = total.max() * factor
        else:
            vmax = vlim
        sample = np.array(narray[0])
        if (np.min(sample) < 0):
            cmap = cm.RdBu_r
            vmin = -vmax
        else:
            cmap = cm.Blues
            vmin = 0
        # get extent
        extent = bio.get_extent(axis=axis, shape=shape, constant=dx)
        #plot
        plt.figure(figsize=(dpix, dpiy))
        for i in range(len(files)):
            ax = plt.gca()
            im = ax.imshow(narray[i]/constants, extent=extent, origin='lower', cmap=cmap, vmax=vmax/constants, vmin=vmin/constants, interpolation='spline36')
            plt.xlabel(labels[0])
            plt.ylabel(labels[1])
            plt.title(title[i])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='3%', pad=0.1)
            plt.colorbar(im, cax=cax)
            # display or save picture
            if (display == True):
                plt.show()
            else:
                s1 = 'figure/'
                plt.savefig(s1+path[i], dpi=300)
                if (i == len(files)-1):
                    plt.close()
                else:
                    plt.clf()
#        return extent
            

