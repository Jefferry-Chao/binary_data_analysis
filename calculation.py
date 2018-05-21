# -*- coding: utf-8 -*-
#*************************************************************************
#***File Name: calculation.py
#***Author: Zhonghai Zhao
#***Mail: zhaozhonghi@126.com 
#***Created Time: 2018年03月25日 星期日 14时39分05秒
#*************************************************************************
class calculation(object):
    '''
    This class contains some functions used to calculate some data.
    '''
    # initialization
    def __init__(self):
        pass
    # 多项式拟合
    def polyfit(self, X, Y, n):
        '''
        This function is poly fit.
        Parameters:
            X         - independent variable.
            Y         - dependent variable.
            n         - index.
        Returns:
            polynomial coefficient,
        Raises:
            KeyError.

        '''
        import numpy as np
        len_x = len(X)
        len_y = len(Y)
        if ((len_x != len_y) or (len_x <2) or (n <= 0)):
            raise ValueError
        else:
            Z = np.polyfit(X, Y, n)
            print np.poly1d(Z)
    def contrast(self, filename='energy.dat', shape=[50, 2200, 2200], axis='x', index=30, iffocusing=True, smooth=1, frac=0.025, b_factor=0.2, iflimit=False, limx=[0, 1.0], limy=[0, 1.0], ifreturn=False, ifdisplay=False,ifadd=False, epsilon=0.5):
        '''
        This function is used to calculate contrast in proton radiography, use as benchmark.
        Parameters:
            filename  - binary file name.
            shape     - array shape.
            axis      - slice axis.
            index     - slice index.
            iffocusing- if focusing case or not.
            smooth    - smooth slice.
            frac      - smooth factor.
            b_factor  - b field factor to critic value.
            iflimit   - if limit axis range artificially.
            limx      - x axis range limiter.
            limy      - y axis range limiter.
            ifreturn  - if return array.
            ifdisplay - if display figure.
            ifadd     - if add all branches.
            epsilon   - small parameter.
        Returns:
            None.
        Raises:
            KeyError.
        '''
        import numpy as np
        from scipy import interpolate
        from statsmodels.nonparametric.smoothers_lowess import lowess
        from matplotlib import pyplot as plt
        from binary_io import binary_io as mysdf
        bio = mysdf()
        # simulation
        n = int(shape[2]/2)
        a = 100.
        b = 300.
        A = 1836.
        W = 14.7
        zs = 1.
        g_factor = 2.0
        magnification = 11.0
        if (smooth != 1):
            narray = []
            for i in range(smooth):
                narray = narray + [bio.polar2line(filename=filename, shape=shape, axis=axis, index=index-smooth/2+i)]
            array = np.sum(narray, 0)/float(smooth)
        else:
            array = bio.polar2line(filename=filename, shape=shape, axis=axis, index=index)
        x1 = np.array(range(n))/magnification/a*g_factor
        filter_array = lowess(array, x1, is_sorted=True, frac=frac, it=0)
        wake_average = np.sum(array[-200:-100])/100.0
        y1 = filter_array[:, 1]/wake_average
        if (iflimit == False):
            limx[1] = np.max(x1)
            limy[1] = 1.2*np.max(y1)
        # theory calculation
        x0 = np.linspace(0, np.max(x1), n)
        if (iffocusing == True):
            b0 = 0.1903*(a/b)*np.sqrt(A*W)/zs
            b0 = b0*b_factor
            mu = 0.533871e-3 * b0 * b/np.sqrt(A*W)
            nv = mu*zs/a*(1e4)
            x2 = x0*np.abs(1 - nv*np.exp(-x0*x0))
            if (b_factor < 1.0):
                y2 = 1.0/np.abs(np.exp(-2*x0*x0)*(nv - np.exp(x0*x0))*(np.exp(x0*x0) - nv*(1 - 2*x0*x0)))
            else:
                y2 = (1.0 + epsilon)/(epsilon + np.abs(np.exp(-2*x0*x0)*(nv - np.exp(x0*x0))*(np.exp(x0*x0) - nv*(1 - 2*x0*x0))))
        else:
            b0 = 0.4262*(a/b)*np.sqrt(A*W)/zs
            b0 = b0*b_factor
            mu = 0.533871e-3 * b0 * b/np.sqrt(A*W)
            nv = mu*zs/a*(1e4)
            x2 = x0*np.abs(1 + nv*np.exp(-x0*x0))
            if (b_factor < 1.0):
                y2 = 1.0/np.abs(np.exp(-2*x0*x0)*(nv + np.exp(x0*x0))*(np.exp(x0*x0) + nv*(1 - 2*x0*x0)))
            else:
                y2 = (1.0 + epsilon)/(epsilon + np.abs(np.exp(-2*x0*x0)*(nv + np.exp(x0*x0))*(np.exp(x0*x0) + nv*(1 - 2*x0*x0))))
        if ((ifadd == True) and (iffocusing == True)):
            n1 = 247
            n2 = 464
            subx1 = x2[n1]
            subx2 = x2[n2]
            # interpolate f1
            f1x1 = x2[0:n1]
            f1y1 = y2[0:n1]
            f1x2 = np.linspace(subx1, 4.0, 1000)
            f1y2 = np.zeros(1000, np.float)
            f1x3 = np.hstack([f1x1, f1x2])
            f1y3 = np.hstack([f1y1, f1y2])
            f1 = interpolate.interp1d(f1x3, f1y3, 'slinear')
            arrayf1 = f1(x0)
            # interpolate f2
            f2x1 = x2[n1:n2]
            f2y1 = y2[n1:n2]
            f2x2 = np.linspace(subx1, 4.0, 1000)
            f2y2 = np.zeros(1000, np.float)
            # reverse
            lx1 = list(f2x1)
            lx1.reverse()
            f2x1 = np.array(lx1)
            f2x1[0] = 0.0
            ly1 = list(f2y1)
            ly1.reverse()
            f2y1 = np.array(ly1)
            f2x3 = np.hstack([f2x1, f2x2])
            f2y3 = np.hstack([f2y1, f2y2])
            f2 = interpolate.interp1d(f2x3, f2y3, 'slinear')
            arrayf2 = f2(x0)
            # interpolate f3
            f3x1 = x2[n2:n-1]
            f3x1[0] =  0.0
            f3x1[-1] =  4.0
            f3y1 = y2[n2:n-1]
            f3 = interpolate.interp1d(f3x1, f3y1, 'slinear')
            arrayf3 = f3(x0)
            y3 = arrayf1 + arrayf2 + arrayf3
        if ((ifadd == True) and (iffocusing == False)):
            n1 = 470
            n2 = 985
            subx1 = x2[n1]
            subx2 = x2[n2]
            # interpolate f1
            f1x1 = x2[0:n1]
            f1y1 = y2[0:n1]
            f1x2 = np.linspace(subx1, 4.0, 1000)
            f1y2 = np.zeros(1000, np.float)
            f1x3 = np.hstack([f1x1, f1x2])
            f1y3 = np.hstack([f1y1, f1y2])
            f1 = interpolate.interp1d(f1x3, f1y3, 'slinear')
            arrayf1 = f1(x0)
            # interpolate f2
            f2x1 = x2[n1:n2]
            f2y1 = y2[n1:n2]
            f2x2 = np.linspace(subx1, 4.0, 500)
            f2y2 = np.zeros(500, np.float)
            f2x3 = np.linspace(0, subx2, 500)
            f2y3 = np.zeros(500, np.float)
            # reverse
            lx1 = list(f2x1)
            lx1.reverse()
            f2x1 = np.array(lx1)
            f2x1[0] = subx2
            f2x1[-1] = subx1
            ly1 = list(f2y1)
            ly1.reverse()
            f2y1 = np.array(ly1)
            f2x4 = np.hstack([f2x3, f2x1, f2x2])
            f2y4 = np.hstack([f2y3, f2y1, f2y2])
            f2 = interpolate.interp1d(f2x4, f2y4, 'slinear')
            arrayf2 = f2(x0)
            # interpolate f3
            f3x1 = x2[n2:n-1]
            f3x1[0] =  subx2
            f3x1[-1] =  4.0
            f3y1 = y2[n2:n-1]
            f3x2 = np.linspace(0, subx2, 1000)
            f3y2 = np.zeros(1000, np.float)
            f3x3 = np.hstack([f3x2, f3x1])
            f3y3 = np.hstack([f3y2, f3y1])
            f3 = interpolate.interp1d(f3x3, f3y3, 'slinear')
            arrayf3 = f3(x0)
            y3 = arrayf1 + arrayf2 + arrayf3
        if (ifreturn == False):
            plt.plot(x1, y1, label='simulation')
            plt.plot(x2, y2, label='exact solution')
            plt.xlabel(r'$ r/R $')
            plt.ylabel(r'$ I/I_0 $')
            plt.legend()
            plt.xlim(xmin=limx[0], xmax=limx[1])
            plt.ylim(ymin=limy[0], ymax=limy[1])
            if (display == True):
                plt.show()
            else:
                path = 'figure/proton_contrast.png'
                plt.savefig(path, dpi=300)
        else:
            if(ifadd == False):
                return (x1, y1, x2, y2)
            else:
                return (x1, y1, x2, y2, x0, y3)
            #return (x0, arrayf1, arrayf2, arrayf3)
        #print b0, mu, nv
    def reconnection_flux(self, namelist):
        '''
        This function is used plot reconnection flux.
            Parameters:
                namelist  - sdf file name list.
            Returns:
                None.
            Raises:
                KeyError.
        '''
        import numpy as np
        import matplotlib.pyplot as plt
        from sdf_class import sdf_class as mysdf
        sc = mysdf()
        # data
        array = sc.line_integrate(namelist)
        x = np.array(namelist) * 0.1
        # line
        x1 = np.linspace(0.8, 0.8, 30)
        y1 = np.linspace(0, 3, 30)
        x2 = np.linspace(2.0, 2.0, 30)
        y2 = np.linspace(0, 3, 30)
        x3 = np.linspace(0, 3, 30)
        y3 = np.linspace(2.6, 2.6, 30)
        # plot
        plt.figure(figsize=(8, 5))
        plt.plot(x, array)
        plt.plot(x1, y1, color='k', linestyle='--')
        plt.plot(x2, y2, color='k', linestyle='--')
        plt.plot(x3, y3, color='k', linestyle='--')
        plt.xlabel(r'$ t/\Omega_i $', fontsize=16)
        plt.ylabel('REconnection Flux', fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlim([0, 3])
        plt.ylim([0, 3])
        plt.show()
    def plot_array(self, x, y, xlabel, ylabel, s=20,  axislimit=[0, 1, 0, 1], font=30):
        '''
        This function is used to plot line figure.
            Parameters:
                x         - x array.
                y         - y array.
                xlabel    - x axis label.
                ylabel    - y axis label.
                s         - marker size.
                axislimit - x and y axis limit.
                font      - font size.
            Returns:
                None.
            Raises:
                KeyError.
        '''
        import numpy as np
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 6))
        plt.scatter(x, y, color='r', marker='s', s=s)
        plt.xlabel(xlabel, fontsize=font)
        plt.ylabel(ylabel, fontsize=font)
        plt.xticks(fontsize=font)
        plt.yticks(fontsize=font)
        plt.xlim(axislimit[:2])
        plt.ylim(axislimit[-2:])
        plt.text(10.5, 27, "b)", size=30)
        #plt.savefig('divergence_angle.png', dpi=120)
        plt.show()
