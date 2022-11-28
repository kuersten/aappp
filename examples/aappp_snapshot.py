
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from os.path import isfile
import sys, getopt
import aappp




plotmargin=0.3
plotarea=5.28
point_size_factor=1.0

def show_snapshot(filename, colorwheel, label, save_filename):
    mysim=aappp.aappp_load(filename)
    x, y, angles, lx, ly= aappp.aappp_get_xythetalxly(mysim)
    N=len(x)
    x=np.array(x)
    y=np.array(y)
    angles=np.array(angles)
    filename='snapshot.png'
    plotsizex=(plotarea*lx/ly)**0.5
    plotsizey=plotsizex*ly/lx
    lxy=(lx*ly)**0.5
    fig, plot1= plt.subplots(figsize=(plotsizex+2*plotmargin, plotsizey+2*plotmargin))
    plt.subplots_adjust(left=plotmargin/(plotsizex+2*plotmargin), right=(plotmargin+plotsizex)/(plotsizex+2*plotmargin), bottom=plotmargin/(plotsizey+2*plotmargin), top=(plotmargin+plotsizey)/(plotsizey+2*plotmargin))
    hsv=plt.get_cmap('hsv')
    cNorm=mpl.colors.Normalize(vmin=0., vmax=1.)
    plt.scatter(x, y, s=point_size_factor*1000./N*plotarea**2/5.28**2, c=(angles+np.pi)/2./np.pi, cmap=hsv, norm=cNorm)
    plt.xlim([0., lx])
    plt.ylim([0., ly])
    #draw color wheel if desired
    if colorwheel:
        colorwheel_axes=fig.add_axes([plotmargin/(2*plotmargin+plotsizex)+1.0*plotsizex/(2*plotmargin+plotsizex)-0.17*(plotsizex*plotsizey)**0.5/(2*plotmargin+plotsizex), plotmargin/(2*plotmargin+plotsizey)+1.0*plotsizey/(2*plotmargin+plotsizey)-0.17*(plotsizex*plotsizey)**0.5/(2*plotmargin+plotsizey), 0.16*(plotsizex*plotsizey)**0.5/(2*plotmargin+plotsizex), 0.16*(plotsizex*plotsizey)**0.5/(2*plotmargin+plotsizey)], projection='polar')
        norm = mpl.colors.Normalize(-np.pi, np.pi)
        quant = 2056
        plt.grid(False)
        cw = mpl.colorbar.ColorbarBase(colorwheel_axes, cmap=mpl.cm.get_cmap('hsv',quant), norm=norm, orientation='horizontal')
        cw.outline.set_visible(False)
        colorwheel_axes.set_axis_off()
        colorwheel_axes.set_rlim([-1,1])
        myrec=plt.Rectangle((lx-0.18*lxy, ly-0.18*lxy), 0.18*lxy, 0.18*lxy, facecolor='w', edgecolor='k', fill=True, figure=fig)
        plot1.add_artist(myrec)
    #write label if desired
    if (label!=False):
        myrec2=plt.Rectangle((lx-0.09*lxy, 0.0), 0.09*lxy, 0.09*lxy, facecolor='w', edgecolor='k', fill=True, figure=fig)
        plot1.add_artist(myrec2)
        plot1.text(lx-0.08*lxy, 0.022*lxy, '$'+label+'$')
    plot1.set_xticks([0.0, lx])
    plot1.set_yticks([0.0, ly])
    plot1.set_yticklabels(['$0$', '$'+str(ly)+'$'], rotation=90)
    plot1.set_xticklabels(['$0$', '$'+str(lx)+'$'])
    if (save_filename==False):
        plt.show()
    else:
        plt.savefig(save_filename, dpi=300)
    plt.close()



label=False
colorwheel=False
save_filename=False
try:
    opts, args=getopt.getopt(sys.argv[1:], "hws:l:f:a:p:")
except getopt.GetoptError:
    print ('create_snapshot.py -f <data file name>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('create_snapshot.py -f <data file name>')
        print('-s <plot file name> option for saving the plot as a file')
        print('-w option to include color wheel')
        print('-l <label> option to include label')
        print('-a <area> option to change the size (area) of the plot, default is 5.28')
        print('-p <point size factor> changing the size of points (default=1.0)')
        sys.exit()
    elif opt=='-f':
        file= arg
    elif opt=='-w':
        colorwheel=True
    elif opt=='-s':
        save_filename=arg
    elif opt=='-l':
        label=arg
    elif opt=='-a':
        plotarea=float(arg)
    elif opt=='-p':
        point_size_factor=float(arg)

show_snapshot(file, colorwheel, label, save_filename)

