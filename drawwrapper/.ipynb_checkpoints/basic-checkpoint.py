import matplotlib
if 'plt' in globals():
    pass
elif 'inline' in matplotlib.get_backend():
    pass
else:
    matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from numpy import histogram

from mymodule import *
from mymodule.linearft import *

class PlotOptions(object):
    def __init__(self,label=20,legend=20,ticks=20,mfc='k',axisscale='normal',saveshow='save',rasvec=False):
        self.label = label
        self.legend = legend
        self.ticks = ticks
        self.mfc = mfc
        self.axisscale = axisscale
        self.rasvec = rasvec
        self.saveshow = saveshow

def color_line_style_generator():
    linestyle = ['-', '--', '-.', ':']
    markerstyle = ['h', '^', 'v', 's', '2', '<', '>', '3', '4', '8', 'p', '*', 'H', '+', ',', '.', 'x', 'o', 'D', 'd', '|', '_']
    mfcstyle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    line_idx = 0
    marker_idx = 0
    mfc_idx = 0
    while True:
        yield mfcstyle[mfc_idx] + linestyle[line_idx] + markerstyle[marker_idx]
        line_idx = (line_idx + 1) % len(linestyle)
        marker_idx = (marker_idx + 1) % len(markerstyle)
        mfc_idx = (mfc_idx + 1) % len(mfcstyle)

def monochrome_style_generator():
    linestyle = ['-', '--', '-.', ':']
    markerstyle = ['h', '^', 'v', 's', '2', '<', '>', '2', '3', '4', '8', 'p', '*', 'H', '+', ',', '.', 'x', 'o', 'D', 'd', '|', '_']
    line_idx = 0
    marker_idx = 0
    while True:
        yield 'k' + linestyle[line_idx] + markerstyle[marker_idx]
        line_idx = (line_idx + 1) % len(linestyle)
        marker_idx = (marker_idx + 1) % len(markerstyle)

def points_style_generator():
    markerstyle = ['h', '^', 'v', 's', '2', '<', '>', '2', '3', '4', '8', 'p', '*', 'H', '+', ',', '.', 'x', 'o', 'D', 'd', '|', '_']
    mfcstyle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']#, 'w']
    marker_idx = 0
    mfc_idx = 0
    while True:
        yield mfcstyle[mfc_idx] + markerstyle[marker_idx]
        marker_idx = (marker_idx + 1) % len(markerstyle)
        mfc_idx = (mfc_idx + 1) % len(mfcstyle)

def monochrome_points_style_generator():
    markerstyle = ['h', '^', 'v', 's', '2', '<', '>', '2', '3', '4', '8', 'p', '*', 'H', '+', ',', '.', 'x', 'o', 'D', 'd', '|', '_']
    marker_idx = 0
    while True:
        yield 'k' + markerstyle[marker_idx]
        marker_idx = (marker_idx + 1) % len(markerstyle)

def white_points_style_generator():
    markerstyle = ['h', '^', 'v', 's', '2', '<', '>', '2', '3', '4', '8', 'p', '*', 'H', '+', ',', '.', 'x', 'o', 'D', 'd', '|', '_']
    marker_idx = 0
    while True:
        yield 'w' + markerstyle[marker_idx]
        marker_idx = (marker_idx + 1) % len(markerstyle)

def SetPlotter(axisscale='normal'):
    if axisscale == 'normal':
        myplot = plt.plot
    elif axisscale == 'log':
        myplot = plt.loglog
    elif axisscale == 'semilogx':
        myplot = plt.semilogx
    elif axisscale == 'semilogy':
        myplot = plt.semilogy
    return myplot

def SetPlotterWithAxe(ax,axisscale='normal'):
    if axisscale == 'normal':
        myplot = ax.plot
    elif axisscale == 'log':
        myplot = ax.loglog
    elif axisscale == 'semilogx':
        myplot = ax.semilogx
    elif axisscale == 'semilogy':
        myplot = ax.semilogy
    return myplot

def DrawSecondAxsisGraph(xdata, ydatas, linetypes, xlabel, ylabels,savefilename, po=PlotOptions()):
# def DrawSecondAxsisGraph(xdata, ydatas, linetypes, xlabel, ylabels,savefilename, fs=PlotOptions(), makerfacecolor='k', axisscale='normal', rasvec=False):
    """
    data : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    """
    # myplotter = SetPlotter(axisscale)
    # myplotter(data, linetype, rasterized=rasvec)
    
    # plt.figure(figsize=(6,3))
    
    fig, ax1 = plt.subplots()
    fig.subplots_adjust(right=0.85)
    ln1 = ax1.plot(xdata, ydatas[0], linetypes[0], label=ylabels[0], mfc=po.mfc)
    ax1.set_xlabel(xlabel,fontsize=po.label)
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel(ylabels[0],fontsize=po.label)

    ax2 = ax1.twinx()
    ln2 = ax2.plot(xdata, ydatas[1], linetypes[1], label=ylabels[1], mfc=po.mfc)
    ax2.set_ylabel(ylabels[1],fontsize=po.label)

    lns = ln1+ln2
    labs = [ l.get_label() for l in lns ]
    
    lg = ax1.legend(lns, labs, loc='best',fontsize=po.legend)
    lg.draw_frame(False)

    # plt.xticks(fontsize=20)
    # plt.yticks(fontsize=20)
    # plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawHistogram(data, nbin, xlabel, ylabel, histlegend, savefilename, pdf=False, po=PlotOptions(20,20,20,'gray','normal',False)):
# def DrawHistogram(data, nbin, xlabel, ylabel, histlegend, savefilename, po=PlotOptions(), axisscale='normal', pdf=False, rasvec=False):
    """
    data : input what you want to draw
    nbin : number of bins 
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    """
    plt.subplots_adjust(bottom=0.15)
    logflag=False
    if po.axisscale=='semilogy':
        logflag = True
    myplotter = SetPlotter(po.axisscale)
    histdata, xbins = histogram(data,bins=nbin,density=pdf)
    dx, offset, xdatafit, xdataplot = histinfo(xbins, nbin)
    if po.axisscale=='log':
        myplotter(xdataplot,histdata,'x',rasterized=po.rasvec,mfc=po.mfc)
    else:
        plt.bar(xdataplot,histdata,width=dx,rasterized=po.rasvec,log=logflag,fc=po.mfc)
    # plt.bar(xdataplot,histdata,width=dx,rasterized=po.rasvec,log=logflag,fc=po.mfc)
    plt.legend(loc='best',fontsize=po.legend)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawHistogramMultidata(data, num, nbin, xlabel, ylabel, histlegend, savefilename, linetypes=False, pdf=False, po=PlotOptions(20,20,20,'gray','normal',False)):
# def DrawHistogram(data, nbin, xlabel, ylabel, histlegend, savefilename, po=PlotOptions(), axisscale='normal', pdf=False, rasvec=False):
    """
    data : input what you want to draw
    num  : number of data
    nbin : number of bins 
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    """
    plt.subplots_adjust(bottom=0.15)
    logflag=False
    if po.axisscale=='semilogy':
        logflag = True
    myplotter = SetPlotter(po.axisscale)
    for i in range(num):
        histdata, xbins = histogram(data[i],bins=nbin,density=pdf)
        dx, offset, xdatafit, xdataplot = histinfo(xbins, nbin)
        if linetypes == False:
            plt.bar(xdataplot,histdata,width=dx,rasterized=po.rasvec,log=logflag,fc=po.mfc)
        else:
            myplotter(xdataplot,histdata,linetypes[i],rasterized=po.rasvec,mfc=po.mfc)
    plt.legend(loc='best',fontsize=po.legend)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawHistDataWithFittingLine(xdataplot, dx, histdata, xdatafit, fitdata, xlabel, ylabel, dataname, fitname, savefilename, po=PlotOptions()):
    """
    data : input what you want to draw
    nbin : number of bins 
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    """
    # myplotter = SetPlotter(po.axisscale)
    plt.subplots_adjust(bottom=0.15)
    logflag=False
    myplotter = SetPlotter(po.axisscale)
    if po.axisscale=='semilogy':
        logflag = True
        plt.bar(xdataplot,histdata,width=dx,label=dataname,rasterized=po.rasvec,log=logflag,fc=po.mfc,alpha=0.5)
    elif po.axisscale=='log':
        myplotter(xdataplot,histdata,label=dataname,rasterized=po.rasvec,mfc=po.mfc,alpha=0.5)
    else:
        plt.bar(xdataplot,histdata,width=dx,label=dataname,rasterized=po.rasvec,log=logflag,fc=po.mfc,alpha=0.5)
    # histdata, xbins = histogram(data,bins=nbin,density=pdf)
    # dx, offset, xdatafit, xdataplot = histinfo(xbins, nbin)
    myplotter(xdataplot+0.5*dx,fitdata,'k-',label=fitname)
    # myplotter(xdataplot,histdata,'kx',label=histlegend,rasterized=po.rasvec, mfc=po.mfc)
    lg = plt.legend(loc='best',fontsize=po.legend)
    lg.draw_frame(False)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.tight_layout()
    if po.saveshow=='save':
        plt.savefig(savefilename)
    elif po.saveshow=='show':
        plt.show()
    plt.close()

def Draw_heatmap(xdata, ydata, nbins, xlabel, ylabel, savefilename, po=PlotOptions(), pdf=False):
# def Draw_heatmap(xdata, ydata, nbins, xlabel, ylabel, savefilename, po=PlotOptions(), pdf=False, rasvec=False):
    """
    xdata : input what you want to draw
    ydata : input what you want to draw
    nbin : number of bins 
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    """
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # ax.hist2d(xdata, ydata, bins=nbins)
    
    heatmap, xedges, yedges, img = ax.hist2d(xdata, ydata, bins=nbins)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    ax.imshow(heatmap, extent=extent)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    ax.set_aspect('auto')
    plt.tight_layout()
    plt.savefig(savefilename)
    ax.clear()
    plt.close()
    
def DrawMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), locstr='best', **kwargs):
# def DrawMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), axisscale='normal', rasvec=False, locstr='best'):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    legens : 
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    locster : set location of label in figure. A default value is 'best'
    """
    # fig=plt.figure(figsize=(10,10))
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(111)
    myplotter = SetPlotter(po.axisscale)
    [ myplotter(xdata, ydata[i], linetypes[i], label=legends[i], rasterized=po.rasvec, mfc=po.mfc, **kwargs) for i in range(yn) ]
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    
    # plt.vlines([0.75, 1.1, 2.5], 0.0, 1.0, linestyles="dashdot")
    
    # def DrawRegion(a,b,fc,ax):
    #     ix = numpy.linspace(a, b)
    #     iy = numpy.ones(len(ix))
    #     verts = [(a, 0)] + list(zip(ix, iy)) + [(b, 0)]
    #     poly = matplotlib.patches.Polygon(verts, facecolor=fc, edgecolor='0.0',alpha=0.2)
    #     ax.add_patch(poly)
    # hight = 1.03
    # [a,b]=(0,0.75)
    # DrawRegion(a,b,'0.2',ax)
    # plt.text(0.3, hight, r"$R_I$", horizontalalignment='center', fontsize=18, backgroundcolor='w')
    # [a,b]=(0.75,1.1)
    # DrawRegion(a,b,'0.4',ax)
    # plt.text(0.92, hight, r"$R_{II}$", horizontalalignment='center', fontsize=18, backgroundcolor='w')
    # [a,b]=(1.1,2.5)
    # DrawRegion(a,b,'0.7',ax)
    # plt.text(1.8, hight, r"$R_{III}$", horizontalalignment='center', fontsize=18, backgroundcolor='w')
    # [a,b]=(2.5,4.5)
    # DrawRegion(a,b,'0.9',ax)
    # plt.text(3.5, hight, r"$R_{IV}$", horizontalalignment='center', fontsize=18, backgroundcolor='w')
    
    # def DrawRegionGradient(center,minx,maxx,mycm,vma=1.0,vmi=None):
    #     def paraf(a,b,x,y):
    #         return -b*(x-a)**2
    #     dx, dy = 0.01, 0.005
    #     x = arange(minx, maxx, dx)
    #     y = arange(0.0, 1.0, dy)
    #     X,Y = numpy.meshgrid(x, y)
        
    #     xmin, xmax, ymin, ymax = numpy.amin(x), numpy.amax(x), numpy.amin(y), numpy.amax(y)
    #     extent = xmin, xmax, ymin, ymax
    #     Z2 = paraf(minx+(maxx-minx)*0.5,2.0/(maxx-minx), X, Y)
    #     im2 = plt.imshow(Z2, cmap=mycm, alpha=0.8, interpolation='bilinear', extent=extent, vmin=vmi, vmax=vma)

    # DrawRegionGradient(0.375,0.3,0.75,plt.cm.Reds)
    # DrawRegionGradient(0.925,0.75,1.1,plt.cm.Blues)
    # DrawRegionGradient(1.75,1.1,2.4,plt.cm.Greens)
    # def DrawRegionGradient2(center,minx,maxx,mycm,vma=1.0,vmi=None):
    #     def paraf(a,b,x,y):
    #         return b*(x-a)**2
    #     dx, dy = 0.01, 0.005
    #     x = arange(minx, maxx, dx)
    #     y = arange(0.0, 1.0, dy)
    #     X,Y = numpy.meshgrid(x, y)
        
    #     xmin, xmax, ymin, ymax = numpy.amin(x), numpy.amax(x), numpy.amin(y), numpy.amax(y)
    #     extent = xmin, xmax, ymin, ymax
    #     Z2 = paraf(minx+(maxx-minx)*0.5,2.0/(maxx-minx), X, Y)
    #     im2 = plt.imshow(Z2, cmap=mycm, alpha=0.8, interpolation='bilinear', extent=extent, vmin=vmi, vmax=vma)
    # DrawRegionGradient2(3.45,2.4,4.5,plt.cm.afmhot,1.0,-1.0)
    # plt.axis([0,4.5,0,1])
    # ax.set_aspect('auto')

    lg = plt.legend(loc=locstr,fontsize=po.legend)
    lg.draw_frame(False)
    # ax=plt.gca()

    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    
    # plt.ylim([0.99,1.01])
    # aspect = 0.9*(ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    # ax.set_aspect(aspect)
    plt.tight_layout()
    
    plt.savefig(savefilename)
    plt.close()

def DrawMultiDataOutLegend(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), figasp=(10,10), locstr='best', **kwargs):
# def DrawMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), axisscale='normal', rasvec=False, locstr='best'):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    legens : 
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    locster : set location of label in figure. A default value is 'best'
    """
    fig=plt.figure(figsize=(figasp[0],figasp[1]))
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(111)
    myplotter = SetPlotter(po.axisscale)
    [ myplotter(xdata, ydata[i], linetypes[i], label=legends[i], rasterized=po.rasvec, mfc=po.mfc, **kwargs) for i in range(yn) ]
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)

    plt.tight_layout() # this command must be above commands to draw legend outside of a main figure
    
    lg = plt.legend(bbox_to_anchor=(1.05,1), loc=locstr,fontsize=po.legend)
    lg.draw_frame(False)
    plt.subplots_adjust(right=0.7)
    plt.savefig(savefilename)
    plt.close()

def DrawBothMultiDataOutLegend(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), figasp=(10,10), locstr='best', **kwargs):
# def DrawMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), axisscale='normal', rasvec=False, locstr='best'):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    legens : 
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    locster : set location of label in figure. A default value is 'best'
    """
    fig=plt.figure(figsize=(figasp[0],figasp[1]))
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(111)
    myplotter = SetPlotter(po.axisscale)
    [ myplotter(xdata[i], ydata[i], linetypes[i], label=legends[i], rasterized=po.rasvec, mfc=po.mfc, **kwargs) for i in range(yn) ]
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)

    plt.tight_layout() # this command must be above commands to draw legend outside of a main figure
    
    lg = plt.legend(bbox_to_anchor=(1.05,1), loc=locstr,fontsize=po.legend)
    lg.draw_frame(False)
    plt.subplots_adjust(right=0.7)
    plt.savefig(savefilename)
    plt.close()

def DrawBothMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), locstr='best', **kwargs):
# def DrawBothMultiData(xdata, yn, ydata, linetypes, xlabel, ylabel, legends, savefilename, po=PlotOptions(), axisscale='normal', rasvec=False, locstr='best'):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    legens : 
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    rasvec : chose raster plot (True) or vector plot (False, default)
    locster : set location of label in figure. A default value is 'best'
    """
    myplotter = SetPlotter(po.axisscale)
    [ myplotter(xdata[i], ydata[i], linetypes[i], label=legends[i], rasterized=po.rasvec, mfc=po.mfc, **kwargs) for i in range(yn) ]
    # plt.xticks(fontsize=25)
    # plt.yticks(fontsize=25)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    lg = plt.legend(loc=locstr,fontsize=po.legend)
    # lg.draw_frame(False)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawSingleDataNolegend(data,linetype,xlabel,ylabel,savefilename, po=PlotOptions()):
# def DrawSingleDataNolegend(data,linetype,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False):
    """
    data : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    """
    myplotter = SetPlotter(po.axisscale)
    myplotter(data, linetype, rasterized=po.rasvec, mfc=po.mfc)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawSingleDataWithlegend(data,linetype,labelstr,xlabel,ylabel,savefilename, po=PlotOptions(), locstr='best'):
# def DrawSingleDataWithlegend(data,linetype,labelstr,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False,locstr='best'):
    """
    data : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    locster : set location of label in figure. A default value is 'best'
    """
    myplotter = SetPlotter(po.axisscale)
    myplotter(data, linetype, label=labelstr, rasterized=po.rasvec, mfc=po.mfc)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    lg = plt.legend(loc=locstr,fontsize=po.legend)
    lg.draw_frame(False)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawXYDataWithlegend(xdata,ydata,linetype,labelstr,xlabel,ylabel,savefilename, po=PlotOptions(), locstr='best'):
# def DrawXYDataWithlegend(xdata,ydata,linetype,labelstr,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False,locstr='best'):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    locster : set location of label in figure. A default value is 'best'
    """
    myplotter = SetPlotter(po.axisscale)
    myplotter(xdata, ydata, linetype, label=labelstr, rasterized=po.rasvec, mfc=po.mfc)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.legend(loc=locstr,fontsize=po.legend)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawXYDataNolegend(xdata,ydata,linetype,xlabel,ylabel,savefilename, po=PlotOptions()):
# def DrawXYDataNolegend(xdata,ydata,linetype,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    locster : set location of label in figure. A default value is 'best'
    """
    myplotter = SetPlotter(po.axisscale)
    myplotter(xdata, ydata, linetype, rasterized=po.rasvec, mfc=po.mfc)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawXYDataNolegendLimitRange(xdata,ydata,xlim,ylim,linetype,xlabel,ylabel,savefilename, po=PlotOptions()):
# def DrawXYDataNolegend(xdata,ydata,linetype,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    locster : set location of label in figure. A default value is 'best'
    """
    myplotter = SetPlotter(po.axisscale)
    myplotter(xdata, ydata, linetype, rasterized=po.rasvec, mfc=po.mfc)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.xlim(xlim)
    plt.ylim(ylim)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

def DrawXYDataNolegendWithErrBar(xdata,ydata,yerrdata,linetype,xlabel,ylabel,savefilename, po=PlotOptions()):
# def DrawXYDataNolegendWithErrBar(xdata,ydata,yerrdata,linetype,xlabel,ylabel,savefilename, po=PlotOptions(), axisscale='normal', rasvec=False):
    """
    xdata, ydata : input what you want to draw
    linetype : define the plot data style. ex) '-', 'o', ...
    xlabel, ylabel : set the x and y axis label respectively
    savefilename : set filename of figure you draw. It must include the file type. ex) .eps, .png, ...
    locster : set location of label in figure. A default value is 'best'
    """
    # myplotter = SetPlotter(axisscale)
    # myplotter(xdata, ydata,yerr, linetype, rasterized=rasvec)
    plt.errorbar(xdata, ydata, yerr=yerrdata, fmt=linetype, rasterized=po.rasvec)
    plt.xticks(fontsize=po.ticks)
    plt.yticks(fontsize=po.ticks)
    plt.xlabel(xlabel,fontsize=po.label)
    plt.ylabel(ylabel,fontsize=po.label)
    plt.tight_layout()
    plt.savefig(savefilename)
    plt.close()

