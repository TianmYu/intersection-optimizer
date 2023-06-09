import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors

TRAVEL_CONSTANT = 0

def EVoptionB(totalTime, crossingTime, EWL, NSL):
    # cross at intersection both ways
    p = NSL/totalTime
    A = 2*crossingTime
    B = 2*crossingTime + EWL
    return A*p + (1/2) * (1-p)*(B**2-A**2)/(B-A)

def EVoptionA(totalTime, crossingTime, EWL, NSL):
    # cross at intersection both ways
    p = NSL/totalTime
    A = 2*crossingTime
    B2 = crossingTime + NSL
    B3 = crossingTime + EWL
    return (1/totalTime)*((B2**2-A**2)/2+(B3**2-A**2)/2+2*A*crossingTime)

def sampleEVFn(evFunction, axisRanges, dim, crossingTime):
    # ranges: range of time to sweep
    # dim: dimension of array to use

    x = np.linspace(axisRanges[0], axisRanges[1], dim)
    y = np.linspace(axisRanges[2], axisRanges[3], dim)

    evArr = np.zeros((dim, dim))
    print(evArr.shape)

    for i, xval in enumerate(x):
        for j, yval in enumerate(y):
            nsl = xval
            ewl = yval

            totalTime = nsl+ewl+2*crossingTime

            ev = evFunction(totalTime, crossingTime, ewl, nsl)

            evArr[j,i] = ev

    X, Y = np.meshgrid(x, y)

    return X, Y, evArr

def plotHeatmapContour(title, X, Y, ev, subplotNum, timeRanges, fig, ax, vmin = None, vmax = None, cmap=None):
    plt.subplot(subplotNum)
    if cmap is not None:
        pc = plt.imshow(ev, norm=colors.CenteredNorm(), cmap=cmap, origin='lower', extent=timeRanges)
    else:
        pc = plt.imshow(ev,  origin='lower', vmin=vmin, vmax=vmax, extent=timeRanges)
    fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    plt.xlabel("ns light on time")
    plt.ylabel("ew light on time")
    plt.title(title)

    contours = plt.contour(X, Y, ev, colors='black')
    plt.clabel(contours, inline=True, fontsize=8)
    return True

def compareEvs(axisRanges, dim, crossingTime, minTime, maxTime):
    # getting for option a
    Xa, Ya, evA = sampleEVFn(EVoptionA, axisRanges, dim, crossingTime)

    # getting for option b
    Xb, Yb, evB = sampleEVFn(EVoptionB, axisRanges, dim, crossingTime)

    diffEv = evA - evB

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1)
    fig.tight_layout(pad=2.0)
    cmap = plt.cm.coolwarm

    plotHeatmapContour("A Mean - B Mean", Xa, Ya, diffEv, 131, axisRanges, fig, ax1, cmap=cmap)
    plotHeatmapContour("A EV", Xa, Ya, evA, 132, axisRanges, fig, ax2, vmin = minTime, vmax = maxTime/1.8)
    plotHeatmapContour("B EV", Xb, Yb, evB, 133, axisRanges, fig, ax3, vmin = minTime, vmax = maxTime/1.8)
    plt.show()
    return True

def main():
    minT = 10
    maxT = 240
    contourLevels = 10
    increment = 5
    numT = int((maxT-minT)/increment)+1
    nsTimes = np.linspace(minT, maxT, numT)
    ewTimes = np.linspace(minT, maxT, numT)

    walkSpeed = 2 # about 7.2 km/h
    crossWidth = 20 # 20 m for a 10s crossing
    deadTime = 6
    minTime = 2*(crossWidth/walkSpeed) # fastest possible path time
    maxTime = 2*(crossWidth/walkSpeed)+2*deadTime+np.maximum(np.amax(nsTimes), np.amax(ewTimes))

    timeRanges = [minT, maxT, minT, maxT]

    compareEvs(timeRanges, numT, crossWidth/walkSpeed, minTime, maxTime)

main()