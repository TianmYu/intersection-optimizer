import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors

PLOT_DEBUG = True

# problem definition: you are travelling north - south
# you can cross east west now, or wait until the 4 way intersection to cross both
# which is optimal and under what conditions is it optimal?
# assumption: constant walk speed, equal width roads

class IntersectionCrossing:
    def __init__(self, walkSpeed, crossWidth, eastWestOnTime, northSouthOnTime, deadTime):
        self.walkSpeed = walkSpeed # in m/s
        self.crossWidth = crossWidth # in m
        self.crossTime = self.crossWidth / self.walkSpeed
        self.ewt = eastWestOnTime # in s
        self.nst = northSouthOnTime
        self.deadTime = deadTime
        self.totalTime = eastWestOnTime + northSouthOnTime + 2*deadTime

        self.rng = np.random.default_rng(42)

        # north south starts at 0, ends at self.nst
        self.eastWestStart = northSouthOnTime+deadTime
        self.eastWestEnd = northSouthOnTime+deadTime+eastWestOnTime

        # order:
        # north south on
        # both off, east west next
        # east west on
        # both off, north south next
        self.states = np.array([0,
                                northSouthOnTime/self.totalTime
                                ,(self.eastWestStart)/self.totalTime
                                ,(self.eastWestEnd)/self.totalTime
                                ,1 ])

    def timeFromState(self, distribArr, state, startingPos):
        # given state type 0, 1, 2, 3 return the time required to cross the intersection for each element in array
        # startingPos: 0=top left (path A), 1=top right (path B, greedy)

        timeArr = distribArr * self.totalTime
        cross0WaitTime = np.zeros(len(distribArr)) # time spent waiting at top left
        cross1WaitTime = np.zeros(len(distribArr))-1 # time spent waiting at top right
        cross2WaitTime = np.zeros(len(distribArr))-1 # time spent waiting at bottom right
        # -1 means the intersection is not taken, and is ignored for statistics later

        if startingPos == 0:
            if state == 0: # we arrive when north/south is on
                # time is time to cross north/south + remaining time before light change + time to cross east west
                remTime = self.eastWestStart-timeArr # remaining time needed to wait before east west turns on

                cross2WaitTime = np.maximum(remTime - self.crossTime, 0)

                crossWaitTime = np.maximum(remTime, self.crossTime) 
                # amount of time needed to wait for light to change
                # we either need to wait after crossing, or the light changes as we cross

                totalTime = crossWaitTime + self.crossTime

            elif state == 1: # we arrive in the dead time when north/south is off and east/west is yet to change
                # first cross east/west. need to wait for it to turn green:
                deadTime = self.eastWestStart-timeArr
                cross0WaitTime = deadTime
                # next wait for north/south to turn green. Time added is time between ew turning green and ns turning green
                # note that ns turns green at totalTime (end of cycle)
                # if the light changes before we can cross, account for it
                crossWaitTime = np.maximum(self.totalTime - self.eastWestStart, self.crossTime) 
                cross1WaitTime = np.maximum((self.totalTime - self.eastWestStart) - self.crossTime,0)
               
                # add deadtime, ew cross time, and ns cross time
                totalTime = deadTime + crossWaitTime + self.crossTime

            elif state == 2: # we arrive when east/west is on
                # time is time to cross east/west + remaining time before light change + time to cross north/south
                remTime = self.totalTime-timeArr # remaining time needed to wait before north/south turns on

                cross1WaitTime = np.maximum(remTime - self.crossTime,0)
                crossWaitTime = np.maximum(remTime, self.crossTime) 
                # amount of time needed to wait for light to change
                # we either need to wait after crossing, or the light changes as we cross

                totalTime = crossWaitTime + self.crossTime

            elif state == 3: # we arrive in the dead time when e/w is off and n/s is yet to change
                # first cross n/s. need to wait for it to turn green:
                deadTime = self.totalTime-timeArr
                cross0WaitTime = deadTime
                # next wait for e/w to turn green. Time added is time between ns turning green and ew turning green
                # if we cant cross ns fast enough, use ns a rate limiting step
                crossWaitTime = np.maximum(self.eastWestStart, self.crossTime)
                cross2WaitTime = np.maximum((self.eastWestStart) - self.crossTime,0)
                # if the light changes before we can cross, account for it
                # sum initial dead time, crossing ew and crossing ns

                totalTime = deadTime + crossWaitTime + self.crossTime
            
        if startingPos == 1:
            remTime = self.totalTime-timeArr # time until the next n/s green light
            remTime[remTime > (self.totalTime-self.nst)] = 0 # checks if we've arrived on a green light
            
            cross1WaitTime = remTime
            totalTime = remTime + 2*self.crossTime # account for crossing time at current intersection, and before at greedy

        return totalTime, cross0WaitTime, cross1WaitTime, cross2WaitTime

    def sampleMonteCarlo(self, numSamples):
        sampleArr = self.rng.random(size = numSamples)
        stateDistrib = [[],[],[],[]]

        for i in range(1, len(self.states)):
            lowerBound = self.states[i-1]
            upperBound = self.states[i]

            stateDistrib[i-1] = sampleArr[np.logical_and(sampleArr>=lowerBound, sampleArr<upperBound)]

        # get times for path A (not greedy) or path B (greedy - cross immediately)
        ATimes = np.zeros(numSamples)
        BTimes = np.zeros(numSamples)

        ACross0s = np.zeros(numSamples)
        ACross1s = np.zeros(numSamples)
        ACross2s = np.zeros(numSamples)
        BCross0s = np.zeros(numSamples)
        BCross1s = np.zeros(numSamples)
        BCross2s = np.zeros(numSamples)

        prevEnd = 0
        for i, distribArr in enumerate(stateDistrib):
            ATimeArr, ACross0, ACross1, ACross2 = self.timeFromState(distribArr, i, 0)
            BTimeArr, BCross0, BCross1, BCross2 = self.timeFromState(distribArr, i, 1)

            currEnd = prevEnd + len(distribArr)

            ATimes[prevEnd:currEnd] = ATimeArr
            BTimes[prevEnd:currEnd] = BTimeArr
            
            ACross0s[prevEnd:currEnd] = ACross0
            ACross1s[prevEnd:currEnd] = ACross1
            ACross2s[prevEnd:currEnd] = ACross2
            BCross0s[prevEnd:currEnd] = BCross0
            BCross1s[prevEnd:currEnd] = BCross1
            BCross2s[prevEnd:currEnd] = BCross2

            prevEnd = currEnd

        return ATimes, BTimes, ACross0s, ACross1s, ACross2s, BCross0s, BCross1s, BCross2s

def contourSetup(axisRanges, valueArr):
    dim = valueArr.shape[0]
    x = np.linspace(axisRanges[0], axisRanges[1], dim)
    y = np.linspace(axisRanges[2], axisRanges[3], dim)

    X, Y = np.meshgrid(x, y)
    return X, Y, valueArr

def plotNonDiffEV(title, subplotNum, valArr, vmin, vmax, extent, fig, ax):
    plt.subplot(subplotNum)
    pc = plt.imshow(valArr,  origin='lower', vmin=vmin, vmax=vmax, extent=extent)
    fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    plt.xlabel("ns light on time")
    plt.ylabel("ew light on time")
    plt.title(title)

    x, y, z = contourSetup(extent, valArr)
    contours = plt.contour(x, y, z, colors='black')
    plt.clabel(contours, inline=True, fontsize=8)

def iterateAndPlot():
    # properties
    minT = 30
    maxT = 60
    increment = 30
    numT = int((maxT-minT)/increment)+1
    nsTimes = np.linspace(minT, maxT, numT)
    ewTimes = np.linspace(minT, maxT, numT)

    walkSpeed = 2 # about 7.2 km/h
    crossWidth = 20 # 20 m for a 10s crossing
    deadTime = 6
    numSamples = 10000
    minTime = 2*(crossWidth/walkSpeed) # fastest possible path time
    maxTime = 2*(crossWidth/walkSpeed)+2*deadTime+np.maximum(np.amax(nsTimes), np.amax(ewTimes))

    # save mean and std dev initially
    # we can also model as 2 uniform PDFs, with 1 at the optimal time, and the other spread across all other times

    expectedDt = np.zeros((numT, numT))

    expectedA = np.zeros((numT, numT))
    stdDevA = np.zeros((numT, numT))
    expectedB = np.zeros((numT, numT))
    stdDevB = np.zeros((numT, numT))

    ABinomial = np.zeros((numT, numT))
    BBinomial = np.zeros((numT, numT))

    AUniformRange = np.zeros((numT, numT))
    BUniformRange = np.zeros((numT, numT))

    total = numT**2
    counter = 0
    for k, nsTime in enumerate(nsTimes):
        for j, ewTime in enumerate(ewTimes):
            crossing = IntersectionCrossing(walkSpeed, crossWidth, nsTime, ewTime, deadTime)
            results = crossing.sampleMonteCarlo(numSamples)
            # results is tuple with A time, B time, A cross wait time at 0, 1, 2, B cross wait time at 0, 1, 2

            optimalSamplesA = np.isclose(results[0], minTime)
            numOptimalA = np.count_nonzero(optimalSamplesA)
            ABinomial[k][j] = numOptimalA/numSamples

            optimalSamplesB = np.isclose(results[1], minTime)
            numOptimalB = np.count_nonzero(optimalSamplesB)
            BBinomial[k][j] = numOptimalB/numSamples

            AUniformRange[k][j] = np.amax(results[0]) - np.amin(results[0])
            BUniformRange[k][j] = np.amax(results[1]) - np.amin(results[1])

            resultsProps = []
            for i, result in enumerate(results):
                mean = np.mean(result)
                stdDev = np.std(result)
                resultsProps += [(mean, stdDev)]

            AMean = resultsProps[0][0]
            AStd = resultsProps[0][1]
            BMean = resultsProps[1][0]
            BStd = resultsProps[1][1]
            Diff = AMean - BMean

            expectedA[k][j] = AMean
            stdDevA[k][j] = AStd
            expectedB[k][j] = BMean
            stdDevB[k][j] = BStd

            expectedDt[k][j] = Diff

            counter += 1
            if counter%100 == 0:
                print("processed: {}/{}".format(counter, total))

            if PLOT_DEBUG:
                ATimes = results[0]
                BTimes = results[1]

                p = np.random.permutation(len(results[0]))
                shuffledResults = []
                for result in results:
                    shuffledResults += [result[p]]

                #for i in range(0,100):
                #   print("time: {}, nw: {}, ne: {}, sw: {}".format(shuffledResults[0][i], shuffledResults[2][i], shuffledResults[3][i], shuffledResults[4][i]))

                print("plotting for ns time: {}, ew time: {}".format(nsTime, ewTime))

                plt.subplot(424)
                n, bins, patches = plt.hist(ATimes, bins=25, alpha=0.75, color="g")
                print("total A size: {}".format(ATimes.shape))
                plt.ylabel('Approach A Times')
                plt.grid(True)
                plt.subplot(421)
                n, bins, patches = plt.hist(results[2][results[2]>=0], bins=25, alpha=0.75)
                print("total NW A size: {}".format(results[2][results[2]>=0].shape))
                plt.ylabel('NW A')
                plt.grid(True)
                plt.subplot(422)
                n, bins, patches = plt.hist(results[3][results[3]>=0], bins=25, alpha=0.75)
                print("total NE A size: {}".format(results[3][results[3]>=0].shape))
                plt.ylabel('NE A')
                plt.grid(True)
                plt.subplot(423)
                n, bins, patches = plt.hist(results[4][results[4]>=0], bins=25, alpha=0.75)
                print("total SW A size: {}".format(results[4][results[4]>=0].shape))
                plt.ylabel('SW A')
                plt.grid(True)
                
                plt.subplot(428)
                n, bins, patches = plt.hist(BTimes, bins=25, alpha=0.75, color="g")
                plt.xlabel('Time')
                plt.ylabel('Approach B Times')
                plt.grid(True)
                plt.subplot(426)
                n, bins, patches = plt.hist(results[6], bins=25, alpha=0.75)
                plt.xlabel('Time')
                plt.ylabel('NE B')
                plt.grid(True)
                plt.show()

    timeRanges = [minT, maxT, minT, maxT]

    #fig, ((ax1, ax2, ax3),(ax4, ax5, ax6),(ax7, ax8, ax9)) = plt.subplots(ncols=3, nrows=3)
    fig, ((ax1, ax2, ax3)) = plt.subplots(ncols=3, nrows=1)
    fig.tight_layout(pad=2.0)
    cmap = plt.cm.coolwarm
    plt.subplot(131)
    pc = plt.imshow(expectedDt, norm=colors.CenteredNorm(), cmap=cmap, origin='lower', extent=timeRanges)
    fig.colorbar(pc, ax=ax1, fraction=0.046, pad=0.04)
    plt.xlabel("ns light on time")
    plt.ylabel("ew light on time")
    plt.title("A Mean - B Mean")
    x, y, z = contourSetup(timeRanges, expectedDt)
    contours = plt.contour(x, y, z, colors='black')
    plt.clabel(contours, inline=True, fontsize=8)

    plotNonDiffEV("A EV", 132, expectedA, minTime, maxTime/1.8, timeRanges, fig, ax2)
    plotNonDiffEV("B EV", 133, expectedB, minTime, maxTime/1.8, timeRanges, fig, ax3)
    plt.show()

    fig, ((ax1, ax2)) = plt.subplots(ncols=2, nrows=1)
    plotNonDiffEV("A Binomal Prob (min time)", 121, ABinomial, 0, 1, timeRanges, fig, ax1)
    plotNonDiffEV("A Uniform Dist Range", 122, AUniformRange, minTime, maxTime, timeRanges, fig, ax2)
    plt.show()

    fig, ((ax1, ax2)) = plt.subplots(ncols=2, nrows=1)
    plotNonDiffEV("B Binomal Prob (min time)", 121, BBinomial, 0, 1, timeRanges, fig, ax1)
    plotNonDiffEV("B Uniform Dist Range", 122, BUniformRange, minTime, maxTime, timeRanges, fig, ax2)
    plt.show()

   

iterateAndPlot()

