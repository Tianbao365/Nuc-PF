import time
import os
import networkx
import numpy
import scipy
from scipy import stats
from scipy import misc


class DeduplicateWorker(BioTasker):

    def __init__(self, pargs, plog ):
        self.pargs   = pargs
        self.log = plog

        if(self.pargs.regex==None):
            self.REGEX = "[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*"
        else:
            self.REGEX = self.pargs.regex
        
        self.optDistance = self.pargs.dist       
        self.qvalue = self.pargs.qvalue
        
        output = self.pargs.output + '/Dedup'
        if(not os.path.exists(self.pargs.output + '/Dedup')):
            tmp = os.mkdir(output)   
    
    def configureTask(self, chromNode):
        taskDatar = TaskDatar()    
        chromDatar = chromNode.data
        fragMapDict = chromDatar.map
        lociDiGraph = networkx.DiGraph()
        for fid in fragMapDict.keys():
            R1 = fragMapDict[fid]['R1'] #to +/-
            R2 = abs(fragMapDict[fid]['R2'])
            if(not lociDiGraph.has_node(R2)):
                lociDiGraph.add_node(R2, toR1s=collections.defaultdict(list))
            lociDiGraph.node[R2]['toR1s'][R1].append(fid)
        nodes = sorted(lociDiGraph.nodes())
        for i in range(0, len(nodes)-1):
            lociDiGraph.add_edge(nodes[i], nodes[i+1])  

        taskDatar._soft_sortedPoints = nodes            
        taskDatar._soft_lociDiGraph = lociDiGraph
        
        delattr(chromDatar, 'map') 
        chromSize = chromDatar.size
        taskDatar._soft_EXTENTION = numpy.ceil(1.0*chromSize/len(nodes))
        return taskDatar         
    
    def getLinkedR1Count(self, R2, lociDiGraph):
        toR1sDict = lociDiGraph.node[R2]['toR1s']
        localPoints = []
        for R1 in map(abs, toR1sDict.keys()):
            fragments = toR1sDict.get(R1,[])
            fragments.extend(toR1sDict.get(-R1,[]))
            size = len(fragments)
            localPoints.append(size)
        return numpy.median(localPoints)
        
    def findDuplicates(self, taskDatar, chrid):
        amplifyFragments, opticalFragments = [], []
        lociDiGraph = taskDatar._soft_lociDiGraph
        sortedPoints = taskDatar._soft_sortedPoints
        for R2 in sortedPoints:
            toR1sDict = lociDiGraph.node[R2]['toR1s']
            localPoints = None # prepare firstly, filled when necessary
            for R1 in map(abs, toR1sDict.keys()):
                fragments = toR1sDict.get(R1,[])
                fragments.extend(toR1sDict.get(-R1,[]))
                size = len(fragments)
                if(size >= 5 ):           
                    if(localPoints is None):
                        localPoints = self.buildLocalPoints(R2, taskDatar)      
                    if(len(localPoints) == 0):
                        lociLambda = 1.0 
                    else:                        
                        lociLambda = numpy.mean(localPoints) 
                    k = numpy.ceil(size - lociLambda)
                    if(k <= 0):
                        continue
                    else:
                        prb = 1.0 - stats.skellam.cdf(k, lociLambda, lociLambda) 
                        if(prb <= self.qvalue):
                            ampids = random.sample(fragments, int(k)-1 ) 
                            for item in ampids:
                                amplifyFragments.append((R2, R1, item, k, size, lociLambda, prb))
                for i in range(0, size):
                    fieldi= self.getLocation(fragments[i])   
                    if(fieldi is None or len(fieldi)!=3):
                        continue
                    for j in range(i+1,size):
                        fieldj= self.getLocation(fragments[j])  
                        if(fieldj is None or len(fieldj)!=3):
                            continue
                        if(fieldi[0]!=fieldj[0]): 
                            continue
                        if(abs(fieldi[1]-fieldj[1]) > self.optDistance): 
                            continue
                        if(abs(fieldi[2]-fieldj[2]) <= self.optDistance):
                            opticalFragments.append( (R2, R1, fragments[j], fragments[i]) )               
        taskDatar._soft_amplifyFragments = amplifyFragments
        taskDatar._soft_opticalFragments = opticalFragments
        return taskDatar
        
    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        #self.log.info("%s: Starting Dedup-task on %s..." % (time.ctime(), chrid) )
        self.log.info("Starting Dedup-task on %s..." % (chrid) )
        taskDatar = self.configureTask(chromNode)

        taskDatar = self.findDuplicates(taskDatar, chrid)
        outfile = self.pargs.output + '/Dedup/dedup_on_' + chrid + '.txt'
        handle.write(headline)
        for item in taskDatar._soft_opticalFragments:
            R2, R1, optfid, matchfid = item
            start, end = min(R2, R1), max(R2, R1)
            linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t' \
                    + optfid + '\t' + matchfid + '\n'
            handle.write(linestr)
        for item in taskDatar._soft_amplifyFragments:
            R2, R1, frag, k, size, lociLambda, prb = item
            start, end = min(R2, R1), max(R2, R1)
            linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t' \
                    + frag + '\t' + str(size) +'\t' + str(lociLambda) + '\t' \
                    + str(k) + '\t' + str(prb) + '\n'
            handle.write(linestr)
        handle.close() 
     
        self.log.info("Dedup-task is finished on %s..." % ( chrid) )  
        self.log.info("Please refer to %s!" % (outfile) )		
     
class PeakScanner(BioTasker):
    
    def __init__(self, pargs, plog):
        self.pargs   = pargs
        self.log = plog 

        # get those parameters
        self.R = self.pargs.rscan
        self.s = self.pargs.minstep
        self.S = self.pargs.maxstep
        self.pvalue = self.pargs.pvalue        
        
        # preparing the output fold
        if(not os.path.exists(self.pargs.output + '/Peak')):
            tmp = os.mkdir(self.pargs.output + '/Peak') 
        self.outDir = self.pargs.output + '/Peak'            
        
    def configureTask(self, chromNode):
        chromDatar = chromNode.data 
        taskDatar = self.configurePairMode(chromDatar)
        return taskDatar
        
        
    def configurePairMode(self, chromDatar): 
        taskDatar = TaskDatar()    
        lociDiGraph = chromDatar.lociDiGraph 
        sortedPoints = chromDatar.sortedPoints
        pRefDiGraph = networkx.DiGraph()
        pHead = 'PHEAD'
        pRefDiGraph.add_node(pHead)
        nRefDiGraph = networkx.DiGraph()
        nHead = 'NHEAD'
        nRefDiGraph.add_node(nHead)
        for R2 in sortedPoints:
            toR1sDict = lociDiGraph.node[R2]['toR1s']
            pcountDict = collections.defaultdict(int)
            ncountDict = collections.defaultdict(int)
            for R1 in toR1sDict.keys():
                if(R1 > 0): 
                    ncountDict[R1] += toR1sDict[R1]
                else: 
                    pcountDict[R1] += toR1sDict[R1]
            pdegree = sum(pcountDict.values())        
            if( pdegree > 0 ):
                pRefDiGraph.add_node(R2, toR1s=pcountDict, depth=pdegree)
                pRefDiGraph.add_edge(pHead, R2)
                pHead = R2
            ndegree = sum(ncountDict.values())    
            if( ndegree > 0 ):
                nRefDiGraph.add_node(R2, toR1s=ncountDict, depth=ndegree)
                nRefDiGraph.add_edge(nHead, R2)
                nHead = R2
        
        pHead = pRefDiGraph.successors('PHEAD')[0]
        pRefDiGraph.remove_node('PHEAD')
        nHead = nRefDiGraph.successors('NHEAD')[0]
        nRefDiGraph.remove_node('NHEAD')
        taskDatar.pRefDiGraph = pRefDiGraph
        taskDatar.pHead = pHead
        taskDatar.nRefDiGraph = nRefDiGraph
        taskDatar.nHead = nHead
        delattr(chromDatar, 'lociDiGraph') 
        delattr(chromDatar, 'sortedPoints')
        return taskDatar        
        
    def peakScanningProc(self, taskDatar, chrid):
        pRefDiGraph = taskDatar.pRefDiGraph 
        nRefDiGraph = taskDatar.nRefDiGraph        
        measure = lambda e: abs(e[1] - e[0])
        pComponents = self.breakingIntoComponents(pRefDiGraph, taskDatar.pHead, measure, k=3.0)
        nComponents = self.breakingIntoComponents(nRefDiGraph, taskDatar.nHead, measure, k=3.0)
        degreepf = lambda group: sum( [pRefDiGraph.node[n]['depth'] for n in group] )
        pComponents = [ item for item in pComponents if degreepf(item) >= self.R ]
        degreenf = lambda group: sum( [ nRefDiGraph.node[n]['depth'] for n in group] )
        nComponents = [ item for item in nComponents if degreenf(item) >= self.R ]
        pPeakList, pCompSerials = [], []
        compid = 0
        for component in pComponents:
            pPeaks = self.processComponent(component, pRefDiGraph)   
            pPeakList.extend(pPeaks)
            pCompSerials.extend([compid for i in range(0, len(pPeaks))])
            compid += 1
        nPeakList, nCompSerials = [], []
        compid = 0
        for component in nComponents:
            nPeaks = self.processComponent(component, nRefDiGraph)
            nPeakList.extend(nPeaks)            
            nCompSerials.extend([compid for i in range(0, len(nPeaks))])
            compid += 1 
        pmergedPeakDict = self.peakMerging(pPeakList, pCompSerials, pRefDiGraph)        
        nmergedPeakDict = self.peakMerging(nPeakList, nCompSerials, nRefDiGraph)
        return pmergedPeakDict, nmergedPeakDict
        
    def breakingIntoComponents(self, refDiGraph, head, measure, k=3.0):
        edgeValList = []
        current = head
        while(True):
            succ = refDiGraph.successors(current)
            if(succ):
                value = measure((current, succ[0]))
                edgeValList.append(value)
                current = succ[0]
            else:
                break  
        obslist, enumber = edgeValList, len(edgeValList)                
        edgeFlagList = [1 for i in range(0, enumber)]
        while(True):   
            valarray = numpy.array(obslist)
            mu, sigma = numpy.mean(valarray), numpy.std(valarray)
            edgeFlagList = [ 0 if val-mu>k*sigma else 1 for val in edgeValList ] 
            obslist = [ val for i,val in enumerate(edgeValList) if edgeFlagList[i] ]
            onumber = len(obslist)
            if(onumber < enumber):
                enumber = onumber
            else:
                break
      
        components = [[head]]  
        current, cursor = head, 0
        while(True):
            succ = refDiGraph.successors(current)
            if(succ):
                if(edgeFlagList[cursor]):
                    components[-1].append(succ[0])
                else:
                    components.append([succ[0]])                
                current = succ[0]
                cursor += 1
            else:
                break              
        return components      
    
    def peakMerging(self, peakList, compSerials, refDiGraph):
        peakMgGraph = networkx.Graph()
        for i in range(0, len(peakList)-1):
            prevPeak, prevCompId = peakList[i], compSerials[i]
            nextPeak, nextCompId = peakList[i+1], compSerials[i+1]
            prevStart, prevEnd = prevPeak.start, prevPeak.end
            nextStart, nextEnd = nextPeak.start, nextPeak.end
            dist = abs(nextStart - prevEnd)
            if(prevCompId==nextCompId): #peak from the same component
                peakMgGraph.add_edge(i, i+1, dist=dist)
            else:
                peakMgGraph.add_nodes_from([i, i+1])            
        measure = lambda e: peakMgGraph.get_edge_data(*e)['dist']
        peakMgGraph = self.outlierEdgeCutting(peakMgGraph, measure, k=3.0)
        
        components = networkx.connected_components(peakMgGraph)
        MergedPeak = collections.namedtuple('MergedPeak', 'chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue')
        mergedPeakDict = collections.defaultdict(None)
        for comp in components:
            comp = sorted(comp)
            chromStart = peakList[comp[0]].start
            chromEnd = peakList[comp[-1]].end #included
            posList, depthList = self.nodeDepthFlowing(chromStart, chromEnd, refDiGraph)
            valList = [1.0*pos*val for pos,val in zip(posList, depthList)] 
            uPoint = int(sum(valList)/sum(depthList)) # the "average" point based on density distribution
            blockCount = len(comp)
            blockSizes = [peakList[index].end - peakList[index].start + 1 for index in comp]
            blockStarts = [ peakList[index].start-chromStart for index in comp ]
            name = str(chromStart) #+ '-' + str(chromEnd) 
            #name = str(uPoint)
            cL, cN = self.getRScanLocalContext(chromStart, chromEnd, refDiGraph)  
            w, r = chromEnd - chromStart + 1, sum(valList) + 1
            score, pValue = self.calculateScanScore(r, w, cL, cN)         
            mgPeak = MergedPeak(chromStart, chromEnd, name, score, blockCount, blockSizes, blockStarts, uPoint, pValue)
            mergedPeakDict[chromStart] = mgPeak      
        return mergedPeakDict    
        
    
    def nodeDepthFlowing(self, start, end, refDiGraph):
        posList, valList = [], []
        pos = start
        while(pos <= end):
            depth = refDiGraph.node[pos]['depth']
            posList.append(pos)
            valList.append(depth)
            succs = refDiGraph.successors(pos)
            if(succs):
                pos = succs[0]
            else:
                break  # to the end          
        return posList, valList            
    
    def calculateScanScore(self, r, w, L, N):
        logx    = numpy.log(w) - numpy.log(L) + (1.0 + 1.0/r) * numpy.log(N)
        score  = r*logx - util.fastLogFactorial(r)
        lda = numpy.exp(score)
        prob = 1.0 - numpy.exp(-lda)
        return score, prob
    
    def peakPairing(self, pmPeakDict, nmPeakDict ):
        pairDiGraph = networkx.DiGraph()
        pPoints, nPoints = sorted(pmPeakDict.keys()), sorted(nmPeakDict.keys())
        pN, nN = len(pPoints), len(nPoints)
        pc, nc   = 0, 0 
        for pPos in pPoints:
            pUPoint = pmPeakDict[pPos].uPoint
            nPos = nPoints[nc]
            nUPoint  = nmPeakDict[nPos].uPoint
            while(nUPoint <= pUPoint and nc < nN-1):
                nc += 1
                nPos = nPoints[nc]
                nUPoint = nmPeakDict[nPos].uPoint
            GAPSIZE = abs(pUPoint - nUPoint)
            pairDiGraph.add_edge( 'P-' + str(pPos), 'N-' + str(nPos), dist=GAPSIZE )   
        for nPos in nPoints:
            nUPoint = nmPeakDict[nPos].uPoint
            pPos = pPoints[pc]
            pUPoint  = pmPeakDict[pPos].uPoint
            if(pUPoint >= nUPoint):
                continue
            while(pUPoint < nUPoint and pc < pN-1):
                pc += 1
                pPos = pPoints[pc]
                pUPoint = pmPeakDict[pPos].uPoint
            pc -= 1
            pPos = pPoints[pc]
            pUPoint = pmPeakDict[pPos].uPoint           
            GAPSIZE = abs(nUPoint - pUPoint)
            pairDiGraph.add_edge('N-' + str(nPos), 'P-' + str(pPos), dist=GAPSIZE )
        ef = lambda e: pairDiGraph.get_edge_data(*e)['dist']        
        pairDiGraph = self.outlierEdgeCutting(pairDiGraph, ef, k=3.5 )         
        return pairDiGraph
        
    def goodPeakPairs(self, peakPairDiGraph):
        PNPairList = set()
        for edge in peakPairDiGraph.edges():
            source, target = edge[0], edge[1]
            if(   source in peakPairDiGraph.successors(target) \
              and target in peakPairDiGraph.successors(source) ):
                pair = (source, target) if source.startswith('P-') else (target, source)
                PNPairList.add(pair)
        return list(PNPairList)     
        
    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        self.log.info("Starting PeakCalling-task on %s..." % (chrid))
        taskDatar = self.configureTask(chromNode)
        pmergedPeakDict, nmergedPeakDict = self.peakScanningProc(taskDatar, chrid)  
        pairDiGraph = self.peakPairing(pmergedPeakDict, nmergedPeakDict)
        pairs = self.goodPeakPairs(pairDiGraph)
        
        unpairNodes = set(pairDiGraph.nodes())
        pNodes = set([pair[0] for pair in pairs])
        nNodes = set([pair[1] for pair in pairs])
        
        unpairNodes = unpairNodes.difference(pNodes)
        unpairNodes = unpairNodes.difference(nNodes)
        outfile = self.outDir + '/peakcall_on_' + chrid + '.paired.bed'
        handle = open(outfile, 'w')
        for pair in pairs:
            mpPeakId = int(pair[0].split('-')[1])
            mpPeak = pmergedPeakDict[mpPeakId]
            linestr = self.buildOneLineStr(chrid, '+', mpPeak, 'pp')
            handle.write(linestr)
            
            mnPeakId = int(pair[1].split('-')[1])
            mnPeak = nmergedPeakDict[mnPeakId]
            linestr = self.buildOneLineStr(chrid, '-', mnPeak, 'np')
            handle.write(linestr)
        handle.close()     
        chromDatar = chromNode.data
        chromDatar.pRefDiGraph = taskDatar.pRefDiGraph
        chromDatar.pHead = taskDatar.pHead
        chromDatar.pmergedPeakDict = pmergedPeakDict
        
        chromDatar.nRefDiGraph = taskDatar.nRefDiGraph
        chromDatar.nHead = taskDatar.nHead
        chromDatar.nmergedPeakDict = nmergedPeakDict
        chromDatar.peakPairs = pairs
        self.log.info("PeakCalling-task is finished on %s..." % (chrid) )

    def buildOneLineStr(self, chrid, strand, peak, flag):
        linestr = chrid + '\t' + str(peak.chromStart) + '\t' + str(peak.chromEnd+1) + '\t'
        if(strand == '+'):
            name = 'PP-' + peak.name
        else:
            name = 'PN-' + peak.name
        score = 0    
        linestr+= name + '\t' + str(score) + '\t' + strand + '\t'
        thickStart, thickEnd = str(peak.chromStart), str(peak.chromEnd+1)
        linestr+= thickStart + '\t' + thickEnd + '\t'
        if(flag == 'pp'):
            itemRgb = "255,0,0"
        elif(flag == 'np'):
            itemRgb = "0,0,255"
        else:
            itemRgb = "0,255,0"
        linestr+= itemRgb + '\t'
        linestr+= str(peak.blockCount) + '\t'
        blockSizes = ','.join(map(str, peak.blockSizes))
        linestr+= blockSizes + '\t'
        blockStarts = ','.join(map(str, peak.blockStarts))
        linestr+= blockStarts + '\n'
        return linestr        


class BorderScanner(BioTasker):
    '''A class implemented for Border-positionate'''
    
    def __init__(self, pargs, plog):
        self.pargs   = pargs
        self.log = plog 
        self.c = self.pargs.chernoff
        self.k = self.pargs.outlier        
        
        # preparing the output fold
        if(not os.path.exists(self.pargs.output + '/Border')):
            tmp = os.mkdir(self.pargs.output + '/Border') 
        self.outDir = self.pargs.output + '/Border'
    
    def configurePairMode(self, chromDatar): 
        taskDatar = TaskDatar()    
        
        pRefDiGraph = chromDatar.pRefDiGraph  # R2
        pHead = chromDatar.pHead
        pmergedPeakDict = chromDatar.pmergedPeakDict
        
        nRefDiGraph = chromDatar.nRefDiGraph
        nHead = chromDatar.nHead
        nmergedPeakDict = chromDatar.nmergedPeakDict
        
        peakPairs = chromDatar.peakPairs
        taskDatar.pRefDiGraph = pRefDiGraph
        taskDatar.pHead = pHead
        taskDatar.pmergedPeakDict = pmergedPeakDict
        taskDatar.nRefDiGraph = nRefDiGraph
        taskDatar.nHead = nHead
        taskDatar.nmergedPeakDict = nmergedPeakDict
        taskDatar.peakPairs = peakPairs
        
        return taskDatar

    def retriveR1SignalsForward(self, startR2, endR2, refR2DiGraph):
        posSignalDict = defaultdict(int)
        current = startR2
        while(current and current <= endR2):
            toR1sDict = refR2DiGraph.node[current]['toR1s']
            for R1 in toR1sDict.keys():
                if(R1 < 0 ): 
                    posSignalDict[abs(R1)] += toR1sDict[R1]
            succs = refR2DiGraph.successors(current)
            if(succs):
                current = succs[0]
            else:
                break
        return posSignalDict  

    def retriveR2BackgroundForward(self, startR2, endR2, refR2DiGraph):
        posSignalDict = defaultdict(int)
        current = startR2
        while(current and current <= endR2):
            degree = refR2DiGraph.node[current]['depth']
            posSignalDict[current] = degree
            succs = refR2DiGraph.successors(current)
            if(succs):
                current = succs[0]
            else:
                break
        return posSignalDict        
    
    def retriveR1SignalsBackward(self, startR2, endR2, refR2DiGraph):
        posSignalDict = defaultdict(int)
        current = endR2
        while(current and current >= startR2):
            toR1sDict = refR2DiGraph.node[current]['toR1s']
            for R1 in toR1sDict.keys():
                if(R1 > 0 ): 
                    posSignalDict[abs(R1)] += toR1sDict[R1]
            preds = refR2DiGraph.predecessors(current)
            if(preds):
                current = preds[0]
            else:
                break
        return posSignalDict
        
    def retriveR2BackgroundBackward(self, startR2, endR2, refR2DiGraph):
        posSignalDict = defaultdict(int)
        current = endR2
        while(current and current >= startR2):
            degree = refR2DiGraph.node[current]['depth']
            posSignalDict[current] = degree
            preds = refR2DiGraph.predecessors(current)
            if(preds):
                current = preds[0]
            else:
                break
        return posSignalDict
        
    def runBorderScanning(self, taskDatar):
        pPeakBorderDict, nPeakBorderDict = {}, {}
        pPeakStrengthDict, nPeakStrengthDict = {}, {}
        for pPeak, nPeak in taskDatar.peakPairs:
            pStart = int(pPeak.split('-')[1])
            nStart = int(nPeak.split('-')[1])
            pEnd = taskDatar.pmergedPeakDict[pStart].chromEnd
            nEnd = taskDatar.nmergedPeakDict[nStart].chromEnd
            
            nR1PosSignalDict = self.retriveR1SignalsForward(pStart, nEnd, taskDatar.pRefDiGraph)
            nR1Strength = sum(nR1PosSignalDict.values())
            nPeakStrengthDict[nPeak] = nR1Strength
            nR2PosBackgroundDict = self.retriveR2BackgroundBackward(pStart, nEnd, taskDatar.nRefDiGraph)
            nR1BorderSites = self.estimateBorder(nR1PosSignalDict, nR2PosBackgroundDict)
            nPeakBorderDict[nPeak] = self.borderMerged(nR1BorderSites)
            pR1PosSignalDict = self.retriveR1SignalsBackward(pStart, nEnd, taskDatar.nRefDiGraph)
            pR1Strength = sum(pR1PosSignalDict.values())
            pPeakStrengthDict[pPeak] = pR1Strength
            pR2PosBackgroundDict = self.retriveR2BackgroundForward(pStart, nEnd, taskDatar.pRefDiGraph)
            pR1BorderSites = self.estimateBorder(pR1PosSignalDict, pR2PosBackgroundDict)   
            pPeakBorderDict[pPeak] = self.borderMerged(pR1BorderSites)

        taskDatar.pPeakBorderDict = pPeakBorderDict
        taskDatar.nPeakBorderDict = nPeakBorderDict
        taskDatar.pPeakStrengthDict = pPeakStrengthDict
        taskDatar.nPeakStrengthDict = nPeakStrengthDict        
        return taskDatar   
    
    
    def estimateBorder(self, posSignalDict, backgroundDict):
        borderSites = []
        
        oPosi = [ pos for pos in posSignalDict.keys() ]
        oPosi = sorted(oPosi)
        oDepth = [ posSignalDict[pos] for pos in oPosi ]
        ou = numpy.mean(oDepth)
        
        df = lambda c: 1.0*c/ou - 1
        pf = lambda d: numpy.power(numpy.exp(d)/numpy.power(1+d, 1+d), ou)
 
        bPosi  = [ pos for pos in backgroundDict.keys() ]
        bPosi = sorted(bPosi)
        bDepth = [ backgroundDict[pos] for pos in bPosi ]
        bu = numpy.mean(bDepth)
        bRugSignals = [ random.uniform(pos-0.5, pos+0.5) for pos, val in backgroundDict.items() for j in range(0, val) ]
        for i in range(0, len(oPosi)):
            pos, depth = oPosi[i], oDepth[i]
            delta = df(depth)
            if(delta > 0):
                p_chernoff = pf(delta)
                if(p_chernoff < self.c): 
                    p_rug = brugsig_kde.integrate_box_1d(pos-0.5, pos+0.5)  
                    prob = 1.0 - stats.binom.cdf(depth, N, p_rug)
                    if(prob < 0.00001):
                        borderSites.append((pos, pos+1, depth, p_chernoff, prob))                        
        return borderSites
    
    def borderMerged(self, borderSites): 
        mergedBorders = []
        for item in borderSites:
            curr_start, curr_end, curr_depth = item[0], item[1], item[2]
            if(len(mergedBorders) == 0):
                mergedBorders.append( item )
                continue
            else:
                    prevBorder = mergedBorders[-1]
                    prev_start, prev_end, prev_depth = prevBorder[0], prevBorder[1], prevBorder[2]
                    if( abs(curr_start - prev_end) <= 3 ):
                        if(curr_depth > prev_depth):
                            mergedBorders[-1] = item 
                    else:
                        mergedBorders.append(item) 
        maxDepth = -numpy.inf
        for border in mergedBorders:
            if(border[2] > maxDepth):
                maxDepth = border[2]
        marked = [] 
        for border in mergedBorders:
            if(1.0*maxDepth/border[2] >= 4.0 or border[2] < 5): 
                marked.append(border)
        mergedBorders = [item for item in mergedBorders if item not in marked]                 
        return mergedBorders  
        
    def borderPairing(self, taskDatar):
        pPeakBorderDict = taskDatar.pPeakBorderDict
        nPeakBorderDict = taskDatar.nPeakBorderDict
        borderPairDiGraph = networkx.DiGraph()
        for pPeak, nPeak in taskDatar.peakPairs:
            pBorderSites = sorted( [ item[0] for item in pPeakBorderDict[pPeak] ] )
            nBorderSites = sorted( [ item[0] for item in nPeakBorderDict[nPeak] ] )
            pN, nN = len(pBorderSites), len(nBorderSites)
            if(max(pN, nN) == 0):
                continue
            pc, nc = 0, 0 
            for pBorder in pBorderSites:
                borderPairDiGraph.add_node('P-' + str(pBorder))
                index = pBorderSites.index(pBorder)
                if(index >= 1):
                    currNode = 'P-' + str(pBorder)
                    prevNode = 'P-' + str(pBorderSites[index-1])
                    paradist = abs(pBorder - pBorderSites[index-1])
                    borderPairDiGraph.add_edge(currNode, prevNode, dist=paradist, type='P-P')
                minDist = numpy.inf
                partner = None
                while(True):
                    if(nN == 0):
                        break
                    nBorder = nBorderSites[nc]
                    dist = nBorder - pBorder  
                    if(dist > 0 and dist < minDist):
                        minDist = dist
                        partner = nBorder
                    else:
                        if(partner):
                            break    
                    if(nc < nN - 1):                  
                        nc += 1
                    else:
                        break                    
                if(nc != 0): 
                    nc -= 1            
                if(partner):   
                    source = 'P-'+str(pBorder)
                    target = 'N-'+str(partner)               
                    borderPairDiGraph.add_edge(source, target, dist=minDist, type='P-N')
            #for the borders on the negative strand        
            for nBorder in nBorderSites:
                borderPairDiGraph.add_node('N-'+str(nBorder))
                index = nBorderSites.index(nBorder)
                if(index >= 1):
                    currNode = 'N-' + str(nBorder)
                    prevNode = 'N-' + str(nBorderSites[index-1])
                    paradist = abs(nBorder - nBorderSites[index-1])
                    borderPairDiGraph.add_edge(currNode, prevNode, dist=paradist, type='N-N')
                minDist = numpy.inf
                partner = None
                while(True):
                    if(pN == 0):
                        break
                    pBorder = pBorderSites[pc]
                    dist = nBorder - pBorder
                    if(dist > 0 and dist < minDist):
                        minDist = dist
                        partner = pBorder
                    else:
                        break  # here is different from the pBorderSites
                    if(pc < pN - 1):
                        pc += 1
                    else:
                        break                    
                if(pc != 0): 
                    pc -= 1            
                if(partner):
                    source = 'N-'+str(nBorder)
                    target = 'P-'+str(partner)           
                    borderPairDiGraph.add_edge(source, target, dist=minDist, type='N-P')
        return borderPairDiGraph                      
     
    def outlierEdgeCutting(self, paramGraph, measure, k=2.0): 
        paramEdges = paramGraph.edges()
        edgeValList = [ measure(edge) for edge in paramEdges ]             
        edgeFlagList = [1 for i in range(0, len(edgeValList))]
        obslist, enumber = edgeValList, len(edgeValList)
        while(True):   
            valarray = numpy.array(obslist)
            mu, sigma = numpy.mean(valarray), numpy.std(valarray)
            edgeFlagList = [ 0 if val-mu > k*sigma else 1 for val in edgeValList ] 
            obslist = [ val for i,val in enumerate(edgeValList) if edgeFlagList[i] ]
            onumber = len(obslist)
            if(onumber < enumber):
                enumber = onumber
            else:
                break
        eraseEdges = [paramEdges[i] for i,flag in enumerate(edgeFlagList) if not flag]
        paramGraph.remove_edges_from(eraseEdges)        
        return paramGraph
    
    
    def postBorderProc(self, taskDatar):
        borderPairDiGraph = self.borderPairing(taskDatar)
        ef = lambda e: borderPairDiGraph.get_edge_data(*e)['dist']        
        borderPairDiGraph = self.outlierEdgeCutting(borderPairDiGraph, ef, k=self.k)
        taskDatar.borderPairDiGraph = borderPairDiGraph
        return taskDatar
        
    def runTask(self, chromNode):
        chrid = chromNode.identifier
        chrid = chrid[0:3] + chrid[3:].upper()
        self.log.info("Starting Border Scanning-task on %s..." % (chrid) )
        taskDatar = self.configureTask(chromNode)
        taskDatar = self.runBorderScanning(taskDatar)
        
        taskDatar = self.postBorderProc(taskDatar)
        borderInfoDict = {}
        for peak in taskDatar.pPeakBorderDict.keys():
            borderSites = taskDatar.pPeakBorderDict[peak]
            peakname = 'P' + peak
            for i in range(0, len(borderSites)):
                 start, end, depth, p_chernoff, prob = borderSites[i]
                 bordername = 'BP-' + str(start)
                 borderInfoDict[bordername] = (start, end, depth, p_chernoff, peakname)
        for peak in taskDatar.nPeakBorderDict.keys():
            borderSites = taskDatar.nPeakBorderDict[peak]
            peakname = 'P' + peak
            for i in range(0, len(borderSites)):
                start, end, depth, p_chernoff, prob = borderSites[i]
                bordername = 'BN-' + str(start)
                borderInfoDict[bordername] = (start, end, depth, p_chernoff, peakname)
        outfile = self.outDir + '/borders_on_' + chrid + '.bed'
        handle = open(outfile, 'w')
        borderPairDiGraph = taskDatar.borderPairDiGraph
        components = networkx.connected_components(borderPairDiGraph.to_undirected())
        for i in range(0, len(components)):
            comp = components[i]
            compid = chrid + '_' + 'Comp_' + str(i+1)
            pnodes = [ node for node in comp if node.__contains__('P-') ]
            nnodes = [ node for node in comp if node.__contains__('N-') ]
            p, n, b = self.getPNBackboneNum(comp, borderPairDiGraph)
            for node in comp:
                bordername = 'B' + node
                start, end, depth, p_chernoff, peakname = borderInfoDict[bordername]
                linestr = chrid + '\t' + str(start) + '\t' + str(end) + '\t'
                linestr+= bordername + '\t' + str(depth) + '\t'
                if(node.startswith('P')):
                    linestr += '+\t'
                else:
                    linestr += '-\t'                
                linestr += str(p_chernoff) + '\t' + peakname + '\t'
                linestr+=  compid + '\t' + str(p) + ',' + str(n) + ',' + str(b) + '\n'
                handle.write(linestr)
        handle.close()
        borderPairDiGraph = taskDatar.borderPairDiGraph
        out = open(self.outDir + '/borderGraphs_on_' + chrid + '.json', 'w')
        json_graph.dump(borderPairDiGraph, out)
        out.close()
        self.log.info("Border Scanning-task has finished on %s..." % ( chrid) ) 

    def getPNBackboneNum(self, comp, hostGraph):
        pnodes = [ node for node in comp if node.__contains__('P-') ]
        nnodes = [ node for node in comp if node.__contains__('N-') ]
        backbone = 0
        for pn in pnodes:
            for nn in nnodes:
                e1 = hostGraph.get_edge_data(pn, nn)
                e2 = hostGraph.get_edge_data(nn, pn)
                if(e1 and e2):
                    backbone += 1
        return len(pnodes), len(nnodes), backbone        


with open('ePENS.py') as f:
    code = compile(f.read(), 'ePENS.py','exec')
    exec(code, None, None)

baseTree = exoBlocker.dataModel.baseTree
treeNodes = baseTree.all_nodes()
branches = [node for node in treeNodes if baseTree.level(node.identifier)==1]
import bll.dedup as dedup
pargs, plog = app.pargs, app.log
dedupTask = dedup.DeduplicateWorker(pargs, plog)
chromNode = branches[0]
dedupTask.runTask(chromNode)

import bll.peak as peak
peakTask = peak.PeakScanner(pargs, plog) 
peakTask.runTask(chromNode)


