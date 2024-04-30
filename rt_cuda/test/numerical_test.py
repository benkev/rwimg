import time, sys
import raytrace_cuda as rt


t = rt.implane([1,10],
                   [-2.0,-2.0,2.0,2.0],
                   (215,0,0),
                   25.0,
                   70e6,
                   mode='Tbr')
                   #mode='TbrIQUV')

before = time.time()

#t.trace(3000)
#t.trace(10, tracePoints=[[20,20],[21,22],[30,20],[25,25]])
#t.trace(10, tracePoints=[[2,2],[2,3],[3,2],[4,4]])
#t.trace(3000, tracePoints=[[45,45],[46,46],[50,50],[51,50]])

#t.trace(3000, tracePoints = [(2,3),['a',47]])
#t.trace(3000, tracePoints = [2,3,4])
t.trace(3000, tpts='All')

print
print 'elapsed time: %8.3f s' % (time.time()-before)

f = open('out.out', 'w')
for i in range(0,len(t.tbr)):
    f.write(str(t.tbr[i]))
