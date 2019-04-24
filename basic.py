"""
 My module : my utility functions -- basic.py
"""
import sys, csv

import SwarmData

FtoS = lambda fnum : str('{:g}'.format(float(fnum)))

# def FtoS(fnum):
#     """
#     Wrapper to convert float to str while removing unnesessarary 0
#     """
#     return str('{:g}'.format(fnum))

def MakeFilenameBody(N,v0,Ca,Cr,la,lr,alpha,init='random'):
    """
    Make file name string to load data. These parameters must be number objects. NOT STRING!!
    """
    if init=='random' or init=='random_Long':
        filename='Num'+str(N)+'v0'+str('{:g}'.format(v0))+'Ca'+str('{:g}'.format(Ca))+'Cr'+str('{:g}'.format(Cr))+'la'+str('{:g}'.format(la))+'lr'+str('{:g}'.format(lr))+'Alpha'+str('{:g}'.format(alpha))
    elif init=='mill':
        filename='Num'+str(N)+'v0'+str('{:g}'.format(v0))+'Ca'+str('{:g}'.format(Ca))+'Cr'+str('{:g}'.format(Cr))+'la'+str('{:g}'.format(la))+'lr'+str('{:g}'.format(lr))+'Alpha'+str('{:g}'.format(alpha))+'Inc0.1'
    else:
        print("Fail to set initial condition, ", init)
        sys.exit(0)
    return filename

def MakeFilenameBodySS(N,v0,Ca,Cr,la,lr,alpha,ss,init='random'):
    """
    Make file name string to load data. These parameters must be number objects. NOT STRING!!
    """
    if init=='random' or init=='random_Long':
        filename='Num'+str(N)+'v0'+str('{:g}'.format(v0))+'Ca'+str('{:g}'.format(Ca))+'Cr'+str('{:g}'.format(Cr))+'la'+str('{:g}'.format(la))+'lr'+str('{:g}'.format(lr))+'Alpha'+str('{:g}'.format(alpha))+'ss'+str('{:g}'.format(ss))
    elif init=='mill':
        filename='Num'+str(N)+'v0'+str('{:g}'.format(v0))+'Ca'+str('{:g}'.format(Ca))+'Cr'+str('{:g}'.format(Cr))+'la'+str('{:g}'.format(la))+'lr'+str('{:g}'.format(lr))+'Alpha'+str('{:g}'.format(alpha))+'ss'+str('{:g}'.format(ss))+'Inc0.1'
    else:
        print("Fail to set initial condition, ", init)
        sys.exit(0)
    return filename

def LoadBinaryState(f, N):
    """
    f : file object you want to load. f is NOT FILENAME
    N : Number of particles
    """
    PhaseDim=4*N
    sys_info=struct.unpack('<ii',f.read(8))
    print(sys_info)
    tmp = np.fromfile(f,dtype='<d')
    tmp = tmp.reshape(sys_info[1]+1,4*sys_info[0]);
    timelen = len(tmp[:,0])
    pos = np.zeros([timelen,2*N])
    vel = np.zeros([timelen,2*N])
    pos[:,0:2*N:2] = tmp[:,0:PhaseDim:4]
    pos[:,1:2*N:2] = tmp[:,1:PhaseDim:4]
    vel[:,0:2*N:2] = tmp[:,2:PhaseDim:4]
    vel[:,1:2*N:2] = tmp[:,3:PhaseDim:4]
    return timelen, pos, vel

def SaveCSV(csvfile, m_data):
    with open(csvfile,'w') as fs:
        writecsv = csv.writer(fs, delimiter='\t', lineterminator='\n')
        writecsv.writerows(m_data)

