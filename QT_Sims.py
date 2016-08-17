# Simulation for Quantum Teleportation
# Last update: 20160524
import numpy as np
import random
import logging  
import logging.handlers  



# QuGates for single qubit
I = np.matrix([[1.,0.],[0.,1.]])
X = np.matrix([[0.,1.],[1.,0.]])
Z = np.matrix([[1.,0.],[0.,-1.]])
nn=1./np.sqrt(2)
H = np.matrix([[nn,nn],[nn,-nn]])

# QuGates for double qubits
CNOT = np.matrix([[1.,0.,0.,0.],\
				  [0.,1.,0.,0.],\
				  [0.,0.,0.,1.],\
				  [0.,0.,1.,0.]])


# Two qubits for constructing EPR pair
phi1 = np.matrix([1,0]).transpose()
phi2 = np.matrix([1,0]).transpose()

# Construct a Logger 
LOG_FILE = 'mylog0.log'  
  
fh = logging.FileHandler(LOG_FILE)
ch = logging.StreamHandler()
fmt = '%(message)s' 
  
formatter = logging.Formatter(fmt) 
fh.setFormatter(formatter)
ch.setFormatter(formatter) 
  
logger = logging.getLogger('MCProcess')  
logger.addHandler(fh)
logger.addHandler(ch)          
logger.setLevel(logging.INFO)  

def DirectProd(mat1,mat2):
	'''
	Obtains the tensor product of mat1 and mat2.
	'''
	# Obtaining matrices' shape:
	(m,n) = mat1.shape
	(k,l) = mat2.shape
	numRow = m*k
	numCol = n*l
	prod = np.zeros([numRow,numCol], dtype=complex)
	for i in range(m):
		for j in range(n):
				tmp = mat1[i,j]*mat2
				for ii in range(k):
					for jj in range(l):
						prod[i*k+ii, j*l+jj] = tmp[ii,jj]

	return np.matrix(prod)

def Telport(q):
	'''
	Formal QT process.
	'''
	# Constructing EPR pair:
	EPR = CNOT*(DirectProd(H,I)*DirectProd(phi1, phi2))
	logger.info("Constructed EPR pair:\n"+str(EPR))
	# Add q to the system, yielding the tri-Qubit:
	triQuBit = DirectProd(q, EPR)
	logger.info("The tri-qubit:\n"+str(triQuBit))

	# Seperate the entangled pair 0 & 1 for Alice's measurement
	node1 = DirectProd(CNOT, I)*triQuBit
	preMeas = DirectProd(DirectProd(H, I),I) * node1
	logger.info("State before Alice's measurement:\n"+str(preMeas))

	# Simulation for collapsing
	[xx] = random.sample([0,1,2,3], 1)
	if xx ==0: [m1,m2] = [0,0]
	elif xx == 1: [m1,m2] = [0,1]
	elif xx == 2: [m1,m2] = [1,0]
	else: [m1,m2] = [1,1]
	# Measurement results:
	[x,y] = [preMeas[2*xx], preMeas[2*xx+1][0]]
	x = 2*x[0,0]
	y = 2*y[0,0]

	# Bob's operation to resotre q
	result =  (Z**m1)*(X**m2)*np.matrix([x,y]).T
	return result

def qubitGenerator():
	'''
	Randomly generates a normalized qubit.
	'''
	r1 = np.random.random()
	r2 = np.sqrt(1-r1**2)
	theta1 = np.random.uniform(0,2*np.pi)
	theta2 = np.random.uniform(0,2*np.pi)
	a = r1*(np.cos(theta1)+np.sin(theta1)*1j)
	b = r2*(np.cos(theta2)+np.sin(theta2)*1j)
	q = np.matrix([a,b]).T
	return q

if __name__=="__main__":
	q = qubitGenerator()
	logger.info("Initial qubit 1 and 2:")
	logger.info(str(phi1.T)+','+str(phi2.T))
	logger.info('Qubit 0 to teleportation\n'+str(q))
	logger.info('\n')
	qq = Telport(q)
	logger.info("Restored qubit 0:\n"+str(qq))


