import networkx as nx
import numpy as np

class SI_Simulation():
    """
    A class to simulate the SI Model.
    ===================================================
    Input: A - Network adjacency matrix (numpy array) or Networkx graph object.
           lam - Contagion parameter.
           gam - Spontaneous infection parameter.
           i0 - Initial infected fraction.
           prop - Propensity calculation method [ALL,TAKE,FANCY,LOOP,SLICE].
    """

    def __init__(self, A, lam=0.1, gam=0.001, i0=0.1, prop='ALL'):

        # Network setup.
        if type(A)==np.ndarray:
            self.A = A
        elif type(A)==nx.classes.graph.Graph:
            self.A = nx.adj_matrix(A).A
        else:
            raise BaseException("Input an adjacency matrix or networkx object only.")

        # Model Parameters.
        self.N = A.shape[0]
        self.lam = lam/self.N
        self.gam = gam
        self.prop = prop

        # Time-keeping.
        self.t = 0
        self.times = [0]
        
        # Node numbers.
        self.I = [int(i0*self.N)]
        self.S = [self.N-self.I[0]]
        
        # Node states.
        self.X = np.array([1]*self.S[0] +[2]*self.I[0]).reshape((self.N,1))
        np.random.shuffle(self.X)
         
        # Initial propensity setup.        
        self.UpdatePropensityALL() 
        
        # Select which propensity scheme to use.
        if self.prop == 'ALL':
            self.UpdatePropensity = self.UpdatePropensityALL 
        elif self.prop == 'FANCY':
            self.UpdatePropensity = self.UpdatePropensityFANCY 
        elif self.prop == 'TAKE':
            self.UpdatePropensity = self.UpdatePropensityTAKE 
        elif self.prop == 'SLICE':
            self.UpdatePropensity = self.UpdatePropensitySLICE
        elif self.prop == 'LOOP':
            self.UpdatePropensity = self.UpdatePropensityLOOP
        else:
            raise BaseException("Please specify a propensity scheme [ALL,TAKE,FANCY,LOOP,SLICE].")
        return None
        
    def UpdatePropensityALL(self, n_nodes=None):
        self.IP = (self.gam + self.lam*self.A.dot(self.X==2))*(self.X==1)
        return None
    
    def UpdatePropensityTAKE(self,n_nodes):   
        self.IP[n_nodes] = (self.gam + self.lam*np.take(self.A,n_nodes, axis=0).dot(self.X==2))*((np.take(self.X,n_nodes, axis=0)==1))
        return None    
    
    def UpdatePropensityFANCY(self,n_nodes):
        self.IP[n_nodes] = (self.gam + self.lam*self.A[n_nodes].dot(self.X==2))*(self.X[n_nodes]==1)
        return None
        
    def UpdatePropensityLOOP(self,n_nodes):        
        for node in n_nodes:
            self.IP[node] = (self.gam + self.lam*self.A[node].dot(self.X==2))*(self.X[node]==1)
        return None
    
    def UpdatePropensitySLICE(self,n_nodes):        
        nmax, nmin = n_nodes.max()+1, n_nodes.min()
        self.IP[nmin:nmax] = (self.gam + self.lam*self.A[nmin:nmax,:].dot(self.X==2))*(self.X[nmin:nmax]==1)
        return None
        
    def RunIteration(self):
        
        # Termination.
        if self.S[-1] == 0:
            self.S = np.array(self.S, dtype=float)
            self.I = np.array(self.I, dtype=float)
            return False

        # 1. Generate random numbers r1,r2 uniformly distributed in (0,1)
        r1 = np.random.rand()
        r2 = np.random.rand()
        
        # 2. Calculate alpha.
        cumsum = self.IP.cumsum()
        self.alpha = cumsum[-1]

        # 3. Compute the time until the next reaction takes place.
        tau = (1.0/self.alpha)*np.log(float(1.0/r1))
        self.t += tau
        self.times.append(self.t)

        # 4. Compute which reaction takes place.
        index = np.searchsorted(cumsum,r2*self.alpha)

        # 5. Update node states. 
        self.X[index%self.N] = 2
        self.S.append(self.S[-1] - 1)
        self.I.append(self.I[-1] + 1)
            
        # 6. Update propensities
        n1 = np.nonzero(self.A[index%self.N])[0]
        n1 = np.append(n1,index%self.N)
        
        self.UpdatePropensity(n1)
        return True

    def RunToConvergence(self):
        running = True
        while running:
            running = self.RunIteration()
        return None
    
    def IntegrateSolution(self):
        from scipy.integrate import odeint
        gam = self.gam
        lam = self.lam*self.N
        k = self.A.sum(axis=0).mean()

        def deriv(x,t):
            xdot = [-x[0]*(gam + lam*x[1]),
                    +x[0]*(gam + lam*x[1])]
            return xdot

        x0 = [float(self.S[0])/self.N,float(self.I[0])/self.N]
        t = np.linspace(0,self.times[-1],100)
        x = odeint(deriv,x0,t)
        self.solution = x
        return None  