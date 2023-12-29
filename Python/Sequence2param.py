import numpy as np;
import matplotlib.pyplot as plt;
import scipy.fft as fft;

class Sequence2param:


    def KLDivergence( P, Q):
        if len(P) != len(Q):
            return None;
        sum =0;

        for i in range(len(P)):
            sum = np.dot(P, np.log2( np.divide(P, Q)));

        return sum;

    def JSDivergence( P, Q):
        M = (P+Q)/2;
        JS = 0.5 * ( KLDivergence(P,M) + KLDivergence(Q, M));
        return JS;

    def __init__ (self, filepath):
        self.AList = np.loadtxt(filepath + "Acceleration.txt");
        self.VList = np.loadtxt(filepath + "LatticeDepth.txt");

        self.AVListIndex = np.loadtxt(filepath + "AVIndex.txt", dtype = int);
        self.MomProb = np.loadtxt(filepath + "MomentumProbability.txt"); # np array with rows containing momentum probabilities for each [a,V] value pair in AVList
        self.a0 = 0.0;
        self.V0 = 10.0;

    def BayesianUpdating(self, totalmeasurements, recordstep, totalrecords):
        PossibleMomentumOutcomes =np.array( [-10+2*i for i in range(0,11)]); # values of momentum in n\hbar k_L
        PossibleOutcomes = range(0,11);
        datamom =np.reshape(self.MomProb, (self.AList.size,self.VList.size,11));

        P_actual=np.array(datamom[50,25,:]); # Fix this hardcoded value

        P_actual=P_actual/np.sum(P_actual);
        P_simulated = P_actual; #No errors

        Runs=totalmeasurements; # How many simulated data do we want
        outcomes = np.random.default_rng().choice(PossibleOutcomes,size=Runs, p = P_simulated);
        unique, frequency = np.unique(outcomes, return_counts = True);

        PaVprior = np.full((self.AList.size, self.VList.size),1)/(self.AList.size*self.VList.size);

        plotPaV=np.array([PaVprior]);
        
        counter =0;
        for m in outcomes:
            for i in range(self.AList.size*self.VList.size):
                indexpair = self.AVListIndex[i];
                MomentumProbabilities = self.MomProb[i];

                PaVprior[indexpair[0], indexpair[1]] *= (MomentumProbabilities[m])
                PaVprior/=np.sum(PaVprior)

            counter+=1;
            if counter % recordstep and counter < totalrecords:
                plotPaV=np.append(plotPaV,[PaVprior], axis=0); 

        return plotPaV;


    def JSDMatrix (self):
        momindices=np.where(self.AVListIndex[:,1]==25)[0];
        momproblist =self.MomProb[momindices];
        indices = self.AVListIndex[momindices];

        accindices = indices[:, 0];
        acc = self.AList[accindices]

        no_of_values = len(accindices);
        JSDivergenceMatrix=np.zeros((no_of_values, no_of_values)) ;

        for i in accindices:
            for j in accindices:
                JSDivergenceMatrix [i][j]= JSDivergence(momproblist[i], momproblist[j])

        return JSDivergenceMatrix;