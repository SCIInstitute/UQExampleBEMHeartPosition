from itertools import chain, combinations
import scipy.io
import os
import matplotlib.pyplot as plt

import numpy as np

from UncertainSCI.distributions import BetaDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion

#place for UQ files
output_dir = "/Users/jess/software/UQExampleBEMHeartPosition/"

# Number of parameters
dimension = 4

# Specifies 1D distribution on [0,1] (alpha=beta=1 ---> uniform)
alpha = 1
beta = 1.
# domain is the range of the hypercube
domain = np.array([[-30, 30], [-30, 30], [-30, 30], [-30, 30]]).T
dist = BetaDistribution(alpha=alpha, beta=beta, dim=dimension, domain=domain)

# # Expressivity setup
order = 3
indices = TotalDegreeSet(dim=dimension, order=order)
pce = PolynomialChaosExpansion(indices, dist)
pce.generate_samples()

#samples = pce.samples
#print(samples)


# solution points (771 nodes x 21 time points)
N = 16191


def run_scirun_model_bash(samples,N):

    scirun_call = "/Users/jess/software/SCIRun/bin_clean/SCIRun/SCIRun_test"
#    scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
    network_file = output_dir+"nets/torso_position_model_all.srn5"
    tmp_dir = output_dir+"tmp/"
    

def run_scirun_model(samples,N):

    # convert samples into transformation matrices, save to disk, load in SCIRun, run model, save solutions to disk, load back to script.
    scirun_call = "/Users/jess/software/SCIRun/bin_clean/SCIRun/SCIRun_test"
#    scirun_call = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"
    network_file = output_dir+"nets/torso_position_model_all.srn5"
    tmp_dir = output_dir+"UQ/tmp/"
    
    

    # experiment files.  It is here for now to make some things easier
    vect_file = output_dir+"data/geom/s_vect.mat"
    heart_pots = output_dir+"Run0001-cappedsock.mat"
    heart_geom = output_dir+"data/geom/capped_sock.mat"
    torso_geom = output_dir+"data/geom/Full_tank_771.mat"
    
    
    script_file_tmp = os.path.join(tmp_dir, "UQ_torso_tmp.py")
    SR_output_file = os.path.join(tmp_dir, "UQ_torso_SR_solutions.txt")
    samples_file = os.path.join(tmp_dir, "UQ_torso_samples.mat")
    
    
    trans_samples = make_transforms(samples, vect_file)
    scipy.io.savemat(samples_file, dict(samples=trans_samples))
    
    s_file=open(script_file_tmp,'w+')
    s_file.write("scirun_load_network('"+network_file+"')\n")
    s_file.write("scirun_set_module_state('ReadField:1','Filename','"+heart_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadField:0','Filename','"+torso_geom+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:0','Filename','"+heart_pots+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:1','Filename','"+samples_file+"')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','FileTypeName','SimpleTextFile')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:1','Filename','n_comp.txt')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','Filename','"+SR_output_file+"')\n")
    s_file.write("scirun_set_module_state('WriteMatrix:0','FileTypeName','SimpleTextFile (*.*)')\n")
    
#    s_file.write("scirun_execute_all()\n")
    s_file.close()
    
#    print( scirun_call+" -0 -S "+script_file_tmp)
    print( "run this script in SCIRun: "+script_file_tmp)
#    output=os.system(scirun_call+" -0  -S "+script_file_tmp)
    input("press enter when SCIRun net is finished.")
    
    print(SR_output_file)
    
    pot_solution = np.loadtxt(SR_output_file)
    
#    os.remove(script_file_tmp)
#    os.remove(output_file)

    return pot_solution.T
    
    
    
    
    
    
def make_transforms(samples, vect_file):

    transformed_samples = np.zeros([samples.shape[0],12])

    m_vect = scipy.io.loadmat(vect_file)

    s_evecs = m_vect['s_evecs']
    cent = m_vect['cent']
    radii = m_vect['radii']

    h_length = radii[0][0]
    
    print(h_length)


    v1 = [0,1,0]
    v2 = [1,0,0]
    v3 = [s_evecs[0][0],s_evecs[1][0], s_evecs[2][0]]


    for ind in range(samples.shape[0]):
        weights = samples[ind,:]
        theta1 = np.arcsin(weights[0]/h_length)
        theta2 = np.arcsin(weights[1]/h_length)
        theta3 = np.radians(weights[2])
        c1 = np.cos(theta1)
        s1 = np.sin(theta1)
        c2 = np.cos(theta2)
        s2 = np.sin(theta2)
        c3 = np.cos(theta3)
        s3 = np.sin(theta3)
        T1 = np.array([ [ 1,0, 0, - cent[0][0]], [ 0,1, 0, - cent[1][0]], [ 0, 0, 1, - cent[2][0]], [0, 0, 0, 1] ])
        T2 = np.array([ [ 1,0, 0, cent[0][0]], [ 0,1, 0, cent[1][0]], [ 0, 0, 1, cent[2][0]], [0, 0, 0, 1]  ])
        Ry = np.array( [ [ c1, 0, s1, 0], [0,1,0,0], [-s1, 0, c1, 0], [0,0,0,1] ])
        Rx = np.array( [ [ 1,0,0,0], [0, c2, -s2, 0], [0, s2, c2, 0], [0,0,0,1] ])
        R3 = np.array( [ [ c3+v3[0]*v3[0]*(1-c3), v3[0]*v3[1]*(1-c3)-v3[2]*s3,  v3[0]*v3[2]*(1-c3)+v3[1]*s3 , 0],
                                  [ v3[1]*v3[0]*(1-c3)+v3[2]*s3,  c3+v3[1]*v3[1]*(1-c3) , v3[2]*v3[1]*(1-c3)-v3[0]*s3 , 0],
                                  [ v3[0]*v3[2]*(1-c3)-v3[1]*s3,  v3[2]*v3[1]*(1-c3)+v3[0]*s3,  c3+v3[2]*v3[2]*(1-c3), 0],
                                  [ 0, 0, 0, 1] ])
        Tz = np.vstack((np.hstack((np.identity(3),np.array([[0],[0],[weights[3]]]))),np.array([[0,0,0,1]])))
        Transform = np.dot(Tz,np.dot(T2,np.dot(Rx,np.dot(Ry,np.dot(R3,T1)))))
        print(Transform)
        transformed_samples[ind, :] = Transform[0:3,:].flatten()
        
    return transformed_samples
    
# run model to get solutions
model_output = run_scirun_model(pce.samples,N)
pce.build(model_output=model_output)

mean = pce.mean()
stdev = pce.stdev()

# Power set of [0, 1, ..., dimension-1]
variable_interactions = list(chain.from_iterable(combinations(range(dimension), r) for r in range(1, dimension+1)))

# "Total sensitivity" is a non-partitive relative sensitivity measure per parameter.
total_sensitivity = pce.total_sensitivity()

# "Global sensitivity" is a partitive relative sensitivity measure per set of parameters.
global_sensitivity = pce.global_sensitivity(variable_interactions)

Q = 3  # Number of quantile bands to plot
dq = 0.5/(Q+1)
q_lower = np.arange(dq, 0.5-1e-7, dq)[::-1]
q_upper = np.arange(0.5 + dq, 1.0-1e-7, dq)
quantile_levels = np.append(np.concatenate((q_lower, q_upper)), 0.5)

quantiles = pce.quantile(quantile_levels, M=int(2e3))

median = median = pce.quantile(0.5, M=int(1e3))

median = median.reshape((771,383))
quantiles = quantiles.reshape((quantiles.shape[0],771,383))

mo_shape = model_output.shape
pots = model_output(reshape,(mo_shape[0],771,383))


UQ_file = os.path.join(output_dir+"UQ_data", "torso_position_UQ_values.mat")

scipy.io.savemat(UQ_file, dict(mean = mean.T, stdev = stdev.T, tot_sensitivity = total_sensitivity.T, glob_sensitivity = global_sensitivity.T, quantiles = quantiles.T, median = median.T))
    

lead_num = 404


print(model_output.shape)
plt.figure()

#plt.plot(model_output[:V, :].T, 'k', alpha=0.8, linewidth=0.2)
plt.plot(model_output[0:383,lead_num],median[lead_num,:], 'b', label='PCE median')

band_mass = 1/(2*(Q+1))

for ind in range(Q):
    alpha = (Q-ind) * 1/Q - (1/(2*Q))
    if ind == 0:
        plt.fill_between(np.array(range(383)),quantiles[ind,lead_num, :], quantiles[Q+ind,lead_num, :],
                         interpolate=True, facecolor='red', alpha=alpha,
                         label='{0:1.2f} probability mass (each band)'.format(band_mass))
    else:
        plt.fill_between(np.array(range(383)),quantiles[ind,lead_num, :], quantiles[Q+ind,lead_num, :], interpolate=True, facecolor='red', alpha=alpha)
