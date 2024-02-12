import tensorflow as tf
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.animation import FFMpegWriter
from PDE_Classes1 import *
from PINN_Modular_Classes1 import *
#from silence_tensorflow import silence_tensorflow
#silence_tensorflow()
import time
import argparse
# import itertools as product
import os
os.environ["KMP_AFFINITY"] = "noverbose"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

#point FFMpeg to the executable
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

# Set data type
DTYPE = 'float64'
tf.keras.backend.set_floatx(DTYPE)

# ----- START FIXED INPUTS -----

# Initialize order parameter for PDE
order = 2

# number of internal collocation points
N = 2**10
# number of boundary and interface points
N_b = 2**6

# Initialize list of points which lie on the system boundaries
domain = [0, 1]

beta = 1

# Declare constant hyperparameters
alpha = 0.2
numEpochs = 2**10
learn_rate = 5e-4#tf.keras.optimizers.schedules.PiecewiseConstantDecay([2**6],[5e-2,1e-3])
schwarz_tol = 1e-3
err_tol = 1e-2

# This parameter is true when you are retrieving a file for a parameter sweep
# Otherwise you are drawing a single test case from the manual test file
retrieve = 1

# ----- END FIXED INPUTS -----

def driver(parameter_file, outdir):

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Initialize dataframe for parameters and results storage
    if retrieve:
        ParameterSweep = pd.read_csv(parameter_file)
    else:
        PS = pd.read_csv(parameter_file)
        pos = 9+PS.loc[0, 'Sub-domains']
        ParameterSweep = PS.drop(PS.iloc[:,pos:len(PS.columns)], axis=1)

    # parameter sweep loop
    for z in range(len(ParameterSweep.index)):

        # Set percentage overlap
        percent_overlap = ParameterSweep.loc[z, 'Percent Overlap']

        # Set the number of subdomains and the desired percentage overlap
        n_subdomains = ParameterSweep.loc[z, 'Sub-domains']

        BC_label = '_WDBC'
        # Choose type of BC enforcement for system BCs (strong (1) or weak (0))
        if ParameterSweep.loc[z, 'System BCs'] == 'weak':
            sysBC = 0
        elif ParameterSweep.loc[z, 'System BCs'] == 'strong':
            sysBC = 1

        # Choose type of BC enforcement for Schwarz BCs (strong (1) or weak (0))
        if ParameterSweep.loc[z, 'Schwarz BCs'] == 'weak':
            schBC = 0
        elif ParameterSweep.loc[z, 'Schwarz BCs'] == 'strong':
            schBC = 1

        # assertion to catch SChwarz boundary enforcement when only domain is used
        if schBC:
            assert n_subdomains > 1, "Schwarz boundaries cannot be strongly enforced because you are using a single-domain model which has no Schwarz boundaries."

        # Store BC enforcement booleans together
        SDBC = [sysBC, schBC]
        if all(SDBC):
            BC_label = '_SDBC_both'
        elif SDBC[0]:
            BC_label = '_SDBC_sys'
        elif SDBC[1]:
            BC_label = '_SDBC_schwarz'

        NN_label = 'PINN'
        # number of snapshots per subdomain used to aid training
        if ParameterSweep.loc[z, 'Snapshots']:
            snap = 2**6
            NN_label = 'NN'
        else:
            snap=0

        # Set subdomains for FOM modeling, if any
        FOM_label = ''
        sub_FOM = np.zeros((n_subdomains,))
        if ParameterSweep.loc[z, 'FOM'] == 'left':
            sub_FOM[0] = 1
            FOM_label = '_FOM_left'
        elif ParameterSweep.loc[z, 'FOM'] == 'right':
            sub_FOM[-1] = 1
            FOM_label = '_FOM_right'

        # Calculate the size of subdomains and overlap to achieve the parameters defined above
        domain_size =  1 / ( n_subdomains*(1 - percent_overlap) + percent_overlap )
        overlap = percent_overlap*domain_size

        # Construct subdomain set
        sub = ()
        for i in range(n_subdomains):
            step = i*(domain_size - overlap)
            sub += ([step, step+domain_size],)

        # Truncate subdomain boundaries to 8 decimal precision
        sub = np.round(sub, 8)

        # Generate boundary points for each subdomain boundary
        X_b_om = [ [ tf.constant(np.tile([i], (N_b, 1)), dtype=DTYPE) for i in sub[j] ] for j in range(len(sub))]

        # Set random seed for reproducible results
        tf.random.set_seed(0)

        # Declare nu based on Peclet number
        Pe = ParameterSweep.loc[z, 'Peclet Number']
        nu = 1/Pe
        # Declare an instance of the PDE class
        pde1 = PDE_1D_Steady_AdvecDiff(nu=nu, order=order)

        # FOM value "true solution" to use as a reference
        x_true = tf.constant(np.linspace(domain[0], domain[1], num=N), shape=(N, 1), dtype=DTYPE)
        u_true = pde1.f(x_true)

        # Set number of FD points and initialize the step size
        n_FD = int(N/2)

        # Set number of hidden layers and nodes per layer
        hl = 2
        nl = 20

        # Initialize tuple to store internal points and models in each subdomain
        X_r_om = ()
        model_om = ()

        # This alernative mode of quasi-random sample point generation may be still be used
    #     # Generate uniform internal points for the whole domain
    #     sampler = stats.qmc.Sobol(1, scramble=False)
    #     points = sampler.random(N)
    #     xwidth = domain[1] - domain[0]
    #     x = xwidth*points + domain[0]

        # Build and store neural networks for applicable subdomains; Store FD class for FOM sub-domains
        for i in range(sub.shape[0]):

            # If subdomain is to be modeled by NN, add random uniform points to subdomain, include sub-domain bounds
            if not sub_FOM[i]:
                xl = sub[i][0]
                xr = sub[i][1]
                temp = np.random.uniform(low=xl, high=xr, size=(N,))
                temp = np.append(temp, sub[i])
    #             temp = (x[(x >= sub[i][0]) & (x <= sub[i][1])] - sub[0][0]) / (sub[-1][1] - sub[0][0])
    #             temp.reshape((len(temp),1))
                X_r_om += ( tf.constant(temp, shape=(temp.shape[0],1), dtype=DTYPE), )

                model_om += ( PINN_Architecture(xl=xl, xr=xr, num_hidden_layers=hl, num_neurons_per_layer=nl), )
                model_om[i].build(input_shape=(None, X_r_om[i].shape[1]))

                continue

            # If subdomain is to be modeled by FD, construct uniform line space for FD on subdomain
            xl = sub[i][0]
            xr = sub[i][1]
            x_FD = np.linspace(xl, xr, num=n_FD)

            # Add line space for subdomain to internal point storage
            X_r_om += ( tf.constant(x_FD, shape=(x_FD.shape[0],1), dtype=DTYPE), )

            model = FD_1D_Steady(x_FD, domain, pde1)

            model_om += ( model, )


        # Initialize schwarz loop operators
        schwarz_conv = 1
        ref_err = 1
        iterCount = 0
        x_schwarz = [tf.constant(np.linspace(s[0], s[1], num=n_FD), shape=(n_FD, 1), dtype=DTYPE) for s in sub]
        u_i_minus1 = [tf.constant(np.random.rand(n_FD,1), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]
        u_i = [tf.constant(np.zeros((n_FD,1)), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]

        # Initialize variables for plotting Schwarz results
        fig = plt.figure(layout='tight')
        #fig.tight_layout(pad=1.0)
        fig.set_size_inches(6, 3)
        plt.xlabel('x', fontsize=14)
        plt.ylabel('u(x)', fontsize=14)
        ref, = plt.plot(x_true, u_true, 'k--')
        subdomain_plots = [plt.plot([], []) for _ in range(n_subdomains)]

        #Set animation parameters
        fps = 6
        metadata = dict(title='SchwarzRecord', author='Will Snyder')
        #writer = PillowWriter(fps=fps, metadata=metadata)
        writer = FFMpegWriter(fps=fps, metadata=metadata)

        # Declare writer to begin recording animations
        mp4file = os.path.join(
            outdir,
            NN_label+BC_label+FOM_label+"_Pe_{:d}_nSub_{:d}_over{:f}.mp4".format(int(beta/nu), n_subdomains, percent_overlap)
        )
        with writer.saving(fig, mp4file, 600):

            # Begin recording time
            start = time.time()

            # Record u_hat at each boundary to stack on each iteration
            u_gamma = np.zeros(sub.shape)

            # Main Schwarz loop
            while (schwarz_conv > schwarz_tol or ref_err > err_tol):

                # initialize error check variables to 0 as they are to be added to during sub-domain iteration
                schwarz_conv = 0
                ref_err = 0

                # Add to Schwarz iteration count
                iterCount += 1

                # Update title for Schwarz iter
                plt.title('Schwarz iteration {:d}; Pe = {:d}'.format(iterCount, Pe), fontsize=14)

                # loop over each model for training
                for s in range(len(model_om)):

                    # Current model domain points
                    X_r = X_r_om[s]
                    # Current model boundary points
                    X_b = X_b_om[s]

                    # Current model
                    model_r = model_om[s]
                    # Adjacent models for interface conditions
                    model_i = [model_om[s-1:s], model_om[s+1:s+2]]

                    # update the current models BCs according to adjacent models and save to model for SDBCs
                    if any(SDBC):
                        if model_i[0]:
                            u_gamma[s][0] = np.interp(sub[s][0], x_schwarz[s-1][:,0], u_i[s-1][:,0])
                        if model_i[1]:
                            u_gamma[s][1] = np.interp(sub[s][1], x_schwarz[s+1][:,0], u_i[s+1][:,0])
                    model_r.u_gamma = u_gamma[s]

                    # Initialize solver
                    p = PINN_Schwarz_Steady(pde1, model_r, model_i, SDBC, X_r, X_b, alpha, snap)

                    # Solve model for current sub-domain
                    p.solve(tf.keras.optimizers.Adam(learning_rate=learn_rate), numEpochs)

                    MSE_tag = "MSE Sub-domain {:d}".format(s+1)

                    # If current model is FD, output FD Error, update current iteration solution and convergence metrics
                    if isinstance(model_r, FD_1D_Steady):
                        print('Model {:d}: '.format(s+1))
                        print('\t'+'Finite Difference error = {:10.8e}'.format(p.err))

                        ParameterSweep.loc[z, MSE_tag] = np.float64(p.err)

                        u_i[s] = model_r(x_schwarz[s])
                        schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s])/tf.math.reduce_euclidean_norm(u_i[s])
                        ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s]))/tf.math.reduce_euclidean_norm(pde1.f(x_schwarz[s]))

                    # If current model is NN output the model loss, update current iteration solution and convergence metrics
                    else:
                        print('Model {:d}: '.format(s+1))
                        print('\t'+'Residual loss = {:10.8e}'.format(p.phi_r))
                        if not schBC:
                            print('\t'+'Interface loss = {:10.8e}'.format(p.phi_i))
                        if not sysBC:
                            print('\t'+'Boundary loss = {:10.8e}'.format(p.phi_b))
                        if snap:
                            print('\t'+'Snapshot loss = {:10.8e}'.format(p.phi_s))
                        print('\t'+'Total loss = {:10.8e}'.format(p.loss))

                        ParameterSweep.loc[z, MSE_tag] = np.float64(p.loss)

                        if any(SDBC):
                            u_i[s] = p.get_u_hat(x_schwarz[s], model_r)
                        else:
                            u_i[s] = model_r(x_schwarz[s])

                        schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s])/tf.math.reduce_euclidean_norm(u_i[s])
                        ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s]))/tf.math.reduce_euclidean_norm(pde1.f(x_schwarz[s]))


                    # update frame for animation
                    subdomain_plots[s][0].set_data(x_schwarz[s][:,0], u_i[s])

                    # Gather frame for animation
                    writer.grab_frame()


                # Calculate the normalized difference between u for the current iteration and u for the previous iteration
                schwarz_conv = schwarz_conv/len(u_i)
                ref_err = ref_err/len(u_i)

                # Update the value of u stored for the previous iteration
                u_i_minus1 = u_i.copy()

                # Save current frame as image if needed for beamer animations
                #fig.savefig("Frames/SchwarzPINN_WDBC_Pe_{:d}_nSub_{:d}-{:d}.pdf".format(int(beta/nu), n_subdomains, iterCount),dpi=600)

                # Output current Schwarz error
                print('\nSchwarz iteration {:d}: Convergence error = {:10.8e}, Reference Error = {:10.8e}'.format(iterCount, schwarz_conv, ref_err), "\n")

                # Cut off simulat at a particular number of Schwarz iterations
                if (iterCount == 100):
                    break

            # End time recording
            end = time.time()

            # Capture 2 more seconds on final frame to pad out video
            for _ in range(fps*2):
                writer.grab_frame()

            # Record results for coupled model
            ParameterSweep.loc[z, 'CPU Time (s)'] = end-start
            ParameterSweep.loc[z, 'Schwarz Iterations'] = iterCount
            ParameterSweep.loc[z, 'Avg L2 Error'] = np.float64(ref_err)

        # Export final results to CSV
        ParameterSweep.to_csv(os.path.join(outdir, "ParamResults.csv"))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='PINN-Schwarz',
        description='Schwarz-based training of PINN-PINN coupling')

    parser.add_argument('parameter_file')
    parser.add_argument('outdir')
    args = parser.parse_args()

    driver(args.parameter_file, args.outdir)