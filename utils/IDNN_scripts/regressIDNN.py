import tensorflow as tf
import numpy as np
import random
import shutil, os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def actFun(x):
        return tf.log(1. + tf.exp(x)) #SoftPlus

def actFunDer(x):
        return 1./(1. + tf.exp(-x)) #Logistic

def model_fn(features, labels, mode, params):

        # Define DNN
        Layers = params['layers']
        n = len(Layers)-2

        # Apply symmetry constraints on y through input function
        x1 = features['x'][:,0:1]
        x2 = features['x'][:,1:]
        x = tf.concat([x1,tf.pow(x2,2)],1) #This applies symmetry only to the order parameter (see lines below as well)
        #x = tf.concat([tf.pow(x1-1./np.sqrt(2),2),tf.pow(x2,2)],1) #This applies symmetry to both

        # Define variables
        W = [tf.Variable(tf.truncated_normal([Layers[0],Layers[1]], stddev=0.1),name='weights0')]
        b = [tf.Variable(0.1*tf.ones([Layers[1]]),name='bias0')]
        for j in range(1,n):
                W += [tf.Variable(tf.truncated_normal([Layers[j],Layers[j+1]], stddev=0.1),name='weights'+str(j))]
                b += [tf.Variable(0.1*tf.ones([Layers[j+1]]),name='bias'+str(j))]
        W += [tf.Variable(tf.truncated_normal([Layers[n],Layers[n+1]], stddev=0.1),name='weights'+str(n))]

        z = tf.matmul(x,W[0]) + b[0]
        X = [actFun(z)]
        Xp1 = [actFunDer(z)*W[0][0]]
        Xp2 = [actFunDer(z)*W[0][1]]
        for j in range(1,n):
                z = tf.matmul(X[j-1],W[j]) + b[j]
                X += [actFun(z)]
                Xp1 += [actFunDer(z)*tf.matmul(Xp1[j-1],W[j])]
                Xp2 += [actFunDer(z)*tf.matmul(Xp2[j-1],W[j])]
        z = tf.matmul(X[n-1],W[n])# + b[n]
        y = z # prediction distribution based on weights/bias for current input
    
        dydx1 = tf.matmul(Xp1[n-1],W[n]) #This applies symmetry only to the order parameter
        #dydx1 = 2.*(x1-1./np.sqrt(2))*tf.matmul(Xp1[n-1],W[n]) #This applies symmetry to both
        dydx2 = 2.*x2*tf.matmul(Xp2[n-1],W[n])

        # Reshape output layer to return predictions
        predictions = tf.concat([dydx1,dydx2,y],1)
    
        # If called by 'predict' function...
        if mode == tf.estimator.ModeKeys.PREDICT:
                return tf.estimator.EstimatorSpec(mode=mode, predictions={'predictions': predictions})

        # Calculate loss using mean squared error
        loss = tf.losses.mean_squared_error(labels, predictions[:,0:2])
        optimizer = tf.train.AdagradOptimizer(learning_rate=params["learning_rate"])
        #optimizer = tf.train.ProximalAdagradOptimizer(learning_rate=params["learning_rate"],l1_regularization_strength=0.04)
        train_op = optimizer.minimize(loss=loss, global_step=tf.train.get_global_step())

        # Else, called by 'evaluate' or 'train' functions
        return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)

if __name__ == '__main__':

        # Read in data
        dataIn = np.genfromtxt('bcc_lattice_monte_carlo.txt',dtype=np.float32)[:,1:]
        features = dataIn[:,0:2]
        labels = dataIn[:,2:4]

        # Scale the output for better training
        labels = 100.*labels

        # Separate into training and validation data
        N = features.shape[0]
        ind = np.arange(N)
        np.random.seed(0) # For repeatedly using the same training/validation split
        np.random.shuffle(ind)
        np.random.seed() # Reset the seed randomly
        x_train = features[ind[:75*N/100]]
        x_valid = features[ind[75*N/100:]]
        y_train = labels[ind[:75*N/100]]
        y_valid = labels[ind[75*N/100:]]

        print x_train.shape, x_valid.shape
        
        # Reset graph directory
        for rnd in range(1):
                restart = False #True
                if restart:
                        shutil.rmtree('graph{}'.format(rnd),ignore_errors=True)
                        open('IDNNtraining{}.txt'.format(rnd),'w').close()

                # Get DNN
                model_params = {"learning_rate": 0.1, "layers": [2,10,10,1]}
                dnn = tf.estimator.Estimator(model_fn=model_fn, params=model_params, model_dir='graph{}'.format(rnd))

                # Fit (train) model
                batch_size=10
                train_input_fn = tf.estimator.inputs.numpy_input_fn(x={"x": x_train},
                                                                    y=y_train,
                                                                    batch_size=batch_size,
                                                                    num_epochs=None,
                                                                    shuffle=True)

                # Validate
                train_loss_fn = tf.estimator.inputs.numpy_input_fn(x={"x": x_train},
                                                                   y=y_train,
                                                                   batch_size=batch_size,
                                                                   num_epochs=1,
                                                                   shuffle=False)
                valid_loss_fn = tf.estimator.inputs.numpy_input_fn(x={"x": x_valid},
                                                                   y=y_valid,
                                                                   batch_size=batch_size,
                                                                   num_epochs=1,
                                                                   shuffle=False)

                # Train
                fout = open('IDNNtraining{}.txt'.format(rnd),'a')
                train_loss = dnn.evaluate(input_fn=train_loss_fn)["loss"]
                valid_loss = dnn.evaluate(input_fn=valid_loss_fn)["loss"]
                print "Iteration: ",0,", Training loss: ",train_loss,", Validation loss: ",valid_loss
                fout.write("Epoch: {}, Iteration: {}, Training loss: {}, Validation loss: {}\n".format(0,0,train_loss,valid_loss))
                
                for i in range(20):
                        dnn.train(input_fn=train_input_fn, steps=10506)
                        train_loss = dnn.evaluate(input_fn=train_loss_fn)["loss"]
                        valid_loss = dnn.evaluate(input_fn=valid_loss_fn)["loss"]
                        print "Epoch: ",i+1,"Iteration: ",(i+1)*10506,", Training loss: ",train_loss,", Validation loss: ",valid_loss
                        fout.write("Epoch: {}, Iteration: {}, Training loss: {}, Validation loss: {}, Learning rate: {}\n".format(i+1,(i+1)*10506,train_loss,valid_loss,model_params["learning_rate"]))
                fout.close()

                # Create test set
                x_min = np.amin(features,axis=0)
                x_max = np.amax(features,axis=0)
                N = 100
                eta_test = np.mgrid[0.:1.:N*1j, 0.:1.:N*1j].astype(np.float32).reshape(2,-1).T
                x_test = np.copy(eta_test)
                x_test[:,0] = (eta_test[:,0]+eta_test[:,1])/np.sqrt(2)
                x_test[:,1] = (eta_test[:,0]-eta_test[:,1])/np.sqrt(2)
    

                test_input_fn = tf.estimator.inputs.numpy_input_fn(x={"x": x_test},
                                                                   num_epochs=1,
                                                                   shuffle=False)
    
                predictions = list(dnn.predict(input_fn=test_input_fn));
                y_test = 0.01*np.array([p['predictions'] for p in predictions])

                n = len(model_params['layers'])-2
                for i in range(0,n):
                        weights = dnn.get_variable_value('weights'+str(i))
                        biases = dnn.get_variable_value('bias'+str(i))
                        np.savetxt('weights_'+str(i)+'.txt',weights,header=str(weights.shape[0])+' '+str(weights.shape[1]))
                        np.savetxt('bias_'+str(i)+'.txt',biases,header=str(biases.shape[0]))
                weights = dnn.get_variable_value('weights'+str(n))
                np.savetxt('weights_'+str(n)+'.txt',weights,header=str(weights.shape[0])+' '+str(weights.shape[1]))

                indices = random.sample(range(features.shape[0]),500) # Plot random selection of data points
                x1 = x_test[:,0].reshape(N,-1)
                x2 = x_test[:,1].reshape(N,-1)
                e1 = (x1+x2)/np.sqrt(2)
                e2 = (x1-x2)/np.sqrt(2)
                y1 = y_test[:,0].reshape(N,-1)
                y2 = y_test[:,1].reshape(N,-1)
                Y = y_test[:,2].reshape(N,-1)

                X1 = features[indices,0]
                X2 = features[indices,1]
                E1 = (X1+X2)/np.sqrt(2)
                E2 = (X1-X2)/np.sqrt(2)
                Y1 = labels[indices,0]
                Y2 = labels[indices,1]

                # Save data for plotting elsewhere
                np.savetxt('dnn_test.txt',np.hstack((x_test,y_test)),header='c eta excess_cp1 excess_cp2 excess_fe')

                # Now, plot here
                
                f = plt.figure(1)
                ax1 = f.add_subplot(111, projection='3d')
                ax1.scatter(X1, X2,
                            0.01*Y1,
                            c='r', marker='o', label='Actual')
                ax1.plot_surface(x1, x2,
                                 y1,
                                 cmap=cm.coolwarm,linewidth=0,rstride=1,cstride=1)
                ax1.set_xlabel('composition')
                ax1.set_ylabel('order parameter')
                ax1.set_title('Excess chem. pot. (comp)')
                
                f.show()
                
                g = plt.figure(2)
                ax2 = g.add_subplot(111, projection='3d')
                ax2.scatter(X1, X2,
                            0.01*Y2,
                            c='r', marker='o', label='Actual')
                ax2.plot_surface(x1, x2,
                                 y2,
                                 cmap=cm.coolwarm,linewidth=0,rstride=1,cstride=1)
                ax2.set_xlabel('composition')
                ax2.set_ylabel('order parameter')
                ax2.set_title('Excess chem. pot. (order par.)')

                g.show()
                
                h = plt.figure(3)
                ax3 = h.add_subplot(111, projection='3d')
                ax3.plot_surface(x1, x2,
                                 Y,
                                 cmap=cm.coolwarm,linewidth=0,rstride=1,cstride=1)
                ax3.set_xlabel('composition')
                ax3.set_ylabel('order parameter')
                ax3.set_title('Free energy')

                h.show()
                
                raw_input()
    
