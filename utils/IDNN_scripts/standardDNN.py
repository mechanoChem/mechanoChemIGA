import tensorflow as tf
import numpy as np
import random
import shutil, os

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
        x = tf.concat([x1,tf.pow(x2,2)],1) #This applies symmetry only to the order parameter (see terms below as well)

        # Define variables
        W = [tf.Variable(tf.truncated_normal([Layers[0],Layers[1]], stddev=0.1),name='weights0')]
        b = [tf.Variable(0.1*tf.ones([Layers[1]]),name='bias0')]
        for j in range(1,n):
                W += [tf.Variable(tf.truncated_normal([Layers[j],Layers[j+1]], stddev=0.1),name='weights'+str(j))]
                b += [tf.Variable(0.1*tf.ones([Layers[j+1]]),name='bias'+str(j))]
        W += [tf.Variable(tf.truncated_normal([Layers[n],Layers[n+1]], stddev=0.1),name='weights'+str(n))]
        b += [tf.Variable(0.1*tf.ones([Layers[n+1]]),name='bias'+str(n))]
        
        z = tf.matmul(x,W[0]) + b[0]
        X = [actFun(z)]
        for j in range(1,n):
                z = tf.matmul(X[j-1],W[j]) + b[j]
                X += [actFun(z)]
        z = tf.matmul(X[n-1],W[n]) + b[n]
        y = z # prediction distribution based on weights/bias for current input
    
        predictions = tf.concat([y[:,0:1],2.*x2*y[:,1:]],1)
        
        # If called by 'predict' function...
        if mode == tf.estimator.ModeKeys.PREDICT:
                return tf.estimator.EstimatorSpec(mode=mode, predictions={'predictions': predictions})

        # Calculate loss using mean squared error
        loss = tf.losses.mean_squared_error(labels, predictions)
        optimizer = tf.train.AdagradOptimizer(learning_rate=params["learning_rate"])
        #optimizer = tf.train.ProximalAdagradOptimizer(learning_rate=params["learning_rate"],l1_regularization_strength=0.04)
        train_op = optimizer.minimize(loss=loss, global_step=tf.train.get_global_step())

        # Else, called by 'evaluate' or 'train' functions
        return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)

if __name__ == '__main__':

        # Read in data
        dataIn = np.genfromtxt('bcc_lattice_monte_carlo.txt',dtype=np.float32)[:,1:]
        features = dataIn[:,0:2]
        labels = np.hstack((dataIn[:,2:3],dataIn[:,3:4]))

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
        
        # Reset graph directory
        for rnd in range(10):
                restart = True
                if restart:
                        shutil.rmtree('stan_graph{}'.format(rnd),ignore_errors=True)
                        open('stan_training{}.txt'.format(rnd),'w').close()

                # Get DNN
                model_params = {"learning_rate": 0.1, "layers": [2,10,10,2]} #100 @0.1 (used for paper - comparison standard DNN w/ IDNN training curve)
                dnn = tf.estimator.Estimator(model_fn=model_fn, params=model_params, model_dir='stan_graph{}'.format(rnd))

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
                fout = open('stan_training{}.txt'.format(rnd),'a')
                train_loss = dnn.evaluate(input_fn=train_loss_fn)["loss"]
                valid_loss = dnn.evaluate(input_fn=valid_loss_fn)["loss"]
                print "Iteration: ",0,", Training loss: ",train_loss,", Validation loss: ",valid_loss
                fout.write("Epoch: {}, Iteration: {}, Training loss: {}, Validation loss: {}\n".format(0,0,train_loss,valid_loss))
                for i in range(20):
                        dnn.train(input_fn=train_input_fn, steps=10506)
                        train_loss = dnn.evaluate(input_fn=train_loss_fn)["loss"]
                        valid_loss = dnn.evaluate(input_fn=valid_loss_fn)["loss"]
                        print "Iteration: ",(i+1)*10506,", Training loss: ",train_loss,", Validation loss: ",valid_loss
                        fout.write("Epoch: {}, Iteration: {}, Training loss: {}, Validation loss: {}\n".format(i+1,(i+1)*10506,train_loss,valid_loss))
                fout.close()
