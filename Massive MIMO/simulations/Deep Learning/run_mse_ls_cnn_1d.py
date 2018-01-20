import sys
import os
# Generate dummy data
import numpy as np
import scipy.io as sio
import time

# Specify the number of receiving antennas at Base Stattion (BS)
M = 70
# Specify the number of cells.
L = 7
# Specify the number of users within each one of the cells.
K = 10
# Specify pilot-sequence length.
N = K
# Specify batch size
batch_size = 100
# Specify number of epochs.
number_of_epochs = 2

# Specify number of neurons in input layer.
input_layer_size = M*2
# Specify number of neurons in 1st hidden layer.
second_layer_size = 0#M*2
# Specify number of neurons in output layer.
output_layer_size = M*2

# Specify the kernel size.
kernel_size = 3

# Specify the optimizer.
optimizer = 'adam'

t = time.localtime()
timestamp = time.strftime('%d%Y%H%M%S', t)
filename = 'mse_ls_cnn_1d_%s.txt' % (timestamp)

#-----------------------------------------------------------------------

M = 70
L = 7
K = 10
N = K
batch_size = 100
number_of_epochs = 100


input_layer_size = M*2
second_layer_size = M*2
output_layer_size = M*2
kernel_size = 1
optimizer = 'sgd'
cmd = 'python deep_learning_mmimo_channel_estimation_loss_mse_ls_cnn_1d.py %d %d %d %d %d %d %d %d %d %d %s %s&' % (M,L,K,N,batch_size,number_of_epochs,input_layer_size,second_layer_size,output_layer_size,kernel_size,filename,optimizer)
os.system(cmd)


input_layer_size = M*2
second_layer_size = M*2
output_layer_size = M*2
kernel_size = 2
optimizer = 'sgd'
cmd = 'python deep_learning_mmimo_channel_estimation_loss_mse_ls_cnn_1d.py %d %d %d %d %d %d %d %d %d %d %s %s&' % (M,L,K,N,batch_size,number_of_epochs,input_layer_size,second_layer_size,output_layer_size,kernel_size,filename,optimizer)
os.system(cmd)


input_layer_size = M*2
second_layer_size = M*2
output_layer_size = M*2
kernel_size = 3
optimizer = 'sgd'
cmd = 'python deep_learning_mmimo_channel_estimation_loss_mse_ls_cnn_1d.py %d %d %d %d %d %d %d %d %d %d %s %s&' % (M,L,K,N,batch_size,number_of_epochs,input_layer_size,second_layer_size,output_layer_size,kernel_size,filename,optimizer)
os.system(cmd)


input_layer_size = M*2
second_layer_size = M*2
output_layer_size = M*2
kernel_size = 4
optimizer = 'sgd'
cmd = 'python deep_learning_mmimo_channel_estimation_loss_mse_ls_cnn_1d.py %d %d %d %d %d %d %d %d %d %d %s %s&' % (M,L,K,N,batch_size,number_of_epochs,input_layer_size,second_layer_size,output_layer_size,kernel_size,filename,optimizer)
os.system(cmd)

