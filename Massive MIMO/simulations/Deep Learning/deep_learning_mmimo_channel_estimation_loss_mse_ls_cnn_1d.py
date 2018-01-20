import sys
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD
from keras.layers import Embedding
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D, Flatten

# Generate dummy data
import numpy as np
import scipy.io as sio

# Specify the number of receiving antennas at Base Stattion (BS)
M = int(sys.argv[1]) #70
# Specify the number of cells.
L = int(sys.argv[2]) #7
# Specify the number of users within each one of the cells.
K = int(sys.argv[3]) #10
# Specify pilot-sequence length.
N = int(sys.argv[4]) #K
# Specify batch size
batch_size = int(sys.argv[5]) #100
# Specify number of epochs.
number_of_epochs = int(sys.argv[6]) #100

# Specify number of neurons in input layer.
input_layer_size = int(sys.argv[7])
# Specify number of neurons in 1st hidden layer.
second_layer_size = int(sys.argv[8])
# Specify number of neurons in output layer.
output_layer_size = int(sys.argv[9])
# Specify the kernel size.
kernel_size = int(sys.argv[10]) #3

# Open file.
filename = sys.argv[11]
fid = open(filename, 'a')

# Load training and test vectors.
data_set = sio.loadmat('data_set_M_70_K_10_SNR_10_static_scenario_1_ls_est_v1.mat')
x_train = data_set['train_data']
x_train = np.expand_dims(x_train, axis=2)
y_train = data_set['train_label']
x_test = data_set['test_data']
x_test = np.expand_dims(x_test, axis=2)
y_test = data_set['test_label']
x_predict = data_set['prediction_data']
x_predict = np.expand_dims(x_predict, axis=2)
y_predict = data_set['prediction_label']
predict_error = data_set['error_prediction']

# Create Model.
model = Sequential()
model.add(Conv1D(input_layer_size, kernel_size, activation='tanh', input_shape=(M*2, 1)))
#model.add(Conv1D(input_layer_size, kernel_size, activation='tanh'))
model.add(MaxPooling1D(kernel_size))
if(second_layer_size > 0):
   model.add(Conv1D(second_layer_size, kernel_size, activation='tanh'))
   #model.add(Conv1D(second_layer_size, kernel_size, activation='tanh'))
   model.add(MaxPooling1D(kernel_size))
model.add(Flatten())
#model.add(GlobalAveragePooling1D())
#model.add(Dropout(0.5))
model.add(Dense(output_layer_size, activation='tanh'))

# Set optimizer parameters.
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
adam = keras.optimizers.Adam()

# Compile model.
model.compile(loss='mean_squared_error', optimizer=adam, metrics=['accuracy'])

# Train model.
train_history = model.fit(x_train, y_train, epochs=number_of_epochs, batch_size=batch_size)

# Test model.
score = model.evaluate(x_test, y_test, batch_size=batch_size)

print('Test loss:', score[0])
print('Test accuracy:', score[1])

# Calculate predictions
if(1):
   predictions = model.predict(x_predict)

   error_prediction_dl = np.zeros(len(y_predict));
   total_avg_error_dl = 0;
   total_avg_error_ls = 0;
   for i in range(len(y_predict)):
      for j in range(2*M):
         #print('expected[%d]: %f - actual[%d]: %f' % (j,y_predict[i][j],j,predictions[i][j]))
         error_prediction_dl[i] = error_prediction_dl[i] + abs(y_predict[i][j] - predictions[i][j]);

      error_prediction_dl[i] = error_prediction_dl[i]/(2*M)
      total_avg_error_dl = total_avg_error_dl + error_prediction_dl[i]
      total_avg_error_ls = total_avg_error_ls + predict_error[0][i]
      #print('prediction error dl[%d]: %f - prediction error ls[%d]: %f' % (i,error_prediction_dl[i],i,predict_error[0][i]))

   print('Avg. DL: %f - Avg. LS: %f' % (total_avg_error_dl/len(y_predict),total_avg_error_ls/len(y_predict)))

# Write result into file.
if(second_layer_size > 0):
   fid.write('%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n' % (input_layer_size,second_layer_size,output_layer_size,kernel_size,train_history.history['loss'][number_of_epochs-1],train_history.history['acc'][number_of_epochs-1],score[0],score[1],total_avg_error_dl/len(y_predict),total_avg_error_ls/len(y_predict)))
else:
   fid.write('%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n' % (input_layer_size,output_layer_size,kernel_size,train_history.history['loss'][number_of_epochs-1],train_history.history['acc'][number_of_epochs-1],score[0],score[1],total_avg_error_dl/len(y_predict),total_avg_error_ls/len(y_predict)))

# Close file.
fid.close()

