import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD

# Generate dummy data
import numpy as np
import scipy.io as sio

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
number_of_epochs = 40

# Load training and test vectors.
data_set = sio.loadmat('data_set_M_70_K_10_SNR_10_static_scenario_1_normv3.mat')
x_train = data_set['train_data']
y_train = data_set['train_label']
x_test = data_set['test_data']
y_test = data_set['test_label']
x_predict = data_set['prediction_data']
y_predict = data_set['prediction_label']

# Create Model.
model = Sequential()
model.add(Dense(M*N*2*2, activation='relu', input_dim=M*N*2))
#model.add(Dropout(0.5))
model.add(Dense(M*N*2, activation='relu'))
#model.add(Dropout(0.5))
model.add(Dense(M*2, activation='linear'))

# Set optimizer parameters.
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
adam = keras.optimizers.Adam()

# Compile model.
model.compile(loss='mean_squared_error', optimizer=sgd, metrics=['accuracy'])

# Train model.
model.fit(x_train, y_train, epochs=number_of_epochs, batch_size=batch_size)

# Test model.
score = model.evaluate(x_test, y_test, batch_size=batch_size)

print('Test loss:', score[0])
print('Test accuracy:', score[1])

# Calculate predictions
predictions = model.predict(x_predict)

for i in range(len(y_predict)):
   for j in range(2*M):
      print('expected: %f - actual: %f' % (y_predict[i][j],predictions[i][j]))


