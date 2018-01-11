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

# Load training and test vectors.

data_set = sio.loadmat('data_set_M_70_K_10_SNR_10_static_scenario_1.mat')
x_train = data_set['train_data']
y_train = data_set['train_label']
x_test = data_set['test_data']
y_test = data_set['test_label']

model = Sequential()
# Dense(M*N*2*4) is a fully-connected layer with M*N*2*4 hidden units.
# in the first layer, you must specify the expected input data shape:
# here, M*N*2-dimensional vectors.
model.add(Dense(M*N*2*2*2, activation='tanh', input_dim=M*N*2))
model.add(Dense(M*N*2*2, activation='tanh'))
model.add(Dense(M*2, activation='linear'))

# Set optimizer parameters.
sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)

# Compile model.
model.compile(loss='mean_squared_error', optimizer=sgd, metrics=['accuracy'])

# Train model.
model.fit(x_train, y_train, epochs=100, batch_size=100)

# Test model.
score = model.evaluate(x_test, y_test, batch_size=100)

print('Test loss:', score[0])
print('Test accuracy:', score[1])

# Calculate predictions
#predictions = model.predict(X)
