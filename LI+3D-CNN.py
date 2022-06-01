'''
LI + 3D Convolutional Neural Network
'''

# To run this code, a listed modules are required
from datetime import datetime
import tensorflow as tf
import pickle
import numpy as np
import pandas as pd
from scipy import interpolate

# Flag for module import
print(datetime.now(), "Successfully imported modules")

# Parameters
# user = "ap2021"
user = "rpe26"
filename_0 = "/home/" + user + "/rds/hpc-work/channel_1_1-100.pkl"
filename_1 = "/home/" + user + "/rds/hpc-work/channel_2_1-100.pkl"
act = 'relu'
snapshots = 200
nx = 256
ny = 128
nz = 160
EPOCHS = 2000
BATCH_SIZE = 20
VAL_SPLIT = 0.2
n_slices = 5  # Only supports 3, 5 and 7

##log which ops are executed on which devices:
physical_devices = tf.config.list_physical_devices('GPU')
# Create a MirroredStrategy.
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
strategy = tf.distribute.MirroredStrategy()


# Configurations of input/output variables are as follows:
#'''
#3D_field: time, nx =256, ny=128, nz=160, components=3
#...

with open(filename_0, 'rb') as f:
    obj = pickle.load(f)
    uvw3D_field = obj
# Flag for first file
print(datetime.now(), "Loaded data from pickle file")
    
with open(filename_1, 'rb') as f:
    obj = pickle.load(f)
    uvw3D_field = np.concatenate((uvw3D_field, obj), axis=0)
# Flag for second file
print(datetime.now(), "Loaded data from second pickle file")


slice_pos = np.zeros(n_slices)
slice_distance = 0.8/ (n_slices - 1)
for i in range(n_slices):
    slice_pos[i] = (nz - 1) * (0.5 + slice_distance * (i - (n_slices - 1) / 2))
slice_pos = slice_pos.astype(int)


training_field_linear = uvw3D_field.copy()
x = slice_pos
x_new = np.arange(0, nz)
y = training_field_linear[:,:,:,x,:]
f = interpolate.interp1d(x, y, kind='linear', axis= 3, bounds_error=False, fill_value=(y[:,:,:,0,:],y[:,:,:,-1,:]))   
training_field_linear[:,:,:,x_new,:] = f(x_new)


# Flag for dataset shape
print(datetime.now(), "Shape of 3D Dataset", uvw3D_field.shape)
print(datetime.now(), "Shape of 3D input dataset", training_field_linear.shape)

### train_test_split

from sklearn.model_selection import train_test_split
print(datetime.now(), "Imported train_test_split")
training_field_linear, training_field_linear_val, uvw3D_field, uvw3D_field_val = train_test_split(training_field_linear, uvw3D_field, test_size=VAL_SPLIT, shuffle=True)
print(datetime.now(), "Data split", VAL_SPLIT)

  
with strategy.scope():
    # Input variables
    input_field = tf.keras.layers.Input(shape=(nx, ny, nz, 3))

    # Network structure
    x = tf.keras.layers.Conv3D(32, (3, 3, 3), activation=act, padding='same')(input_field)
    x = tf.keras.layers.Conv3D(32, (3, 3, 3), activation=act, padding='same')(x)
    x = tf.keras.layers.MaxPooling3D((2, 2, 4))(x) ##(2, 2, 8)
    x = tf.keras.layers.Conv3D(16, (3, 3, 3), activation=act, padding='same')(x)
    x = tf.keras.layers.Conv3D(16, (3, 3, 3), activation=act, padding='same')(x)
    x = tf.keras.layers.UpSampling3D((2, 2, 4))(x) ##(2, 2, 8)
    x = tf.keras.layers.Conv3D(32, (3, 3, 3), activation=act, padding='same')(x)
    x = tf.keras.layers.Conv3D(32, (3, 3, 3), activation=act, padding='same')(x)
    x_final = tf.keras.layers.Conv3D(3, (3, 3, 3), activation='linear', padding='same')(x)
    # -------- #
    ##model = tf.keras.models.Model(input_field, x_final)
    model = tf.keras.models.Model(input_field, x_final)
    model.compile(optimizer='adam', loss='mse')

# Flag for compiling model
print(datetime.now(), "NN Model compiled")

##################

model_cb = tf.keras.callbacks.ModelCheckpoint('/home/' + user + '/rds/hpc-work/Test_Model_Checkpoint_2.hdf5',monitor='val_loss', save_best_only=True, verbose=1)
early_cb = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=20, verbose=1)
cb = [model_cb, early_cb]
history = model.fit(training_field_linear, uvw3D_field, epochs=EPOCHS, batch_size=BATCH_SIZE, verbose=1, callbacks=cb, shuffle=True, validation_data=(training_field_linear_val,uvw3D_field_val))
df_results = pd.DataFrame(history.history)
df_results['epoch'] = history.epoch
df_results.to_csv(path_or_buf='/home/' + user + '/rds/hpc-work/Test_Model_Results_2.csv', index=False)
    
model.save("/home/" + user + "/rds/hpc-work/Test_Model_2")

# Flag for model.save
print(datetime.now(), "Successfully saved trained model")