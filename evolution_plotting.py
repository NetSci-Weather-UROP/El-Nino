import numpy as np

with open('bulk_temp.npy', 'rb') as f:
    results = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

plot_data = results[1959:1972, :, :, :]
c_average = np.mean(results, axis = 3)[1959:1972, :, 0]
print(np.shape(c_avrage))