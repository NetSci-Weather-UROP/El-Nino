import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

with open('bulk_temp_1949_2016.npy', 'rb') as f:
    results = np.load(f)
    lat = np.load(f)
    lon = np.load(f)

results = results[1949:2017, :, :, :]
print('Shape of distilled results: ', np.shape(results))

n_y = np.count_nonzero(results[:, :, 1, :], axis = 2)
print('Shape of N_y array: ', np.shape(n_y))
n_y_flat = n_y.flatten()
print('Shape of flattened N_y array: ', np.shape(n_y_flat))
n_y_moving_average = moving_average(n_y_flat, 3)
print('Shape of moving average (3-month) N_y array: ', 
       np.shape(n_y_moving_average))
n_y_moving_average_adjusted = np.concatenate((n_y_flat[0:2], 
                                              n_y_moving_average))
print('Shape of adjusted moving average N_y array: ',
       np.shape(n_y_moving_average_adjusted))

c_theta = np.sum(abs(results[:, :, 0, :]), axis = 2)
print('Shape of C_theta array: ', np.shape(c_theta))
c_y = c_theta / n_y
print('Shape of C_y array: ', np.shape(c_y))
c_y_flat = c_y.flatten()
print('Shape of flattened C_y array: ', np.shape(c_y_flat))
c_y_moving_average = moving_average(c_y_flat, 3)
print('Shape of moving average (3-month) C_y array: ',
       np.shape(c_y_moving_average))
c_y_moving_average_adjusted = np.concatenate((c_y_flat[0:2],
                                              c_y_moving_average))
print('Shape of adjusted moving average C_y array: ',
       np.shape(c_y_moving_average_adjusted))

year_labels = np.empty([68, 12])
years = np.arange(1949, 2017)
for i in range(12):
    year_labels[:, i] = years + i/12
year_labels = year_labels.flatten()

fig = plt.figure()

gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
ax0 = plt.subplot(gs[0])
ax0.set_ylabel('$N^y$', fontsize = 16)
plot0, = ax0.plot(year_labels, n_y_moving_average_adjusted, c = 'purple')
plt.xticks(np.arange(1950, 2016, 1))


plt.xlim([1950, 2015])
ax1 = plt.subplot(gs[1], sharex = ax0)
ax1.set_ylabel('$C^y$', fontsize = 16)
plot1, = ax1.plot(year_labels, c_y_moving_average_adjusted, c = 'navy')
plt.setp(ax0.get_xticklabels(), visible=False)
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
plt.subplots_adjust(hspace=.0)

#ax0.axvspan(1971.8, 1972.2, alpha=0.5, color='red')
#ax1.axvspan(1971.8, 1972.2, alpha=0.5, color='red')

ax0.xaxis.grid(True)
ax1.xaxis.grid(True)
ax1.set_xlabel('Year', fontsize = 16)
plt.xticks(rotation=45)
n = 5  # Keeps every 7th label
[l.set_visible(False) for (i,l) in enumerate(ax1.xaxis.get_ticklabels()) if i % n != 0]
fig.align_ylabels([ax0, ax1])
fig.set_size_inches(16, 6)
fig.savefig('test2png.png', dpi=60)
plt.show()
