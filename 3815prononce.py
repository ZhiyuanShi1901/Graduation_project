import os
import obspy
import numpy as np
import sys
import math

import matplotlib.pyplot as plt
sys.path.append(r'G:\Pythonworks\graduate\noise')
from read_mseed_myself import read_mseed_myself_waves

# %matplotlib notebook
from scipy.fftpack import fft, fftshift  # 进行fft
from matplotlib.animation import FuncAnimation
import matplotlib.animation

path = 'changednameshallowfile'

x_axis_data = []
y_axis_data = []

for i in range(1, 181):

    filepath = path + '/' + f'{i}' +'.SAC'
    seismogram = obspy.read(filepath)
    trace = seismogram[0]
    print('read path successfully')
    # 去均值
    trace.detrend('demean')
    # 去线性趋势
    trace.detrend('linear')
    # 波形尖灭
    trace.taper(max_percentage=0.05, max_length=5.)
    # trace.spectrogram()

    # 数据截取
    trace.data = trace.data[0:3600*6]
    trace.stats.npts = len(trace.data)

    fftresult = fft(trace.data)

    # print(type(fftresult))
    # print(fftresult)
    N = len(trace.data)

    fft_amp0 = np.array(np.abs(fftresult) / N * 2)  # 计算双边谱
    N_2 = int(N/2) # 单边谱数量
    fft_amp1 = fft_amp0[0:N_2] # 计算单边谱
    fft_amp1 = np.concatenate((fft_amp1, [i]))  # 先将p_变成list形式进行拼接，注意输入为一个tuple:

    list1 = np.array(range(0, int(N/2)))
    sample_freq = trace.stats.sampling_rate
    freq1 = sample_freq*list1/N        # 单边谱的频率轴
    freq1 = np.concatenate((freq1, [1]))
    x_axis_data.append(freq1)
    y_axis_data.append(fft_amp1)
    print(trace.stats.station)


number = len(freq1)
fig,ax = plt.subplots(num=1)
line = plt.plot([])
plt.xlim(0.00381, 0.00390)
leg = [i for i in range(1, 181)]
def update(yaxis):
    ax.clear()
    ax.plot(freq1, yaxis)
    plt.xlim(0.003, 0.004)
    plt.ylim(0, 4)
    plt.legend([f'Epicentral distance = {yaxis[number-1]} degree'], bbox_to_anchor=(1,1))
    plt.title('3.815mHz Spectrum.')
    return line

print(leg)
ani = FuncAnimation(fig, update, frames=y_axis_data, interval=200)
ani.save('3_815mHz.gif', writer='imagemagick')
plt.show()





