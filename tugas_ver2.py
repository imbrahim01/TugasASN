import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from streamlit_option_menu import option_menu
import math
import streamlit as st 

column_names = ['ECG']
data=pd.read_csv('ECG5minutes.txt',delimiter="\t", names=column_names)
data["sample interval"] = np.arange(len(data))
data["elapsed time"] = (data["sample interval"])*(1/200)
x=data["elapsed time"]
y=data["ECG" ] - (sum(data["ECG" ]/len(data["ECG"]))) #agar turun ke baseline

fs=int(round(1/(data.iloc[1,2]-data.iloc[0,2])))
jumlahdata = int(np.size(x))

#fc_lpf =input("FREQUENCY CUT-OFF FOR LOWPASS FILTER :")
fc_lpf = 11
fc_lpf=float(fc_lpf)

lpf_ecg = np.zeros(jumlahdata) 
for n in range(3):
    lpf_ecg[-n] = lpf_ecg[0]
    y[-n] = y[0]

# Coefficients of LPF
T = 1 / fs
w = 2 * math.pi * fc_lpf
a0 = w**2
a1 = 2*(w**2)
b1 = (8/(T**2)) - (2*(w**2))
C0 = ((4/(T**2)) - ((2 * (math.sqrt(2))*w)/T) + (w**2))
C1 = ((4/(T**2)) + ((2 * (math.sqrt(2))*w)/T) + (w**2))

# BUTTERWORTH LOWPASS FILTER EQUATION
for n in range(jumlahdata):
    lpf_ecg[n] = ((b1 * lpf_ecg[n-1]) - (C0 * lpf_ecg[n-2]) + (a0 * y[n]) + (a1 * y[n-1]) + (a0 * y[n-2])) / C1

#fc_hpf =input("FREQUENCY CUT-OFF FOR HIGHPASS FILTER :")
fc_hpf = 5
fc_hpf=float(fc_hpf)

hpf_ecg = np.zeros(np.size(lpf_ecg))
for n in range(3):
    hpf_ecg[-n] = hpf_ecg[0]

# Coefficients of LPF
T = 1 / fs
w = 2 * math.pi * fc_hpf
e0 = 4*T
e1 = 8*T
e2 = 4*T
d0 = ((2*(w**2)*(T**2))-8)
d1 = (((w**2)*(T**2)) - (2 * (math.sqrt(2)) * T * w)+4)
d2 = (((w**2)*(T**2)) + (2 * (math.sqrt(2)) * T* w)+4)
# BUTTERWORTH LOWPASS FILTER EQUATION
for n in range(np.size(lpf_ecg)):
    hpf_ecg[n] = ((e0 * lpf_ecg[n]) - (e1 * lpf_ecg[n-1]) +(e2 * lpf_ecg[n-2])-(d0 * hpf_ecg[n-1])- (d1 * hpf_ecg[n-2]))/ d2

drv=np.zeros(np.size(hpf_ecg))
for n in range (np.size(hpf_ecg)-2):
   drv[n]= (1/8)*(-(hpf_ecg[n-2]) - (2*hpf_ecg[n-1]) + (2*hpf_ecg[n+1]) + (hpf_ecg[n+2]))

# SQUARING PROCEDURE METHOD
sqr=np. zeros(np.size(drv) )
for n in range (np.size(drv)):
  sqr[n]=(drv[n])**2

w=10
mav=np. zeros (np. size(sqr))
for n in range (np.size(sqr)) :
  for i in range (w):
   mav[n]=mav[n]+sqr [n-i]
mav [n]=mav[n]/w

tinggi=0
tinggi=np.zeros(np.size(mav))
for n in range (np. size(mav) -1): 
    if (tinggi < mav[n]) .all():
      tinggi [n]=mav[n]

thr=tinggi*0.5
thrqrs=np.zeros(np.size(mav))
for n in range (np. size(mav)-1):
  if (mav[n] >= thr).all():
     thrqrs [n]=1
  elif (mav[n]<thr).all():
    thrqrs [n]=0

ptp = 0
waktu = np.zeros(np.size(thrqrs))
selisih = np.zeros(np.size(thrqrs))

for n in range(np.size(thrqrs) - 1):
    if thrqrs[n] < thrqrs[n + 1]:
        waktu[ptp] = n / fs;
        selisih[ptp] = waktu[ptp] - waktu[ptp - 1]
        ptp += 1

ptp = ptp - 1

j = 0
peak = np.zeros(np.size(thrqrs))
for n in range(np.size(thrqrs)-1):
    if thrqrs[n] == 1 and thrqrs[n-1] == 0:
        peak[j] = n
        j += 1

temp = 0
interval = np.zeros(np.size(thrqrs))
BPM = np.zeros(np.size(thrqrs))

for n in range(ptp):
    interval[n] = (peak[n] - peak[n-1]) * (1/fs)
    BPM[n] = 60 / interval[n]
    temp = temp+BPM[n]
    rata = temp / (n - 1)

RR_SDNN=0
for n in range (ptp):
   RR_SDNN += (((selisih[n])-(60/rata))**2)

SDNN = math.sqrt (RR_SDNN/ (ptp-1))

RR_RMSSD=0
for n in range (ptp):
   RR_RMSSD += ((selisih[n+1]-selisih[n])**2)
RMSSD =  math. sqrt (RR_RMSSD/(ptp-1))

# FIND NN50 ALGORITHM
NN50 = 0

for n in range (ptp): 
    if (abs(selisih[n+1]-selisih[n])>0.05):
      NN50 +=1
pNN50 = (NN50/ (ptp-1)) *100 

dif = 0
for n in range (ptp):
  dif += abs(selisih[n]-selisih[n+1])
RRdif = dif/(ptp-1)

RR_SDSD = 0
for n in range (ptp):
  RR_SDSD += (((abs(selisih[n]-selisih[n+1]))-RRdif)**2)
SDSD = math.sqrt(RR_SDSD/(ptp-2))

bpm_rr = np.zeros(ptp)
for n in range (ptp):
  bpm_rr[n] = 60/selisih[n]
  if bpm_rr [n]>100:
    bpm_rr[n]=rata
bpm_rr
n = np. arange(0,ptp,1,dtype=int)


with st.sidebar:
    selected = option_menu("TUGAS 1", ["PAGE 1", "PAGE 2", "PAGE 3","PAGE 4","PAGE 5"], default_index=0)

if selected == "PAGE 1":
    st.title('Data & Graphic Input')
    st.header("Data Input")
    st.write(data)

    fig_data = go.Figure(data=go.Scatter(x=x[0:2000], y=y[0:2000], mode='lines'))
    fig_data.update_layout(
        title="Time Signal",
        xaxis_title="Elapsed Time",
        yaxis_title="Amplitude (mV)",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True)
    )
    st.header("Graphic Input")
    st.plotly_chart(fig_data)
