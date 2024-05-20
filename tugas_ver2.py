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

fc_lpf = 13
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




with st.sidebar:
    selected = option_menu("TUGAS 1", ["Home", "PAGE 1", "PAGE 2","PAGE 3","PAGE 4"], default_index=0)

if selected == "Home":
   st.title('Project ASN Kelompok 1')
   st.subheader("Anggota kelompok")
   new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">Farhan Majid Ibrahim - 5023211049</p>'
   st.markdown(new_title, unsafe_allow_html=True)
   new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">Nayla Pramudhita Putri Pertama - 5023211012</p>'
   st.markdown(new_title, unsafe_allow_html=True)
   new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">Mohammad Rayhan Amirul Haq Siregar - 5023211045</p>'
   st.markdown(new_title, unsafe_allow_html=True)
   new_title = '<p style="font-family:Georgia; color: black; font-size: 15px;">Reynard Prastya Savero - 5023211042</p>'
   st.markdown(new_title, unsafe_allow_html=True)


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
if selected == "PAGE 2":
     new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Nilai FS</p>'
     st.markdown(new_title, unsafe_allow_html=True)
     st.write(fs)
     new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Jumlah Semua Data</p>'
     st.markdown(new_title, unsafe_allow_html=True)
     st.write(jumlahdata)
     st.header("LPF")
     fc_lpf = st.number_input("FREQUENCY CUT-OFF FOR LOWPASS FILTER :", value=0)
    
     fig_LPF = go.Figure(data=go.Scatter(x=x[0:2000], y=lpf_ecg[0:1000], mode='lines'))
     fig_LPF.update_layout(
        title="LPF",
        xaxis_title="Sequence (n)",
        yaxis_title="Amplitude",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True)

     )
     st.plotly_chart(fig_LPF)

     
