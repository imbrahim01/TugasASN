import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from streamlit_option_menu import option_menu
import math
import streamlit as st 
from streamlit_lottie import st_lottie
import requests
from st_click_detector import click_detector

########
column_names = ['ECG']
data=pd.read_csv('ECG5minutes.txt',delimiter="\t", names=column_names)
data["sample interval"] = np.arange(len(data))
data["elapsed time"] = (data["sample interval"])*(1/200)
x=data["elapsed time"]
y=data["ECG" ] - (sum(data["ECG" ]/len(data["ECG"]))) #agar turun ke baseline

fs=int(round(1/(data.iloc[1,2]-data.iloc[0,2])))
jumlahdata = int(np.size(x))

#LPF
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

#HPF
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

#DERIVATIVE
drv=np.zeros(np.size(hpf_ecg))
for n in range (np.size(hpf_ecg)-2):
   drv[n]= (1/8)*(-(hpf_ecg[n-2]) - (2*hpf_ecg[n-1]) + (2*hpf_ecg[n+1]) + (hpf_ecg[n+2]))

# SQUARING PROCEDURE METHOD
sqr=np. zeros(np.size(drv) )
for n in range (np.size(drv)):
  sqr[n]=(drv[n])**2
# MAV PROCEDURE METHOD
w = 10 
mav = np.zeros(np.size(sqr))
for n in range(np.size(sqr)):
    for i in range(w):
        mav[n] = mav[n] + sqr[n - i]
    mav[n] = mav[n] / w
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
      
#NUMBERS OF R TO R CALCULATIONS
ptp = 0
waktu = np.zeros(np.size(thrqrs))
selisih = np.zeros(np.size(thrqrs))

for n in range(np.size(thrqrs) - 1):
    if thrqrs[n] < thrqrs[n + 1]:
        waktu[ptp] = n / fs;
        selisih[ptp] = waktu[ptp] - waktu[ptp - 1]
        ptp += 1
ptp = ptp - 1

#CALCULATION OF THE AMOUNT OF R
j = 0
peak = np.zeros(np.size(thrqrs))
for n in range(np.size(thrqrs)-1):
    if thrqrs[n] == 1 and thrqrs[n-1] == 0:
        peak[j] = n
        j += 1



#BPM CALCULATIONS':
temp = 0
interval = np.zeros(np.size(thrqrs))
BPM = np.zeros(np.size(thrqrs))

for n in range(ptp):
    interval[n] = (peak[n] - peak[n-1]) * (1/fs)
    BPM[n] = 60 / interval[n]
    temp = temp+BPM[n]
    rata = temp / (n - 1)
    
#TIME DOMAIN
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

n = np. arange(0,ptp,1,dtype=int)




with st.sidebar:
    selected = option_menu("TUGAS 1", ["Home","Ecyclopedia", "Data & Graphic", "Filter","Method","HRV Analysis"], icons=['house', 'book',None,None,None,None], menu_icon="cast", default_index=1)

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
   st_lottie("https://lottie.host/50914bbf-8aa3-4ac1-8ab7-d7d7882d51d5/QVzC4aV82R.json", height=400, width=400)

if selected == "Ecyclopedia":
     # Main title
    st.markdown("<h1 style='text-align: center; color: red;'>🫀ECYCLOPEDIA</h1>", unsafe_allow_html=True)
     # Subtitle
    new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">1. Apa yang dimaksud HRV?</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: justify;">HRV secara sederhana adalah ukuran variasi waktu antara setiap detak jantung. Variasi ini dikendalikan oleh bagian primitif dari sistem saraf yang disebut sistem saraf otonom (ANS). Sistem ini bekerja di belakang layar, secara otomatis mengatur detak jantung, tekanan darah, pernapasan, dan pencernaan di antara tugas-tugas utama lainnya. ANS dibagi lagi menjadi dua komponen besar: sistem saraf simpatis dan parasimpatis, yang juga dikenal sebagai mekanisme fight-or-flight dan respons relaksasi.</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">2. Bagaimana cara kerja HRV?</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: Justify;">Jantung Anda berdetak dengan kecepatan tertentu setiap saat. Denyut tersebut berubah tergantung pada apa yang sedang Anda lakukan saat itu. Denyut jantung yang lebih lambat terjadi ketika Anda sedang beristirahat atau santai, dan denyut yang lebih cepat terjadi ketika Anda sedang aktif, stres, atau ketika Anda dalam bahaya. Terdapat variabilitas dalam detak jantung Anda berdasarkan kebutuhan tubuh dan pola pernapasan Anda. Obat-obatan tertentu dan perangkat medis - seperti alat pacu jantung - juga dapat memengaruhi variabilitas detak jantung Anda. Variabilitas detak jantung Anda juga cenderung menurun secara normal seiring bertambahnya usia.</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    
  

# HTML content with the new YouTube video embedded
    content = """

     <iframe id='Video 1' width='560' height='315' src='https://www.youtube.com/embed/MUhtAXPvVnE' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe>
     """

# Display the HTML content
    st.markdown(content, unsafe_allow_html=True)
    st.link_button("Go to video", "https://youtu.be/MUhtAXPvVnE?si=rvYo04B8FCIcPT3I")

 




if selected == "Data & Graphic":
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
    new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Nilai FS</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    st.write(fs)
    new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Jumlah Semua Data</p>'
    st.markdown(new_title, unsafe_allow_html=True)
    st.write(jumlahdata)
if selected == "Filter":
     st.header("LPF")

     fig_LPF = go.Figure(data=go.Scatter(x=x[0:2000], y=lpf_ecg[0:1000], mode='lines'))
     fig_LPF.update_layout(
        title="LPF",
        xaxis_title="Sequence (n)",
        yaxis_title="Amplitude",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True)

     )
     st.plotly_chart(fig_LPF)
    
     st.header("HPF")
     fig_HPF = go.Figure(data=go.Scatter(x=x[0:2000], y=hpf_ecg[0:1000], mode='lines'))
     fig_HPF.update_layout(
        title="HPF",
        xaxis_title="Sequence (n)",
        yaxis_title="Amplitude",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True)
     )
     st.plotly_chart(fig_HPF)
if selected == "Method":

    optimizer_options = ['','Derivative', 'Squaring',"Moving Average",'Thresholding']
    selected_optimizer = st.selectbox('pilih metode', optimizer_options)
    
    if selected_optimizer == 'Derivative':
        
        fig_DRV = go.Figure(data=go.Scatter(x=x[9:1000], y=drv[0:1000], mode='lines'))
        fig_DRV.update_layout(
            title="DERIVATIVE",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
        )
        st.header("DERIVATIVE")
        st.plotly_chart(fig_DRV)
    elif selected_optimizer == 'Squaring':

        fig_sqr = go.Figure(data=go.Scatter(x=x[0:1000], y=sqr[0:1000], mode='lines'))
        fig_sqr.update_layout(
            title="SQUARING",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
        )
        st.header("SQUARING")
        st.plotly_chart(fig_sqr)

    elif selected_optimizer == 'Moving Average':
        fig_mav = go.Figure(data=go.Scatter(x=x[0:1000], y=mav[0:1000], mode='lines'))
        fig_mav.update_layout(
            title="MAV",
            xaxis_title="Time",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
        )
        st.header("MAV")
        st.plotly_chart(fig_mav)
    elif selected_optimizer == 'Thresholding':
        fig = go.Figure(data=go.Scatter(x=x[0:4000], y=y[0:4000], mode='lines'))
        fig.update_layout(
        title="RAW SIGNAL",
        xaxis_title="Elapsed Time",
        yaxis_title="Amplitude (mV)",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True)
        )
        st.subheader("THRESHOLDING")
        st.plotly_chart(fig)
    
        fig = go.Figure(data=go.Scatter(x=x[0:4000], y=thrqrs[0:4000], mode='lines'))
        fig.update_layout(
        title="SIGNAL THRESHOLD",
        xaxis_title="Time",
        yaxis_title="Amplitude",
        xaxis=dict(showline=True,showgrid=True),
        yaxis=dict(showline=True,showgrid=True)  
         )
        st.plotly_chart(fig)
        
    optimizer_options = ['NUMBERS OF R TO R CALCULATIONS', 'CALCULATION OF THE AMOUNT OF R',"BPM CALCULATIONS"]
    selected_optimizer = st.selectbox('Calculation of HR', optimizer_options)
                                      
    if selected_optimizer == 'NUMBERS OF R TO R CALCULATIONS':
        st.write(ptp)
    elif selected_optimizer == 'CALCULATION OF THE AMOUNT OF R':
        st.write(j)
    elif selected_optimizer == 'BPM CALCULATIONS':
        st.write(rata)
        
if selected == "HRV Analysis":
    sub_selected = st.sidebar.radio(
        "Pilih Metode HRV Analysis",
        ["Time Domain Analysis", "Risiko Penyakit Jantung", "Gejala Penyakit Jantung"],
        index=0
    )

    if sub_selected == 'Time Domain Analysis':
        optimizer_options1 = ['SDNN', 'RMSSD', "pNN50", "SDSD"]
        selected_optimizer1 = st.selectbox('Time-domain analysis', optimizer_options1)

        if selected_optimizer1 == 'SDNN':
            st.write(SDNN)
        elif selected_optimizer1 == 'RMSSD':
            st.write(RMSSD)
        elif selected_optimizer1 == 'pNN50':
            st.write(pNN50)
        elif selected_optimizer1 == 'SDSD':
            st.write(SDSD)

        fig_Tachogram = go.Figure(data=go.Scatter(x=n, y=bpm_rr, mode='lines'))
        fig_Tachogram.update_layout(
            title="TACHOGRAM",
            xaxis_title="n",
            yaxis_title="BPM",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
        )
        st.plotly_chart(fig_Tachogram)

        fig_histogram = go.Figure(data=go.Histogram(x=bpm_rr, nbinsx=ptp))

        fig_histogram.update_layout(
            title="Histogram Interval RR",
            xaxis_title="Interval RR",
            yaxis_title="Banyak Data",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            bargap=0.2,  # Optional: Adjusts the gap between bars
            bargroupgap=0.1,  # Optional: Adjusts the gap between groups
        )

        st.plotly_chart(fig_histogram)



    


         

        







     
