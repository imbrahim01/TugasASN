import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objs as go
from streamlit_option_menu import option_menu
import math
import streamlit as st 
from streamlit_lottie import st_lottie
from st_click_detector import click_detector
import plotly.express as px


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

def fourier_transform(signal):
    N = len(signal)
    fft_result = np.zeros(N, dtype=complex)
    for k in range(N):
        for n in range(N):
            fft_result[k] += signal[n] * np.exp(-2j * np.pi * k * n / N)
    return fft_result

def calculate_frequency(N, sampling_rate):
    return np.arange(N) * sampling_rate / N
    
bpm_rr_baseline = bpm_rr - 70
# Ambil subset data dari 0 sampai 49
n_subset = n[0:50]
bpm_rr_baseline_subset = bpm_rr_baseline[0:50]
M = len(bpm_rr_baseline_subset) - 1
hamming_window = np.zeros(M+1)
for i in range(M+1):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i / M)
bpm_rr_baseline_windowed = bpm_rr_baseline_subset * hamming_window
fft_result = fourier_transform(bpm_rr_baseline_windowed)
sampling_rate = 1
fft_freq = calculate_frequency(len(bpm_rr_baseline_windowed), sampling_rate)
half_point = len(fft_freq) // 2
fft_freq_half = fft_freq[:half_point]
fft_result_half = fft_result[:half_point]

# Ambil subset data dari 50 sampai 100
n_subset1 = n[50:100]
bpm_rr_baseline_subset1 = bpm_rr_baseline[50:100]
M = len(bpm_rr_baseline_subset1) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed1 = bpm_rr_baseline_subset1 * hamming_window
fft_result1 = fourier_transform(bpm_rr_baseline_windowed1)
fft_freq1 = calculate_frequency(len(bpm_rr_baseline_windowed1), sampling_rate)
half_point1 = len(fft_freq1) // 2
fft_freq_half1 = fft_freq1[:half_point1]
fft_result_half1 = fft_result1[:half_point1]

# Ambil subset data dari 101 sampai 151
n_subset2 = n[101:151]
bpm_rr_baseline_subset2 = bpm_rr_baseline[101:151]
M = len(bpm_rr_baseline_subset2) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed2 = bpm_rr_baseline_subset2 * hamming_window
fft_result2 = fourier_transform(bpm_rr_baseline_windowed2)
fft_freq2 = calculate_frequency(len(bpm_rr_baseline_windowed2), sampling_rate)
half_point2 = len(fft_freq2) // 2
fft_freq_half2 = fft_freq2[:half_point2]
fft_result_half2 = fft_result2[:half_point2]

# Ambil subset data dari 151 sampai 201
n_subset3 = n[151:201]
bpm_rr_baseline_subset3 = bpm_rr_baseline[151:201]
M = len(bpm_rr_baseline_subset3) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed3 = bpm_rr_baseline_subset3 * hamming_window
fft_result3 = fourier_transform(bpm_rr_baseline_windowed3)
fft_freq3 = calculate_frequency(len(bpm_rr_baseline_windowed3), sampling_rate)
half_point3 = len(fft_freq3) // 2
fft_freq_half3 = fft_freq3[:half_point3]
fft_result_half3 = fft_result3[:half_point3]

# Ambil subset data dari 201 sampai 251
n_subset4 = n[201:251]
bpm_rr_baseline_subset4 = bpm_rr_baseline[201:251]
M = len(bpm_rr_baseline_subset4) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed4 = bpm_rr_baseline_subset4 * hamming_window
fft_result4 = fourier_transform(bpm_rr_baseline_windowed4)
fft_freq4 = calculate_frequency(len(bpm_rr_baseline_windowed4), sampling_rate)
half_point4 = len(fft_freq4) // 2
fft_freq_half4 = fft_freq4[:half_point4]
fft_result_half4 = fft_result4[:half_point4]

# Ambil subset data dari 251 sampai 301
n_subset5 = n[251:301]
bpm_rr_baseline_subset5 = bpm_rr_baseline[251:301]
M = len(bpm_rr_baseline_subset5) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed5 = bpm_rr_baseline_subset5 * hamming_window
fft_result5 = fourier_transform(bpm_rr_baseline_windowed5)
fft_freq5 = calculate_frequency(len(bpm_rr_baseline_windowed5), sampling_rate)
half_point5 = len(fft_freq5) // 2
fft_freq_half5 = fft_freq5[:half_point5]
fft_result_half5 = fft_result5[:half_point5]

# Ambil subset data dari 301 sampai 351
n_subset6 = n[301:351]
bpm_rr_baseline_subset6 = bpm_rr_baseline[301:351]
M = len(bpm_rr_baseline_subset6) -1
hamming_window = np.zeros(M+1)
for i in range(M):
    hamming_window[i] = 0.54 - 0.46 * np.cos(2 * np.pi * i /M )
bpm_rr_baseline_windowed6 = bpm_rr_baseline_subset6 * hamming_window
fft_result6 = fourier_transform(bpm_rr_baseline_windowed6)
sampling_rate = 1
fft_freq6 = calculate_frequency(len(bpm_rr_baseline_windowed6), sampling_rate)
half_point6 = len(fft_freq6) // 2
fft_freq_half6 = fft_freq6[:half_point6]
fft_result_half6 = fft_result6[:half_point6]

# FFT Total dari 0 sampai 351
FFT_TOTAL = (fft_result + fft_result1 +  fft_result2 + fft_result3 + fft_result4 + fft_result5 + fft_result6) / 7
FFT_FREQ_TOTAL = (fft_freq + fft_freq1 + fft_freq2 +  fft_freq3 +  fft_freq4 +  fft_freq5 +  fft_freq6 )/ 7
half_point_total = len(FFT_FREQ_TOTAL) // 2
fft_freq_total = FFT_FREQ_TOTAL[:half_point_total]
fft_result_total = FFT_TOTAL[:half_point_total]

def manual_interpolation(x, xp, fp):
    return np.interp(x, xp, fp)
x_vlf = np.linspace(0.003, 0.04, 99)
x_lf = np.linspace(0.04, 0.15, 99)
x_hf = np.linspace(0.15, 0.4, 99)
y_vlf = manual_interpolation(x_vlf, fft_freq_total, np.abs(fft_result_total))
y_lf = manual_interpolation(x_lf, fft_freq_total, np.abs(fft_result_total))
y_hf = manual_interpolation(x_hf, fft_freq_total, np.abs(fft_result_total))


def trapezoidal_rule(y, x):
    return np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2)
TP = trapezoidal_rule(np.abs(fft_result_total), fft_freq_total)
VLF = trapezoidal_rule(y_vlf, x_vlf)
LF = trapezoidal_rule(y_lf, x_lf)
HF = trapezoidal_rule(y_hf, x_hf)

total_power = VLF + LF + HF
# Hitung LF dan HF yang dinormalisasi
LF_norm = LF / (total_power - VLF)
HF_norm = HF / (total_power - VLF)
LF_HF = LF / HF




def dirac(x): 
    if (x == 0) :
        dirac_delta = 1
    else:
        dirac_delta = 0
    
    return dirac_delta
    return result


h = []
g = []
n_list = []

for n in range(-2, 2):
    n_list.append(n)
    temp_h = 1/8 * (dirac(n-1) + 3*dirac(n) + 3*dirac(n+1) + dirac(n+2))
    h.append(temp_h)
    temp_g = -2 * (dirac(n) - dirac(n+1))
    g.append(temp_g)

# Hw = []
# Gw = []

Hw = np.zeros(20000)
Gw = np.zeros(20000)
i_list = []
for i in range (0,fs+1):
    i_list.append(i)
    reG = 0
    imG = 0
    reH = 0
    imH = 0
    for k in range(-2,2):
      reG = reG + g[k+abs(-2)]*np.cos(k*2*np.pi*i/fs)
      imG = imG - g[k+abs(-2)]*np.sin(k*2*np.pi*i/fs)
      reH = reH + h[k+abs(-2)]*np.cos(k*2*np.pi*i/fs)
      imH = imH - h[k+abs(-2)]*np.sin(k*2*np.pi*i/fs)
    temp_Hw = np.sqrt((reH**2)+(imH**2))
    temp_Gw = np.sqrt((reG**2)+(imG**2))


    Hw[i] = temp_Hw
    Gw[i] = temp_Gw
    
i_list = i_list[0:round(fs/2)+1]

Q = np.zeros((9, round(fs/2) + 1))

# Generate the i_list and fill Q with the desired values
i_list = []
for i in range(0, round(fs/2) + 1):
    i_list.append(i)
    Q[1][i] = Gw[i]
    Q[2][i] = Gw[2*i] * Hw[i]
    Q[3][i] = Gw[4*i] * Hw[2*i] * Hw[i]
    Q[4][i] = Gw[8*i] * Hw[4*i] * Hw[2*i] * Hw[i]
    Q[5][i] = Gw[16*i] * Hw[8*i] * Hw[4*i] * Hw[2*i] * Hw[i]
    Q[6][i] = Gw[32*i] * Hw[16*i] * Hw[8*i] * Hw[4*i] * Hw[2*i] * Hw[i]
    Q[7][i] = Gw[64*i] * Hw[32*i] * Hw[16*i] * Hw[8*i] * Hw[4*i] * Hw[2*i] * Hw[i]
    Q[8][i] = Gw[128*i] * Hw[64*i] * Hw[32*i] * Hw[16*i] * Hw[8*i] * Hw[4*i] * Hw[2*i] * Hw[i]

traces = []

T1= round (2**(1-1))-1
T2 = round(2** (2-1)) - 1
T3 = round(2** (3-1)) - 1
T4 = round(2**(4-1)) - 1
T5 = round(2**(5-1))- 1
Delay1= T5-T1
Delay2= T5-T2
Delay3= T5-T3
Delay4= T5-T4
Delay5= T5-T5

ecg=y

min_n = 0 * fs
max_n = 4 * fs 


def process_ecg(min_n, max_n, ecg, g, h):
    w2fm = np.zeros((5, max_n - min_n + 1))
    s2fm = np.zeros((5, max_n - min_n + 1))

    for n in range(min_n, max_n + 1):
        for j in range(1, 6):
            w2fm[j-1, n - min_n] = 0
            s2fm[j-1, n - min_n] = 0
            for k in range(-1, 3):
                index = round(n - 2**(j-1) * k)
                if 0 <= index < len(ecg):  # Ensure the index is within bounds
                    w2fm[j-1, n - min_n] += g[k+1] * ecg[index]  # g[k+1] to match Pascal's array index starting from -1
                    s2fm[j-1, n - min_n] += h[k+1] * ecg[index]  # h[k+1] to match Pascal's array index starting from -1

    return w2fm, s2fm

# Compute w2fm and s2fm
w2fm, s2fm = process_ecg(min_n, max_n, ecg, g, h)

# Prepare data for plotting
n_values = np.arange(min_n, max_n + 1)
w2fm_values = [w2fm[i, :] for i in range(5)]  # Equivalent to w2fm[1,n] to w2fm[5,n] in original code (0-based index)
s2fm_values = [s2fm[i, :] for i in range(5)]  




with st.sidebar:
    selected = option_menu("TUGAS 1", ["Home","Encyclopedia", "Signal Processing","HRV Analysis","DWT"], default_index=0)

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

if selected == "Encyclopedia":
     # Main title
    st.markdown("<h1 style='text-align: center; color: red;'>🫀ENCYCLOPEDIA</h1>", unsafe_allow_html=True)
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
if selected == "Signal Processing":
    selected1 = option_menu(None, ["Information","Data & Graphic", "Filter","Method & Calculation"], 
    menu_icon="cast", default_index=0, orientation="horizontal")
    
    if selected1 == 'Information':
        st.title("d")
    elif selected1 == 'Data & Graphic':
        st.title('Data & Graphic Input')
        st.header("Data Input")
        st.write(data)

        # Create the figure with Plotly
        fig = go.Figure(data=go.Scatter(x=x[0:2000], y=y[0:2000], mode='lines'))
        fig.update_layout(
            title="Original Signal",
            xaxis_title="Elapsed Time",
            yaxis_title="Amplitude (mV)",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            
        )
        
        # Display the figure in Streamlit
        st.header("Graphic Input")
        st.plotly_chart(fig) 
        
        
        new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Nilai FS</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        st.write(fs)
        new_title = '<p style="font-family:Georgia; color: black; font-size: 20px;">Jumlah Semua Data</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        st.write(jumlahdata)
    elif selected1 == 'Filter':
        st.header("LPF")

        fig_LPF = go.Figure(data=go.Scatter(x=x[0:2000], y=lpf_ecg[0:1000], mode='lines'))
        fig_LPF.update_layout(
            title="LPF",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'

         )
        st.plotly_chart(fig_LPF)
    
        st.header("HPF")
        fig_HPF = go.Figure(data=go.Scatter(x=x[0:2000], y=hpf_ecg[0:1000], mode='lines'))
        fig_HPF.update_layout(
            title="HPF",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
         )
        st.plotly_chart(fig_HPF)
    elif selected1 == 'Method & Calculation':
     optimizer_options = ['', 'Derivative', 'Squaring', 'Moving Average', 'Thresholding','Calculation']
     selected_optimizer = st.selectbox('Method & Calculation', optimizer_options)

     if selected_optimizer == 'Derivative':
        fig_DRV = go.Figure(data=go.Scatter(x=x[9:1000], y=drv[0:1000], mode='lines'))
        fig_DRV.update_layout(
            title="DERIVATIVE",
            xaxis_title="Sequence (n)",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
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
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
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
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
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
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.subheader("THRESHOLDING")
        st.plotly_chart(fig)

        fig = go.Figure(data=go.Scatter(x=x[0:4000], y=thrqrs[0:4000], mode='lines'))
        fig.update_layout(
            title="SIGNAL THRESHOLD",
            xaxis_title="Time",
            yaxis_title="Amplitude",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True),
            template='plotly_dark'
        )
        st.plotly_chart(fig)
     elif selected_optimizer == 'Calculation':
    # Define the data for the table
        data = {
            "Calculation of HR": ["NUMBERS OF R TO R CALCULATIONS", "CALCULATION OF THE AMOUNT OF R", "BPM CALCULATIONS"],
            "Hasil": [ptp, j, rata]
        }
        df = pd.DataFrame(data)

        # Create the table using Plotly
        fig = go.Figure(data=[go.Table(
            columnwidth=[80, 20],  # Set column width
            header=dict(values=list(df.columns),
                    fill_color='red',  # Change header color to red
                    align='left',
                    line_color='black',
                    height=30),  # Set header height
        cells=dict(values=[df["Calculation of HR"], df["Hasil"]],
                   fill_color='white',  # Change cell color to white
                   align='left',
                   line_color='black',
                   height=25,  # Set cell height
                   font_size=12,  # Set font size
                   ),
    )])

    # Set layout to adjust the table size
        fig.update_layout(
            width=700,
            height=200,
            margin=dict(l=10, r=10, t=10, b=10)
        )

    # Display the table
        st.plotly_chart(fig) 

if selected == "HRV Analysis":
    sub_selected = st.sidebar.radio(
        "Pilih Metode HRV Analysis",
        ["Time Domain Analysis", "Frequency Domain analysis", "Non Liniear analysis"],
        index=0
    )

    if sub_selected == 'Time Domain Analysis':
        new_title = '<p style="font-family:Georgia; color:black; font-size: 25px; text-align: center;">Time Domain Analysis</p>'
        st.markdown(new_title, unsafe_allow_html=True)
        selected2 = option_menu(None, ["Result", "Information"], 
            menu_icon="cast", default_index=0, orientation="horizontal")
        if selected2 == "Result":

            new_title = '<p style="font-family:Georgia; color: black; font-size: 18px;">Statistical measures</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            data = {
            "Statistical measures": ["SDNN", "RMSSD", "pNN50","SDSD"],
            "Hasil": [SDNN,RMSSD, pNN50, SDSD]
            }
            df = pd.DataFrame(data)

            fig = go.Figure(data=[go.Table(
            columnwidth=[80, 20],  # Set column width
            header=dict(values=list(df.columns),
                fill_color='red',  # Ubah warna header menjadi merah
                align='left',
                line_color='darkslategray',
                height=30),  # Set header height
            cells=dict(values=[df["Statistical measures"], df["Hasil"]],
               fill_color='white',  # Ubah warna sel menjadi merah
               align='left',
               line_color='darkslategray',
               height=25,  # Set cell height
               font_size=12,  # Set font size
               ),
            )])

       
            fig.update_layout(
             width=800,
            height=200,
            margin=dict(l=10, r=10, t=10, b=10)
            )

        
            st.plotly_chart(fig)


        if selected2 == "Information":
            new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">SDNN</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: justify;">SDNN adalah ukuran yang digunakan untuk mengukur variasi dalam interval RR. SDNN menghitung standar deviasi dari semua interval RR dalam suatu rekaman. SDNN digunakan sebagai indikator aktivitas sistem saraf otonom dan respons tubuh terhadap stres. SDNN juga digunakan untuk menganalisis bagaimana tubuh menanggapi stres dan bagaimana sistem saraf otonom bekerja untuk mengatur respons tubuh terhadap stres.</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">SDSD</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: Justify;">SDSD adalah ukuran yang digunakan untuk mengukur variasi dalam perbedaan antara interval RR yang berdekatan. SDSD menghitung standar deviasi dari perbedaan antara interval RR yang berdekatan. SDSD digunakan sebagai indikator aktivitas sistem saraf otonom dan respons tubuh terhadap stres. SDSD juga digunakan untuk menganalisis bagaimana tubuh menanggapi stres dan bagaimana sistem saraf otonom bekerja untuk mengatur respons tubuh terhadap stres.</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">RMSSD </p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: Justify;">RMSSD adalah ukuran statistik yang digunakan untuk mengukur variasi dalam interval RR. RMSSD menghitung perbedaan antara interval RR berurutan dan kemudian menghitung akar rata-rata kuadrat dari perbedaan tersebut. RMSSD digunakan sebagai indikator aktivitas sistem saraf otonom, khususnya cabang parasympathetic. RMSSD juga digunakan sebagai dasar untuk menghitung skor HRV, yang memberikan informasi tentang kemampuan tubuh dalam menanggapi stres.</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:blue; font-size: 23px; text-align: left;">PNN50</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            new_title = '<p style="font-family:Georgia; color:black; font-size: 20px; text-align: Justify;">PNN50 adalah ukuran yang digunakan untuk mengukur persentase interval RR yang berbeda lebih dari 50 ms. PNN50 digunakan sebagai indikator aktivitas sistem saraf otonom, khususnya cabang parasympathetic. PNN50 juga digunakan sebagai indikator stres dan keseimbangan sistem saraf otonom. PNN50 dapat digunakan untuk menganalisis bagaimana tubuh menanggapi stres dan bagaimana sistem saraf otonom bekerja untuk mengatur respons tubuh terhadap stres.</p>'
            st.markdown(new_title, unsafe_allow_html=True)

   
    
    elif sub_selected == 'Frequency Domain analysis':
        selected3 = option_menu(None, ["Baseline", "Segmentation","Spektrum"], 
            menu_icon="cast", default_index=0, orientation="horizontal")
        if selected3 == "Baseline":
             # Plotting dengan Plotly
            n = np.arange(0, ptp, 1, dtype=int)
            fig = go.Figure(data=go.Scatter(x=n, y=bpm_rr_baseline, mode='lines'))
            fig.update_layout(
            title="TACHOGRAM",
            xaxis_title="n",
            yaxis_title="BPM",
            xaxis=dict(showline=True, showgrid=True),
            yaxis=dict(showline=True, showgrid=True)
            )
            st.plotly_chart(fig)
        elif selected3 == 'Segmentation':
             optimizer_options3 = ['', 'Data 0-50', 'Data 50-100', 'Data 101-151', 'Data 151-201','Data 201-251','Data 251-301','Data 301-351','FFT TOTAL']
             selected_optimizer3 = st.selectbox('Segmentation', optimizer_options3)
             if selected_optimizer3 == 'Data 0-50':
                fig = go.Figure(data=go.Scatter(x=n_subset, y=bpm_rr_baseline_subset, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 0-49)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                fig1 = go.Figure(data=go.Scatter(x=n_subset, y=bpm_rr_baseline_windowed, mode='lines'))
                fig1.update_layout(
                title="TACHOGRAM (Data 0-49) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig1)
        
               # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half, y=np.abs(fft_result_half),mode="lines"))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 50-100':
                # Plotting dengan Plotly untuk subset data 50:100
                fig = go.Figure(data=go.Scatter(x=n_subset1, y=bpm_rr_baseline_subset1, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 50-101)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset1, y=bpm_rr_baseline_windowed1, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 50-100) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)

                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half1, y=np.abs(fft_result_half1), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 50:100",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 101-151':
        
                # Plotting dengan Plotly untuk subset data 101:151
                fig = go.Figure(data=go.Scatter(x=n_subset2, y=bpm_rr_baseline_subset2, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 101-151)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset2, y=bpm_rr_baseline_windowed2, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 101-151) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)

                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half2, y=np.abs(fft_result_half2), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 101:151",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 151-201':
                # Plotting dengan Plotly untuk subset data 151:201
                fig = go.Figure(data=go.Scatter(x=n_subset3, y=bpm_rr_baseline_subset3, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 151-201)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
                
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset3, y=bpm_rr_baseline_windowed3, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 151-201) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half3, y=np.abs(fft_result_half3), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 151:201",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 201-251':
        
                # Plotting dengan Plotly untuk subset data 201:251
                fig = go.Figure(data=go.Scatter(x=n_subset4, y=bpm_rr_baseline_subset4, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 201-251)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
                
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset4, y=bpm_rr_baseline_windowed4, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 201-251) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half4, y=np.abs(fft_result_half4), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 201:251",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 251-301':

                # Plotting dengan Plotly untuk subset data 251:301
                fig = go.Figure(data=go.Scatter(x=n_subset5, y=bpm_rr_baseline_subset5, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 251-301)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
                
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset5, y=bpm_rr_baseline_windowed5, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 251-301) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half5, y=np.abs(fft_result_half5), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 251:301",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'Data 301-351':
                # Plotting dengan Plotly untuk subset data 301:351
                fig = go.Figure(data=go.Scatter(x=n_subset6, y=bpm_rr_baseline_subset6, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 301-351)",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
                
                #  windowing Hamming
                fig = go.Figure(data=go.Scatter(x=n_subset6, y=bpm_rr_baseline_windowed6, mode='lines'))
                fig.update_layout(
                title="TACHOGRAM (Data 301-351) with Hamming Window",
                xaxis_title="n",
                yaxis_title="BPM",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig)
        
                # Membuat grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_half6, y=np.abs(fft_result_half6), mode='lines'))
                fig_fft.update_layout(
                title="FFT of TACHOGRAM 301:351",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
             if selected_optimizer3 == 'FFT TOTAL':
        
                # Membuat Total grafik FFT
                fig_fft = go.Figure(data=go.Scatter(x=fft_freq_total, y=np.abs(fft_result_total), mode='lines'))
                fig_fft.update_layout(
                title="FFT tOTAL Of TACHOGRAM",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Magnitude",
                xaxis=dict(showline=True, showgrid=True),
                yaxis=dict(showline=True, showgrid=True)
                )
                st.plotly_chart(fig_fft)
        if selected3 == "Spektrum":
                fig = go.Figure()
        
                fig.add_trace(go.Scatter(
                x=x_vlf,
                y=y_vlf,
                fill='tozeroy',
                fillcolor='rgba(166, 81, 216, 0.2)',
                line=dict(color='rgba(166, 81, 216, 0.5)'),
                name='VLF'
                ))
        
              
                fig.add_trace(go.Scatter(
                x=x_lf,
                y=y_lf,
                fill='tozeroy',
                fillcolor='rgba(81, 166, 216, 0.2)',
                line=dict(color='rgba(81, 166, 216, 0.5)'),
                name='LF'
                ))
        
             
                fig.add_trace(go.Scatter(
                x=x_hf,
                y=y_hf,
                fill='tozeroy',
                fillcolor='rgba(216, 166, 81, 0.2)',
                line=dict(color='rgba(216, 166, 81, 0.5)'),
                name='HF'
                ))
        
             
                fig.update_layout(
                title="FFT Spectrum (Welch's periodogram)",
                xaxis_title="Frequency (Hz)",
                yaxis_title="Density",
                xaxis=dict(range=[0, 0.5]),
                yaxis=dict(range=[0, max(np.abs(fft_result_total))]),
                legend=dict(x=0.8, y=0.95)
               )
                st.plotly_chart(fig)
                
                data = {
                "Metric": ["Total Power (TP)", "VLF", "LF", "HF", "LF/HF"],
                "Value": [total_power, VLF, LF_norm, HF_norm, LF_HF]
                }
                df = pd.DataFrame(data)
                fig = go.Figure(data=[go.Table(
                  header=dict(values=list(df.columns),
                        fill_color='paleturquoise',
                        align='left'),
                  cells=dict(values=[df.Metric, df.Value],
                       fill_color='lavender',
                       align='left'))
                   ])
                st.plotly_chart(fig)
                
                categories = ['Total Power (TP)', 'VLF', 'LF', 'HF']
                values = [total_power, VLF, LF_norm, HF_norm]
        
                fig = go.Figure()
        
                fig.add_trace(go.Bar(
                x=categories,
                y=values,
                marker_color=['blue', 'orange', 'green', 'red']
                ))
        
               
                fig.update_layout(
                title='Bar Series dari VLF, LF, HF',
                xaxis_title='Kategori',
                yaxis_title='Nilai'
                )
                st.plotly_chart(fig)
        
                def determine_category(LF_norm, HF_norm, LF_HF):
                    if LF_norm < 0.2 and HF_norm < 0.2:
                        return 1  # Low - Low
                    elif LF_norm >= 0.2 and LF_norm <= 0.6 and HF_norm < 0.2:
                        return 2  # Normal - Low
                    elif LF_norm > 0.6 and HF_norm < 0.2:
                        return 3  # High - Low
                    elif LF_norm < 0.2 and HF_norm >= 0.2 and HF_norm <= 0.6:
                        return 4  # Low - Normal
                    elif LF_norm >= 0.2 and LF_norm <= 0.6 and HF_norm >= 0.2 and HF_norm <= 0.6:
                        return 5  # Normal - Normal
                    elif LF_norm > 0.6 and HF_norm >= 0.2 and HF_norm <= 0.6:
                        return 6  # High - Normal
                    elif LF_norm < 0.2 and HF_norm > 0.6:
                        return 7  # Low - High
                    elif LF_norm >= 0.2 and LF_norm <= 0.6 and HF_norm > 0.6:
                        return 8  # Normal - High
                    elif LF_norm > 0.6 and HF_norm > 0.6:
                        return 9  # High - High
                    else:
                        return 0  # Undefined
                
                
                st.title("Autonomic Balance Diagram")
                
                category = determine_category(LF_norm, HF_norm, LF_HF)
                st.write("Category:", category)
                
                
                data = [
                    [7, 8, 9],
                    [4, 5, 6],
                    [1, 2, 3]
                 ]
                
                coordinates = {
                    1: (2, 0),
                    2: (2, 1),
                    3: (2, 2),
                    4: (1, 0),
                    5: (1, 1),
                    6: (1, 2),
                    7: (0, 0),
                    8: (0, 1),
                    9: (0, 2)
                  }
        # Create heatmap with Plotly Express
                fig = px.imshow(data, labels=dict(x="Sympathetic Level", y="Parasympathetic Level"), x=["Low", "Normal", "High"], y=["High", "Normal", "Low"])
        
        # Mark category on the heatmap
                coord = coordinates.get(category, None)
                if coord:
                     fig.add_shape(
                         type="circle",
                         xref="x",
                         yref="y",
                         x0=coord[1],
                         y0=coord[0],
                         x1=coord[1] + 0.5,  
                         y1=coord[0] + 0.5,  
                        line_color="black"
                    )
        
        
        # Add annotations for numbers
                annotations = []
                for i, row in enumerate(data):
                    for j, val in enumerate(row):
                        annotations.append(dict(
                        x=j, y=i, text=str(val), showarrow=False,
                        font=dict(color="black", size=16)
                        ))
        
                fig.update_layout(
                title="Autonomic Balance Diagram",
                annotations=annotations
                )
                fig.update_xaxes(ticks="outside", tickvals=[0, 1, 2])
                fig.update_yaxes(ticks="outside", tickvals=[0, 1, 2])
        
        # Display heatmap in Streamlit
                st.plotly_chart(fig)
                        
    elif sub_selected == 'Non Liniear analysis':
        selected3 = option_menu(None, ["Detrended Fluctuation Analysis", "Poincare Plot Analysis", "Sample Entropy"], 
            menu_icon="cast", default_index=0, orientation="horizontal")
        if selected3 == "Detrended Fluctuation Analysis":
             new_title = '<p style="font-family:Georgia; color:black; font-size: 18px;">Detrended Fluctuation Analysis</p>'
             st.markdown(new_title, unsafe_allow_html=True)
             fig = go.Figure()
             fig.add_trace(go.Scatter(x=list(range(1, len(dfa_ydata) + 1)), y=dfa_ydata, mode='lines+markers'))
             fig.update_layout(title='Detrended Fluctuation Analysis (DFA)',
                               xaxis_title='n',
                               yaxis_title='Fluctuation')
             st.plotly_chart(fig)
        elif selected3 == "Poincare Plot Analysis":
            new_title = '<p style="font-family:Georgia; color:black; font-size: 18px;">Poincare Plot Analysis</p>'
            st.markdown(new_title, unsafe_allow_html=True)
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=poincare_x, y=poincare_y, mode='markers'))
            fig.update_layout(title='Poincare Plot',
                              xaxis_title='RR(n)',
                              yaxis_title='RR(n+1)')
            st.plotly_chart(fig)
            data = {
                "Measures": ["SD1", "SD2"],
                "Hasil": [sd1, sd2]
            }
            df = pd.DataFrame(data)

            fig = go.Figure(data=[go.Table(
                columnwidth=[80, 20],
                header=dict(values=list(df.columns),
                            fill_color='red',
                            align='left',
                            line_color='darkslategray',
                            height=30),
                cells=dict(values=[df["Measures"], df["Hasil"]],
                           fill_color='white',
                           align='left',
                           line_color='darkslategray',
                           height=25,
                           font_size=12,
                           ),
            )])

            fig.update_layout(
                width=400,
                height=200,
                margin=dict(l=10, r=10, t=10, b=10)
            )

            st.plotly_chart(fig)
        elif selected3 == "Sample Entropy":
             new_title = '<p style="font-family:Georgia; color:black; font-size: 18px;">Sample Entropy</p>'
             st.markdown(new_title, unsafe_allow_html=True)
             new_title = f'<p style="font-family:Georgia; color:black; font-size: 18px;">{sample_entropy}</p>'
             st.markdown(new_title, unsafe_allow_html=True)


if selected == "DWT":
        sub_selected1 = st.sidebar.radio(
        "",
        ["Filter Coef", "Mallat", "Filter Bank "],
        index=0
    )
        if sub_selected1 == 'Filter Coef':
            optimizer_options4 = [' ', 'h(n)', 'g(n)', 'Hw', 'Gw','Qj(f)']
            selected_optimizer4 = st.selectbox('Penurunan inverse fourier transform', optimizer_options4)
            if selected_optimizer4 == 'h(n)':
                fig = go.Figure(data=[go.Bar(x=n_list, y=h)])
                fig.update_layout(title='h(n) Plot', xaxis_title='n', yaxis_title='g(n)',template='plotly_dark')
                st.plotly_chart(fig)
            if selected_optimizer4 == 'g(n)':
                fig = go.Figure(data=[go.Bar(x=n_list, y=g)])
                fig.update_layout(title='g(n) Plot', xaxis_title='n', yaxis_title='g(n)',template='plotly_dark')
                st.plotly_chart(fig)
            if selected_optimizer4 == 'Hw':
                fig = go.Figure(data=go.Scatter(x=i_list, y=Hw[:len(i_list)]))
                fig.update_layout(title='Hw Plot', xaxis_title='i', yaxis_title='Gw',template='plotly_dark')
                st.plotly_chart(fig)
            if selected_optimizer4 == 'Gw':
                fig = go.Figure(data=go.Scatter(x=i_list, y=Gw[:len(i_list)]))
                fig.update_layout(title='Gw Plot', xaxis_title='i', yaxis_title='Gw',template='plotly_dark')
                st.plotly_chart(fig)
            if selected_optimizer4 == 'Qj(f)':
                traces = []
                for i in range(1, 9):
                 trace = go.Scatter(x=i_list, y=Q[i], mode='lines', name=f'Q[{i}]')
                 traces.append(trace)
                
                
                 layout = go.Layout(title='Qj (f)',
                                   xaxis=dict(title=''),
                                   yaxis=dict(title=''),
                                   template='plotly_dark'
                )
                
                
                 fig = go.Figure(data=traces, layout=layout)
                
                
                 st.plotly_chart(fig)
        if sub_selected1 == 'Mallat':
            optimizer_options5 = ['', 'Delay', 'w2fm','s2fm','gabungan']
            selected_optimizer5 = st.selectbox('', optimizer_options5)
            if selected_optimizer5 == 'Delay':
                data = {
                    "": ["T1", "T2", "T3","T4","T5"],
                    "Hasil": [T1, T2, T3,T4,T5]
                }
                df = pd.DataFrame(data)
                
                # Buat tabel menggunakan Plotly
                fig = go.Figure(data=[go.Table(
                    columnwidth=[80, 20],  # Set column width
                    header=dict(values=list(df.columns),
                                fill_color='red',  # Ubah warna header menjadi merah
                                align='left',
                                line_color='darkslategray',
                                height=30),  # Set header height
                    cells=dict(values=[df[""], df["Hasil"]],
                               fill_color='white',  # Ubah warna sel menjadi merah
                               align='left',
                               line_color='darkslategray',
                               height=25,  # Set cell height
                               font_size=12,  # Set font size
                               ),
                )])
                
                # Set layout to adjust the table size
                fig.update_layout(
                    width=800,
                    height=200,
                    margin=dict(l=10, r=10, t=10, b=10)
                )
                
                # Tampilkan tabel
                st.plotly_chart(fig)
                
                data = {
                    "": ["Delay1","Delay2","Delay3","Delay4","Delay5"],
                    "Hasil": [Delay1,Delay2,Delay3,Delay4,Delay5]
                }
                df = pd.DataFrame(data)
                
                # Buat tabel menggunakan Plotly
                fig = go.Figure(data=[go.Table(
                    columnwidth=[80, 20],  # Set column width
                    header=dict(values=list(df.columns),
                                fill_color='red',  # Ubah warna header menjadi merah
                                align='left',
                                line_color='darkslategray',
                                height=30),  # Set header height
                    cells=dict(values=[df[""], df["Hasil"]],
                               fill_color='white',  # Ubah warna sel menjadi merah
                               align='left',
                               line_color='darkslategray',
                               height=25,  # Set cell height
                               font_size=12,  # Set font size
                               ),
                )])
                
                # Set layout to adjust the table size
                fig.update_layout(
                    width=800,
                    height=200,
                    margin=dict(l=10, r=10, t=10, b=10)
                )
                
                # Tampilkan tabel
                st.plotly_chart(fig)
            
            if selected_optimizer5 == 'w2fm':
                # Function to create and display a plot for a given series
                def create_plot(n_values, series, index, series_name):
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=n_values, y=series, mode='lines', name=f'{series_name}[{index+1},n]'))
                    fig.update_layout(
                        title=f'{series_name}({index+1})f'
                    )
                    st.plotly_chart(fig)
                
                # Create and show separate plots for each w2fm series
                for i in range(5):
                    create_plot(n_values, w2fm_values[i], i, 'w2fm')
                    
            if selected_optimizer5 == 's2fm':
                # Function to create and display a plot for a given series
                def create_plot(n_values, series, index, series_name):
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=n_values, y=series, mode='lines', name=f'{series_name}[{index+1},n]'))
                    fig.update_layout(
                        title=f'{series_name}({index+1})f'
                    )
                    st.plotly_chart(fig)
                
                # Create and show separate plots for each w2fm series
                for i in range(5):
                    create_plot(n_values, s2fm_values[i], i, 's2fm')
            if selected_optimizer5 == 'gabungan':
                n_values = np.arange(min_n, max_n + 1)
                for i in range(0, 5):
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=n_values, y=w2fm[i], mode='lines', name=f'w2fm {i+1}'))
                    fig.add_trace(go.Scatter(x=n_values, y=s2fm[i], mode='lines', name=f's2fm {i+1}'))
                    fig.update_layout(
                        title=f'w2fm and s2fm ({i+1})f',
                        template='plotly_dark'
                    )
                    st.plotly_chart(fig)




        


    
    



        
    






 


        
        






        





        
        
    




    


         
