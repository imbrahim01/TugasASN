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




with st.sidebar:
    selected = option_menu("TUGAS 1", ["Home", "PAGE 1", "PAGE 3","PAGE 4","PAGE 5"], default_index=0)

if selected == "Home":
   st.title('Project ASN Kelompok 1')
   st.subheader("Anggota kelompok")
   st.text("Farhan Majid Ibrahim - 5023211049")
   st.text("Nayla Pramudhita Putri Pertama - 5023211012")
   st.text("Mohammad Rayhan Amirul Haq Siregar - 5023211045")
   new_title = '<p style="font-family:Georgia; color:#FF0000; font-size: 20px;">Reynard Prastya Savero - 5023211042</p>'
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
