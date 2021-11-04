from io import StringIO
from Bio import SeqIO
import pandas as pd
import streamlit as st
from PIL import Image
import FeatureExtractor as fe
from keras.layers import LSTM
from keras.models import Sequential
from keras.layers.core import Dense, Dropout
import numpy as np

icon = Image.open('fav.png')
st.set_page_config(page_title='ORI-Deep', page_icon = icon)

def seqValidator(seq):
    for i in range(len(seq)):
        if (seq[i] != 'A' and seq[i] != 'G' and seq[i] != 'C' and seq[i] != 'T' and seq[i] != 'a' and seq[i] != 'g' and seq[i] != 'c' and seq[i] != 't'):
            return False
    return True

def createModel():
    model = Sequential()
    model.add(LSTM(100, input_shape=(822,1)))
    model.add(Dropout(0.5))
    model.add(Dense(512, activation='relu'))
    model.add(Dense(256, activation='relu'))
    model.add(Dense(128, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    model.load_weights("model.h5")
    return model

final_df = pd.DataFrame(columns=['Sequence ID', 'Sequence','Indices','Label'])
seq = ""
len_seq = 0
image = Image.open('WebPic.png')
st.subheader("""ORI-Deep""")
st.image(image, use_column_width=True)
st.sidebar.subheader(("Input Sequence(s) of Length 300 or greater (FASTA FORMAT ONLY)"))
fasta_string  = st.sidebar.text_area("Sequence Input", height=200)
            
st.subheader("Click the Example Button for Sample Data")

if st.button('Example'):
    st.info("ORI Sequences")
    st.code(">Seq1\nTGGAGGTACCTGCTCTCCACGCAGTACAGCATCTGGCTGTGTGGTGGCTGCAGTCCAGTGGGTCTGAGAGGGGTGGAGGGGGAAGGTGGAGAAGGAGAGGAGGGAGGAGGGGAGGGAGGGAAAGGTGAGAGGGGAAGACGAAGGGAAGGAGAAGAGGTGGGGGAGGGGGCAGGGGGGCAGCATTGGCCAGGTAAGGGCGTGTCCACTGTGAGGCACCGTGGGATTCCTGGGCGAGGACTGTGAGACACTGGTGAGCGTGTGTGCAGGGCTCTTTCGGGCTTCCTCTGTTTTTAATGACAT", language="markdown")
    st.code(">Seq2\nCTGGTGCCCCTTGCTTCCCCCATTCTGGAGCCACATATATACCCTACTCTTCTTCAGAACCTGTCCTCAGTGAAGGTTCCTGACCCTCAGGGATACATCCATGACTAATGGCCTTATTTGTTGCAAGAAAAGAAGTACTCCAATGGAAAGAGATCTTTGAGGCCATCTCCCTGAGTATCCCCATGTTATAGTTGAGTAATTTGAAGGTTGGGAGTCTCAAGCCTGAGAGGAGAATTGCATGTGGTGATTCTGACTACAGACCCTGGAGGCAGACTGCCTGGCTTGGAGTCCTGGCTCATC", language="markdown")
    st.info("Non-ORI Sequences")
    st.code(">Seq1\nTGTTGGGCACTGTGCTGGGGCATTTGGATATGATTTTACTTCAAGCCTGGGTAGTAGCAATGATTCCATTTATCGGCACTGGTAATGCGGCCCAGAGAGGTGAGGTGACTTGCCCCAGGCTACACAGCTGCTTTGCCTGACTGCACAGCCCAAGTCTTCCCATCTTCCTTGGCTGCCTCCCTTTGGCCTTTGGCCCCCCACTGAAGCCAGGTTTAGAGCCGCACAGGCCTCTTGCCACCCTGGCGGCTGTCGGGGGAAGCGGTTAGCGGTGATGGACCGTCTAGGACCCTTGGGCCGCCC", language="markdown")
    st.code(">Seq2\nCTCAGGACCTAGCATTGTCTCTTTCTTTGTTCTTAATAGGTCTCCTCAACCACATTGCAAATGCCATGAAGGCAGGATTGAGTCTATCTTGATCACATGTGTACTCCTCGTTCCTAGAACAGCACCTCTAGGTGCTTAATAAAGCACGCAGTCCAAACACAGTTGTGTTATGAGTGATGGGACAAATGAAGAGCTCCCAGCTTGGAGGATGATGTTGGAGGACATTAGAGGAGGATGGCAAGCACTCAGTCAACCAGCCTATGGGTACAGGGGTGGTTTAAGAAGTCATGGAGGAGGGGC", language="markdown")


if st.sidebar.button("SUBMIT"):
    if(fasta_string==""):
        st.info("Please input the sequence first. Length should be 300 or greater.")
    fasta_io = StringIO(fasta_string) 
    records = SeqIO.parse(fasta_io, "fasta") 
    for rec in records:
        seq_id = str(rec.id)
        seq=str(rec.seq)
        if(seqValidator(seq)):
            len_seq = len(seq)
            if (len_seq < 300):
                st.info("Please input the sequence again. Length should be 300 or greater. Currently length is " + str(len_seq))
            elif (len_seq == 300):
                df_temp = pd.DataFrame([[seq_id, seq,'Complete(1-300)','None']], columns=['Sequence ID', 'Sequence','Indices','Label'] )
                final_df = pd.concat([final_df,df_temp], ignore_index=True)
            else:
                n_seqs = len_seq - 300
                for i in range(n_seqs + 1):
                    sub_seq = seq[i: i+300]
                    df_temp = pd.DataFrame([[seq_id, sub_seq,str(i+1)+'-'+str(i+300),'None']], columns=['Sequence ID', 'Sequence','Indices','Label'] )
                    final_df = pd.concat([final_df,df_temp], ignore_index=True)
        else:
            st.info("Sequence with Sequence ID: " + str(seq_id) + " is invalid, containing letters other than A,G,C,T.")
    fasta_io.close()
    if(final_df.shape[0]!=0):
        model = createModel()
        for iter in range(final_df.shape[0]):
            temp_seq =  final_df.iloc[iter, 1]
            fv_array = fe.extractFeature(temp_seq)
            score = model.predict(fv_array)
            pred_label = np.round_(score, decimals=0, out=None)
            if(pred_label==1):
                pred_label="ORI"
            else:
                pred_label="Non-ORI"
            final_df.iloc[iter, 3] = str(pred_label)

    st.dataframe(final_df)