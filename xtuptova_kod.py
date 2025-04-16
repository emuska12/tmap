# Assignment 2 Analysis and visualization of biological data
# Subject: BIOINF_B FIIT STU 2025
# GENE-X scientist: Ema Tuptová
# Date: 14. April 2025

############################################################################################################

# we found out that some biological samples from the metheorite react weirdly when in contact with H20
# hypoothesis: hydrofobic and hydrophilic patterns in the samples function strucrtures
# for a amino sqn 
#   - count hydropatic values of each aminoacids
#   - plot hydrophobicity scale
#   - mark hydrophobic regions (under treshold)
# if we know where the hydrophobic regions are, we can predict the structure and shape of the protein --> based on shape we can predict function
# hydrophilic (negative values) on the outside, hydrophobic (positive values) aminoacis will be in the inside of the protein (hiding from water) or near the oily membrane (tunnel through membrane where ions can pass in or out) 
# so we need to find regions where aminoacids have high hydrophobicity values --> Kyte-Doolittle scale
# if long hydrophobic region (around 20 aminoacids with high hydrophobicity score) it can be part of the protein that forms tunnel

############################################################################################################

import streamlit as st # UI
from Bio import SeqIO # parse FASTA
from io import StringIO # convert to correct format
import matplotlib.pyplot as plt # plot the graph
from datetime import datetime # timestamp

def calculate_hydro(sqn, N):
    hydro_scores = []
    for acid in sqn: # assign value to each amino acid
        hydro_scores.append(kd_scale[acid])

    smoothed = []
    # for each position, where we can start the window count mean of window
    for i in range(0, len(hydro_scores) - N + 1):
        total = 0
        for j in range(i, i + N):
            total += hydro_scores[j]
        smoothed.append(total / N)
    
    plot_hydro(smoothed, sampleID)

    detect_transmembrane_regions(smoothed, sample, sampleID, len(sample))


def plot_hydro(scores, sequence_id):
    graph = plt.figure()
    # relative positions of aminoacids in the sequence (0-1)
    positions = []
    for x in range(len(scores)):
        positions.append(x / (len(scores) - 1))
    
    # visualize score historgam
    plt.plot(positions, scores, color="blue", linewidth=0.8, label="KD skóre") 
    plt.plot([0, 1], [0, 0], color="black", linestyle="--", label="Neutralita")
    plt.plot([0, 1], [threshold, threshold], color="red", label="Threshold")
    plt.title("Hydrofilný profil - " + sequence_id)
    plt.xlabel("Relatívna pozícia aminokyselín v sekvencii")
    plt.ylabel("Hydropatická hodnota")
    plt.grid(True)

    # filling the area under the curve when score is above threshold = possible transmembrane region
    # color area between fill_y_points and threshold
    fill_y_points = []
    for score in scores:
        if score >= threshold: # use just a score as border
            fill_y_points.append(score) 
        else:
            fill_y_points.append(threshold) # use a threshold so the area below the threshold doesnt fill
    plt.fill_between(positions, fill_y_points, threshold, color="green", alpha=0.3, label="Potenciálne transmembránové oblasti")

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    st.pyplot(graph)
    plt.close(graph)


def detect_transmembrane_regions(scores, sequence, sequence_id, total_length):
    found_regions = []
    position = 0
    while position < len(scores):
        if scores[position] > threshold: # if the score is above threshold we can start looking for a region
            start = position # save the start index of the possible region
            while position < len(scores) and scores[position] > threshold:
                position = position + 1 # find the end of the region where the score is still above threshold
            end = position
            region_size = end - start + window_size - 1
            if region_size >= region_length:  # consider it a region only if it is long enough
                region_start = start
                region_end = min(region_start + region_size, len(sequence))
                sequence_piece = sequence[region_start:region_end]
                region_info = {
                    'start': region_start + 1,
                    'end': region_end,
                    'sequence': sequence_piece
                }
                found_regions.append(region_info)
        else:
            position += 1
    
    output = []
    output.append("#=======================================")
    output.append(f"# Program: tmap GENE-X")
    output.append(f"# Čas merania: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}\n")
    output.append(f"# Sekvencia: {sequence_id}\tod: 1\tdo: {total_length}")
    output.append(f"# HitCount: {len(found_regions)}\n")
    if found_regions:
        output.append(f"\tŠtart\tKoniec\tRegión")
        for region in found_regions:
            output.append(f"     {region['start']:>5}     {region['end']:>5}        {region['sequence']}")
    else:
        output.append("# Transmembránové oblasti nenájdené.")
    output.append("#=======================================")
    st.text_area("output", "\n".join(output), height=500)


# https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
kd_scale = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
    'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

st.title("🔬 Analýza a vizualizácia biologických sekvencií")
st.markdown("### Ema Tuptová")

# load the sequences as FASTA file or type it into textare
# used the code from previous assignment
sample = st.text_area("Sekvencia", ">sp|Q92952|KCNN1_HUMAN Small conductance calcium-activated potassium channel protein 1 OS=Homo sapiens OX=9606 GN=KCNN1 PE=1 SV=2\nMNSHSYNGSVGRPLGSGPGALGRDPPDPEAGHPPQPPHSPGLQVVVAKSEPARPSPGSPRGQPQDQDDDEDDEEDEAGRQRASGKPSNVGHRLGHRRALFEKRKRLSDYALIFGMFGIVVMVTETELSWGVYTKESLYSFALKCLISLSTAILLGLVVLYHAREIQLFMVDNGADDWRIAMTCERVFLISLELAVCAIHPVPGHYRFTWTARLAFTYAPSVAEADVDVLLSIPMFLRLYLLGRVMLLHSKIFTDASSRSIGALNKITFNTRFVMKTLMTICPGTVLLVFSISSWIIAAWTVRVCERYHDKQEVTSNFLGAMWLISITFLSIGYGDMVPHTYCGKGVCLLTGIMGAGCTALVVAVVARKLELTKAEKHVHNFMMDTQLTKRVKNAAANVLRETWLIYKHTRLVKKPDQARVRKHQRKFLQAIHQAQKLRSVKIEQGKLNDQANTLTDLAKTQTVMYDLVSELHAQHEELEARLATLESRLDALGASLQALPGLIAQAIRPPPPPLPPRPGPGPQDQAARSSPCRWTPVAPSDCG")
sample_file = st.file_uploader("Vybrať súbor", type=["fasta", "txt"])
# added components for user to play with settings
window_size = st.slider("Veľkosť sliding window", min_value=1, max_value=20, value=9)
threshold = st.slider("Threshold hodnota", min_value=-2.0, max_value=2.0, value=0.5, step=0.1)
region_length = st.slider("Min. počet aminokyselín tvoriaci region", min_value=1, max_value=50, value=18)

if sample_file:
    try:
        file_content = sample_file.getvalue().decode("utf-8")
        file_handle = StringIO(file_content)
        record = next(SeqIO.parse(file_handle, "fasta"))
        sample = str(record.seq)
        sampleID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri načítaní súboru: {e}")
        sample = ""
else:
    try:
        file_handle = StringIO(sample)
        record = next(SeqIO.parse(file_handle, "fasta"))
        sample = str(record.seq)
        sampleID = str(record.id)
    except Exception as e:
        st.error(f"Chyba pri spracovaní sekvencie: {e}")
        sample = ""

if st.button("Analyzovať"):
    if sample:
        st.title("📊 Výsledky")
        scores = calculate_hydro(sample, window_size)
        