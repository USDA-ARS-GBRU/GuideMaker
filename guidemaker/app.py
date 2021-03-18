import os
import subprocess
import base64
import streamlit as st
import pandas as pd
import pathlib
import uuid
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4
import shutil
import pandas as pd
import numpy as np
import altair as alt
from PIL import Image
import guidemaker




######

@contextmanager
def genome_connect(db_bytes):
    """Write input genome to local disk and clean after using
    """
    fp = Path(str(uuid4()))
    fp.write_bytes(db_bytes.getvalue())
    conn = str(fp)
    try:
        yield conn
    finally:
        fp.unlink()

@st.cache(suppress_st_warning=True)
def run_command(args):
    """Run command, transfer stdout/stderr back into Streamlit and manage error"""
    #st.info(f"Running:: '{' '.join(args)}'")
    result = subprocess.run(args, capture_output=True, text=True)
    try:
        result.check_returncode()
        #st.info(result.stdout)
        #st.text(result.stderr)    
    except subprocess.CalledProcessError as e:
        st.error(result.stderr)
        raise e

def get_binary_file_downloader_html(bin_file, file_label='File'):
    """Binary file downloader in html format
    """
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{file_label}</a>'
    return href

def read_markdown_file(markdown_file):
    """Read markdown file
    """
    return Path(markdown_file).read_text()

def GuideMakerPlot(df):
    """Returns guidemaker plot describing PAM targets
    """
    source = df
    brush = alt.selection(type='interval', encodings=['x'])
    binNum = int(round(source['Feature end'].max()/200,0))
    display_info = source.columns.tolist()

    # Feature density
    densityF = alt.Chart(source).transform_density(
    'Feature start',
    as_=['Feature start', 'Feature Density'],
    extent=[1, source['Feature end'].max()],
    bandwidth=binNum,
    ).mark_area(color='black',opacity=0.6).encode(
    x = alt.X('Feature start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
    y='Feature Density:Q',
    ).properties(height=50)

    # gRNA density
    densityG = alt.Chart(source).transform_density(
    'Guide start',
    as_=['Guide start', 'Guide Density'],
    extent=[1, source['Feature end'].max()],
    bandwidth=binNum,
    ).mark_area(color='pink',opacity=0.6).encode(
    x = alt.X('Guide start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
    y='Guide Density:Q',
    ).properties(height=50).add_selection(brush)

    # locus tag
    locus = alt.Chart(source).mark_bar(cornerRadiusTopLeft=3,cornerRadiusTopRight=3).encode(
    x='count(locus_tag):Q',
    y = alt.Y('locus_tag', axis=alt.Axis( title='Locus')),
    color='PAM:N',
    tooltip=display_info
    ).transform_filter(
    brush
    ).interactive().properties(height=500)

    gmplot = densityF & densityG & locus
    return( gmplot)


def main(arglist: list=None):
    """Run web App
    """
    subheader = "Globally design gRNAs for any CRISPR-Cas system in any small genome 🦠🧬"
    st.markdown(f'<h1 style="color: #003087">GuideMaker</h1>',unsafe_allow_html=True)
    st.markdown(f'<strong style="font-family:Hoefler Text;font-size: 18px;color: #FA4616">{subheader}</strong>',unsafe_allow_html=True)

    st.markdown("---")

    sessionID = str(uuid.uuid1())
    #st.write(sessionID)
    logfilename = sessionID+"_log.txt"


    # Create a downloads directory within the streamlit static asset directory
    # and write output files to it. Binary downloader as html file  has limited
    # file size for conversion. So, we need to write it to local folder.
    STREAMLIT_STATIC_PATH = pathlib.Path(st.__path__[0]) / 'static'
    DOWNLOADS_PATH = (STREAMLIT_STATIC_PATH / "downloads")
    if not DOWNLOADS_PATH.is_dir():
        DOWNLOADS_PATH.mkdir()

    # Define input paramters and widgets
    genome = st.sidebar.file_uploader("Upload a Genome file [ gbk, gbk.gz ]", type=["gbk","gz"])
    pam = st.sidebar.text_input("Input PAM Motif [ E.g. NGG ] ", "NGG")
    restriction_enzyme_list = st.sidebar.text_input("Restriction Enzymes list [ E.g. NGRT ] ", "NGRT")
    pam_orientation = st.sidebar.selectbox("PAM Orientation [ Options: 3prime, 5prime ]", ("3prime","5prime"))
    guidelength = st.sidebar.number_input('Guidelength [ Options: 10 - 27 ]', 10, 27)
    lsr = st.sidebar.number_input('Length of seed region[ Options: 0 - 27 ]', 0, 27)
    dist = st.sidebar.number_input('Hamming Distance [Options: 0 - 5 ]', 0, 5)
    before = st.sidebar.number_input('Before [Options: 1 - 500 ]', 1, 500)
    into = st.sidebar.number_input('Into [Options: 1 - 500 ]', 1, 500)
    threads = st.sidebar.number_input('Number of Threads [ Options: 2, 4, 6, 8]', 2, 8, step = 2)

    # name_cols = st.beta_columns(3)
    # First_name = name_cols[0].text_input("Name")
    # Middle_name = name_cols[1].text_input("Middle Name")
    # Last_name = name_cols[2].text_input("Last Name")

    # name2_cols = st.beta_columns(3)
    # First_name2 = name2_cols[0].text_input("Name2")
    # Middle_name2 = name2_cols[1].text_input("Middle Name2")
    # Last_name2 = name2_cols[2].text_input("Last Name2")

    # st.write(Last_name2)

    if genome:
        with genome_connect(genome) as conn:
            #st.write("Connection object:", conn)
            args = ["guidemaker",
            "-i",conn,
            "-p",pam,
            "--guidelength",str(guidelength),
            "--pam_orientation",pam_orientation,
            "--lsr",str(lsr),
            "--dist",str(dist),
            "--outdir", sessionID ,
            "--log", logfilename,
            "--into",str(into),
            "--before",str(before),
            "--restriction_enzyme_list", restriction_enzyme_list,
            "--threads",str(threads)]
            if(st.sidebar.button("SUBMIT")):
                #st.markdown("""🏃🏃🏃🏃🏃🏃🏃🏃""")
                run_command(args)

    if os.path.exists(sessionID):
        
        source = pd.read_csv(os.path.join("./", sessionID,'targets.csv'))
        source = pd.read_csv(os.path.join("./", sessionID,'targets.csv'),low_memory=False)

        accession_list = list(set(source['Accession']))
        for accession in accession_list:
            accession_df = source[source["Accession"] == accession]
            accession_info = f"**Accession:** {accession}"
            st.markdown(accession_info)
            st.write(GuideMakerPlot(accession_df))

        # Targets
        target_tab = "✅ [Target Data](downloads/targets.csv)"
        targets = pd.read_csv(os.path.join("./", sessionID,'targets.csv'),low_memory=False)
        targets.to_csv(str(DOWNLOADS_PATH / "targets.csv"), index=False)

        # Controls
        control_tab = "✅ [Control Data](downloads/controls.csv)"
        controls = pd.read_csv(os.path.join("./", sessionID,'controls.csv'),low_memory=False)
        controls.to_csv(str(DOWNLOADS_PATH / "controls.csv"), index=False)

        # logs
        with st.beta_expander("Results"):
            st.write(target_tab)
            st.write(control_tab)
            st.write(get_binary_file_downloader_html(logfilename, '✅ Log File'), unsafe_allow_html=True)
        
    else:
        pass

    # Parameters Dictionary 
    image = Image.open(guidemaker.APP_PARAMETER_IMG)
    optionals = st.beta_expander("Parameter Dictionary", False)
    optionals.image(image, caption='GuideMaker Parameters',use_column_width=True)

    with st.beta_expander("Designing Experiments with GuideMaker Results"):
        intro_markdown = read_markdown_file(guidemaker.APP_EXPERIMENT_FILE)
        st.markdown(intro_markdown, unsafe_allow_html=True)
        

    st.markdown("""
    ##### API documentation 📖

    API documentation for the module can be found [here](https://guidemaker.org)


    ##### License information ©️

    *As a work of the United State Government Department of Agriculture - Agricultural Research Service (USDA-ARS)
    this software is available under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication*

    """)

    ## Check if the output dir exist, if yes, delete so that previous results are not display
    try:
        shutil.rmtree(sessionID, ignore_errors=True)
        os.remove(logfilename)
    except FileNotFoundError as e:
        pass

if __name__ == "__main__":
    main()