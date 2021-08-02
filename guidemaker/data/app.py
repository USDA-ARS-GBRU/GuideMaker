"""guidemaker/app.py A configuration template fo a streamlit web application version of Guidemaker."""

import os
import subprocess
import base64
import pathlib
from pathlib import Path
import uuid
from uuid import uuid4
from contextlib import contextmanager
import shutil

import pandas as pd
import altair as alt
from PIL import Image
import streamlit as st
from streamlit_tags import st_tags_sidebar
from copy import deepcopy

import guidemaker

apptitle = 'GuideMaker'
st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.sidebar.markdown("## Select Parameters to Design gRNAs")



# @contextmanager
# def genome_connect(db_bytes):
#     """Write input genome to local disk and clean after using."""
#     fp = Path(str(uuid4()))
#     with open("input.gbk", 'wb') as fp:
#         for i in db_bytes:
#             if guidemaker.is_gzip(i):
#                 with zip.open(i, 'rt') as f:
#                     fp = bytearray(f)
#                     file.write(fp)
#             else:
#                 with open(i, 'r') as f:
#                     fp = bytearray(f)
#                     file.write(fp)
#     conn = str(fp)
#     try:
#         yield conn
#     finally:
#         pass


@contextmanager
def genome_connect(db_bytes):
    """Write input genome to local disk and clean after using."""
    fp = Path(str(uuid4()))
    with open('input.gbk', 'wb') as file:
        for i in db_bytes:
            fp = bytearray(i)
            file.write(fp)
    conn = str(fp)
    try:
        yield conn
    finally:
        pass

@st.cache(suppress_st_warning=True)
def run_command(args):
    """Run command, transfer stdout/stderr back into Streamlit and manage error."""
    st.info(f"Running:: '{' '.join(args)}'")
    result = subprocess.run(args, capture_output=True, text=True)
    try:
        result.check_returncode()
        # st.info(result.stdout)
        # st.text(result.stderr)
    except subprocess.CalledProcessError as e:
        st.error(result.stderr)
        raise e


def get_binary_file_downloader_html(bin_file, file_label='File'):
    """Binary file downloader in html format."""
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{file_label}</a>'
    return href


def read_markdown_file(markdown_file):
    """Read markdown file."""
    return Path(markdown_file).read_text()


def guidemakerplot(df):
    """Returns guidemaker plot describing PAM targets."""
    source = df
    brush = alt.selection(type='interval', encodings=['x'])
    binnum = int(round(source['Feature end'].max() / 200, 0))
    display_info = source.columns.tolist()

    # Feature density
    densityF = alt.Chart(source).transform_density(
    'Feature start',
    as_=['Feature start', 'Feature Density'],
    extent=[1, source['Feature end'].max()],
    bandwidth=binnum,
    ).mark_area(color='black', opacity=0.6).encode(
    x=alt.X('Feature start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
    y='Feature Density:Q',
    ).properties(height=50)

    # gRNA density
    densityg = alt.Chart(source).transform_density(
    'Guide start',
    as_=['Guide start', 'Guide Density'],
    extent=[1, source['Feature end'].max()],
    bandwidth=binnum,
    ).mark_area(color='pink', opacity=0.6).encode(
    x=alt.X('Guide start', axis=alt.Axis(title='Genome Coordinates (bp)', tickCount=5)),
    y='Guide Density:Q',
    ).properties(height=50).add_selection(brush)

    # locus tag
    locus = alt.Chart(source).mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
    x='count(locus_tag):Q',
    y=alt.Y('locus_tag', axis=alt.Axis(title='Locus')),
    color='PAM:N',
    tooltip=display_info
    ).transform_filter(
    brush
    ).interactive().properties(height=500)

    gmplot = densityF & densityg & locus
    return gmplot




def main(arglist: list = None):
    """Run web App."""
    header = "GuideMaker"
    subheader = "Software to design CRISPR-Cas guide RNA pools in non-model genomes 🦠 🧬"
    st.markdown(f'<strong style="font-family:Hoefler Text;font-size: 36px;color: #0021A5">{header}</strong>',
                unsafe_allow_html=True)
    st.markdown(
        f'<strong style="font-family:Hoefler Text;font-size: 18px;color: #FA4616">{subheader}</strong>', unsafe_allow_html=True)

    st.markdown("---")

    sessionID = str(uuid.uuid1())
    # st.write(sessionID)
    logfilename = sessionID + "_log.txt"

    # Create a downloads directory within the streamlit static asset directory
    # and write output files to it. Binary downloader as html file  has limited
    # file size for conversion. So, we need to write it to local folder.
    STREAMLIT_STATIC_PATH = pathlib.Path(st.__path__[0]) / 'static'
    DOWNLOADS_PATH = (STREAMLIT_STATIC_PATH / "downloads")
    if not DOWNLOADS_PATH.is_dir():
        DOWNLOADS_PATH.mkdir()

    # Define input parameters and widgets

    multiple_files = st.sidebar.file_uploader("Upload one or more Genome file [ .gbk, .gbk.gz]", type=[".gbk", ".gz"], accept_multiple_files=True)
    genome = list( map(lambda x: x.getvalue(), multiple_files))

    DemoGenome = st.sidebar.selectbox("OR Use Demo GBK",['Carsonella_ruddii.gbk.gz','Pseudomonas_aeruginosa.gbk.gz'])
    demo = open(DemoGenome,"rb")
    #st.write("You selected this option ",xx)

  

    pam = st.sidebar.text_input("Input PAM Motif [ E.g. NGG ] ", "NGG")
    restriction_enzyme_list = st_tags_sidebar(label = 'Restriction Enzymes[e.g. NGRT]:',
                                                                  text  = 'Enter to add more', 
                                                                  value = ['NGRT'])
    #restriction_enzyme_list= st.sidebar.text_input("Restriction Enzymes list [ E.g. NGRT ] ", "NGRT")
    pam_orientation = st.sidebar.selectbox(
        "PAM Orientation [ Options: 3prime, 5prime ]", ("3prime", "5prime"))
    guidelength = st.sidebar.number_input('Guidelength [ Options: 10 - 27 ]', 10, 27, value=20)
    lsr = st.sidebar.number_input('Length of seed region[ Options: 0 - 27 ]', 0, 27, value=10)
    dtype = st.sidebar.selectbox(
        "Type of Edit Distance [ Options: hamming, leven ]", ("hamming", "leven"))
    dist = st.sidebar.number_input('Hamming Distance [Options: 0 - 5 ]', 0, 5, value=2)
    before = st.sidebar.number_input('Before [Options: 1 - 500 ]', 1, 500, value=100, step=50)
    into = st.sidebar.number_input('Into [Options: 1 - 500 ]', 1, 500, value=200, step=50)
    knum = st.sidebar.number_input('Similar Guides[Options: 2 - 20 ]', 2, 20, value=3)
    controls = st.sidebar.number_input('Control RNAs', 1, 1000, value=10, step=100)
    #threads = st.sidebar.number_input('Threads [ Options: 2, 4, 6, 8]', 2, 8, step=2)
    filter_by_locus = st_tags_sidebar(label = 'Filter by Locus Tag [e.g. NGRT]:', text  = 'Enter to add more')
                                                            
    filter_by_gene = st_tags_sidebar(label = 'Filter by Gene Name [e.g. NGRT]:', text  = 'Enter to add more')
                                                                  
    filter_by_genome_cordinate_start = st.sidebar.number_input('Genome Coordinate Start Position', 0, 500000000, value=0, step=1000)
    filter_by_genome_cordinate_end = st.sidebar.number_input('Genome Coordinate End Position', 0, 500000000, value=0, step=1000)


    # st.write(genome)
    # st.write(type(genome))
  

    if demo and "input.gbk":
        with genome_connect(demo):
            #st.write("Connection object:", conn)
            args = ["guidemaker",
            "-i", "input.gbk",
            "-p", pam,
            "--guidelength", str(guidelength),
            "--pam_orientation", pam_orientation,
            "--lsr", str(lsr),
            "--dtype", str(dtype),
            "--dist", str(dist),
            "--outdir", sessionID,
            "--log", logfilename,
            "--into", str(into),
            "--before", str(before),
            "--knum", str(knum),
            "--controls", str(controls),
            "--threads", str(2),
            "--restriction_enzyme_list"]
            scriptorun = args + restriction_enzyme_list

    if genome and "input.gbk":
        with genome_connect(genome):
            #st.write("Connection object:", conn)
            args = ["guidemaker",
            "-i", "input.gbk",
            "-p", pam,
            "--guidelength", str(guidelength),
            "--pam_orientation", pam_orientation,
            "--lsr", str(lsr),
             "--dtype", str(dtype),
            "--dist", str(dist),
            "--outdir", sessionID,
            "--log", logfilename,
            "--into", str(into),
            "--before", str(before),
            "--knum", str(knum),
            "--controls", str(controls),
            "--threads", str(2),
            "--restriction_enzyme_list"]
            scriptorun = args + restriction_enzyme_list

    if(st.sidebar.button("SUBMIT")):
        # st.markdown("""🏃🏃🏃🏃🏃🏃🏃🏃""")
        run_command(scriptorun)
    

    if os.path.exists(sessionID):

        #source = pd.read_csv(os.path.join("./", sessionID,'targets.csv'))
        source = pd.read_csv(os.path.join("./", sessionID, 'targets.csv'), low_memory=False)

        accession_list = list(set(source['Accession']))
        for accession in accession_list:
            accession_df = source[source["Accession"] == accession]
            accession_info = f"**Accession:** {accession}"
            st.markdown(accession_info)
            st.write(guidemakerplot(accession_df))


        # Targets
        target_tab = "✅ [Target Data](downloads/targets.csv)"
        targets = pd.read_csv(os.path.join("./", sessionID, 'targets.csv'), low_memory=False)
        targets.to_csv(str(DOWNLOADS_PATH / "targets.csv"), index=False)

        # Controls
        control_tab = "✅ [Control Data](downloads/controls.csv)"
        controls = pd.read_csv(os.path.join("./", sessionID, 'controls.csv'), low_memory=False)
        controls.to_csv(str(DOWNLOADS_PATH / "controls.csv"), index=False)

        # logs
        with st.beta_expander("Results"):
            st.write(target_tab)
            st.write(control_tab)
            st.write(get_binary_file_downloader_html(
                logfilename, '✅ Log File'), unsafe_allow_html=True)
    else:
        pass

    # Parameters Dictionary
    image = Image.open(guidemaker.APP_PARAMETER_IMG)
    optionals = st.beta_expander("Parameter Dictionary", False)
    optionals.image(image, caption='GuideMaker Parameters', use_column_width=True)

    with st.beta_expander("Designing Experiments with GuideMaker Results"):
        intro_markdown = read_markdown_file(guidemaker.APP_EXPERIMENT_FILE)
        st.markdown(intro_markdown, unsafe_allow_html=True)

    st.markdown("""
    ##### API documentation 📖

    API documentation for the module can be found [here](https://github.com/USDA-ARS-GBRU/GuideMaker)


    ##### License information ©️

    *Guidemaker was created by the United States Department of Agriculture - Agricultural Research Service (USDA-ARS). As a work of the United States Government this software is available under the CC0 1.0 Universal Public Domain Dedication (CC0 1.0)*

    """)
    

    # Check if the output dir exist, if yes, delete so that previous results are not display
    try:
        shutil.rmtree(sessionID, ignore_errors=True)
        os.remove(logfilename)
        os.remove('input.gbk')
    except FileNotFoundError as e:
        pass

if __name__ == "__main__":
    main()