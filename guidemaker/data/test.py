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

import guidemaker











multiple_files = st.sidebar.file_uploader("Upload one or more Genome file [ .gbk, .gz]", type=["gbk", ".gz"], accept_multiple_files=True)
genome = list( map(lambda x: x.getvalue(), multiple_files))
st.write(genome)
st.write(multiple_files)
