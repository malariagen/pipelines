FROM python:3.7.2

RUN pip install --upgrade pip

COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY sample_select_variants.py .
COPY test_sample_select_variants.py .
COPY fixture ./fixture
