FROM ubuntu:bionic

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y python3.8
RUN apt-get -y install python3-pip
RUN pip3 install numpy==1.19.5 pandas==1.1.5 torch==1.10.2 scikit-learn==0.24.2

RUN mkdir selected_feature

COPY test_submission.save /usr/local/bin/
COPY Run_model.py /usr/local/bin/
COPY selected_feature/* selected_feature

RUN chmod a+x /usr/local/bin/Run_model.py

ENTRYPOINT ["python3", "/usr/local/bin/Run_model.py"]