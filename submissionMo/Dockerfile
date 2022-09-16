FROM python:3.8-slim-buster

RUN pip install numpy torch Counter joblib scikit-learn==1.0.2 pandas

RUN mkdir selected_feature

COPY Mtd_waspreterm.save /usr/local/bin/
COPY pty_waspreterm.save /usr/local/bin/
COPY txy_waspreterm.save /usr/local/bin/
COPY L2logistic_waspreterm.save /usr/local/bin/

COPY run_model.py /usr/local/bin/.
COPY selected_feature/* selected_feature

RUN chmod a+x /usr/local/bin/run_model.py

ENTRYPOINT ["python3", "/usr/local/bin/run_model.py"]