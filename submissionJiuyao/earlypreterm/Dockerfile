FROM python:3.8-slim-buster

RUN pip install numpy torch Counter joblib scikit-learn==1.0.2 pandas tsai

COPY learner_earlypreterm.pkl /usr/local/bin/

COPY run_model_earlypreterm.py /usr/local/bin/.

RUN chmod a+x /usr/local/bin/run_model_earlypreterm.py

ENTRYPOINT ["python3", "/usr/local/bin/run_model_earlypreterm.py"]