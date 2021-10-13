FROM python:3.8-slim
WORKDIR /app 
COPY . . 
RUN pip install dash   
RUN pip install dash-html-components                                         
RUN pip install dash-core-components                                     
RUN pip install plotly
RUN pip install pandas
RUN pip install requests
RUN pip install numpy
RUN pip install matplotlib
RUN pip freeze > requirements.txt
EXPOSE 8050 
CMD ["python", "app.py"]