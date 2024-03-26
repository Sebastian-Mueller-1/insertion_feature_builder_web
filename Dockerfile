FROM python:3.10-slim
RUN apt-get update && apt-get install -y \
    bedtools \
    gcc \
    g++ \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir -r requirements.txt
CMD ["python", "server.py"]
