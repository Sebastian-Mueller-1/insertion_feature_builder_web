FROM python:3.10-slim
RUN apt-get update && apt-get install -y \
    bedtools \
    gcc \
    g++ \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*
ENV FLASK_ENV=development
ENV FLASK_DEBUG=1
ENV PYTHONUNBUFFERED=1
WORKDIR /app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install gunicorn
COPY . /app
CMD ["gunicorn", "--timeout", "6000", "-w", "4", "-b", "0.0.0.0:5000", "server:app"]
