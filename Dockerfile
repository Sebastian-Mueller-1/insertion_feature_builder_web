# Use a specific version of the Python image for consistency
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    bedtools \
    gcc \
    g++ \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Set environment variables
# Set to "production" when deploying to a live environment
ENV FLASK_ENV=development
ENV FLASK_DEBUG=1
# This ensures Python output is set straight to the terminal without being first buffered
ENV PYTHONUNBUFFERED=1

# Set the working directory inside the container
WORKDIR /app

# Copy the application's requirements file and install dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Gunicorn will be added to the requirements.txt
# Add Gunicorn, and any other additional requirements you might have
RUN pip install gunicorn

# Copy the rest of the application's code
COPY . /app

# Start Gunicorn
# Adjust the number of workers and threads as per your application's requirement
# "server:app" should match the naming convention of your Flask application script and app instance
CMD ["gunicorn", "--timeout", "6000", "-w", "4", "-b", "0.0.0.0:5000", "server:app"]
