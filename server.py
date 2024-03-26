import os
import zipfile
from io import BytesIO

from flask import Flask, render_template, request, send_file
from werkzeug.utils import secure_filename
from waitress import serve

from insertion_feature_builder_v7 import process_data

app = Flask(__name__)
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)


@app.route("/", methods=["GET", "POST"])
def upload_and_process():
    if request.method == "POST":
        # List to hold uploaded file paths
        file_paths = []

        # Process each uploaded file
        for i in range(1, 6):
            file = request.files.get(f"file{i}")
            if file:
                filename = secure_filename(file.filename)
                save_path = os.path.join(UPLOAD_FOLDER, filename)
                file.save(save_path)
                file_paths.append(save_path)

        # Ensure we received all necessary files
        if len(file_paths) < 5:
            return "Error: Not all files were uploaded.", 400

        # Call the processing function from insertion_feature_builder_v7.py
        output_df, multiple_warnings_df, remove_df = process_data(*file_paths)

        # Prepare the DataFrames for download (example with zipping them together)
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for df, name in zip(
                [output_df, multiple_warnings_df, remove_df],
                ["output.csv", "warnings.csv", "remove.csv"],
            ):
                df_buffer = BytesIO()
                df.to_csv(df_buffer, index=False)
                df_buffer.seek(0)
                zip_file.writestr(name, df_buffer.getvalue())

        # Reset buffer's cursor to the beginning
        zip_buffer.seek(0)
        return send_file(
            zip_buffer,
            mimetype="application/zip",
            as_attachment=True,
            attachment_filename="processed_files.zip",
        )

    return render_template("upload.html")


if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=8000)
