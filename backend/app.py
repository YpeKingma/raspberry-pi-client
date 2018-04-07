import os
from flask import Flask, send_from_directory

app = Flask(__name__, static_folder='../frontend/build')

# Serve React App
@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve(path):
    if(path == ""):
        return send_from_directory('../frontend/build', 'index.html')
    else:
        if(os.path.exists("../frontend/build/" + path)):
            return send_from_directory('../frontend/build', path)
        else:
            return send_from_directory('../frontend/build', 'index.html')

# Serve API backend
@app.route('/api/identity')
def identity():
    return 'I am Pi number 5'
