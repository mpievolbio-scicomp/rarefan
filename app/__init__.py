from flask import Flask
from flask_dropzone import Dropzone
from flask_wtf import CSRFProtect
from .config import Config
import os
import logging
logging.basicConfig(level=logging.DEBUG)

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')

app.testing = False
app.debug = False

# email
app.config.from_object(Config)
app.config['UPLOAD_DIR'] = upload_dir

dropzone = Dropzone(app)
csrf = CSRFProtect(app)
from app import views, routes
