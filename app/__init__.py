# app/__init__.py

from flask import Flask
from .bp import auto_bp
import os

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')
app.register_blueprint(auto_bp, url_prefix='/files')

from app import views, routes 

app.testing = False
app.debug = False

app.config['SECRET_KEY'] = 'meq348vyojdc9p42micniorq93eakg'
app.config['UPLOAD_DIR'] = upload_dir