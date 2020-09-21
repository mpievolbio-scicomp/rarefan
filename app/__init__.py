# app/__init__.py

from flask import Flask
import os

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')

# from werkzeug.debug import DebuggedApplication

from app import views, routes 

# app.wsgi_app = DebuggedApplication(app.wsgi_app, True)
app.testing = False
app.debug = False

app.config['SECRET_KEY'] = 'meq348vyojdc9p42micniorq93eakg'
app.config['UPLOAD_DIR'] = os.path.join(app.static_folder, 'uploads')