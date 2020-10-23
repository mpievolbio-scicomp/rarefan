from flask import Flask
import os

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')

from app import views, routes 

app.testing = False
app.debug = False

app.config['SECRET_KEY'] = 'meq348vyojdc9p42micniorq93eakg'
app.config['UPLOAD_DIR'] = upload_dir