# app/__init__.py

from flask import Flask

app = Flask(__name__, instance_relative_config=True)

from app import views, routes

app.testing = True
app.debug = True
app.config['SECRET_KEY'] = 'meq348vyojdc9p42micniorq93eakg'
app.config['UPLOAD_DIR'] = '/tmp'