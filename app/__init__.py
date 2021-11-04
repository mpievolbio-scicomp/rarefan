from flask import Flask
from flask_wtf import CSRFProtect
from flask_mongoengine import MongoEngine, MongoEngineSessionInterface
from .config import Config
import os
from redis import Redis
import rq

import logging
from flask_debugtoolbar import DebugToolbarExtension
from flask_mail import Mail

logging.basicConfig(level=logging.DEBUG)

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')

app.testing = app.debug = True
# app.testing = app.debug = False

# # email
app.config.from_object(Config)
app.config['UPLOAD_DIR'] = upload_dir

db = MongoEngine()
app.session_interface = MongoEngineSessionInterface(db)
db.init_app(app)

mail = Mail(app)

app.redis = Redis.from_url(app.config['REDIS_URL'])
app.queue = rq.Queue('rarefan', connection=app.redis)

csrf = CSRFProtect(app)
from app import views, routes

