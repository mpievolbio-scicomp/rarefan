import os
import rq
from .config import Config
from flask import Flask
from flask_wtf.csrf import CSRFProtect
from flask_mail import Mail
from flask_mongoengine import MongoEngine, MongoEngineSessionInterface
from redis import Redis
import datetime
import logging
import datetime

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')

app.testing = app.debug = False

# email
app.config.from_object(Config)
app.config['UPLOAD_DIR'] = upload_dir

# Logging
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(module)s: %(message)s')
handler = logging.FileHandler(filename="/tmp/rarefan.log")
handler.setFormatter(formatter)

# Configure root logger (this config will trickle down to all module loggers.)
app.logger.addHandler(handler)

app.logger.setLevel(logging.INFO)
if app.config['DEBUG']:
    app.logger.setLevel(logging.DEBUG)

if app.debug:
    app.logger.debug("****************** Debug mode is active ******************")


# csrf = CSRFProtect()
# csrf.init_app(app)

db = MongoEngine()
# app.session_interface = MongoEngineSessionInterface(db)
db.init_app(app)

mail = Mail(app)

app.redis = Redis.from_url(app.config['REDIS_URL'])
app.queue = rq.Queue('rarefan', connection=app.redis)


# Has to be the last import to avoid cyclic dependencies.
from app import views, routes

