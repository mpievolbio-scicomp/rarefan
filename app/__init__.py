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

logger = logging.getLogger('rarefan')
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(module)s: %(message)s')
timestamp = datetime.datetime.now().strftime(format="%Y%m%d-%H%M%S")
handler = logging.FileHandler("/tmp/rarefan.log".format(timestamp))
handler.setFormatter(formatter)
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

app = Flask(__name__, instance_relative_config=True, static_url_path='/static')
upload_dir = os.path.join(app.static_folder, 'uploads')

app.testing = app.debug = True

# email
app.config.from_object(Config)
app.config['UPLOAD_DIR'] = upload_dir

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
