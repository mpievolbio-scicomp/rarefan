# config.py
import os

class Config(object):
    SECRET_KEY = 'meq348vyojdc9p42micniorq93eakg'
    MONGODB_SETTINGS = {
        'db': 'rarefan',
        'host': 'localhost',
        'port': 27017,
        'username': 'rarefan',
        'password': '!rarefan$',
        'authentication_source': 'admin',
    }
    REDIS_URL = os.environ.get("REDIS_URL") or 'redis://'
    MAIL_SERVER = 'zimbra.evolbio.mpg.de'
    MAIL_USERNAME='rarefan@evolbio.mpg.de'
    MAIL_PASSWORD='7SaaZv34Xw5isyu'

    MAIL_USE_TLS=True
    MAIL_USE_SSL=False
    MAIL_PORT=25

    MAIL_DEBUG=False
    DEFAULT_MAIL_SENDER='rarefan@evolbio.mpg.de'
