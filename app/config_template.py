# config.py
import os

class Config(object):
    SECRET_KEY = ''
    MONGODB_SETTINGS = {
        'db': '',
        'host': '',
        'port': ,
        'username': '',
        'password': '',
        'authentication_source': 'admin',
    }
    REDIS_URL = os.environ.get("REDIS_URL") or 'redis://'
    MAIL_SERVER = ''
    MAIL_USERNAME=''
    MAIL_PASSWORD=''

    MAIL_USE_TLS=True
    MAIL_USE_SSL=False
    MAIL_PORT=25

    MAIL_DEBUG=False
    DEFAULT_MAIL_SENDER=''
