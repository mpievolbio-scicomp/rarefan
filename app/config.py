# config.py
import os

class Config(object):
    SECRET_KEY = 'meq348vyojdc9p42micniorq93eakg',
    MONGODB_SETTINGS = {
        'db': 'rarefan',
        'host': 'localhost',
        'port': 27017
    }
