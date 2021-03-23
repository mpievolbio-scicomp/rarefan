# config.py
import os

class Config(object):
    SECRET_KEY = 'meq348vyojdc9p42micniorq93eakg'
 # Flask-Dropzone config:
    DROPZONE_MAX_FILE_SIZE = 100000
    DROPZONE_MAX_FILES = 10000
    DROPZONE_ENABLE_CSRF = True  # enable CSRF protection
    DROPZONE_UPLOAD_MULTIPLE = True
    DROPZONE_PARALLEL_UPLOADS = 50000
    DROPZONE_ALLOWED_FILE_CUSTOM = True
    DROPZONE_TIMEOUT = 3600000
    DROPZONE_ALLOWED_FILE_TYPE = '.faa, .fa, .fasta, .fas, .fasta, .fna, .fn, .nwk'