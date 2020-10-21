from flask import Blueprint
from flask_autoindex import AutoIndexBlueprint
import os
auto_bp = Blueprint('auto_bp', __name__)
AutoIndexBlueprint(auto_bp, browse_root=os.path.join(os.path.dirname(__file__), 'static', 'uploads'))