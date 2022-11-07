import re
from jinja2 import evalcontextfilter
from markupsafe import Markup, escape


@evalcontextfilter
def nl2br(eval_ctx, value):
    """Converts newlines in text to HTML-tags"""
    result = "<br>".join(re.split(r'(?:\r\n|\r|\n)', escape(value)))

    if eval_ctx.autoescape:
        result = Markup(result)
    return result