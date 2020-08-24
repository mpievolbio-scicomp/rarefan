from flask import render_template
from app.views import LoginForm
from app import app

@app.route('/')
def index():
    return render_template("index.html")


@app.route('/submit')
def submit():
    form = LoginForm()
    return render_template('submit.html', title='Submit', form=form)
